import math
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


BASES = ("A", "C", "G", "T")
BASE_TO_CH = {b: i for i, b in enumerate(BASES)}
ABIF_DIR_ENTRY_STRUCT = struct.Struct(">4sIHHIIII")


def parse_fastq(path):
    """
    Simple FASTQ parser for standard 4-line FASTQ.
    Assumes Sanger-style PHRED+33 unless you change phred_offset below.
    """
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not (seq and plus and qual):
                raise ValueError("Truncated FASTQ record")

            header = header.rstrip("\n\r")
            seq = seq.rstrip("\n\r")
            plus = plus.rstrip("\n\r")
            qual = qual.rstrip("\n\r")

            if not header.startswith("@"):
                raise ValueError(f"Bad FASTQ header: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Bad FASTQ separator: {plus!r}")
            if len(seq) != len(qual):
                raise ValueError(
                    f"Sequence/quality length mismatch: {len(seq)} vs {len(qual)}"
                )

            yield header[1:].split()[0], seq.upper(), qual


def phred_from_ascii(qual, phred_offset=33):
    return np.array([max(0, ord(c) - phred_offset) for c in qual], dtype=np.int16)


def q_to_amplitude(q):
    """
    Map PHRED score to peak amplitude.
    High quality -> taller/narrower peaks.
    Keeps amplitudes in a nice synthetic range.
    """
    q = np.asarray(q, dtype=np.float64)
    return 60.0 + 4.0 * np.clip(q, 0, 45)


def q_to_sigma(q):
    """
    Low quality -> broader peaks.
    """
    q = np.asarray(q, dtype=np.float64)
    return 5.0 - 2.0 * np.clip(q, 0, 45) / 45.0


def q_to_crosstalk(q):
    """
    Low quality -> more bleed into other channels.
    """
    q = np.asarray(q, dtype=np.float64)
    return 0.02 + 0.18 * (1.0 - np.clip(q, 0, 45) / 45.0)


def add_gaussian(trace, center, sigma, amplitude):
    """
    Add a Gaussian peak to a 1D trace.
    """
    radius = int(max(8, math.ceil(4 * sigma)))
    left = max(0, center - radius)
    right = min(len(trace), center + radius + 1)
    x = np.arange(left, right, dtype=np.float64)
    trace[left:right] += amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def deterministic_unit_wave(index, phase=0.0):
    """
    Smooth deterministic value in [0, 1] from an integer index.
    """
    x = float(index)
    w = (
        math.sin(0.37 * x + phase)
        + 0.6 * math.sin(0.11 * x + 1.7 + phase)
        + 0.4 * math.sin(0.07 * x + 2.9 + phase)
    )
    return 0.5 + 0.5 * (w / 2.0)


def synthesize_chromatogram(seq, q_scores, samples_per_base=12):
    """
    Build approximate 4-channel traces from sequence + PHRED qualities.

    Returns:
        traces: np.ndarray shape (4, n_points)
        peak_locs: np.ndarray shape (len(seq),)
        base_calls: str
        q_scores: np.ndarray
    """
    n_bases = len(seq)
    if n_bases == 0:
        raise ValueError("Empty sequence")

    q_scores = np.asarray(q_scores, dtype=np.int16)
    amps = q_to_amplitude(q_scores)
    sigmas = q_to_sigma(q_scores)
    bleed = q_to_crosstalk(q_scores)

    # Deterministic spacing variation: realistic but non-random.
    idx = np.arange(n_bases, dtype=np.float64)
    step_wave = (
        0.55 * np.sin(0.19 * idx + 0.8)
        + 0.25 * np.sin(0.047 * idx + 2.1)
        + 0.20 * np.sin(0.011 * idx + 1.3)
    )
    steps = np.maximum(8.0, samples_per_base + step_wave)
    peak_locs = np.cumsum(steps).astype(np.int32)
    peak_locs -= peak_locs[0]
    n_points = int(peak_locs[-1] + 8 * samples_per_base)

    traces = np.zeros((4, n_points), dtype=np.float32)

    # Deterministic background floor and drift per channel.
    x = np.arange(n_points, dtype=np.float64)
    for ch in range(4):
        baseline = (
            1.8
            + 0.55 * np.sin(0.017 * x + 0.7 * ch)
            + 0.25 * np.sin(0.051 * x + 1.2 + 0.4 * ch)
            + 0.12 * np.sin(0.13 * x + 2.3 + 0.9 * ch)
        )
        traces[ch] += baseline.astype(np.float32)

    for i, base in enumerate(seq):
        center = int(peak_locs[i])
        sigma = float(sigmas[i])
        amp = float(amps[i])
        crosstalk = float(bleed[i])

        if base in BASE_TO_CH:
            main_ch = BASE_TO_CH[base]
            add_gaussian(traces[main_ch], center, sigma, amp)

            # Bleed into non-called channels
            for ch in range(4):
                if ch != main_ch:
                    v = deterministic_unit_wave(i * 7 + ch, phase=0.4 * (main_ch + 1))
                    frac = crosstalk * (0.72 + 0.56 * v)
                    sigma_scale = 1.04 + 0.18 * deterministic_unit_wave(i * 5 + 3 * ch, phase=1.1)
                    add_gaussian(traces[ch], center, sigma * sigma_scale, amp * frac)

        else:
            # For N or ambiguous bases, spread weak signal across all channels
            for ch in range(4):
                add_gaussian(traces[ch], center, sigma * 1.2, amp * 0.28)

    # Clip to nonnegative, convert to integer-like signal.
    traces = np.clip(traces, 0, None)
    traces = np.rint(traces).astype(np.uint16)

    return traces, peak_locs, seq, q_scores


def save_synthetic_bundle(out_prefix, traces, peak_locs, seq, q_scores):
    """
    Save everything an ABIF writer would need.
    """
    np.savez_compressed(
        f"{out_prefix}.npz",
        trace_A=traces[0],
        trace_C=traces[1],
        trace_G=traces[2],
        trace_T=traces[3],
        peak_locations=peak_locs,
        base_calls=np.array(list(seq), dtype="U1"),
        q_scores=q_scores,
        channel_order=np.array(["A", "C", "G", "T"], dtype="U1"),
    )


def plot_chromatogram(out_png, traces, peak_locs, seq, max_bases=120):
    if plt is None:
        return

    # Plot only the first part by default so the PNG stays readable.
    max_idx = min(len(seq), max_bases)
    last_x = int(peak_locs[max_idx - 1] + 40)
    x = np.arange(last_x)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, traces[0][:last_x], label="A")
    ax.plot(x, traces[1][:last_x], label="C")
    ax.plot(x, traces[2][:last_x], label="G")
    ax.plot(x, traces[3][:last_x], label="T")

    for i in range(max_idx):
        px = int(peak_locs[i])
        ax.text(px, max(traces[:, px]) + 8, seq[i], ha="center", va="bottom", fontsize=8)

    ax.set_xlabel("synthetic scan index")
    ax.set_ylabel("intensity")
    ax.set_title("Synthetic chromatogram from FASTQ")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_real_ab1(out_ab1_path, traces, peak_locs, seq, q_scores):
    """
    Write an ABIF-encoded AB1 file with synthetic traces and calls.
    """

    @dataclass(frozen=True)
    class AbifEntry:
        tag: str
        tag_number: int
        elem_type: int
        elem_size: int
        num_elem: int
        payload: bytes

    def pack_entry(tag, tag_number, elem_type, elem_size, num_elem, payload, data_value):
        if len(tag) != 4:
            raise ValueError(f"ABIF tag must be 4 chars, got {tag!r}")
        data_size = len(payload)
        return ABIF_DIR_ENTRY_STRUCT.pack(
            tag.encode("ascii"),
            int(tag_number),
            int(elem_type),
            int(elem_size),
            int(num_elem),
            int(data_size),
            int(data_value),
            0,
        )

    def to_pstring(text):
        raw = text.encode("ascii", "replace")
        if len(raw) > 255:
            raw = raw[:255]
        return bytes([len(raw)]) + raw

    traces = np.asarray(traces, dtype=np.uint16)
    peak_locs = np.asarray(peak_locs)
    q_scores = np.asarray(q_scores)

    if traces.ndim != 2 or traces.shape[0] != 4:
        raise ValueError(f"Expected traces shape (4, n), got {traces.shape}")
    if len(seq) == 0:
        raise ValueError("Sequence is empty")
    if len(seq) != len(peak_locs) or len(seq) != len(q_scores):
        raise ValueError("seq, peak_locs, and q_scores must have equal lengths")
    if np.any(peak_locs < 0):
        raise ValueError("peak_locs must be nonnegative")
    if len(seq) > 65535:
        raise ValueError("Sequence length exceeds 16-bit AB1 PLOC capacity (65535 bases)")

    max_ploc = int(np.max(peak_locs))
    max_scan = max(max_ploc, int(traces.shape[1] - 1))
    if max_scan > 65535:
        # AB1 PLOC entries are 16-bit; compress the scan axis if synthetic data is wider.
        scale = 65535.0 / float(max_scan)
        old_n = traces.shape[1]
        new_n = max(2, int(np.floor((old_n - 1) * scale)) + 1)

        old_x = np.arange(old_n, dtype=np.float64)
        new_x = np.arange(new_n, dtype=np.float64) / scale
        new_traces = np.empty((4, new_n), dtype=np.uint16)
        for ch in range(4):
            interp = np.interp(new_x, old_x, traces[ch].astype(np.float64))
            new_traces[ch] = np.rint(np.clip(interp, 0, 65535)).astype(np.uint16)

        peak_locs = np.rint(peak_locs.astype(np.float64) * scale).astype(np.int32)
        peak_locs = np.maximum.accumulate(np.clip(peak_locs, 0, 65535))
        traces = new_traces

    n_points = int(traces.shape[1])
    n_bases = int(len(seq))
    seq_bytes = seq.encode("ascii")
    # ABI per-base confidence is byte-valued; clamp to uint8 range.
    q_bytes = np.clip(q_scores, 0, 255).astype(np.uint8, copy=False).tobytes()

    peak_payload = peak_locs.astype(">u2", copy=False).tobytes()
    base_payload = seq_bytes

    # ABI channel order is defined by FWO_ (here: GATC).
    channel_payloads = {
        "G": traces[BASE_TO_CH["G"]].astype(">u2", copy=False).tobytes(),
        "A": traces[BASE_TO_CH["A"]].astype(">u2", copy=False).tobytes(),
        "T": traces[BASE_TO_CH["T"]].astype(">u2", copy=False).tobytes(),
        "C": traces[BASE_TO_CH["C"]].astype(">u2", copy=False).tobytes(),
    }

    sn_values = []
    for ch in range(4):
        tr = traces[ch].astype(np.float64)
        baseline = np.percentile(tr, 20.0)
        noise_window = tr[tr <= np.percentile(tr, 30.0)]
        noise = float(np.std(noise_window)) if len(noise_window) else 1.0
        signal = float(np.percentile(tr, 99.5) - baseline)
        sn = int(np.clip(round(signal / max(noise, 1.0)), 1, 32767))
        sn_values.append(sn)
    sn_payload = np.asarray(sn_values, dtype=">u2").tobytes()

    pdmf_payload = to_pstring("KB_3500_POP7_BDTv3.mob")
    smpl_payload = to_pstring(Path(out_ab1_path).stem)
    cmnt_payload = to_pstring("Generated from FASTQ by fastq_to_ab1.py")

    entries = [
        AbifEntry("DATA", 1, 4, 2, n_points, channel_payloads["G"]),
        AbifEntry("DATA", 2, 4, 2, n_points, channel_payloads["A"]),
        AbifEntry("DATA", 3, 4, 2, n_points, channel_payloads["T"]),
        AbifEntry("DATA", 4, 4, 2, n_points, channel_payloads["C"]),
        AbifEntry("DATA", 9, 4, 2, n_points, channel_payloads["G"]),
        AbifEntry("DATA", 10, 4, 2, n_points, channel_payloads["A"]),
        AbifEntry("DATA", 11, 4, 2, n_points, channel_payloads["T"]),
        AbifEntry("DATA", 12, 4, 2, n_points, channel_payloads["C"]),
        AbifEntry("FWO_", 1, 2, 1, 4, b"GATC"),
        AbifEntry("LANE", 1, 4, 2, 1, struct.pack(">H", 1)),
        AbifEntry("PBAS", 1, 2, 1, n_bases, base_payload),
        AbifEntry("PBAS", 2, 2, 1, n_bases, base_payload),
        AbifEntry("PCON", 1, 2, 1, n_bases, q_bytes),
        AbifEntry("PCON", 2, 2, 1, n_bases, q_bytes),
        AbifEntry("PDMF", 1, 18, 1, len(pdmf_payload), pdmf_payload),
        AbifEntry("PDMF", 2, 18, 1, len(pdmf_payload), pdmf_payload),
        AbifEntry("PLOC", 1, 4, 2, n_bases, peak_payload),
        AbifEntry("PLOC", 2, 4, 2, n_bases, peak_payload),
        AbifEntry("S/N%", 1, 4, 2, 4, sn_payload),
        AbifEntry("SMPL", 1, 18, 1, len(smpl_payload), smpl_payload),
        AbifEntry("CMNT", 1, 18, 1, len(cmnt_payload), cmnt_payload),
    ]

    data_blocks = []
    data_start = 128
    cursor = data_start
    directory_chunks = []

    for entry in entries:
        if entry.num_elem * entry.elem_size != len(entry.payload):
            raise ValueError(
                f"{entry.tag}{entry.tag_number}: num_elem*elem_size does not match payload size"
            )

        if len(entry.payload) <= 4:
            data_value = int.from_bytes(entry.payload.ljust(4, b"\x00"), "big")
        else:
            data_value = cursor
            data_blocks.append((cursor, entry.payload))
            cursor += len(entry.payload)

        directory_chunks.append(
            pack_entry(
                entry.tag,
                entry.tag_number,
                entry.elem_type,
                entry.elem_size,
                entry.num_elem,
                entry.payload,
                data_value,
            )
        )

    directory_payload = b"".join(directory_chunks)
    directory_offset = cursor

    root_entry = ABIF_DIR_ENTRY_STRUCT.pack(
        b"tdir",
        1,
        1023,
        28,
        len(entries),
        len(directory_payload),
        directory_offset,
        0,
    )

    out_ab1_path = Path(out_ab1_path)
    out_ab1_path.parent.mkdir(parents=True, exist_ok=True)

    with out_ab1_path.open("wb") as f:
        f.write(b"ABIF")
        f.write(struct.pack(">H", 101))
        f.write(root_entry)

        if f.tell() > data_start:
            raise RuntimeError("ABIF header overflowed past data start offset")
        f.write(b"\x00" * (data_start - f.tell()))

        for offset, payload in data_blocks:
            if f.tell() != offset:
                raise RuntimeError("ABIF data offset mismatch while writing payload blocks")
            f.write(payload)

        if f.tell() != directory_offset:
            raise RuntimeError("ABIF directory offset mismatch")
        f.write(directory_payload)

    # Validate core ABIF invariants immediately after writing.
    blob = out_ab1_path.read_bytes()
    if blob[:4] != b"ABIF":
        raise RuntimeError("Wrote file without ABIF signature")
    root = ABIF_DIR_ENTRY_STRUCT.unpack(blob[6:34])
    if root[0] != b"tdir" or root[4] != len(entries) or root[5] != len(directory_payload):
        raise RuntimeError("ABIF root directory entry is invalid")
    if root[6] + root[5] != len(blob):
        raise RuntimeError("ABIF directory does not end at file end as expected")


def convert_fastq_to_synthetic_ab1_inputs(fastq_path, out_dir, max_reads=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for idx, (name, seq, qual) in enumerate(parse_fastq(fastq_path), start=1):
        q_scores = phred_from_ascii(qual, phred_offset=33)
        traces, peak_locs, seq, q_scores = synthesize_chromatogram(seq, q_scores, samples_per_base=12)

        safe_name = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)
        prefix = out_dir / safe_name

        save_synthetic_bundle(prefix, traces, peak_locs, seq, q_scores)
        plot_chromatogram(f"{prefix}.png", traces, peak_locs, seq)

        write_real_ab1(f"{prefix}.ab1", traces, peak_locs, seq, q_scores)

        print(f"Wrote {prefix}.npz")
        print(f"Wrote {prefix}.ab1")
        if plt is not None:
            print(f"Wrote {prefix}.png")

        if max_reads is not None and idx >= max_reads:
            break


if __name__ == "__main__":
    # Example:
    #   python synth_trace.py
    convert_fastq_to_synthetic_ab1_inputs(
        fastq_path="./data/barcode01.final.fastq",
        out_dir="synthetic_trace_output",
        max_reads=10,
    )
