#!/usr/bin/env python3
"""
Build customer-facing WPS order folders from an EPI2ME export.

This workflow:
1. discovers per-barcode EPI2ME files
2. maps barcodes to customer metadata
3. realigns raw unmapped BAM reads to the final consensus FASTA
4. generates per-base CSVs, QC plots, and synthetic AB1 files
5. writes a customer-ready package grouped by order number
6. renders a 2-page PDF report with an Alta-style layout
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import zipfile
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pysam

from align_bam_pipeline import run_pipeline
from fastq_to_ab1 import phred_from_ascii, synthesize_chromatogram, write_real_ab1
from generate_report import (
    MIN_MULTIMER_ALIGNMENT_FRACTION,
    MIN_MULTIMER_MAPQ,
    count_fasta_records,
    generate_report_data,
    read_first_fasta_record,
)


PACKAGE_SUBDIRS = {
    "ab1": "CHROMATOGRAM_FILES_ab1",
    "fasta": "FASTA_FILES",
    "gbk": "GENBANK_FILES",
    "per_base": "PER_BASE_BREAKDOWN",
    "qc": "QC REPORTS",
}

HEADER_ALIASES = {
    "barcode": {"barcode", "barcodeid", "barcodenumber", "barcode_name"},
    "sample_name": {"samplename", "sample", "sampleid", "plasmidname", "name"},
    "sample_id": {"id", "sampleidentifier", "samplecode"},
    "order_number": {"ordernumber", "order", "orderno", "ordernum"},
    "serial_number": {"sn", "serialnumber", "samplenumber", "sample_no"},
    "order_date": {"orderdate", "dateordered"},
    "report_date": {"reportdate"},
    "run_date": {"rundate", "sequencingdate"},
    "wps_date": {"wpsdate"},
    "received_date": {"receiveddate"},
    "customer_name": {"customer", "customername"},
    "concentration": {"concngul", "concentration", "concentrationngul"},
    "size": {"size", "plasmidsize"},
    "note": {"note", "notes", "comments"},
}

THEME = {
    "title": "#343C67",
    "heading": "#003B73",
    "teal": "#0D8686",
    "table_header": "#C6D9F1",
    "table_body": "#F3FBFB",
    "table_edge": "#D9D9D9",
    "purple": "#51459A",
    "green": "#44D8A4",
    "cyan": "#12A9D4",
    "text": "#23313f",
    "muted": "#5b6b77",
    "line": "#d5dde5",
}

# With qscore now taken directly from BAM base qualities, use a modest ONT-style cutoff.
LOW_CONFIDENCE_QSCORE = 12
SIZE_MISMATCH_TOLERANCE_FRACTION = 0.10
SIZE_MISMATCH_TOLERANCE_BP = 100


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("_") or "sample"


def normalize_header(header: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", header.lower())


def normalize_barcode(value: str | int | None) -> str | None:
    """Normalize barcode values like 3, 03, 3.0, barcode3, or barcode03."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    explicit = re.fullmatch(r"barcode\s*#?\s*0*(\d+)", text, flags=re.IGNORECASE)
    if explicit:
        return f"barcode{int(explicit.group(1)):02d}"
    numeric = re.fullmatch(r"0*(\d+)(?:\.0+)?", text)
    if numeric:
        return f"barcode{int(numeric.group(1)):02d}"
    return None


def normalize_excel_number_text(value: str | int | float | None) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    try:
        number = Decimal(text)
    except InvalidOperation:
        return text
    if number == number.to_integral_value():
        return str(int(number))
    return text


def normalize_metadata_date(value: str | int | float | None) -> str:
    text = normalize_excel_number_text(value)
    if not text:
        return ""
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        return text
    try:
        serial = Decimal(text)
    except InvalidOperation:
        return text
    if serial != serial.to_integral_value():
        return text
    days = int(serial)
    if not 20000 <= days <= 60000:
        return text
    # Excel's 1900 date system, including its historical leap-year offset.
    return (dt.date(1899, 12, 30) + dt.timedelta(days=days)).isoformat()


def normalize_order_number(value: str | int | float | None) -> str:
    text = normalize_excel_number_text(value)
    digits = re.fullmatch(r"\d+", text)
    if digits:
        return text
    return text.strip()


def parse_expected_size_bp(value: str | int | float | None) -> int | None:
    text = normalize_excel_number_text(value)
    if not text or text.strip().lower() in {"unk", "unknown", "na", "n/a", "none"}:
        return None
    cleaned = text.strip().lower().replace(",", "")
    match = re.search(r"(\d+(?:\.\d+)?)\s*(kb|kilobase|kilobases|k)?\b", cleaned)
    if not match:
        return None
    number = float(match.group(1))
    unit = match.group(2)
    if unit in {"kb", "kilobase", "kilobases", "k"}:
        number *= 1000
    return int(round(number))


def lengths_match(observed_bp: int, expected_bp: int) -> bool:
    allowed = max(SIZE_MISMATCH_TOLERANCE_BP, expected_bp * SIZE_MISMATCH_TOLERANCE_FRACTION)
    return abs(observed_bp - expected_bp) <= allowed


def build_wps_sample_name(meta: dict[str, str], row_number: int | None) -> str | None:
    sample_name = meta.get("sample_name")
    sample_id = meta.get("sample_id")
    if not sample_name:
        return None
    if sample_id:
        prefix = f"{row_number:03d}_" if row_number is not None else ""
        return f"{prefix}{sample_id}_{sample_name}"
    return sample_name


def load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    path = "xl/sharedStrings.xml"
    if path not in zf.namelist():
        return []
    root = ET.fromstring(zf.read(path))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    strings = []
    for si in root.findall("x:si", ns):
        parts = []
        for node in si.iter():
            if node.tag.endswith("}t") and node.text:
                parts.append(node.text)
        strings.append("".join(parts))
    return strings


def column_index_from_ref(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    if not letters:
        return 0
    total = 0
    for ch in letters.group(0):
        total = total * 26 + (ord(ch) - ord("A") + 1)
    return total - 1


def read_xlsx_rows(path: Path) -> list[dict[str, str]]:
    with zipfile.ZipFile(path) as zf:
        sheet_name = "xl/worksheets/sheet1.xml"
        if sheet_name not in zf.namelist():
            raise ValueError(f"Could not find {sheet_name} in {path}")
        shared = load_shared_strings(zf)
        root = ET.fromstring(zf.read(sheet_name))

    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    raw_rows: list[list[str]] = []
    for row in root.findall(".//x:sheetData/x:row", ns):
        cells: dict[int, str] = {}
        max_col = -1
        for cell in row.findall("x:c", ns):
            ref = cell.attrib.get("r", "A1")
            col_idx = column_index_from_ref(ref)
            max_col = max(max_col, col_idx)
            cell_type = cell.attrib.get("t")
            value = ""
            if cell_type == "inlineStr":
                node = cell.find("x:is/x:t", ns)
                value = node.text if node is not None and node.text is not None else ""
            else:
                node = cell.find("x:v", ns)
                raw = node.text if node is not None and node.text is not None else ""
                if cell_type == "s" and raw:
                    value = shared[int(raw)]
                else:
                    value = raw
            cells[col_idx] = value
        if max_col < 0:
            continue
        raw_rows.append([cells.get(idx, "") for idx in range(max_col + 1)])

    header = None
    records: list[dict[str, str]] = []
    for row in raw_rows:
        if header is None:
            if any(cell.strip() for cell in row):
                header = row
            continue
        if not any(cell.strip() for cell in row):
            continue
        row = row + [""] * max(0, len(header) - len(row))
        records.append({header[idx]: row[idx] for idx in range(len(header))})
    return records


def read_table_rows(path: Path) -> list[dict[str, str]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    if suffix == ".tsv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    if suffix == ".xlsx":
        return read_xlsx_rows(path)
    raise ValueError(f"Unsupported metadata format: {path}")


def canonicalize_metadata_row(row: dict[str, str]) -> dict[str, str]:
    normalized = {normalize_header(key): (value or "").strip() for key, value in row.items()}
    result: dict[str, str] = {}
    for canonical, aliases in HEADER_ALIASES.items():
        for alias in aliases:
            if alias in normalized and normalized[alias]:
                result[canonical] = normalized[alias]
                break
    if "barcode" in result:
        barcode = normalize_barcode(result["barcode"])
        if barcode is None:
            raise ValueError(f"Could not parse barcode value: {result['barcode']!r}")
        result["barcode"] = barcode
    if "order_number" in result:
        result["order_number"] = normalize_order_number(result["order_number"])
    if "sample_id" in result:
        result["sample_id"] = normalize_excel_number_text(result["sample_id"])
    for date_key in ("order_date", "report_date", "run_date", "received_date", "wps_date"):
        if date_key in result:
            result[date_key] = normalize_metadata_date(result[date_key])
    if "wps_date" in result:
        for date_key in ("order_date", "report_date", "run_date"):
            result.setdefault(date_key, result["wps_date"])
    return result


def resolve_metadata_path(path: Path) -> Path:
    if path.is_dir():
        candidates = sorted(
            item
            for item in path.iterdir()
            if item.is_file()
            and item.suffix.lower() in {".xlsx", ".csv", ".tsv"}
            and not item.name.startswith("~$")
        )
        if not candidates:
            raise ValueError(f"No metadata .xlsx, .csv, or .tsv file found in {path}")
        if len(candidates) > 1:
            names = ", ".join(item.name for item in candidates)
            raise ValueError(f"Multiple metadata files found in {path}: {names}")
        return candidates[0].resolve()
    return path.resolve()


def load_metadata_lookup(path: Path) -> dict[str, dict[str, str]]:
    metadata_path = resolve_metadata_path(path)
    lookup = {}
    barcode_rows: dict[str, int] = {}
    for row_number, row in enumerate(read_table_rows(metadata_path), start=1):
        meta = canonicalize_metadata_row(row)
        wps_sample_name = build_wps_sample_name(meta, row_number)
        if wps_sample_name:
            meta["sample_name"] = wps_sample_name
        barcode = meta.get("barcode")
        if barcode:
            if barcode in lookup:
                raise ValueError(
                    f"Duplicate metadata barcode {barcode}: rows {barcode_rows[barcode]} and {row_number}"
                )
            barcode_rows[barcode] = row_number
            lookup[barcode] = meta
    return lookup


def group_packaged_by_order(packaged: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    grouped: dict[str, dict[str, object]] = {}
    for item in packaged:
        order_number = item["order_number"]
        order = grouped.setdefault(
            order_number,
            {
                "order_dir": item["order_dir"],
                "sample_count": 0,
                "samples": [],
                "reports": [],
            },
        )
        order["sample_count"] = int(order["sample_count"]) + 1
        order["samples"].append(item["sample_name"])
        order["reports"].append(item["pdf"])
    return grouped


def add_discovered_file(
    records: dict[str, dict[str, Path]],
    discovery_errors: list[dict[str, str]],
    barcode: str,
    key: str,
    path: Path,
) -> None:
    existing = records[barcode].get(key)
    if existing is not None:
        discovery_errors.append(
            {
                "barcode": barcode,
                "reason": f"multiple {key} files found: {existing}, {path}",
            }
        )
    else:
        records[barcode][key] = path


def discover_files_for_key(
    records: dict[str, dict[str, Path]],
    discovery_errors: list[dict[str, str]],
    key: str,
    directory: Path,
    pattern: re.Pattern,
) -> None:
    if not directory.exists():
        raise ValueError(f"{key} directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"{key} path is not a directory: {directory}")
    for path in directory.iterdir():
        if not path.is_file():
            continue
        match = pattern.match(path.name)
        if not match:
            continue
        barcode = normalize_barcode(match.group(1))
        add_discovered_file(records, discovery_errors, barcode, key, path)


def infer_bam_barcode(path: Path) -> tuple[str | None, str | None]:
    filename_barcodes = sorted(
        {
            normalize_barcode(match)
            for match in re.findall(r"barcode\d+", path.name, flags=re.IGNORECASE)
        }
    )
    filename_barcodes = [barcode for barcode in filename_barcodes if barcode is not None]
    folder_barcodes = [
        normalize_barcode(parent.name)
        for parent in path.parents
        if re.fullmatch(r"barcode\d+", parent.name, re.IGNORECASE)
    ]
    folder_barcode = folder_barcodes[0] if folder_barcodes else None

    if len(filename_barcodes) > 1:
        if folder_barcode and folder_barcode in filename_barcodes:
            return folder_barcode, None
        return None, f"BAM filename contains multiple barcode tokens {filename_barcodes}: {path}"

    filename_barcode = filename_barcodes[0] if filename_barcodes else None
    if filename_barcode and folder_barcode and filename_barcode != folder_barcode:
        return None, f"BAM filename barcode {filename_barcode} does not match folder barcode {folder_barcode}: {path}"
    return folder_barcode or filename_barcode, None


def discover_bam_files(
    records: dict[str, dict[str, Path]],
    discovery_errors: list[dict[str, str]],
    directory: Path,
) -> None:
    if not directory.exists():
        raise ValueError(f"bam directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"bam path is not a directory: {directory}")

    for path in directory.rglob("*.bam"):
        if not path.is_file():
            continue
        barcode, error = infer_bam_barcode(path)
        if error:
            discovery_errors.append({"barcode": "unknown", "reason": error})
            continue
        if barcode is None:
            continue
        add_discovered_file(records, discovery_errors, barcode, "bam", path)


def discover_input_records(
    fasta_dir: Path,
    genbank_dir: Path,
    bam_dir: Path,
    fastq_dir: Path | None = None,
    maf_dir: Path | None = None,
) -> tuple[dict[str, dict[str, Path]], list[dict[str, str]]]:
    records: dict[str, dict[str, Path]] = defaultdict(dict)
    discovery_errors: list[dict[str, str]] = []
    discover_files_for_key(
        records,
        discovery_errors,
        "fasta",
        fasta_dir,
        re.compile(r"^(barcode\d+)\.final\.(?:fasta|fa)$", re.IGNORECASE),
    )
    discover_files_for_key(
        records,
        discovery_errors,
        "gbk",
        genbank_dir,
        re.compile(r"^(barcode\d+)\.annotations\.gbk$", re.IGNORECASE),
    )
    discover_bam_files(records, discovery_errors, bam_dir)
    if fastq_dir is not None:
        discover_files_for_key(
            records,
            discovery_errors,
            "fastq",
            fastq_dir,
            re.compile(r"^(barcode\d+)\.final\.(?:fastq|fq)$", re.IGNORECASE),
        )
    if maf_dir is not None:
        discover_files_for_key(
            records,
            discovery_errors,
            "maf",
            maf_dir,
            re.compile(r"^(barcode\d+)\.assembly\.maf$", re.IGNORECASE),
        )
    return dict(records), discovery_errors


def sample_stem_for_barcode(barcode: str, metadata: dict[str, str]) -> str:
    return slugify(metadata.get("sample_name") or barcode)


def find_sample_stem_collisions(
    records: dict[str, dict[str, Path]],
    metadata_lookup: dict[str, dict[str, str]],
    requested: set[str] | None,
    invalid_barcodes: set[str],
) -> dict[str, str]:
    stems: dict[str, list[str]] = defaultdict(list)
    for barcode in records:
        if requested is not None and barcode not in requested:
            continue
        if barcode in invalid_barcodes:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            continue
        stem = sample_stem_for_barcode(barcode, metadata_lookup.get(barcode, {}))
        stems[stem].append(barcode)
    collisions = {}
    for stem, barcodes in stems.items():
        if len(barcodes) > 1:
            joined = ", ".join(sorted(barcodes))
            for barcode in barcodes:
                collisions[barcode] = f"sample name collision after slugify: {stem!r} used by {joined}"
    return collisions


def write_renamed_fasta(src_fasta: Path, dst_fasta: Path, sequence_name: str) -> dict[str, str | int]:
    record = read_first_fasta_record(src_fasta)
    sequence = "".join(ch if ord(ch) < 128 else "N" for ch in record["sequence"])
    dst_fasta.parent.mkdir(parents=True, exist_ok=True)
    with dst_fasta.open("w", encoding="ascii") as handle:
        handle.write(f">{sequence_name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")
    return {"name": sequence_name, "sequence": sequence, "length_bp": len(sequence)}


def rewrite_genbank_locus(src_gbk: Path, dst_gbk: Path, locus_name: str) -> None:
    dst_gbk.parent.mkdir(parents=True, exist_ok=True)
    lines = src_gbk.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    safe_locus = re.sub(r"[^A-Za-z0-9_.-]+", "_", locus_name).strip("_")[:16] or "record"
    out_lines = []
    for idx, line in enumerate(lines):
        if idx == 0 and line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 3:
                bp = parts[2]
                remainder = ""
                if "bp" in line:
                    remainder = line[line.index("bp") + 2 :]
                line = f"LOCUS       {safe_locus:<16}{bp:>11} bp{remainder}"
            else:
                line = f"LOCUS       {safe_locus}"
        out_lines.append(line)
    dst_gbk.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def parse_fastq_record(path: Path) -> tuple[str, str, str]:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
        seq = handle.readline().rstrip("\n\r")
        handle.readline()
        qual = handle.readline().rstrip("\n\r")
    if not header.startswith("@") or len(seq) != len(qual):
        raise ValueError(f"Could not parse FASTQ record from {path}")
    return header[1:], seq.upper(), qual


def generate_ab1_files(
    fasta_path: Path,
    fastq_path: Path | None,
    output_dir: Path,
    sample_stem: str,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    record = read_first_fasta_record(fasta_path)
    if fastq_path and fastq_path.exists():
        _, fastq_seq, fastq_qual = parse_fastq_record(fastq_path)
        if fastq_seq == record["sequence"]:
            q_scores = phred_from_ascii(fastq_qual, phred_offset=33)
        else:
            q_scores = None
    else:
        q_scores = None

    if q_scores is None:
        q_scores = [30] * len(record["sequence"])

    traces, peak_locs, seq, q_scores = synthesize_chromatogram(
        record["sequence"],
        q_scores,
        samples_per_base=12,
    )

    for stale_name in (
        f"{sample_stem}_contig_trace.ab1",
        f"{sample_stem}_contig_trace1.ab1",
        f"{sample_stem}_contig_trace2.ab1",
    ):
        stale_path = output_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()

    out_ab1 = output_dir / f"{sample_stem}_contig_trace.ab1"
    write_real_ab1(out_ab1, traces, peak_locs, seq, q_scores)
    return [out_ab1]


def read_lengths_from_bam(bam_path: Path) -> list[int]:
    lengths = []
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            lengths.append(read.query_length or 0)
    return lengths


def plot_pdf_coverage_map(per_base_csv: Path, low_conf_csv: Path, out_path: Path) -> Path:
    positions = []
    depths = []
    depth_by_pos = {}
    with per_base_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            pos = int(row["pos"])
            depth = int(row["depth"])
            positions.append(pos)
            depths.append(depth)
            depth_by_pos[pos] = depth

    low_positions = []
    with low_conf_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            low_positions.append(int(row["pos"]))
    low_depths = [depth_by_pos[pos] for pos in low_positions if pos in depth_by_pos]

    fig, ax = plt.subplots(figsize=(7.6, 3.2))
    if positions:
        ax.plot(positions, depths, color="#4f9da6", linewidth=1.3)
        if low_positions:
            ax.scatter(low_positions, low_depths, marker="x", color="#e67e22", s=16, linewidths=0.8)
        ax.set_xlim(left=0, right=max(positions))
    else:
        ax.text(0.5, 0.5, "No coverage data available", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(left=0, right=1)
        ax.set_ylim(bottom=0, top=1)
    ax.set_xlabel("Base Position")
    ax.set_ylabel("Depth")
    ax.grid(alpha=0.12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def plot_read_length_vs_bases(raw_bam: Path, aligned_bam: Path, contig_length: int, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapped_names = set()
    with pysam.AlignmentFile(aligned_bam, "rb") as bam:
        seen_aligned_names = set()
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in seen_aligned_names:
                raise ValueError(f"Duplicate primary read name in aligned BAM: {read_name!r}")
            seen_aligned_names.add(read_name)
            if read.is_unmapped:
                continue
            alignment_fraction = ((read.query_alignment_length or 0) / read.query_length) if read.query_length else 0.0
            if read.mapping_quality >= MIN_MULTIMER_MAPQ and alignment_fraction >= MIN_MULTIMER_ALIGNMENT_FRACTION:
                mapped_names.add(read_name)

    mapped_lengths = []
    other_lengths = []
    with pysam.AlignmentFile(raw_bam, "rb", check_sq=False) as bam:
        seen_raw_names = set()
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in seen_raw_names:
                raise ValueError(f"Duplicate primary read name in raw BAM: {read_name!r}")
            seen_raw_names.add(read_name)
            qlen = read.query_length or 0
            if qlen <= 0:
                continue
            if read_name in mapped_names:
                mapped_lengths.append(qlen)
            else:
                other_lengths.append(qlen)

    all_lengths = mapped_lengths + other_lengths
    fig, ax = plt.subplots(figsize=(7.6, 3.7))
    if all_lengths:
        max_len = max(all_lengths)
        bin_size = max(250, int(math.ceil(max_len / 24 / 50.0)) * 50)
        start = int(math.floor(min(all_lengths) / bin_size) * bin_size)
        stop = int(math.ceil(max_len / bin_size) * bin_size) + bin_size
        bins = np.arange(start, stop + bin_size, bin_size)
        centers = bins[:-1] + bin_size / 2
        mapped_counts, _ = np.histogram(mapped_lengths, bins=bins)
        other_counts, _ = np.histogram(other_lengths, bins=bins)

        band_left = contig_length * 0.9
        band_right = contig_length * 1.1
        ymax = max((mapped_counts + other_counts).max(), 1)
        ax.axvspan(band_left, band_right, color="#d9d9d9", alpha=0.6, lw=0)
        ax.text((band_left + band_right) / 2, ymax * 0.96, "monomer", ha="center", va="bottom", fontsize=9, fontweight="bold")

        width = bin_size * 0.78
        ax.bar(centers, mapped_counts, width=width, color=THEME["purple"], label="Eligible mapped reads")
        ax.bar(centers, other_counts, width=width, bottom=mapped_counts, color=THEME["green"], label="Other reads")
        ax.set_xlim(left=max(0, start - bin_size * 0.25), right=stop)
        ax.set_ylim(bottom=0)
    else:
        ax.text(0.5, 0.5, "No read-length data available", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(left=0, right=1)
        ax.set_ylim(bottom=0, top=1)
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Read Count")
    ax.legend(loc="upper right", frameon=True, facecolor="white")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def find_existing_path(paths: Iterable[Path | None]) -> Path | None:
    for path in paths:
        if path is not None and path.exists():
            return path
    return None


def validate_length_consistency(metadata: dict[str, str], fasta_length_bp: int, report_summary: dict) -> None:
    expected_size_bp = parse_expected_size_bp(metadata.get("size"))
    if expected_size_bp is not None and not lengths_match(fasta_length_bp, expected_size_bp):
        raise ValueError(
            "metadata Size does not match FASTA contig length: "
            f"metadata={expected_size_bp:,} bp, FASTA={fasta_length_bp:,} bp"
        )

    genbank_summary = report_summary.get("genbank_summary") or {}
    genbank_length_bp = genbank_summary.get("length_bp")
    if genbank_length_bp is not None and int(genbank_length_bp) != int(fasta_length_bp):
        raise ValueError(
            "GenBank LOCUS length does not match FASTA contig length: "
            f"GenBank={int(genbank_length_bp):,} bp, FASTA={fasta_length_bp:,} bp"
        )


def add_logo(fig, logos: list[Path]) -> None:
    for logo in logos[:1]:
        if not logo.exists():
            continue
        ax = fig.add_axes([0.055, 0.94, 0.16, 0.055])
        ax.set_axis_off()
        ax.imshow(mpimg.imread(logo))


def draw_report_title(fig, logos: list[Path]) -> None:
    add_logo(fig, logos)
    header_ax = fig.add_axes([0.05, 0.905, 0.90, 0.03])
    header_ax.set_axis_off()
    header_ax.add_patch(Rectangle((0, 0.45), 1, 0.1, color=THEME["teal"], lw=0))

    ax = fig.add_axes([0.06, 0.855, 0.88, 0.038])
    ax.set_axis_off()
    ax.text(
        0.5,
        0.5,
        "Whole Plasmid Sequencing Report",
        ha="center",
        va="center",
        fontsize=18.5,
        fontweight="bold",
        color=THEME["title"],
        family="DejaVu Serif",
    )


def draw_section_heading(fig, y: float, title: str) -> None:
    ax = fig.add_axes([0.06, y, 0.88, 0.04])
    ax.set_axis_off()
    ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=13.5, fontweight="bold", color=THEME["heading"])


def draw_table(fig, bbox, headers, values, col_widths=None) -> None:
    ax = fig.add_axes(bbox)
    ax.set_axis_off()
    table = ax.table(
        cellText=[values],
        colLabels=headers,
        cellLoc="center",
        colLoc="center",
        loc="center",
        bbox=[0, 0, 1, 1],
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.1)
    table.scale(1, 1.95)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(THEME["table_edge"])
        cell.set_linewidth(0.8)
        if row == 0:
            cell.set_facecolor(THEME["table_header"])
            cell.set_text_props(weight="bold", color="#333333", fontsize=8.6)
        else:
            cell.set_facecolor(THEME["table_body"])
            cell.set_text_props(color="#333333", fontsize=9.1)


def draw_image(fig, bounds, image_path: Path) -> None:
    ax = fig.add_axes(bounds)
    ax.set_axis_off()
    ax.imshow(mpimg.imread(image_path))


def draw_footer(fig) -> None:
    footer_ax = fig.add_axes([0.06, 0.03, 0.88, 0.035])
    footer_ax.set_axis_off()
    footer_ax.text(
        0.5,
        0.5,
        "Alta Biotech, LLC  |  2115 N Scranton St Ste 3040B, Aurora CO 80045  |  Tel: 720-640-9400  |  Support@altabiotech.com  |  www.altabiotech.com",
        ha="center",
        va="center",
        fontsize=7.5,
        color=THEME["muted"],
    )


def report_date_value(metadata: dict[str, str]) -> str:
    return (
        metadata.get("report_date")
        or metadata.get("run_date")
        or metadata.get("order_date")
        or metadata.get("received_date")
        or ""
    )


def ecoli_contamination_pct(report_summary: dict) -> float | None:
    contamination = report_summary.get("contamination") or {}
    if contamination.get("ecoli_genomic_contamination_pct") is not None:
        return contamination["ecoli_genomic_contamination_pct"]
    return None


def format_percent(value: float | None) -> str:
    if value is None:
        return "N/A"
    return f"{value:.2f}%"


def format_yes_no(value: bool | None) -> str:
    if value is None:
        return "N/A"
    return "Yes" if value else "No"


def validate_percent_value(name: str, value: object, allow_none: bool = False) -> float | None:
    if value is None and allow_none:
        return None
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} is not a finite percentage: {value!r}")
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"{name} is outside 0..100%: {value}")
    return value


def multimer_pdf_values(assembly: dict) -> list[float | None]:
    required = [
        ("monomer_pct", "Monomer"),
        ("dimer_pct", "Dimer"),
        ("trimer_pct", "Trimer"),
        ("tetramer_pct", "Tetramer"),
        ("unclassified_multimer_read_pct", "Unclassified"),
    ]
    missing = [key for key, _label in required if key not in assembly]
    if missing:
        raise ValueError(f"Missing multimer fields: {', '.join(missing)}")

    if not assembly.get("multimer_calculated"):
        return [None for _key, _label in required]

    values = [validate_percent_value(label, assembly[key]) for key, label in required]
    total = sum(value for value in values if value is not None)
    if not 99.0 <= total <= 101.0:
        raise ValueError(f"Multimer percentages should sum to ~100%, got {total:.2f}%")
    return values


def render_pdf_report(
    report_summary: dict,
    metadata: dict[str, str],
    output_pdf: Path,
    sample_stem: str,
    coverage_png: Path,
    read_length_bases_png: Path,
    feature_map_png: Path | None,
    logos: list[Path],
) -> None:
    del feature_map_png
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    contig = report_summary["contig"]
    assembly = report_summary["assembly_status"]
    coverage = report_summary["coverage"]
    contamination = ecoli_contamination_pct(report_summary)
    multimer_values = multimer_pdf_values(assembly)
    with PdfPages(output_pdf) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor("white")
        draw_report_title(fig, logos)

        draw_section_heading(fig, 0.785, "Sample & Order Information")
        draw_table(
            fig,
            [0.07, 0.69, 0.86, 0.08],
            ["Sample Name", "Order\nNumber", "Report\nDate"],
            [
                metadata.get("sample_name") or sample_stem,
                metadata.get("order_number", "UNKNOWN"),
                report_date_value(metadata),
            ],
            col_widths=[0.50, 0.25, 0.25],
        )

        draw_section_heading(fig, 0.615, "Assembly Summary")
        draw_table(
            fig,
            [0.07, 0.52, 0.86, 0.08],
            ["Contig Length\n(bp)", "Bases Mapped", "Reads Mapped", "E. coli DNA %", "Is Circular?"],
            [
                f"{contig['length_bp']:,}",
                f"{assembly.get('bases_mapped', 0):,} ({assembly.get('bases_mapped_pct', 0):.2f}%)",
                f"{assembly.get('reads_mapped', 0):,} ({assembly.get('reads_mapped_pct', 0):.2f}%)",
                format_percent(contamination),
                format_yes_no(contig.get("is_circular")),
            ],
        )

        draw_section_heading(fig, 0.445, "Nanopore Performance")
        draw_table(
            fig,
            [0.07, 0.35, 0.86, 0.08],
            ["Mean Read\nDepth", "Min/Max Depth", "Mean\nCoverage", "Low Confidence\nCount", "Single Contig?"],
            [
                f"{round(coverage.get('mean_depth', 0)):,}",
                f"{coverage.get('min_depth', 0):,} / {coverage.get('max_depth', 0):,}",
                f"{round(coverage.get('mean_depth', 0)):,}x",
                f"{coverage.get('low_confidence_count', 0):,}",
                format_yes_no(assembly.get("single_contig")),
            ],
        )

        draw_section_heading(fig, 0.275, "Multimer Analysis")
        draw_table(
            fig,
            [0.07, 0.18, 0.86, 0.08],
            ["Monomer", "Dimer", "Trimer", "Tetramer", "Unclassified"],
            [format_percent(value) for value in multimer_values],
            col_widths=[0.20, 0.20, 0.20, 0.20, 0.20],
        )

        pdf.savefig(fig)
        plt.close(fig)

        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor("white")
        draw_report_title(fig, logos)

        draw_section_heading(fig, 0.785, "Coverage & Distribution Maps")

        cov_head = fig.add_axes([0.08, 0.735, 0.84, 0.04])
        cov_head.set_axis_off()
        cov_head.text(0.0, 0.7, "Coverage Map", ha="left", va="center", fontsize=11.5, fontweight="bold", color=THEME["heading"])
        cov_head.text(
            0.0,
            0.15,
            'low confidence positions are marked with orange "X"',
            ha="left",
            va="center",
            fontsize=9,
            color=THEME["muted"],
            style="italic",
        )
        draw_image(fig, [0.08, 0.47, 0.84, 0.23], coverage_png)

        dist_head = fig.add_axes([0.08, 0.39, 0.84, 0.04])
        dist_head.set_axis_off()
        dist_head.text(0.0, 0.5, "Read Length Distribution", ha="left", va="center", fontsize=11.5, fontweight="bold", color=THEME["heading"])
        draw_image(fig, [0.08, 0.095, 0.84, 0.255], read_length_bases_png)
        draw_footer(fig)

        pdf.savefig(fig)
        plt.close(fig)


def collect_logos(paths: Iterable[str]) -> list[Path]:
    logos = []
    for item in paths:
        path = Path(item)
        if path.exists():
            logos.append(path)
    return logos


def path_is_inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def remove_output_path(path: Path, output_root: Path) -> None:
    if not path_is_inside(path, output_root) or not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def remove_empty_package_dirs(order_dir: Path) -> None:
    for subdir in PACKAGE_SUBDIRS.values():
        path = order_dir / subdir
        if path.exists() and path.is_dir():
            try:
                path.rmdir()
            except OSError:
                pass
    try:
        order_dir.rmdir()
    except OSError:
        pass


def cleanup_previous_barcode_output(output_root: Path, barcode: str, sample_stem: str) -> None:
    work_root = output_root / "_work"
    touched_order_dirs: set[Path] = set()
    if work_root.exists():
        for summary_path in work_root.glob("*/package_summary.json"):
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if summary.get("barcode") != barcode:
                continue
            paths = summary.get("paths") or {}
            order_dir = paths.get("order_dir")
            if order_dir:
                touched_order_dirs.add(Path(order_dir))
            for value in paths.values():
                if isinstance(value, list):
                    for item in value:
                        remove_output_path(Path(item), output_root)
                elif value and Path(value).name != "":
                    path = Path(value)
                    if path.is_file() or path.suffix:
                        remove_output_path(path, output_root)
            remove_output_path(summary_path.parent, output_root)

    current_work_dir = work_root / sample_stem
    current_summary = current_work_dir / "package_summary.json"
    if current_summary.exists():
        try:
            summary = json.loads(current_summary.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            summary = {}
        existing_barcode = summary.get("barcode")
        if existing_barcode and existing_barcode != barcode:
            raise ValueError(
                f"existing work directory {current_work_dir} belongs to {existing_barcode}, not {barcode}"
            )
    remove_output_path(current_work_dir, output_root)

    for order_dir in touched_order_dirs:
        remove_empty_package_dirs(order_dir)


def package_sample(
    barcode: str,
    record: dict[str, Path],
    metadata: dict[str, str],
    output_root: Path,
    logos: list[Path],
    threads: int = 1,
    sort_memory: str = "768M",
    keep_intermediates: bool = False,
    allow_aligned_input: bool = False,
) -> dict[str, str]:
    sample_name = metadata.get("sample_name") or barcode
    sample_stem = sample_stem_for_barcode(barcode, metadata)
    sequence_name = f"{sample_stem}_contig"
    order_number = metadata.get("order_number")
    if not order_number:
        raise ValueError("metadata is missing order_number; refusing to package under WPS Data_Order #UNKNOWN")
    fasta_record_count = count_fasta_records(record["fasta"])
    if fasta_record_count != 1:
        raise ValueError(f"expected exactly one FASTA record for {barcode}, found {fasta_record_count}")

    order_dir = output_root / f"WPS Data_Order #{order_number}"
    cleanup_previous_barcode_output(output_root, barcode, sample_stem)
    package_dirs = {name: order_dir / subdir for name, subdir in PACKAGE_SUBDIRS.items()}
    for path in package_dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    work_dir = output_root / "_work" / sample_stem
    work_dir.mkdir(parents=True, exist_ok=True)

    fasta_out = package_dirs["fasta"] / f"{sample_stem}_contig.fa"
    gbk_out = package_dirs["gbk"] / f"{sample_stem}_contig.gbk"
    renamed = write_renamed_fasta(record["fasta"], fasta_out, sequence_name)
    rewrite_genbank_locus(record["gbk"], gbk_out, sequence_name)

    alignment_dir = work_dir / "alignment"
    alignment_result = run_pipeline(
        record["bam"],
        fasta_out,
        alignment_dir,
        minimap2_preset="map-ont",
        threads=threads,
        sort_memory=sort_memory,
        keep_intermediates=keep_intermediates,
        allow_aligned_input=allow_aligned_input,
    )
    aligned_bam = Path(alignment_result["sorted_bam"])
    fasta_index = Path(f"{fasta_out}.fai")
    if fasta_index.exists():
        fasta_index.unlink()

    report_dir = work_dir / "report"
    report_summary = generate_report_data(
        aligned_bam=aligned_bam,
        contig_fasta=fasta_out,
        out_dir=report_dir,
        reference_fasta=fasta_out,
        maf_path=record.get("maf"),
        gbk_path=gbk_out,
        sample_name=sample_stem,
        low_confidence_qscore=LOW_CONFIDENCE_QSCORE,
    )
    validate_length_consistency(metadata, renamed["length_bp"], report_summary)

    per_base_src = Path(report_summary["outputs"]["per_base_details_csv"])
    low_conf_src = Path(report_summary["outputs"]["low_confidence_bases_csv"])
    coverage_png = plot_pdf_coverage_map(per_base_src, low_conf_src, work_dir / "coverage_map_pdf.png")
    feature_map_value = report_summary["outputs"].get("feature_map_png")
    feature_map_png = find_existing_path([Path(feature_map_value)]) if feature_map_value else None

    per_base_dst = package_dirs["per_base"] / f"{sample_stem}_contig_per_base_details.csv"
    low_conf_dst = package_dirs["per_base"] / f"{sample_stem}_contig_low_confidence_bases.csv"
    shutil.copyfile(per_base_src, per_base_dst)
    shutil.copyfile(low_conf_src, low_conf_dst)

    ab1_paths = generate_ab1_files(fasta_out, record.get("fastq"), package_dirs["ab1"], sample_stem)

    bases_plot = plot_read_length_vs_bases(
        record["bam"],
        aligned_bam,
        renamed["length_bp"],
        work_dir / "read_length_vs_bases.png",
    )

    pdf_out = package_dirs["qc"] / f"{sample_stem}_report.pdf"
    render_pdf_report(
        report_summary=report_summary,
        metadata={**metadata, "barcode": barcode},
        output_pdf=pdf_out,
        sample_stem=sample_stem,
        coverage_png=coverage_png,
        read_length_bases_png=bases_plot,
        feature_map_png=feature_map_png,
        logos=logos,
    )

    summary_out = work_dir / "package_summary.json"
    summary_out.write_text(
        json.dumps(
            {
                "barcode": barcode,
                "sample_name": sample_stem,
                "order_number": order_number,
                "paths": {
                    "order_dir": str(order_dir),
                    "pdf": str(pdf_out),
                    "fasta": str(fasta_out),
                    "gbk": str(gbk_out),
                    "ab1": [str(path) for path in ab1_paths],
                    "per_base_details": str(per_base_dst),
                    "low_confidence": str(low_conf_dst),
                    "aligned_bam": str(aligned_bam),
                },
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    return {
        "barcode": barcode,
        "sample_name": sample_stem,
        "order_number": order_number,
        "order_dir": str(order_dir),
        "pdf": str(pdf_out),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta-dir",
        required=True,
        help="Directory containing barcodeXX.final.fasta or barcodeXX.final.fa files.",
    )
    parser.add_argument(
        "--genbank-dir",
        required=True,
        help="Directory containing barcodeXX.annotations.gbk files.",
    )
    parser.add_argument(
        "--bam-dir",
        required=True,
        help="Directory containing raw/unmapped BAM files. BAMs may be directly inside it or inside barcodeXX subfolders.",
    )
    parser.add_argument(
        "--fastq-dir",
        default=None,
        help="Optional directory containing barcodeXX.final.fastq or barcodeXX.final.fq files.",
    )
    parser.add_argument(
        "--maf-dir",
        default=None,
        help="Optional directory containing barcodeXX.assembly.maf files.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Required metadata CSV, TSV, XLSX, or directory containing exactly one metadata file.",
    )
    parser.add_argument(
        "--output-dir",
        default="customer_packages",
        help="Output directory for WPS order folders. Default: ./customer_packages",
    )
    parser.add_argument(
        "--barcodes",
        nargs="*",
        default=None,
        help="Optional barcode filter, for example: barcode01 barcode02 or 1 2",
    )
    parser.add_argument(
        "--logo",
        action="append",
        default=[],
        help="Optional logo image(s) to place in the PDF header. Pass up to two times.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for minimap2 and samtools alignment steps. Default: 1.",
    )
    parser.add_argument(
        "--sort-memory",
        default="768M",
        help="Memory per samtools sort thread, for example 768M or 2G. Default: 768M.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep alignment intermediates reads.fastq, aligned.sam, and aligned.unsorted.bam.",
    )
    parser.add_argument(
        "--allow-aligned-input",
        action="store_true",
        help="Allow BAMs that already contain mapped primary reads.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    metadata_path = resolve_metadata_path(Path(args.metadata).resolve())
    metadata_lookup = load_metadata_lookup(metadata_path)
    input_dirs = {
        "fasta_dir": str(Path(args.fasta_dir).resolve()),
        "genbank_dir": str(Path(args.genbank_dir).resolve()),
        "bam_dir": str(Path(args.bam_dir).resolve()),
        "fastq_dir": str(Path(args.fastq_dir).resolve()) if args.fastq_dir else None,
        "maf_dir": str(Path(args.maf_dir).resolve()) if args.maf_dir else None,
    }
    records, discovery_errors = discover_input_records(
        fasta_dir=Path(args.fasta_dir).resolve(),
        genbank_dir=Path(args.genbank_dir).resolve(),
        bam_dir=Path(args.bam_dir).resolve(),
        fastq_dir=Path(args.fastq_dir).resolve() if args.fastq_dir else None,
        maf_dir=Path(args.maf_dir).resolve() if args.maf_dir else None,
    )

    requested = None
    if args.barcodes:
        requested = {normalize_barcode(item) for item in args.barcodes}

    logos = collect_logos(args.logo)
    packaged = []
    skipped = list(discovery_errors)
    invalid_barcodes = {item["barcode"] for item in discovery_errors}
    considered_barcodes = set(records)
    if requested is not None:
        considered_barcodes &= requested

    sample_stem_collisions = find_sample_stem_collisions(records, metadata_lookup, requested, invalid_barcodes)
    for barcode, reason in sorted(sample_stem_collisions.items()):
        skipped.append({"barcode": barcode, "reason": reason})
    invalid_barcodes.update(sample_stem_collisions)

    for barcode in sorted(considered_barcodes):
        if barcode in invalid_barcodes:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            skipped.append({"barcode": barcode, "reason": "barcode not present in metadata"})

    for barcode in sorted(set(metadata_lookup) - set(records)):
        if requested is not None and barcode not in requested:
            continue
        skipped.append({"barcode": barcode, "reason": "metadata row has no matching EPI2ME files"})

    for barcode in sorted(records):
        if requested is not None and barcode not in requested:
            continue
        if barcode in invalid_barcodes:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            continue
        record = records[barcode]
        missing = [key for key in ("fasta", "gbk", "bam") if key not in record]
        if missing:
            skipped.append({"barcode": barcode, "reason": f"missing required files: {', '.join(missing)}"})
            continue
        metadata = metadata_lookup.get(barcode, {})
        try:
            packaged.append(
                package_sample(
                    barcode,
                    record,
                    metadata,
                    output_dir,
                    logos,
                    threads=args.threads,
                    sort_memory=args.sort_memory,
                    keep_intermediates=args.keep_intermediates,
                    allow_aligned_input=args.allow_aligned_input,
                )
            )
        except subprocess.CalledProcessError as exc:
            skipped.append({"barcode": barcode, "reason": f"command failed ({exc.returncode}): {' '.join(map(str, exc.cmd))}"})
        except Exception as exc:
            skipped.append({"barcode": barcode, "reason": str(exc)})

    grouped_orders = group_packaged_by_order(packaged)

    summary = {
        "input_dirs": input_dirs,
        "metadata": str(metadata_path) if metadata_path else None,
        "output_dir": str(output_dir),
        "packaged": packaged,
        "orders": grouped_orders,
        "skipped": skipped,
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(f"run_summary_json: {summary_path}")
    print(f"packaged_count: {len(packaged)}")
    print(f"order_count: {len(grouped_orders)}")
    print(f"skipped_count: {len(skipped)}")
    for order_number, order in sorted(grouped_orders.items()):
        print(f"order: {order_number} -> {order['order_dir']} ({order['sample_count']} sample(s))")
    for item in packaged:
        print(f"packaged: {item['barcode']} -> {item['order_dir']}")
    for item in skipped:
        print(f"skipped: {item['barcode']} ({item['reason']})")


if __name__ == "__main__":
    main()
