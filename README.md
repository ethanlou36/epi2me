# EPI2ME to WPS Report Workflow

This folder turns an EPI2ME/ONT plasmid sequencing export into customer-facing
WPS order folders with one PDF report per sample.

If several samples have the same `Order #` in the metadata sheet, their reports
go into the same `WPS Data_Order #...` folder. Samples with different order
numbers go into separate folders.

The main instructions below are for Windows using Ubuntu/WSL. Run the commands
in the Ubuntu terminal, not in PowerShell, Command Prompt, or Anaconda Prompt.

## 1. Open Ubuntu and Activate Python

1. Open **Ubuntu** from the Windows Start menu.

2. Move into the project folder:

```bash
cd /mnt/c/Users/altab/epi2me
```

That Ubuntu path points to this Windows folder:

```text
C:\Users\altab\epi2me
```

3. Activate the Python virtual environment named `.venv`:

```bash
source .venv/bin/activate
```

After activation, the prompt should usually show `(.venv)` at the beginning.
While `.venv` is active, use `python`, not `python3`, to run the script:

```bash
python --version
minimap2 --version
samtools --version
```

If `source .venv/bin/activate` says the file does not exist, the environment
has not been set up in this folder yet. Set it up once with:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip minimap2 samtools
python3 -m venv .venv
source .venv/bin/activate
python -m pip install pysam numpy matplotlib
```

## 2. Understand Windows Paths in Ubuntu

Ubuntu can read files from the Windows `C:` drive, but the path looks different.

Windows path:

```text
C:\Users\altab\epi2me\runs\Run_2026_04_29
```

Ubuntu path:

```text
/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29
```

Rules:

- `C:\` becomes `/mnt/c/`
- Backslashes `\` become forward slashes `/`
- Wrap paths in quotes if they contain spaces
- Do not use raw `C:\...` paths inside Ubuntu commands

Example:

```bash
--metadata "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/WPS Working Sheet.xlsx"
```

## 3. Prepare the Input Folder

Start with the EPI2ME output files for a sequencing run. The FASTA, GenBank,
BAM, FASTQ, and MAF files do not need to be in the same folder. A clean layout
looks like this:

```text
C:\WPS_Runs\Run_2026_04_29\
  fasta_files\
    barcode01.final.fasta
    barcode02.final.fasta

  genbank_files\
    barcode01.annotations.gbk
    barcode02.annotations.gbk

  bam_files\
    barcode01\
      FBD...bam
    barcode02\
      FBD...bam

  fastq_files\                  optional
    barcode01.final.fastq
    barcode02.final.fastq

  maf_files\                    optional
    barcode01.assembly.maf
    barcode02.assembly.maf

  WPS_Working_Sheet_2026_04_29.xlsx
```

Each barcode must have:

- `barcodeXX.final.fasta`
- `barcodeXX.annotations.gbk`
- one raw/unmapped `.bam` file. The BAM can either have the barcode in the
  filename, or it can be inside a barcode folder such as `bam_files\barcode01\`.

Optional files:

- `barcodeXX.final.fastq`
- `barcodeXX.assembly.maf`

The metadata sheet must contain the same barcode numbers as the data files. For
example, a row with `Barcode #` equal to `1` matches `barcode01`.
The barcode does not need leading zeroes in the sheet: `3`, `03`, `3.0`, and
`barcode3` all match files named `barcode03...`.

You do not need to rename the metadata file. The report command takes the exact
metadata path with `--metadata`.

## 4. Run the Report Generator

In Ubuntu, make sure you are in the project folder and `.venv` is active:

```bash
cd /mnt/c/Users/altab/epi2me
source .venv/bin/activate
```

Then run the report command. This example processes barcode 1 and barcode 2:

```bash
python epi2me_to_final_package.py \
  --fasta-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/fasta_files" \
  --genbank-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/genbank_files" \
  --bam-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/bam_files" \
  --fastq-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/fastq_files" \
  --maf-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/maf_files" \
  --metadata "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/WPS Working Sheet.xlsx" \
  --output-dir "/mnt/c/Users/altab/epi2me/runs/Run_2026_04_29/output" \
  --barcodes 1 2 \
  --threads 4 \
  --sort-memory 1G
```

What the command parts mean:

- `python epi2me_to_final_package.py` starts the report generator.
- `--fasta-dir`, `--genbank-dir`, and `--bam-dir` point to required input folders.
- `--fastq-dir` and `--maf-dir` are optional, but include them when those files exist.
- `--metadata` points to the WPS working sheet Excel/CSV/TSV file.
- `--output-dir` is where the finished customer package will be written.
- `--barcodes 1 2` limits the run to barcode01 and barcode02. Omit this option to process every barcode found.
- `--threads 4` makes alignment faster.
- `--sort-memory 1G` gives `samtools sort` more memory.

It is okay to reuse the same `--output-dir`. If the same barcode is run again,
the script removes the previous files for that barcode and writes fresh ones.
Reports for other barcodes in the same order folder are left alone.

Most runs should not use `--keep-intermediates` or `--allow-aligned-input`.
Those are debugging/override options.

## 5. Find the Finished Reports

After the run, look inside the output folder:

```text
C:\WPS_Runs\Run_2026_04_29\output\
  WPS Data_Order #145011068\
    QC REPORTS\
      001_A39569_MYO2A-KAN_report.pdf
    FASTA_FILES\
    GENBANK_FILES\
    CHROMATOGRAM_FILES_ab1\
    PER_BASE_BREAKDOWN\

  run_summary.json
```

Open the PDF files in `QC REPORTS`.

`run_summary.json` tells you which barcodes were packaged and which were
skipped. Always check it after a run.

## 6. Try the Included Example

This repository includes a small `barcode01` example. From this code folder:

```bash
python epi2me_to_final_package.py \
  --fasta-dir "example_data/epi2me_export" \
  --genbank-dir "example_data/epi2me_export" \
  --bam-dir "example_data/epi2me_export" \
  --fastq-dir "example_data/epi2me_export" \
  --maf-dir "example_data/epi2me_export" \
  --metadata "example_data/barcode01_wps_working_sheet.csv" \
  --output-dir "example_data/output"
```

The example report will appear under:

```text
example_data\output\WPS Data_Order #145011068\QC REPORTS\
```

## macOS Option

If you are using a Mac, do not use WSL. Use Conda or Homebrew.

Conda is the most self-contained option:

```bash
conda create -n wps-report -c conda-forge -c bioconda python=3.11 pysam numpy matplotlib minimap2 samtools
conda activate wps-report
```

Homebrew also works:

```bash
brew install python minimap2 samtools
python3 -m pip install pysam numpy matplotlib
```

The `brew install` or `conda create` command installs `minimap2` and
`samtools`. The `python3 -m pip install ...` command only installs Python
packages.

Run the script on macOS with Unix-style paths:

```bash
python3 epi2me_to_final_package.py \
  --fasta-dir /Users/yourname/WPS_Runs/Run_2026_04_29/fasta_files \
  --genbank-dir /Users/yourname/WPS_Runs/Run_2026_04_29/genbank_files \
  --bam-dir /Users/yourname/WPS_Runs/Run_2026_04_29/bam_files \
  --fastq-dir /Users/yourname/WPS_Runs/Run_2026_04_29/fastq_files \
  --maf-dir /Users/yourname/WPS_Runs/Run_2026_04_29/maf_files \
  --metadata /Users/yourname/WPS_Runs/Run_2026_04_29/WPS_Working_Sheet_2026_04_29.xlsx \
  --output-dir /Users/yourname/WPS_Runs/Run_2026_04_29/output \
  --threads 4 \
  --sort-memory 1G
```

## Installing Ubuntu/WSL

If Ubuntu is not installed yet, install WSL from PowerShell as Administrator:

```powershell
wsl --install
```

Then open **Ubuntu** from the Start menu and follow the main instructions at the
top of this README.

## Troubleshooting

If a barcode is skipped, check `run_summary.json`.

Common causes:

- The barcode exists in EPI2ME files but not in the WPS sheet.
- The WPS sheet has a barcode row but the matching EPI2ME files are missing.
- More than one BAM file was found for the same barcode.
- The BAM is already aligned instead of raw/unmapped.
- Two samples produce the same filename after cleanup.

If you intentionally need to debug intermediate alignment files, add:

```bat
--keep-intermediates
```

If you intentionally need to realign a BAM that already contains mapped reads,
add:

```bat
--allow-aligned-input
```

Use that only when you are sure the BAM is supposed to be realigned.

## Files in This Folder

- `epi2me_to_final_package.py`: main command to run
- `align_bam_pipeline.py`: converts raw BAM reads to FASTQ and aligns them
- `bam_to_per_base_data.py`: creates per-base support CSVs
- `generate_report.py`: computes QC metrics and plots
- `fastq_to_ab1.py`: creates one synthetic AB1 chromatogram per sample
- `example_data\`: small test data and example metadata
