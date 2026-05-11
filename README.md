# EPI2ME to WPS Report Workflow

This folder turns an EPI2ME/ONT plasmid sequencing export into customer-facing
WPS order folders with one PDF report per sample.

If several samples have the same `Order #` in the metadata sheet, their reports
go into the same `WPS Data_Order #...` folder. Samples with different order
numbers go into separate folders.

These instructions assume you are using Windows.

## 1. Install Required Software

The easiest Windows setup is Miniconda or Anaconda. It installs Python and the
bioinformatics tools in one environment.

1. Install Miniconda for Windows:
   <https://docs.conda.io/en/latest/miniconda.html>

2. Open **Anaconda Prompt** from the Start menu.

3. Create the environment:

```bat
conda create -n wps-report -c conda-forge -c bioconda python=3.11 pysam numpy matplotlib minimap2 samtools
```

4. Activate it:

```bat
conda activate wps-report
```

5. Check that everything is available:

```bat
python --version
minimap2 --version
samtools --version
```

If `minimap2` or `samtools` does not work on native Windows, use the WSL option
near the bottom of this README. WSL is often the most reliable Windows setup for
bioinformatics command-line tools.

## 2. Prepare the Input Folder

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
    FBD...barcode01...bam
    FBD...barcode02...bam

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
- one raw/unmapped `.bam` file with that barcode in the filename

Optional files:

- `barcodeXX.final.fastq`
- `barcodeXX.assembly.maf`

The metadata sheet must contain the same barcode numbers as the data files. For
example, a row with `Barcode #` equal to `1` matches `barcode01`.
The barcode does not need leading zeroes in the sheet: `3`, `03`, `3.0`, and
`barcode3` all match files named `barcode03...`.

You do not need to rename the metadata file. The report command takes the exact
metadata path with `--metadata`.

## 3. Run the Report Generator on Windows

Open **Anaconda Prompt**.

Move into this code folder. Example:

```bat
cd C:\Users\YourName\final_epi2me
```

Run the report generator:

```bat
python epi2me_to_final_package.py ^
  --fasta-dir "C:\WPS_Runs\Run_2026_04_29\fasta_files" ^
  --genbank-dir "C:\WPS_Runs\Run_2026_04_29\genbank_files" ^
  --bam-dir "C:\WPS_Runs\Run_2026_04_29\bam_files" ^
  --metadata "C:\WPS_Runs\Run_2026_04_29\WPS_Working_Sheet_2026_04_29.xlsx" ^
  --output-dir "C:\WPS_Runs\Run_2026_04_29\output"
```

For faster alignment on a larger run, use more threads:

```bat
python epi2me_to_final_package.py ^
  --fasta-dir "C:\WPS_Runs\Run_2026_04_29\fasta_files" ^
  --genbank-dir "C:\WPS_Runs\Run_2026_04_29\genbank_files" ^
  --bam-dir "C:\WPS_Runs\Run_2026_04_29\bam_files" ^
  --metadata "C:\WPS_Runs\Run_2026_04_29\WPS_Working_Sheet_2026_04_29.xlsx" ^
  --output-dir "C:\WPS_Runs\Run_2026_04_29\output" ^
  --threads 4 ^
  --sort-memory 1G
```

The metadata file can also be a CSV or TSV if it has the expected columns:

```bat
python epi2me_to_final_package.py ^
  --fasta-dir "C:\WPS_Runs\Run_2026_04_29\fasta_files" ^
  --genbank-dir "C:\WPS_Runs\Run_2026_04_29\genbank_files" ^
  --bam-dir "C:\WPS_Runs\Run_2026_04_29\bam_files" ^
  --metadata "C:\WPS_Runs\Run_2026_04_29\metadata.csv" ^
  --output-dir "C:\WPS_Runs\Run_2026_04_29\output"
```

PowerShell users can run the same command on one line, or use backticks instead
of `^` for line continuation.

Example command with every supported option:

```bat
python epi2me_to_final_package.py ^
  --fasta-dir "C:\WPS_Runs\Run_2026_04_29\fasta_files" ^
  --genbank-dir "C:\WPS_Runs\Run_2026_04_29\genbank_files" ^
  --bam-dir "C:\WPS_Runs\Run_2026_04_29\bam_files" ^
  --fastq-dir "C:\WPS_Runs\Run_2026_04_29\fastq_files" ^
  --maf-dir "C:\WPS_Runs\Run_2026_04_29\maf_files" ^
  --metadata "C:\WPS_Runs\Run_2026_04_29\WPS_Working_Sheet_2026_04_29.xlsx" ^
  --output-dir "C:\WPS_Runs\Run_2026_04_29\output" ^
  --barcodes barcode01 barcode02 ^
  --logo "C:\WPS_Runs\alta_logo.png" ^
  --threads 4 ^
  --sort-memory 1G ^
  --keep-intermediates ^
  --allow-aligned-input
```

Most runs should not use `--keep-intermediates` or `--allow-aligned-input`.
Those are debugging/override options.

## 4. Find the Finished Reports

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

## 5. Try the Included Example

This repository includes a small `barcode01` example. From this code folder:

```bat
python epi2me_to_final_package.py ^
  --fasta-dir "example_data\epi2me_export" ^
  --genbank-dir "example_data\epi2me_export" ^
  --bam-dir "example_data\epi2me_export" ^
  --fastq-dir "example_data\epi2me_export" ^
  --maf-dir "example_data\epi2me_export" ^
  --metadata "example_data\barcode01_wps_working_sheet.csv" ^
  --output-dir "example_data\output"
```

The example report will appear under:

```text
example_data\output\WPS Data_Order #145011068\QC REPORTS\
```

## WSL Option

If native Windows has trouble installing or running `minimap2` or `samtools`,
use Windows Subsystem for Linux.

1. Install WSL from PowerShell as Administrator:

```powershell
wsl --install
```

2. Open Ubuntu from the Start menu.

3. Install dependencies:

```bash
sudo apt update
sudo apt install python3 python3-pip minimap2 samtools
python3 -m pip install pysam numpy matplotlib
```

4. Go to the Windows folder from WSL. A Windows path like:

```text
C:\Users\YourName\final_epi2me
```

appears in WSL as:

```text
/mnt/c/Users/YourName/final_epi2me
```

5. Run the script:

```bash
cd /mnt/c/Users/YourName/final_epi2me
python3 epi2me_to_final_package.py \
  --fasta-dir /mnt/c/WPS_Runs/Run_2026_04_29/fasta_files \
  --genbank-dir /mnt/c/WPS_Runs/Run_2026_04_29/genbank_files \
  --bam-dir /mnt/c/WPS_Runs/Run_2026_04_29/bam_files \
  --fastq-dir /mnt/c/WPS_Runs/Run_2026_04_29/fastq_files \
  --maf-dir /mnt/c/WPS_Runs/Run_2026_04_29/maf_files \
  --metadata /mnt/c/WPS_Runs/Run_2026_04_29/WPS_Working_Sheet_2026_04_29.xlsx \
  --output-dir /mnt/c/WPS_Runs/Run_2026_04_29/output
```

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
