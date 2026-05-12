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

The report command assumes all input run folders live under:

```text
/mnt/c/WPS data/
```

That is the Ubuntu path for:

```text
C:\WPS data\
```

## 3. Prepare the Input Folder

Start with the EPI2ME output files for a sequencing run. Put the FASTA,
GenBank, BAM, metadata, and any optional FASTQ/MAF files in one folder under
`C:\WPS data\`. A clean layout looks like this:

```text
C:\WPS data\Run_2026_04_29\
  barcode01.final.fasta
  barcode02.final.fasta
  barcode01.annotations.gbk
  barcode02.annotations.gbk
  FBD...barcode01...bam
  barcode02\
    FBD...bam
  barcode01.final.fastq         optional
  barcode02.final.fastq         optional
  barcode01.assembly.maf        optional
  barcode02.assembly.maf        optional
  WPS_Working_Sheet_2026_04_29.xlsx
```

Each barcode must have:

- `barcodeXX.final.fasta`
- `barcodeXX.annotations.gbk`
- one raw/unmapped `.bam` file whose filename contains the barcode, such as
  `FBD...barcode01...bam`, or one raw/unmapped `.bam` file inside a barcode
  folder such as `barcode02\FBD...bam`.

Optional files:

- `barcodeXX.final.fastq`
- `barcodeXX.assembly.maf`

The metadata sheet must contain the same barcode numbers as the data files. For
example, a row with `Barcode #` equal to `1` matches `barcode01`.
The barcode does not need leading zeroes in the sheet: `3`, `03`, `3.0`, and
`barcode3` all match files named `barcode03...`.

You do not need to rename the metadata file. The input folder must contain
exactly one metadata `.xlsx`, `.csv`, or `.tsv` file.

## 4. Run the Report Generator

In Ubuntu, make sure you are in the project folder and `.venv` is active:

```bash
cd /mnt/c/Users/altab/epi2me
source .venv/bin/activate
```

Minimum command using the folder name and barcode filter:

```bash
python epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29" \
  --barcodes 1 2
```

Then run the fuller report command if you want to choose an output folder or
alignment settings. Optional FASTQ/MAF files are discovered from the same input
folder under `C:\WPS data\` when they are present. This example processes
barcode 1 and barcode 2:

```bash
python epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29" \
  --output-dir "/mnt/c/WPS data/Run_2026_04_29/output" \
  --barcodes 1 2 \
  --threads 4 \
  --sort-memory 1G
```

What the command parts mean:

- `python epi2me_to_final_package.py` starts the report generator.
- `--folder-name` names the folder under `/mnt/c/WPS data/` containing all run input files.
- `--output-dir` is where the finished customer package will be written. If
  omitted, the output goes into `C:\WPS data\<folder-name>\output\`.
- `--barcodes 1 2` limits the run to barcode01 and barcode02. Omit this option to process every barcode found.
- `--multimer-denominator classified-reads` reports monomer/dimer/trimer/tetramer percentages only among reads that were close enough to 1x/2x/3x/4x plasmid length to classify. This is the default.
- `--multimer-denominator all-eligible-reads` includes eligible mapped reads that were not classifiable and adds an `Unclassified` column to the multimer table.
- `--threads 4` makes alignment faster.
- `--sort-memory 1G` gives `samtools sort` more memory.

It is okay to reuse the same `--output-dir`. If the same barcode is run again,
the script removes the previous files for that barcode and writes fresh ones.
Reports for other barcodes in the same order folder are left alone.
The default `output` folder is ignored during input discovery, so rerunning the
same folder will not treat generated alignment files as new input BAMs.

Most runs should not use `--keep-intermediates` or `--allow-aligned-input`.
Those are debugging/override options.

## 5. Find the Finished Reports

After the run, look inside the output folder:

```text
C:\WPS data\Run_2026_04_29\output\
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

Copy or place the example folder at `C:\WPS data\epi2me_export`, then run:

```bash
python epi2me_to_final_package.py \
  --folder-name "epi2me_export" \
  --output-dir "example_data/output"
```

The example report will appear under:

```text
example_data\output\WPS Data_Order #145011068\QC REPORTS\
```

## macOS Development Option

If you are testing or developing on a Mac, use Conda or Homebrew.

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

Routine report runs are configured for the WSL input root `/mnt/c/WPS data/`.
For those runs, use Ubuntu/WSL and pass only the folder name:

```bash
python3 epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29" \
  --output-dir "/mnt/c/WPS data/Run_2026_04_29/output" \
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
