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

2. Copy and paste this command block:

```bash
cd /mnt/c/Users/altab/plasmid_report
git pull
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

## 2. Input Folder Location

The report command assumes all input run folders live under:

```text
/mnt/c/WPS data/
```

## 3. Prepare the Input Folder

Start with the EPI2ME output files for a sequencing run. Put the FASTA,
GenBank, BAM, FASTQ, metadata, and any optional MAF files in one folder under
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
  barcode01.final.fastq
  barcode02.final.fastq
  barcode01.assembly.maf        optional
  barcode02.assembly.maf        optional
  WPS_Working_Sheet_2026_04_29.xlsx
```

Each barcode must have:

- `barcodeXX.final.fasta`
- `barcodeXX.annotations.gbk`
- `barcodeXX.final.fastq`
- one raw/unmapped `.bam` file whose filename contains the barcode, such as
  `FBD...barcode01...bam`, or one raw/unmapped `.bam` file inside a barcode
  folder such as `barcode02\FBD...bam`.

Optional files:

- `barcodeXX.assembly.maf`

If a FASTQ is missing, the run still completes and the AB1 is generated from the
FASTA with default quality scores, but the sample is reported with a warning.

The metadata sheet must contain the same barcode numbers as the data files. For
example, a row with `Barcode #` equal to `1` matches `barcode01`.
The barcode does not need leading zeroes in the sheet: `3`, `03`, `3.0`, and
`barcode3` all match files named `barcode03...`.

The metadata file should be named `WPS Working Sheet` or something very similar,
such as `WPS_Working_Sheet_2026_04_29.xlsx`. The input folder must contain
exactly one matching metadata `.xlsx`, `.csv`, or `.tsv` file.

## 4. Mixed-Contig Samples

Mixed or contaminated samples can produce more than one contig for the same
barcode and metadata row. The script handles this by generating one report and
one output file set per contig.

The clearest input style is to include an explicit contig suffix in each
contig-specific filename:

```text
C:\WPS data\Run_2026_04_29\
  barcode01.contig001.final.fasta
  barcode01.contig002.final.fasta
  barcode01.contig001.annotations.gbk
  barcode01.contig002.annotations.gbk
  FBD...barcode01...bam
  barcode01.final.fastq
  WPS_Working_Sheet_2026_04_29.xlsx
```

This produces separate outputs such as `sample_contig001.fa`,
`sample_contig002.fa`, and one report PDF per contig. A single barcode-level
BAM or FASTQ is reused for each contig. Contig-specific BAM, FASTQ, or MAF files
are also supported when their filenames include the same `contig001`,
`contig002`, etc. suffix.

The script also has a looser fallback for files that do not include contig
suffixes. If multiple `final` FASTA files are found for the same barcode, it
prints that it thinks the sample has multiple contigs and assigns `contig001`,
`contig002`, etc. in sorted filename order. If there are also multiple
unlabelled GenBank, BAM, FASTQ, or MAF files with the same barcode, they are
paired to those inferred contigs by the same sorted order when the counts match.

Explicit contig suffixes are still safer when possible, because sorted filename
order is only a fallback.

## 5. Run Report

Example command:

```bash
python3 epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29"
```

This processes every barcode found in `C:\WPS data\Run_2026_04_29\` and writes
the output to `C:\WPS data\Run_2026_04_29\output\`.

### Optional Commands

Use these options only when you need to choose specific barcodes, change the
output folder, or adjust alignment settings:

```bash
python3 epi2me_to_final_package.py \
  --folder-name "Run_2026_04_29" \
  --output-dir "/mnt/c/WPS data/Run_2026_04_29/output" \
  --barcodes 1 2 \
  --threads 4 \
  --sort-memory 1G
```

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

## Troubleshooting

If a barcode is skipped, check `run_summary.json`.

Common causes:

- The barcode exists in EPI2ME files but not in the WPS sheet.
- The WPS sheet has a barcode row but the matching EPI2ME files are missing.
- More than one unpaired file of the same type was found for the same barcode.
  Use `contig001`, `contig002`, etc. in filenames for mixed-contig samples.
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
