# Cryptosporidium MLST Nanopore Pipeline

This repository contains an automated pipeline for Multi-Locus Sequence Typing (MLST) of *Cryptosporidium* based on Oxford Nanopore sequencing of PCR-amplified MLST loci. The pipeline is designed to be used within **iRODS** as an automated workflow, but can also be executed manually from the command line.

The pipeline assumes that each barcode directory corresponds to one sample, and it produces per-sample MLST consensus sequences, quality-control reports, and coverage plots.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
  - [Manual execution](#manual-execution)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
- [Detailed Pipeline Description](#detailed-pipeline-description)
  - [Top-level driver script](#top-level-driver-script)
  - [Per-barcode processing](#per-barcode-processing)
- [Coverage Plot Script](#coverage-plot-script)
- [Mpileup Parse Script](#mpileup-parse-script)
- [MLST Reference](#mlst-reference)
---

## Overview

### Context

- **Goal**: MLST analysis of *Cryptosporidium* using amplicon-based Nanopore sequencing.
- **Input**: Nanopore FASTQ files (gzip-compressed) for multiple barcodes; each barcode represents one sample.
- **Core tools**: `minimap2`, `samtools`, `bcftools`, `pandas`, `matplotlib`.
- **Environment**: Conda/mamba environment defined in `envs/mlst_env.conda.yml`.
- **Automation**: Compatible with iRODS, allowing automated triggering based on data arrival.

### High-level workflow

For each barcode (sample):

1. Concatenate all FASTQ files for that barcode.
2. Map reads to a multi-gene MLST reference.
3. Generate an initial consensus sequence using `bcftools`.
4. Remap reads to the initial consensus.
5. Generate a pileup on the consensus and parse base support metrics.
6. Construct a refined final consensus, trimming low-support bases and ends.
7. Identify “dubious” positions (low fraction or low coverage) for QC.
8. Produce an alignment depth plot.
9. Clean up intermediate files.

---

## Repository Structure

```text
.
├── run_mlst_pipeline.sh           # Top-level driver script
├── envs
│   └── mlst_env.conda.yml         # Conda/mamba environment specification
└── src
    ├── get_mlst_consensus_for_iRODS.sh  # Main per-barcode MLST pipeline
    ├── fasta_seq_mlst.fasta             # MLST reference sequences
    ├── parse_mpileup.py                 # Generates base-support statistics
    └── create_depth_plot.py             # Coverage plot generator

```

---

## Usage

### Manual execution

To run the full pipeline manually:
```bash
bash run_mlst_pipeline.sh /path/to/INPUTDIR /path/to/OUTPUTDIR
```
#### Arguments


- `INPUTDIR`  
  Directory containing one subdirectory per barcode. Each barcode directory holds all `.fastq.gz` files for that sample:

```text
  INPUTDIR/  
  ├── barcode01/  
  │   ├── file1.fastq.gz  
  │   └── file2.fastq.gz  
  ├── barcode02/  
  │   └── ...  
  └── ...

```
- `OUTPUTDIR`  
  Directory where output files (consensus FASTA, QC files, plots, flags) will be written. Must exist before running the script.


If either `INPUTDIR` or `OUTPUTDIR` is missing or does not exist, the script will terminate with an error.

### Inputs

- Nanopore **FASTQ** files (`*.fastq.gz`) of PCR-amplified MLST genes.
- One sample per barcode (one directory per barcode in `INPUTDIR`).
- MLST reference FASTA: `src/fasta_seq_mlst.fasta`.

### Outputs

For a barcode directory named `barcode01`, the pipeline will create:

- **Final consensus FASTA**  
  `OUTPUTDIR/barcode01_consensus_final.fasta`

- **Dubious bases report** (QC)  
  `OUTPUTDIR/barcode01_consensus_final_dubious_bases.txt`  
  Contains positions with low consensus fraction or low coverage that may warrant manual review.

- **Coverage plot PDF**  
  `OUTPUTDIR/barcode01_depth_plot.pdf`  
  Depth of coverage across the MLST loci.

If no reads map for a barcode, a flag file is created instead:

- `OUTPUTDIR/no_reads_mapped_for_barcode_barcode01.txt`

Most intermediate files are deleted automatically to keep the output directory clean.

---

## Detailed Pipeline Description

### Top-level driver script

File: `run_mlst_pipeline.sh`

Core steps:

1. **Set working directory**

   The script determines its own location and changes the current directory to that location:
```bash
   DIR="\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" > /dev/null 2>&1 && pwd )"  
   cd "\${DIR}"
```
2. **Check input arguments**
```bash 
   if [ ! -z "\${1}" ] || [ ! -z "\${2}" ] ; then  
       INPUTDIR=\$1  
       OUTPUTDIR=\$2  
   else  
       echo "Error: no inputdir or outputdir given."  
       exit 1  
   fi  

   if [ ! -d "\${INPUTDIR}" ] || [ ! -d "\${OUTPUTDIR}" ] ; then  
       echo "Error: inputdir \$INPUTDIR or output dir \$OUTPUTDIR not exists"  
       exit 2  
   fi
```
3. **Set up environment**
```bash
   mamba env update -f envs/mlst_env.conda.yml -q  
   source activate mlst_env
```
4. **Run the main MLST pipeline**
```bash
   bash src/get_mlst_consensus_for_iRODS.sh "\${INPUTDIR}" "\${OUTPUTDIR}"
```
---

### Per-barcode processing

File: `src/get_mlst_consensus_for_iRODS.sh`

The script is called with:
```bash
bash src/get_mlst_consensus_for_iRODS.sh INPUTDIR OUTPUTDIR
```
#### 1. Initial setup
```bash
input_dir=\$1  
output_dir=\$2  

reference="src/fasta_seq_mlst.fasta"
```
- `reference` is the MLST gene reference FASTA file.

#### 2. Iterate over barcode directories
```
for barcode in \$input_dir/*; do  
    if [ -d "\$barcode" ]; then  
        echo "Processing directory: \$barcode"  
        ...  
    fi  
done
```
Each subdirectory under `input_dir` is treated as a separate barcode/sample.

#### 3. Concatenate FASTQ files per barcode
```bash
cat \$barcode/*.fastq.gz > \$output_dir/\$(basename \$barcode)\_raw.fastq.gz

- All `.fastq.gz` files in `barcode` are concatenated into a single file:
  - e.g. `OUTPUTDIR/barcode01_raw.fastq.gz`
```
#### 4. Map reads to MLST reference and sort BAM
```bash
minimap2 -ax map-ont \$reference \  
    \$output_dir/\$(basename \$barcode)\_raw.fastq.gz | \  
    samtools view -u | \  
    samtools sort -o \$output_dir/\$(basename \$barcode)\_sorted.bam
```
- `minimap2 -ax map-ont`: align Nanopore reads to the MLST gene reference.
- `samtools view -u`: convert SAM to BAM.
- `samtools sort`: sort BAM and write:
  - `OUTPUTDIR/barcode01_sorted.bam`

#### 5. Check if any reads mapped
```bash
if [ "\$(samtools view -c -F 4 \$output_dir/\$(basename \$barcode)\_sorted.bam)" -gt 0 ]; then  
    # at least one mapped read  
    ...  
else  
    rm -f "\$output_dir/\$(basename \$barcode)\_raw.fastq.gz"  
    rm -f "\$output_dir/\$(basename \$barcode)\_sorted.bam"  
    rm -f "\$output_dir/\$(basename \$barcode)\_sorted.bam.bai"  
    touch "\$output_dir/no_reads_mapped_for_barcode_\$(basename \$barcode).txt"  
fi
```
- `samtools view -c -F 4` counts reads that are **not** marked as unmapped.
- If count is zero:
  - Remove temporary files.
  - Write a flag file indicating no reads mapped for this barcode.

#### 6. First consensus construction

If reads are present:

1. **Variant calling with `bcftools`**
```bash
   bcftools mpileup -Ou --fasta-ref \$reference \  
       \$output_dir/\$(basename \$barcode)\_sorted.bam | \  
       bcftools call -c -Ob -o \$output_dir/\$(basename \$barcode)\_call.bcf.gz  

   bcftools index \$output_dir/\$(basename \$barcode)\_call.bcf.gz
```
   - `bcftools mpileup`: generate pileup from the BAM using the reference.
   - `bcftools call -c`: perform variant calling (consensus-based calling) and output BCF.
   - `bcftools index`: index the BCF for efficient access.

2. **Generate consensus FASTA**
```bash
   bcftools consensus -f \$reference \  
       -o \$output_dir/\$(basename \$barcode)\_consensus.fasta \  
       \$output_dir/\$(basename \$barcode)\_call.bcf.gz
```
   - Creates a consensus FASTA (`_consensus.fasta`) based on the variants and the original MLST reference.

#### 7. Remap reads to the initial consensus
```bash
minimap2 -ax map-ont \  
    \$output_dir/\$(basename \$barcode)\_consensus.fasta \  
    \$output_dir/\$(basename \$barcode)\_raw.fastq.gz \  
    | samtools view -u - \  
    | samtools sort -o \$output_dir/\$(basename \$barcode)\_consensus_mapped.bam  

samtools index \$output_dir/\$(basename \$barcode)\_consensus_mapped.bam
```
- Align the original reads to this consensus sequence.
- Produce a sorted and indexed BAM (`_consensus_mapped.bam`).

#### 8. Generate pileup on the consensus
```bash
samtools mpileup -aa --reference \  
    \$output_dir/\$(basename \$barcode)\_consensus.fasta \  
    \$output_dir/\$(basename \$barcode)\_consensus_mapped.bam \  
    > \$output_dir/\$(basename \$barcode)\_consensus_mapped.pileup
```
- `-aa` ensures positions with zero coverage are also included.
- Output: `_consensus_mapped.pileup` containing per-position information on the consensus.

#### 9. Parse mpileup and create base support file
```bash
python3 src/parse_mpileup.py \  
    "\$output_dir/\$(basename \$barcode)\_consensus_mapped.pileup" \  
    -o "\$output_dir/\$(basename \$barcode)\_consensus_base_support.txt"
```
- The external script `parse_mpileup.py` parses the pileup and writes a tab-separated table:
  - `OUTPUTDIR/barcode01_consensus_base_support.txt`
- Each row likely contains:
  - `chrom`, `pos`, `coverage`, `fraction`, `cons_base`, etc.

#### 10. Build final consensus and mark dubious bases

A small in-line Python script refines the consensus:
```bash
python3 <<EOF  
import pandas as pd  

base_support = pd.read_csv("\$output_dir/\$(basename \$barcode)\_consensus_base_support.txt", sep="\t")  

file_out = open("\$output_dir/\$(basename \$barcode)\_consensus_final.fasta", "w")  
for gene in base_support["chrom"].unique():  
    print(f">{gene}", file=file_out)  
    print("".join(list(base_support.loc[  
        (base_support["chrom"] == gene) &  
        # Remove all bases that are not supported by reads  
        ~((base_support["fraction"] < 0.01) & (base_support["coverage"] > 0)),  
        # Remove the first and last 10 bases  
        "cons_base"  
    ]))[10:-10], file=file_out)  
file_out.close()  

# Identify dubious variants and low coverage bases  
dubious_variants = base_support.loc[  
    (  
        # fraction < 80% but >= 1%  
        ((base_support["fraction"] < 0.8) & (base_support["fraction"] >= 0.01)) |  
        # coverage < 20 (excluding bases with coverage == 0)  
        ((base_support["coverage"] < 20) & (base_support["coverage"] != 0))  
    ) &  
    # Exclude first 10 bases and gene-specific tail regions  
    ~(  
        (base_support["pos"] < 11) |  
        (((base_support["chrom"] == "CP003930") & (base_support["pos"] > 441)) |  
         ((base_support["chrom"] == "CP0028230") & (base_support["pos"] > 466)) |  
         ((base_support["chrom"] == "CP0031960") & (base_support["pos"] > 577)) |  
         ((base_support["chrom"] == "CP0021750") & (base_support["pos"] > 497)) |  
         ((base_support["chrom"] == "CP0024650") & (base_support["pos"] > 516)) |  
         ((base_support["chrom"] == "CP0012400") & (base_support["pos"] > 493)) |  
         ((base_support["chrom"] == "CP007000")  & (base_support["pos"] > 548)) |  
         ((base_support["chrom"] == "CP001010")  & (base_support["pos"] > 543)))  
    )  
]  

dubious_variants.to_csv(  
    "\$output_dir/\$(basename \$barcode)\_consensus_final_dubious_bases.txt",  
    sep="\t",  
    index=None  
)  
EOF
```
This script performs:

- **Final consensus generation**:
  - For each MLST gene (`chrom`):
    - Remove bases where:
      - `fraction < 0.01` and `coverage > 0` (unsupported bases).
    - Remove the first and last 10 bases of each gene (potential primer regions or low-quality tails).
  - Write the resulting sequence as:
    - `OUTPUTDIR/barcode01_consensus_final.fasta`

- **Dubious bases detection**:
  - Select positions where:
    - `fraction < 0.8` but `>= 0.01`, or
    - `coverage < 20` and `coverage != 0`.
  - Exclude:
    - First 10 bases (`pos < 11`).
    - Gene-specific tail positions beyond known lengths for each MLST gene.
  - Output a tab-separated table:
    - `OUTPUTDIR/barcode01_consensus_final_dubious_bases.txt`
  - This file can be used for downstream QC and manual inspection.

#### 11. Generate depth plot
```bash
samtools depth -a \$output_dir/\$(basename \$barcode)\_sorted.bam \  
    > \$output_dir/\$(basename \$barcode)\_depth.txt  

python src/create_depth_plot.py \  
    \$output_dir/\$(basename \$barcode)\_depth.txt \  
    \$output_dir/\$(basename \$barcode)\_depth_plot.pdf
```
- `samtools depth -a`:
  - Produces per-position coverage across the MLST reference for the original mapped BAM:
    - `OUTPUTDIR/barcode01_depth.txt`
- `create_depth_plot.py`:
  - Uses `matplotlib` to generate a PDF with coverage vs. position per gene:
    - `OUTPUTDIR/barcode01_depth_plot.pdf`

#### 12. Cleanup of intermediate files

For barcodes with mapped reads:
```bash
rm \$output_dir/\$(basename \$barcode)\_depth.txt  
rm \$output_dir/\$(basename \$barcode)\_call.bcf.gz  
rm \$output_dir/\$(basename \$barcode)\_call.bcf.gz.csi  
rm \$output_dir/\$(basename \$barcode)\*.bam  
rm \$output_dir/\$(basename \$barcode)\*.bam.bai  
rm \$output_dir/\$(basename \$barcode)\*.pileup  
rm \$output_dir/\$(basename \$barcode)\_raw.fastq.gz  
rm \$output_dir/\$(basename \$barcode)\_consensus.fasta  
rm \$output_dir/\$(basename \$barcode)\_consensus.fasta.fai
```
- Removes intermediate BAM, BCF, pileup, raw FASTQ, and first consensus FASTA files.
- Leaves only:
  - `*_consensus_final.fasta`
  - `*_consensus_final_dubious_bases.txt`
  - `*_depth_plot.pdf`

---

## Coverage Plot Script

File: `src/create_depth_plot.py`

This script generates a coverage (depth) plot per gene and saves it as a PDF.

### Usage

python src/create_depth_plot.py input_depth.txt output_plot.pdf

Where:

- `input_depth.txt` is the output from `samtools depth -a`, with columns:
  - `gene` (reference sequence ID)
  - `position`
  - `depth`
- `output_plot.pdf` is the desired path for the PDF output.

### Internal logic

1. Read the TSV file with `pandas`.
2. Rename columns to: `["gene", "position", "depth"]`.
3. Use `matplotlib` to:
   - Group data by `gene`.
   - Plot depth vs. position for each gene in a different color (using `viridis` colormap with 8 colors).
   - Add vertical dashed lines for gene lengths (hard-coded as `[553, 503, 507, 526, 476, 587, 451, 558]`).
4. Add labels and legend, then save the figure into a PDF using `PdfPages`.

The script prints a message indicating the output PDF path when finished.

---

## Mpileup Parse Script

File: `src/parse_mpileup.py`

This script parses a `samtools mpileup` file (generated on the consensus-mapped BAM) into a per-position base-support table. The output is used by the inline Python code to refine the consensus and flag dubious positions.

### Usage
```bash
python src/parse_mpileup.py input.pileup -o output_base_support.txt
```
- `input.pileup`: produced by `samtools mpileup -aa --reference <consensus.fasta> <consensus_mapped.bam>`
- `output_base_support.txt`: tab-separated table used in the subsequent consensus-refinement step  

### Logic and output

For each pileup line (`chrom`, `pos`, `ref_base`, `coverage`, `pile_bases`), the script:

- Parses the mpileup string:
  - Treats `.` and `,` as matches to the consensus base.
  - Counts explicit bases `A/C/G/T`.
  - Skips mpileup markers (`^` + quality, `$`) and inserted bases (`+<len><seq>`).
  - Counts deletions (`-<len><seq>`) as `DEL` events.
- Computes:
  - `support`: number of reads matching the consensus base.
  - `fraction`: `support / coverage` (0 if coverage is 0).
  - `OTHER`: residual coverage not accounted for by `A/C/G/T/DEL`.

The output file has columns:
```
chrom, pos, cons_base, coverage, A, C, G, T, DEL, OTHER, fraction
```
This `*_consensus_base_support.txt` file is then used by the inline Python code to:

- filter bases with very low support (fraction < 0.01),
- trim 10 bases from each gene’s ends,
- build the final consensus FASTA,
- and select dubious bases for QC (`*_consensus_final_dubious_bases.txt`).

---

## MLST Reference

File: `src/fasta_seq_mlst.fasta`

Contains MLST loci for *Cryptosporidium* with gene identifiers such as:

- `CP003930`
- `CP0028230`
- `CP0031960`
- `CP0021750`
- `CP0024650`
- `CP0012400`
- `CP007000`
- `CP001010`

These sequences are used as the mapping reference and as the base for consensus construction.

---

