#!/bin/bash

# Paths to raw FASTQ
FASTQ1="../data/raw/SRR36145750_1.fastq.gz"
FASTQ2="../data/raw/SRR36145750_2.fastq.gz"

# Output folder
OUT="../results/QC"

# Run FastQC
fastqc -o $OUT $FASTQ1 $FASTQ2

# Summarize all reports
multiqc -o $OUT $OUT
