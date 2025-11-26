#!/bin/bash

RAW1="../data/raw/SRR36145750_1.fastq.gz"
RAW2="../data/raw/SRR36145750_2.fastq.gz"
OUT1="../data/processed/SRR36145750_trim_1.fastq.gz"
OUT2="../data/processed/SRR36145750_trim_2.fastq.gz"

FWD_PRIMER="CCTACGGGNGGCWGCAG"
REV_PRIMER="GACTACHVGGGTATCTAATCC"

cutadapt -g $FWD_PRIMER -G $REV_PRIMER -q 20,20 -m 100 -o $OUT1 -p $OUT2 $RAW1 $RAW2
