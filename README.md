# Processing and analysis of plant-associated microbial communities from Illumina 16S rRNA sequencing

## Project Goal
This project processes and analyzes plant-associated microbial communities using Illumina 16S rRNA sequencing data. The goal is to generate high-quality ASVs (Amplicon Sequence Variants), assign taxonomy, and visualize microbial diversity.

## Dataset
- **Source:** NCBI SRA  
- **Accession:** [SRX31186795](https://www.ncbi.nlm.nih.gov/sra/SRX31186795)  
- **Description:** Paired-end Illumina reads of plant-associated microbial communities.

## Workflow
1. **Quality Control:** FastQC & MultiQC on raw FASTQ files using Bash scripts
2. **Trimming:** Remove primers and low-quality bases using Cutadapt (`trim.sh`)
3. **DADA2 Analysis (R):**
   - Filter and dereplicate reads
   - Learn error rates
   - Infer ASVs
   - Merge paired reads and remove chimeras
   - Assign taxonomy using SILVA DADA2-formatted databases
4. **Diversity Analysis & Visualization (R):**
   - Alpha diversity (Shannon, Simpson, Observed)
   - Taxonomic composition barplots (Phylum, Genus)
   - Optional beta diversity and heatmaps
5. **Results:** ASV tables, taxonomy tables, and plots saved in the `results/` folder

## Software & Tools
- **Bash:** FastQC, MultiQC, Cutadapt
- **R:** DADA2, phyloseq, ggplot2

## Folder Structure
Illumina_Project/
├── data/
│   ├── raw/            # raw FASTQ files
│   └── processed/      # trimmed FASTQ files
├── results/
│   ├── ASV_tables/     # ASV and taxonomy tables
│   ├── QC/             # QC reports
│   ├── diversity/      # alpha/beta diversity plots
│   └── plots/          # taxonomic composition plots
├── scripts/
│   ├── bash/           # QC & trimming scripts
│   │   ├── qc.sh
│   │   └── trim.sh
│   └── R/              # DADA2 pipeline script
│       └── dada2_pipeline.R
└── README.md

## How to Reproduce
1. Run QC and trim Bash scripts in `scripts/bash/`  
   ```bash
   ./qc.sh
   ./trim.sh```
2. Run the DADA2 R script in scripts/R/
	- source("scripts/R/dada2_pipeline.R")
3. Results (tables & plots) will be automatically saved in results/

Notes

Download the SILVA DADA2-formatted databases and place them in data/:

silva_nr99_v138.2_toGenus_trainset.fa.gz

silva_v138.2_assignSpecies.fa.gz

Beta diversity and heatmaps may not generate if only one sample is present.

All scripts are designed to run in the provided folder structure; do not move files around.
