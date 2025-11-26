#Step 1: Install and load required R packages
# Install if not already installed
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq"))

# Load packages
library(dada2)
library(phyloseq)
library(ggplot2)

#Step 2: Set project paths (relative paths)

# Path to trimmed FASTQ files
path <- "data/processed"

# Path to save ASV and taxonomy tables
asv_file <- "results/ASV_tables/ASV_table.csv"
tax_file <- "results/ASV_tables/Taxonomy_table.csv"

#Step 3: List FASTQ files

fnFs <- sort(list.files(path, pattern="_trim_1.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_trim_2.fastq.gz", full.names=TRUE))

fnFs  # check forward reads
fnRs  # check reverse reads

#Step 4: Filter and trim reads

# Create output folders for filtered reads
dir.create(file.path(path, "filtF"), showWarnings = FALSE)
dir.create(file.path(path, "filtR"), showWarnings = FALSE)

filtFs <- file.path(path, "filtF", basename(fnFs))
filtRs <- file.path(path, "filtR", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)


#Step 5: Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)



#Step 6: Dereplicate and infer ASVs
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#Step 7: Merge paired reads and remove chimeras
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


#Step 8: Assign taxonomy
# Download SILVA database beforehand (v138)
#silva_db <- "data/SILVA_SSU_r138_train_set.fa.gz"
# SILVA DADA2-formatted databases
silva_genus <- "data/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species <- "data/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

#tax <- assignTaxonomy(seqtab_nochim, silva_genus, multithread=TRUE)

#tax_species <- addSpecies(tax, silva_species)

# Genus-level taxonomy
tax <- assignTaxonomy(seqtab_nochim, "data/silva_nr99_v138.2_toGenus_trainset.fa.gz",
                      multithread = TRUE)

# Species-level assignment
tax_species <- addSpecies(tax, "data/silva_v138.2_assignSpecies.fa.gz",
                          allowMultiple = TRUE)


#Step 9: Save ASV and taxonomy tables
dir.create("results/ASV_tables", showWarnings = FALSE)
write.csv(seqtab_nochim, asv_file)
write.csv(tax, tax_file)


# Step 10: Optional - Create Phyloseq object
# -------------------------------
sample_names <- gsub("_trim_1.fastq.gz","",basename(fnFs))
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
               tax_table(tax_species))

sample_names(ps) <- sample_names

# -------------------------------
# Step 11: Quick visualization (example)
# -------------------------------
# Bar plot of top 10 genera
top_taxa <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
ps_top <- prune_taxa(top_taxa, ps)
plot_bar(ps_top, fill="Genus") + theme(axis.text.x = element_text(angle=90,hjust=1))
ggsave("results/plots/barplot_top10genera.png", width=8, height=6)
# -------------------------------
# End of DADA2 pipeline
# ===============================

# Step 12: Alpha Diversity Metrics

plot_richness(ps, measures = c("Shannon", "Simpson", "Observed")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# this generates three alpha diversity plots
ggsave("results/plots/alpha_diversity.png", width=8, height=6)



# STEP 13 — Taxonomic Composition (beautiful barplots)
# Barplot at Phylum level:

ps_phylum <- tax_glom(ps, taxrank = "Phylum")

plot_bar(ps_phylum, fill="Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("results/plots/phylum_level_barplot.png", width=10, height=6)

#Barplot at Genus level:
ps_genus <- tax_glom(ps, taxrank = "Genus")

plot_bar(ps_genus, fill="Genus") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggsave("results/plots/genus_barplot.png", width=10, height=6)


