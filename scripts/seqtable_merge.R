#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read seqtables

print("Reading seq tables")

files <- list.files(path = "./QC/clustered/", pattern = "*_seqtable.txt", full.names = TRUE)

# Reduce seqtables to a single table

print("Reducing seqtables into single table")

seqtable.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "sequence") %>%
  mutate(id = row_number() - 1) %>%
  select(id, everything()) %>%
  mutate_if(is.numeric, as.integer)

# Write count table

print("Writing count table")

dir.create(path = "results", showWarnings = FALSE)
write_tsv(seqtable.all, path = "./results/seqtable.all", col_names=TRUE)

# Write	tabular	fasta (tab2fx)

print("Writing tabular fasta")

seqs <- tibble(`sequence` = seqtable.all$sequence)
seqs.df <- seqs %>%
	mutate(id = row_number() - 1) %>%
	select(id, everything()) %>%
	mutate_if(is.numeric, as.integer)
write_tsv(seqs.df, "./results/seqtable.tab2fx", col_names = FALSE)

# Save session information
workingDirectory <- getwd()
savePath <- paste(workingDirectory, "/R_session_info/", sep = "")
dir.create(path = savePath, showWarnings = FALSE)
saveFile <- file(paste(savePath, "seqtable_merge_R_session_info.txt", sep = ""))
writeLines(capture.output(sessionInfo()), saveFile)
close(saveFile)
