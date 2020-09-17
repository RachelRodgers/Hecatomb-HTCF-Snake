#!/usr/bin/env Rscript

#----- Load or install required packages -----#
source("./scripts/hecatomb_helpers.R")

requiredPackages <- "tidyverse"

for (package in requiredPackages) {
  TryInstall(package)
  
  if (!require(package, character.only = )) {
    stop(cat("Problem loading R package", package, "\n"),
         call. = FALSE)
    }
  }

#----- Collect seqtables -----#
print("seqtable_merge: Reading sequencing tables.\n")

files <- list.files(path = "./QC/step_8/clustered/", 
                    pattern = "*_seqtable.txt", full.names = TRUE)

#----- Reduce seqtables to a single table -----#
print("seqtable merge: Reducing sample sequencing tables into single table.\n")

seqtable.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "sequence") %>%
  mutate(id = row_number() - 1) %>% # changing to base 0 for downstream compatability
  select(id, everything()) %>%
  mutate_if(is.numeric, as.integer)

#----- Write count table -----#
print("seqtable merge: Writing count table.\n")

dir.create(path = "./results", showWarnings = FALSE)

write_tsv(seqtable.all, path = "./results/seqtable.all", col_names=TRUE)

#----- Write	tabular	fasta (tab2fx) -----#
print("seqtable merge: Writing tabular fasta.\n")

seqs <- tibble(`sequence` = seqtable.all$sequence)

seqs.df <- seqs %>%
	mutate(id = row_number() - 1) %>% # changing to base 0 for downstream compatability
	select(id, everything()) %>%
	mutate_if(is.numeric, as.integer)

write_tsv(seqs.df, "./results/seqtable.tab2fx", col_names = FALSE)

#----- Save session information -----#
print("seqtable merge: Saving session info (retain for debugging).\n")

workingDirectory <- getwd()
savePath <- paste(workingDirectory, "./R_session_info/", sep = "")
dir.create(path = savePath, showWarnings = FALSE)
saveFile <- file(paste(savePath, "seqtable_merge_R_session_info.txt", sep = ""))
writeLines(capture.output(sessionInfo()), saveFile)
close(saveFile)
