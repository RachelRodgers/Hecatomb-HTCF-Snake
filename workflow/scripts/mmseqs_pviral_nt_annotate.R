#!/usr/bin/env Rscript

# mmseqs_pviral_nt_annotate.R

#----- Load or install required packages -----#

source("./workflow/scripts/snakemake_helpers/snakemake_helpers.R")

options(warn = -1) # suppress warning messages for clarity

requiredPackages <- c("tidyverse", "taxonomizr")

for (package in requiredPackages) {
  TryInstall(package)
  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    stop(cat("FATAL: Problem loading R package:", package, 
             "(mmseqs_pviral_nt_annotate.R)\n\n"),
         call. = FALSE)
  }
}

#----- Read in path to accessionTaxa.sql -----#

ncbiAccPath <- readLines(con = "./results/taxonomizr_ncbi_accession_path.txt", n = 1,
                         warn = FALSE)

#----- Annotate mmseqs_nt_out/resultDB.firsthit.m8 Sequences -----#

# Read *.m8 file with best hit NCBI accession number
m8 <- read_tsv("./results/results/mmseqs_nt_out/resultDB.firsthit.m8", 
               col_names = FALSE)

colnames(m8) <- c("query","target","pident","alnlen","mismatch","gapopen",
                  "qstart","qend","tstart","tend","evalue","bits")

# Create NCBI accession number vector
acc <- m8$target

# Convert accessions to taxonomic IDs
ids <- accessionToTaxa(acc, ncbiAccPath)

# Get taxonomic lineage
ncbi_tax <- as_tibble(getTaxonomy(ids, ncbiAccPath))

# Bind m8 file to lineage
seqids <- select(m8, "query")
mmseqs_pviral_nt_lineage <- cbind(seqids, ncbi_tax)

# Write results to table
dir.create(path = "./results/results/mmseqs_nt_checked_out/", 
           showWarnings=FALSE)

write_tsv(mmseqs_pviral_nt_lineage, 
          path = "./results/results/mmseqs_nt_checked_out/mmseqs_pviral_nt_lineage.tsv")

#----- Save session information -----#

print("mmseqs_pviral_nt_annotate: Saving session info (retain for debugging).\n")

workingDirectory <- getwd()
savePath <- paste(workingDirectory, "/results/R_session_info/", sep = "")
dir.create(path = savePath, showWarnings = FALSE)
saveFile <- file(paste(savePath, 
                       "mmseqs_pviral_nt_annotate_session_info.txt", sep = ""))
writeLines(capture.output(sessionInfo()), saveFile)
close(saveFile)
