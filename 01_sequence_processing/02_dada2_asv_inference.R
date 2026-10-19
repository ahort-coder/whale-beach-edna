############################################################
# LOAD LIBRARIES
############################################################
library(dada2)
library(ShortRead)
library(Biostrings)
library(tibble)
library(dplyr)
library(readr)
library(tidyr)

############################################################
# SET PATHS
############################################################
base_path <- "~/Documents/DADA2_output"
dir.create(base_path, showWarnings = FALSE)

fastq_path <- "~/Desktop"

# SILVA references
silva_file <- "~/Desktop/silva_nr99_v138.1_train_set.fa.gz"
species_file <- "~/Desktop/silva_species_assignment_v138.1.fa.gz"

############################################################
# DEFINE SAMPLES  
############################################################
samples <- list(
  list(
    name = "WB1_25_26",
    R1 = "trimmed_WB1_25_26_S5_L001_R1_001.fastq.gz",
    R2 = "trimmed_WB1_25_26_S5_L001_R2_001.fastq.gz"
  )
)
############################################################
# PROCESS EACH SAMPLE
############################################################
for (s in samples) {
  
  sample_name <- s$name
  fnF <- file.path(fastq_path, s$R1)
  fnR <- file.path(fastq_path, s$R2)
  
  message("Processing sample: ", sample_name)
  
  # FILTERING
  filt_path <- file.path(base_path, "filtered")
  dir.create(filt_path, showWarnings = FALSE)
  
  filtF <- file.path(filt_path, paste0(sample_name, "_R1_filt.fastq.gz"))
  filtR <- file.path(filt_path, paste0(sample_name, "_R2_filt.fastq.gz"))
  
  filterAndTrim(fnF, filtF, fnR, filtR,
                truncLen = c(300,280),
                maxEE = c(2,2),
                truncQ = 2,
                maxN = 0,
                rm.phix = TRUE,
                compress = TRUE,
                multithread = TRUE)
  
  # LEARN ERRORS
  errF <- learnErrors(filtF, multithread = TRUE)
  errR <- learnErrors(filtR, multithread = TRUE)
  
  # DEREPLICATION
  derepF <- derepFastq(filtF)
  derepR <- derepFastq(filtR)
  
  # DADA2
  dadaF <- dada(derepF, err = errF, multithread = TRUE)
  dadaR <- dada(derepR, err = errR, multithread = TRUE)
  
  # MERGE
  mergers <- mergePairs(dadaF, derepF, dadaR, derepR)
  
  # SEQ TABLE
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, file.path(base_path, paste0(sample_name, "_seqtab_raw.rds")))
  
  # REMOVE CHIMERAS
  seqtab.nochim <- removeBimeraDenovo(seqtab,
                                      method = "pooled",
                                      multithread = TRUE)
  saveRDS(seqtab.nochim, file.path(base_path, paste0(sample_name, "_seqtab_nochim.rds")))
  
  # TAXONOMY + BOOTSTRAPS
  taxa <- assignTaxonomy(seqtab.nochim,
                         refFasta = silva_file,
                         outputBootstraps = TRUE,
                         multithread = TRUE)
  
  # addSpecies only modifies the TAX table
  if (file.exists(species_file)) {
    taxa$tax <- addSpecies(taxa$tax, refFasta = species_file)
  }
  
  saveRDS(taxa, file.path(base_path, paste0(sample_name, "_taxonomy_assignments.rds")))
  
  # FINAL ASV TABLE
  seqtab_df <- as.data.frame(t(seqtab.nochim)) %>%
    rownames_to_column("sequence")
  colnames(seqtab_df)[2] <- sample_name
  
  taxa_df <- as.data.frame(taxa$tax) %>%
    rownames_to_column("sequence")
  
  boot_df <- as.data.frame(taxa$boot) %>%
    rownames_to_column("sequence") %>%
    dplyr::rename(
      Kingdom_boot = Kingdom,
      Phylum_boot = Phylum,
      Class_boot = Class,
      Order_boot = Order,
      Family_boot = Family,
      Genus_boot = Genus
    )
  
  tax_full <- left_join(taxa_df, boot_df, by = "sequence")
  
  asv_table <- left_join(seqtab_df, tax_full, by = "sequence") %>%
    mutate(
      asv_number = row_number(),
      asv_code = paste0("ASV_", asv_number)
    ) %>%
    select(
      asv_number,
      sequence,
      !!sample_name,
      asv_code,
      Kingdom, Phylum, Class, Order, Family, Genus,
      Kingdom_boot, Phylum_boot, Class_boot, Order_boot, Family_boot, Genus_boot
    )
  
  write_tsv(asv_table, file.path(base_path, paste0(sample_name, "_DADA2_ASV_full_table.tsv")))
  
  message("Finished processing ", sample_name)
}

message("All samples processed. ASV tables saved in ", base_path)