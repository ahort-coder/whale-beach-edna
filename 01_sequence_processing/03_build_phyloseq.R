#to make Phyloseq files

# -----------------------------
# Load libraries
# -----------------------------
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)      # for replace_na()
library(tibble)     # for column_to_rownames()
library(phyloseq)
library(writexl)    # optional if you want to export Excel

# -----------------------------
# Paths
# -----------------------------
input_dir <- "~/Desktop/WB2509"
output_excel <- "~/Desktop/WB2509_all_samples_ASV_table.xlsx"

# -----------------------------
# List TSV files
# -----------------------------
files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# -----------------------------
# Function to read one TSV and format
# -----------------------------
read_asv_file <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Extract sample name from filename: WB2501_10_11 → 10-11
  depth <- str_extract(basename(file), "(?<=WB2509_)\\d+_\\d+")
  sample_name <- str_replace(depth, "_", "-")
  
  # Identify the count column (the one starting with WB2501_)
  count_col <- grep("^WB2509_", colnames(df), value = TRUE)
  
  # Ensure counts are numeric
  df <- df %>%
    mutate(across(all_of(count_col), ~as.numeric(trimws(.x))))
  
  # Keep sequence + counts + taxonomy + ASV info
  df %>%
    select(sequence, all_of(count_col),
           Kingdom:Species, Kingdom_boot:Genus_boot,
           asv_number, asv_code) %>%
    rename(!!sample_name := all_of(count_col))
}

# -----------------------------
# Read all TSV files into list
# -----------------------------
asv_list <- map(files, read_asv_file)

# -----------------------------
# Merge all samples by sequence + taxonomy
# -----------------------------
merged_asv <- reduce(asv_list, full_join,
                     by = c("sequence",
                            "Kingdom","Phylum","Class","Order","Family","Genus","Species",
                            "Kingdom_boot","Phylum_boot","Class_boot","Order_boot",
                            "Family_boot","Genus_boot",
                            "asv_number","asv_code"))

# -----------------------------
# Identify sample count columns
# -----------------------------
sample_cols <- grep("^\\d+-\\d+$", colnames(merged_asv), value = TRUE)

# -----------------------------
# Replace NA counts with 0
# -----------------------------
merged_asv <- merged_asv %>%
  mutate(across(all_of(sample_cols), ~replace_na(.x, 0)))

# -----------------------------
# Collapse duplicate sequences by summing counts
# -----------------------------
asv_df_unique <- merged_asv %>%
  group_by(sequence) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
  # Keep first taxonomy entry for each sequence
  left_join(
    merged_asv %>%
      select(sequence, Kingdom:Species, Kingdom_boot:Genus_boot, asv_number, asv_code) %>%
      distinct(sequence, .keep_all = TRUE),
    by = "sequence"
  )

# -----------------------------
# Create OTU table
# -----------------------------
otu_mat <- asv_df_unique %>%
  select(all_of(sample_cols)) %>%
  as.matrix()
rownames(otu_mat) <- asv_df_unique$sequence
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

# -----------------------------
# Create Taxonomy table
# -----------------------------
tax_mat <- asv_df_unique %>%
  select(sequence, Kingdom:Species) %>%
  column_to_rownames("sequence") %>%
  as.matrix()
tax <- tax_table(tax_mat)

# -----------------------------
# Create Sample metadata
# -----------------------------
sample_df <- data.frame(sample_id = sample_cols) %>%
  mutate(
    depth_start = as.numeric(sub("-.*", "", sample_id)),
    depth_end   = as.numeric(sub(".*-", "", sample_id))
  ) %>%
  column_to_rownames("sample_id")
samp <- sample_data(sample_df)

# -----------------------------
# Build phyloseq object
# -----------------------------
ps <- phyloseq(otu, tax, samp)

# -----------------------------
# Optional: save ASV table to Excel
# -----------------------------
write_xlsx(asv_df_unique, output_excel)

# -----------------------------
# Optional: save phyloseq object
# -----------------------------
saveRDS(ps, "~/Desktop/WB2509_phyloseq_collapsed.rds")

# -----------------------------
# Check counts and number of taxa
# -----------------------------
sample_sums(ps)   # numeric counts
ntaxa(ps)         # should match unique sequences