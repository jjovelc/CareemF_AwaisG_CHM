# Load required libraries
library(phyloseq)
library(tidyverse)
library(ape)
library(Maaslin2)

# Set working directory
setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/phyloseq")

# Read the Kraken2 results table
infile <- "all_samples_chicken_240620_edited_renamed-columns_metadata_kraken2.tsv"
kraken_table <- read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Filter features with low mean abundance
kraken_table <- kraken_table[rowMeans(kraken_table) >= 10, ]

# Read the metadata
metadata <- read.table("metadata_kraken2.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Ensure that 'housing' and 'age' are factors
metadata$housing <- factor(metadata$housing)
metadata$age <- factor(metadata$age)

# Create taxonomy table
tax_mat <- tibble(
  taxonomy = rownames(kraken_table),
  OTU = paste0("OTU", seq_len(nrow(kraken_table)))
) %>%
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~ gsub("^.*__", "", .))) %>%
  column_to_rownames("OTU")

# Replace empty strings with NA
tax_mat[tax_mat == ""] <- NA

TAX <- tax_table(as.matrix(tax_mat))

# Create OTU table with matching row names
otu_mat <- as.matrix(kraken_table)
rownames(otu_mat) <- rownames(tax_mat)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Create sample data
sampledata <- sample_data(metadata)

# Create phylogenetic tree (optional)
tree <- ape::rtree(ntaxa(OTU), rooted = TRUE, tip.label = taxa_names(OTU))

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, sampledata, tree)

# Subset the phyloseq object based on the 'group' variable
physeq_litter <- subset_samples(physeq, group == "litter")
physeq_dust <- subset_samples(physeq, group == "dust")

# Prune taxa that are not present in the subsets
phyL <- prune_taxa(taxa_sums(physeq_litter) > 0, physeq_litter)
phyD <- prune_taxa(taxa_sums(physeq_dust) > 0, physeq_dust)

# Perform MaAsLin 2 analysis separately for phyL and phyD
##########
run_maaslin2 <- function(physeq_obj, analysis_name, tax_level = "Genus") {
  # Load necessary libraries within the function
  library(phyloseq)
  library(dplyr)
  library(Maaslin2)
  
  # 1. Extract and Prepare Counts Data
  counts_data <- as.data.frame(otu_table(physeq_obj))
  
  # Transpose if taxa are rows
  if (taxa_are_rows(physeq_obj)) {
    counts_data <- t(counts_data)
  }
  
  # Ensure counts_data is a data frame
  counts_data <- as.data.frame(counts_data)
  
  # 2. Extract and Process Taxonomy Table
  tax_table_df <- as.data.frame(tax_table(physeq_obj))
  
  # Function to determine the taxa name based on the desired taxonomic level
  get_taxa_name <- function(tax_row, desired_level) {
    levels_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain")
    desired_index <- match(desired_level, levels_order)
    for (i in desired_index:length(levels_order)) {
      level <- levels_order[i]
      if (!is.na(tax_row[[level]]) && tax_row[[level]] != "") {
        return(tax_row[[level]])
      }
    }
    return("Unknown")
  }
  
  # Apply the function to create a TaxaName column
  tax_table_df$TaxaName <- apply(tax_table_df, 1, function(row) get_taxa_name(row, tax_level))
  
  # Ensure taxa names are unique by appending OTU IDs
  tax_table_df$TaxaName <- paste0(tax_table_df$TaxaName, "_", rownames(tax_table_df))
  
  # 3. Map OTU IDs to Taxa Names in Counts Data
  otu_ids <- colnames(counts_data)
  taxa_names <- tax_table_df$TaxaName[match(otu_ids, rownames(tax_table_df))]
  
  # Handle unmatched taxa
  if (any(is.na(taxa_names))) {
    warning("Some OTU IDs could not be matched to taxa names. They will be labeled as 'Unknown_OTUID'.")
    taxa_names[is.na(taxa_names)] <- paste0("Unknown_", otu_ids[is.na(taxa_names)])
  }
  
  # Assign taxa names as column names
  colnames(counts_data) <- taxa_names
  
  # Remove taxa with all zero counts
  counts_data <- counts_data %>% dplyr::select_if(~ !all(. == 0))
  
  # 4. Extract and Clean Metadata
  # Convert sample_data to a pure data frame
  metadata <- as.data.frame(sample_data(physeq_obj))
  
  # Remove any residual classes by coercing to data frame
  metadata <- data.frame(metadata, stringsAsFactors = FALSE)
  
  # Trim whitespace from column names
  colnames(metadata) <- trimws(colnames(metadata))
  
  # Trim whitespace from character columns in metadata
  metadata <- metadata %>%
    mutate(across(where(is.character), ~ trimws(.)))
  
  # Additionally, trim whitespace in factor levels if factors
  if ("age" %in% colnames(metadata)) {
    if (is.factor(metadata$age)) {
      levels(metadata$age) <- trimws(levels(metadata$age))
    }
  }
  
  if ("housing" %in% colnames(metadata)) {
    if (is.factor(metadata$housing)) {
      levels(metadata$housing) <- trimws(levels(metadata$housing))
    }
  }
  
  # 5. Ensure Sample ID Matching
  counts_samples <- rownames(counts_data)
  metadata_samples <- rownames(metadata)
  common_samples <- intersect(counts_samples, metadata_samples)
  
  # Subset counts_data and metadata to common samples
  counts_data <- counts_data[common_samples, , drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]
  
  # 6. Convert Counts to Relative Abundances
  relative_abundance <- sweep(counts_data, 1, rowSums(counts_data), "/")
  
  # Replace NA or NaN values resulting from division by zero
  relative_abundance[is.na(relative_abundance)] <- 0
  
  # Ensure input_data is a data frame with samples as rows and features as columns
  input_data <- as.data.frame(relative_abundance)
  
  # Remove any taxa with all zero counts after intersection
  input_data <- input_data %>% dplyr::select_if(~ !all(. == 0))
  
  # 7. Format Metadata Variables
  # Convert 'housing' and 'age' to factors
  if ("housing" %in% colnames(metadata)) {
    metadata$housing <- factor(metadata$housing)
  } else {
    stop("'housing' variable not found in metadata.")
  }
  
  if ("age" %in% colnames(metadata)) {
    metadata$age <- factor(metadata$age)
  } else {
    stop("'age' variable not found in metadata.")
  }
  
  # 8. Handle Special Characters in Taxa Names
  special_chars <- grepl("[^a-zA-Z0-9_]", colnames(input_data))
  if (any(special_chars)) {
    warning("Some taxa names contain special characters. Replacing them with underscores.")
    colnames(input_data) <- gsub("[^a-zA-Z0-9_]", "_", colnames(input_data))
  }
  
  # 9. Define Fixed Effects and Output Directory
  fixed_effects <- c("housing", "age")
  output_dir <- paste0("Maaslin2_results_", analysis_name)
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    message(paste("Created output directory:", output_dir))
  }
  
  # 10. Run MaAsLin 2
  fit_data <- Maaslin2(
    input_data = input_data,
    input_metadata = metadata,
    output = output_dir,
    fixed_effects = fixed_effects,
    normalization = "NONE",  # Data is already normalized to relative abundances
    transform = "AST",       # Arc-sine square root transformation
    plot_heatmap = TRUE,
    plot_scatter = TRUE
  )
  
  # 11. Print Significant Results
  significant_results <- fit_data$results %>% filter(qval < 0.05)
  
  if (nrow(significant_results) > 0) {
    print(paste("Significant results for", analysis_name))
    print(significant_results)
  } else {
    print(paste("No significant results found for", analysis_name, "at q < 0.05"))
  }
}

##########



# Run MaAsLin 2 for phyL (litter samples)
run_maaslin2(phyL, "Litter")

# Run MaAsLin 2 for phyD (dust samples)
run_maaslin2(phyD, "Dust")
