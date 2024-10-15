# Load required libraries
library(phyloseq)
library(tidyverse)
library(ape)
library(pals)
library(randomcoloR)

# Set working directory
setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/phyloseq")

# Read the Kraken2 results table
infile <- "all_samples_chicken_240620_edited_renamed-columns_metadata_kraken2.tsv"
kraken_table <- read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Filter features with low mean abundance
kraken_table <- kraken_table[rowMeans(kraken_table) >= 10, ]

# Read the metadata
metadata <- read.table("metadata_kraken2.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Ensure that 'housing' is a factor with levels 'cage' and 'floor'
metadata$housing <- factor(metadata$housing, levels = c("cage", "floor"))

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

# List of taxonomic levels
tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Function to generate plots for a phyloseq object
generate_plots <- function(physeq_obj, physeq_name) {
  for (tax_level in tax_levels) {
    # Agglomerate data to the desired taxonomic level
    physeq_tax <- tax_glom(physeq_obj, taxrank = tax_level, NArm = TRUE)
    
    # Calculate relative abundances
    physeq_rel <- transform_sample_counts(physeq_tax, function(x) x / sum(x))
    
    # Melt the phyloseq object for ggplot2
    phylo_melt <- psmelt(physeq_rel)
    
    # Replace NA in taxonomic level with "Unknown"
    phylo_melt[[tax_level]][is.na(phylo_melt[[tax_level]])] <- "Unknown"
    
    # Identify the top 50 taxa at this level across all samples
    top_taxa <- phylo_melt %>%
      group_by(!!sym(tax_level)) %>%
      summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
      arrange(desc(Abundance)) %>%
      slice_head(n = 50) %>%  # Changed from 20 to 50
      pull(!!sym(tax_level))
    
    # Assign taxa to 'xOthers' if they are not in the top 50
    phylo_melt[[tax_level]] <- ifelse(phylo_melt[[tax_level]] %in% top_taxa,
                                      phylo_melt[[tax_level]],
                                      "xOthers")
    
    # Recalculate relative abundances after grouping 'xOthers'
    phylo_melt <- phylo_melt %>%
      group_by(Sample, housing, !!sym(tax_level)) %>%
      summarize(Abundance = sum(Abundance), .groups = 'drop')
    
    # Order samples by 'housing' (cage first, then floor)
    sample_order <- sample_data(physeq_obj) %>%
      data.frame() %>%
      arrange(housing) %>%
      rownames()
    
    phylo_melt$Sample <- factor(phylo_melt$Sample, levels = sample_order)
    
    # Generate color palette
    taxa_levels <- unique(phylo_melt[[tax_level]])
    num_taxa <- length(taxa_levels)
    
    # Handle palette generation based on the number of taxa
    if (num_taxa == 1) {
      # Only one taxon (possibly "xOthers")
      if ("xOthers" %in% taxa_levels) {
        palette <- c("xOthers" = "darkgray")
      } else {
        palette <- setNames(distinctColorPalette(1), taxa_levels)
      }
    } else {
      # More than one taxon
      num_colors_needed <- num_taxa - ifelse("xOthers" %in% taxa_levels, 1, 0)
      # Generate a palette with num_colors_needed colors
      palette_colors <- distinctColorPalette(num_colors_needed)
      # Assign gray color to 'xOthers'
      if ("xOthers" %in% taxa_levels) {
        palette <- setNames(c(palette_colors, "gray"), c(setdiff(taxa_levels, "xOthers"), "xOthers"))
      } else {
        palette <- setNames(palette_colors, taxa_levels)
      }
    }
    
    # Plot stacked bar plot for each sample
    p_sample <- ggplot(phylo_melt, aes(x = Sample, y = Abundance, fill = !!sym(tax_level))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = palette) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = paste("Top 50", tax_level, "per Sample in", physeq_name),
           x = "Sample", y = "Relative Abundance", fill = tax_level)
    
    # Save the sample plot
    ggsave(filename = paste0("StackedBar_", physeq_name, "_", tax_level, "_PerSample.png"),
           plot = p_sample, width = 12, height = 8)
    
    # Summarize the data by 'housing'
    phylo_housing <- phylo_melt %>%
      group_by(housing, !!sym(tax_level)) %>%
      summarize(Abundance = mean(Abundance), .groups = 'drop')
    
    # Plot stacked bar plot for each housing type
    p_housing <- ggplot(phylo_housing, aes(x = housing, y = Abundance, fill = !!sym(tax_level))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = palette) +
      theme_bw() +
      labs(title = paste("Top 50", tax_level, "per Housing Type in", physeq_name),
           x = "Housing", y = "Relative Abundance", fill = tax_level)
    
    # Save the housing plot
    ggsave(filename = paste0("StackedBar_", physeq_name, "_", tax_level, "_PerHousing.png"),
           plot = p_housing, width = 8, height = 6)
  }
}



# Generate plots for phyL (litter samples)
generate_plots(phyL, "Litter")

# Generate plots for phyD (dust samples)
generate_plots(phyD, "Dust")