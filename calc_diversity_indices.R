# Load required libraries
library(phyloseq)
library(tidyverse)
library(ape)
library(pals)

setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/taxa/phyloseq")

infile <- "all_samples_chicken_240620_edited_renamed-columns_metadata_kraken2.tsv"

# Read the Kraken2 results table
kraken_table <- read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
kraken_table <- kraken_table[rowMeans(kraken_table) >= 10, ]

# Read the metadata
metadata <- read.table("metadata_kraken2.txt", header = TRUE, sep = "\t", row.names = 1)

# Create taxonomy table
tax_mat <- tibble(
  taxonomy = rownames(kraken_table), 
  OTU = paste0("OTU", seq_len(nrow(kraken_table)))
) %>%
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(Domain:Species, ~gsub("^.*__", "", .))) %>%
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

# Create phylogenetic tree (optional, using a simple method)
tree <- ape::rtree(ntaxa(OTU), rooted = TRUE, tip.label = taxa_names(OTU))

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, sampledata, tree)

# Print summary of the phyloseq object
print(physeq)

# Subset the phyloseq object based on the 'group' variable
physeq_litter <- subset_samples(physeq, group == "litter")
physeq_dust <- subset_samples(physeq, group == "dust")

# Print summaries of the new phyloseq objects
print("Summary of litter samples:")
print(physeq_litter)

print("Summary of dust samples:")
print(physeq_dust)

# Prune taxa that are not present in the subsets
phyL <- prune_taxa(taxa_sums(physeq_litter) > 0, physeq_litter)
phyD <- prune_taxa(taxa_sums(physeq_dust) > 0, physeq_dust)


########## ALPHA-DIVERSITY ANALYSIS

# Define a function to calculate U test (Wilcoxon rank-sum test) for each alpha diversity measure
perform_u_test <- function(physeq, measure, grouping_var) {
  # Extract alpha diversity data
  alpha_data <- estimate_richness(physeq, measures = measure)
  
  # Merge alpha diversity data with sample data
  sample_data_df <- data.frame(sample_data(physeq))
  alpha_data <- cbind(sample_data_df, alpha_data)
  
  # Perform the U test (Wilcoxon rank-sum test) between groups for each measure
  test_result <- wilcox.test(alpha_data[[measure]] ~ alpha_data[[grouping_var]])
  
  return(test_result$p.value)
}

# Define alpha diversity measures
alpha_meas <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")

# Perform U test for each measure and correct the p-values for litter and dust
p_values_litter <- sapply(alpha_meas, function(measure) perform_u_test(phyL, measure, "housing"))
p_values_dust <- sapply(alpha_meas, function(measure) perform_u_test(phyD, measure, "housing"))

# Apply multiple testing correction (e.g., Benjamini-Hochberg)
corrected_p_values_litter <- p.adjust(p_values_litter, method = "BH")
corrected_p_values_dust <- p.adjust(p_values_dust, method = "BH")

# Function to add p-values to the title of each plot
add_pvalue_to_title <- function(p, corrected_p_values) {
  # Iterate over the alpha diversity measures and update each plot title
  p <- p + facet_wrap(~variable, scales = "free", ncol = 6, labeller = labeller(
    .cols = function(x) paste(x, "\np.adj =", signif(corrected_p_values[alpha_meas == x], 3))))
  return(p)
}

# Plot richness for litter samples and add p-values to title
p_litter <- plot_richness(phyL, "housing", "age", measures = alpha_meas) +
  geom_boxplot(aes(x = housing, y = value, color = NULL), alpha = 0.1) +
  theme_bw()

p_litter <- add_pvalue_to_title(p_litter, corrected_p_values_litter)

# Adjust title size
p_litter <- p_litter + theme(strip.text = element_text(size = 10)) # Adjust font size

# Print and save the litter plot
print(p_litter)
ggsave("litter_alpha_diversity_with_pvalues.png", plot = p_litter, width = 10, height = 6)

# Plot richness for dust samples and add p-values to title
p_dust <- plot_richness(phyD, "housing", "age", measures = alpha_meas) +
  geom_boxplot(aes(x = housing, y = value, color = NULL), alpha = 0.1) +
  theme_bw()

p_dust <- add_pvalue_to_title(p_dust, corrected_p_values_dust)

# Adjust title size
p_dust <- p_dust + theme(strip.text = element_text(size = 10)) # Adjust font size

# Print and save the dust plot
print(p_dust)
ggsave("dust_alpha_diversity_with_pvalues.png", plot = p_dust, width = 10, height = 6)


##########
### BETA-DIVERSITY ANALYSIS
# Make Principal coordinate analysis for litter
# Bray-Curtis distance for litter
library(vegan)

# Bray-Curtis distance for litter
bray_dist_litter <- phyloseq::distance(phyL, method = "bray")

# Run the PERMANOVA test
permanova_litter <- adonis2(bray_dist_litter ~ housing, data = as(sample_data(phyL), "data.frame"), permutations = 999)
pvalue_litter <- permanova_litter$`Pr(>F)`[1]  # Extract the p-value for the 'housing' term

print("PERMANOVA for Litter Samples")
print(permanova_litter)

# Plot PCoA for litter and add p-value in the subtitle
pcoa_res_litter <- ordinate(phyL, method = "PCoA", distance = bray_dist_litter)
pcoa_plot_litter <- plot_ordination(phyL, pcoa_res_litter, color = "housing") +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCoA of Phylogenetic Composition - Litter Samples",
       subtitle = paste("PERMANOVA p-value: ", format(pvalue_litter, digits = 3)),
       x = paste("PC1 - ", round(pcoa_res_litter$values$Relative_eig[1] * 100, 2), "%", sep = ""),
       y = paste("PC2 - ", round(pcoa_res_litter$values$Relative_eig[2] * 100, 2), "%", sep = ""))

print(pcoa_plot_litter)
ggsave("pcoa_plot_litter.png", pcoa_plot_litter, width = 8, height = 6)

# Bray-Curtis distance for dust
bray_dist_dust <- phyloseq::distance(phyD, method = "bray")

# Run the PERMANOVA test
permanova_dust <- adonis2(bray_dist_dust ~ housing, data = as(sample_data(phyD), "data.frame"), permutations = 999)
pvalue_dust <- permanova_dust$`Pr(>F)`[1]  # Extract the p-value for the 'housing' term

print("PERMANOVA for Dust Samples")
print(permanova_dust)

# Plot PCoA for dust and add p-value in the subtitle
pcoa_res_dust <- ordinate(phyD, method = "PCoA", distance = bray_dist_dust)
pcoa_plot_dust <- plot_ordination(phyD, pcoa_res_dust, color = "housing") +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCoA of Phylogenetic Composition - Dust Samples",
       subtitle = paste("PERMANOVA p-value: ", format(pvalue_dust, digits = 3)),
       x = paste("PC1 - ", round(pcoa_res_dust$values$Relative_eig[1] * 100, 2), "%", sep = ""),
       y = paste("PC2 - ", round(pcoa_res_dust$values$Relative_eig[2] * 100, 2), "%", sep = ""))

print(pcoa_plot_dust)
ggsave("pcoa_plot_dust.png", pcoa_plot_dust, width = 8, height = 6)

# Save the phyloseq object
saveRDS(physeq, "phyloseq_object.rds")

# Save the subsetted phyloseq objects
saveRDS(physeq_litter, "physeq_litter.rds")
saveRDS(physeq_dust, "physeq_dust.rds")
