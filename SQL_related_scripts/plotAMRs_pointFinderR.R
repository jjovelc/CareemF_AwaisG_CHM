library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
library(patchwork)

# Set working directory
setwd('/Users/juanjovel/jj/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database')

# Read the TSV file
data <- read.table("pointMutation_results_report.tsv", row.names = 1, header = TRUE, sep = '\t')
plot_outfile <- "pointMutation_combined_plot.png"

# Read metadata
metadata <- read.table("metadata.tsv", row.names = 1, header = TRUE, sep = '\t')
  
# Filter out the specific columns
is_not_Box1_Vial_58_53412 <- colnames(data) != "Box1_Vial_58_53412"
#data_filtered <- data[, is_not_Box1_Vial_58_53412]



# Make metadata and amrs tables IDs the same
metadata$isolate <- gsub("-", "_", metadata$isolate)
data_filtered <- data[, metadata$isolate]
#colnames(data_filtered) <- metadata$isolate

# Check if all isolates in data_filtered have corresponding metadata
isolates_not_in_metadata <- setdiff(colnames(data_filtered), metadata$isolate)
if (length(isolates_not_in_metadata) > 0) {
  stop("The following isolates in data_filtered do not have corresponding metadata: ", paste(isolates_not_in_metadata, collapse = ", "))
}

# Transpose the data
data_t <- t(data_filtered)

# Melt the data for ggplot
melted_data <- melt(as.matrix(data_t))
colnames(melted_data) <- c('Sample', 'AMR', 'Count')

# Check for any discrepancies between Sample and isolate
samples_not_in_metadata <- setdiff(melted_data$Sample, metadata$isolate)
if (length(samples_not_in_metadata) > 0) {
  warning("The following samples in melted_data do not have corresponding metadata: ", paste(samples_not_in_metadata, collapse = ", "))
}

# Merge metadata with melted data using 'isolate'
melted_data <- merge(melted_data, metadata, by.x = "Sample", by.y = "isolate", all.x = TRUE)

# Check the merged data
head(melted_data)

# Sort the samples by 'source'
melted_data <- melted_data %>% arrange(source, Sample)

# Ensure the Sample factor levels are in the sorted order
melted_data$Sample <- factor(melted_data$Sample, levels = unique(melted_data$Sample))

# Create a color palette for the 'source' using RColorBrewer
unique_sources <- unique(melted_data$source)
source_colors <- setNames(brewer.pal(length(unique_sources), "Dark2"), unique_sources)

# Create a color palette for the 'mash_group' using RColorBrewer
unique_mash_groups <- unique(melted_data$mash_group)
mash_group_colors <- setNames(brewer.pal(length(unique_mash_groups), "Set3"), unique_mash_groups)

# Create the histogram for AMR abundance per sample
abundance_data <- melted_data %>% group_by(Sample) %>% summarise(Abundance = sum(Count))
abundance_data$Sample <- factor(abundance_data$Sample, levels = levels(melted_data$Sample))
histogram <- ggplot(abundance_data, aes(x = Sample, y = Abundance)) +
  geom_bar(stat = "identity", fill = "dodgerblue4") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text for the histogram
    axis.ticks.x = element_blank(),  # Remove x-axis ticks for the histogram
    axis.text.y = element_text(size = 24)
  ) + 
  labs(x = NULL, y = "", title = "")

# Create the heatmap with a solid color for counts and add a gray grid
heatmap <- ggplot() +
  geom_tile(data = melted_data, aes(x = Sample, y = AMR, fill = as.factor(Count)), color = "gray") +
  scale_fill_manual(values = c("0" = "gray100", "1" = "dodgerblue1", "2" = "red", "3" = "limegreen", "4" = "black"), name = "Count") +
  new_scale_fill() +
  
  geom_tile(data = unique(melted_data[, c("Sample", "mash_group")]), 
            aes(x = Sample, y = -0.5, fill = mash_group), height = 0.5) +
  scale_fill_manual(values = mash_group_colors, name = "Mash group", guide = guide_legend(order = 2)) +
  new_scale_fill() +
  
  geom_tile(data = unique(melted_data[, c("Sample", "source")]), 
            aes(x = Sample, y = -1.1, fill = source), height = 0.5) +
  scale_fill_manual(values = source_colors, name = "Source", guide = guide_legend(order = 1)) +
  
  
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_text(size = 28),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right",
    axis.title.y = element_text(size = 32),
    legend.text = element_text(size = 32),  # Increase legend text size
    legend.title = element_text(size = 32)  # Increase legend title size
  ) +
  labs(x = NULL, y = "Point mutation in AMR Genes", title = "")

# Combine the histogram and heatmap using patchwork
combined_plot <- histogram / heatmap + plot_layout(heights = c(2, 10))

# Display the combined plot
print(combined_plot)

# Optionally, save the plot
ggsave(plot_outfile, combined_plot, width = 30, height = 20, dpi = 300)