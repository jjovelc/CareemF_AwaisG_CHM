library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
library(patchwork)
library(ggnewscale)

# Set working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database')

# Read the AMR metadata TSV file we generated
amr_metadata <- read.table("amrc_resFinder_metadata.tsv", row.names = 1, header = TRUE, sep = '\t')

# Read the original metadata file
metadata <- read.table("metadata.tsv", row.names = 1, header = TRUE, sep = '\t')

plot_outfile <- "amrc_resFinder_metadata_plot.png"

# Replace hyphens with underscores in metadata isolate names
metadata$isolate <- gsub("-", ".", metadata$isolate)

# Ensure metadata columns match metadata isolates
amr_metadata_filtered <- amr_metadata[, metadata$isolate]

# Melt the AMR metadata
melted_data <- melt(as.matrix(amr_metadata_filtered))
colnames(melted_data) <- c('AMR', 'Sample', 'Source')

# Merge with original metadata
melted_data <- merge(melted_data, metadata, by.x = "Sample", by.y = "isolate", all.x = TRUE)

# Sort the samples by 'source'
melted_data <- melted_data %>% arrange(source, Sample)

# Ensure the Sample factor levels are in the sorted order
melted_data$Sample <- factor(melted_data$Sample, levels = unique(melted_data$Sample))

# Create color palettes
source_colors <- c(
  "AMR" = "dodgerblue1",
  "resFinder" = "limegreen",
  "AMR|resFinder" = "red",
  "none" = "white"
)

unique_mash_groups <- unique(melted_data$mash_group)
mash_group_colors <- setNames(brewer.pal(length(unique_mash_groups), "Set3"), unique_mash_groups)

unique_sources <- unique(melted_data$source)
sample_source_colors <- setNames(brewer.pal(length(unique_sources), "Dark2"), unique_sources)

# Create the heatmap
heatmap <- ggplot() +
  geom_tile(data = melted_data, aes(x = Sample, y = AMR, fill = Source), color = "gray") +
  scale_fill_manual(values = source_colors, name = "AMR Source") +
  new_scale_fill() +
  
  geom_tile(data = unique(melted_data[, c("Sample", "mash_group")]), 
            aes(x = Sample, y = -1, fill = mash_group), height = 3) +
  scale_fill_manual(values = mash_group_colors, name = "Mash group", guide = guide_legend(order = 2)) +
  new_scale_fill() +
  
  geom_tile(data = unique(melted_data[, c("Sample", "source")]), 
            aes(x = Sample, y = -4.5, fill = source), height = 3) +
  scale_fill_manual(values = sample_source_colors, name = "Sample Source", guide = guide_legend(order = 1)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right",
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) + 
  labs(x = "Samples", y = "Resistance genes", title = "AMR Gene Sources Across Samples") +
  theme(axis.title.x = element_text(size = 20),  # Change x-axis title font size
        axis.title.y = element_text(size = 20)) +
        theme(plot.title = element_text(size = 24)) 

# Display the plot
print(heatmap)

# Save the plot
ggsave(plot_outfile, heatmap, width = 30, height = 20, dpi = 300)
