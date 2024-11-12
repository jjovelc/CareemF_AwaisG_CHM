library(ggplot2)
library(tidyverse)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/figures')
data <- read.table("data_for_dotplot.tsv", sep = '\t', row.names = 1, header = TRUE)
metadata <- read.table("metadata.tsv", sep = '\t', row.names = 1, header = TRUE)

# Filter and preprocess data
row_means <- rowMeans(data)
data_10up <- data[row_means >= 10, ]
col_sums <- colSums(data_10up)
data_10up_cpm <- (data_10up / col_sums) * 1000000  
data_10up_cpm <- data_10up_cpm + 1  # Add pseudocount

# Function to calculate mean and SD statistics by housing type
calculate_housing_stats <- function(data, cage_cols, floor_cols) {
  results <- data.frame(
    species = rownames(data),
    cage_mean = rowMeans(data[, cage_cols], na.rm = TRUE),
    cage_sd = apply(data[, cage_cols], 1, sd, na.rm = TRUE),
    floor_mean = rowMeans(data[, floor_cols], na.rm = TRUE),
    floor_sd = apply(data[, floor_cols], 1, sd, na.rm = TRUE)
  )
  return(results)
}

# Function to process and plot each group (litter or dust)
process_group <- function(data, metadata, prefix, bar_color_cage, bar_color_floor) {
  # Subset metadata and data for the specific group (litter or dust)
  meta_group <- metadata[metadata$group == prefix, ]
  data_cpm_group <- data[, rownames(meta_group)]
  
  # Define cage and floor samples based on metadata
  cage_samples <- rownames(meta_group)[meta_group$housing == "cage"]
  floor_samples <- rownames(meta_group)[meta_group$housing == "floor"]
  
  # Calculate statistics
  stats_results <- calculate_housing_stats(data_cpm_group, cage_samples, floor_samples)
  
  # Reshape data from wide to long format and add pseudo-count for plotting
  plot_data <- stats_results %>%
    tidyr::pivot_longer(
      cols = c(cage_mean, floor_mean),
      names_to = "housing",
      values_to = "abundance"
    ) %>%
    mutate(
      abundance = abundance + 1,  # Add pseudo-count
      sd = ifelse(housing == "cage_mean", cage_sd, floor_sd),
      housing = ifelse(housing == "cage_mean", "Cage", "Floor")
    )
  
  # Plot the data
  plot <- ggplot(plot_data, aes(x = abundance, y = species, fill = housing)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      alpha = 0.8,
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        xmin = pmax(abundance - sd, 1),
        xmax = abundance + sd
      ),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    scale_x_log10(
      labels = scales::label_number(accuracy = 0.1),
      breaks = scales::trans_breaks("log10", function(x) 10^x)
    ) +
    scale_fill_manual(
      values = c("Cage" = bar_color_cage, "Floor" = bar_color_floor),
      name = "Housing Type"
    ) +
    labs(
      x = "Abundance + 1 (log scale)",
      y = "Species",
      title = paste("Species Abundance by Housing Type in", prefix, "Group")
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 14),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # Save the plot
  filename <- paste0(prefix, "_taxa_horizontal_barplot_cage_vs_floor.png")
  ggsave(filename = filename, plot = plot, width = 12, height = 15, units = "in", dpi = 300)
}


# Run the function for both litter and dust groups
# change colors if you wish
process_group(data_10up_cpm, metadata, "litter" , "orange", "darkgreen")
process_group(data_10up_cpm, metadata, "dust", "red", "dodgerblue1")
