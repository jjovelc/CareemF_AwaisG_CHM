# Load necessary libraries
library(ggplot2)
library(reshape2)
library(cowplot)  # For combining plots

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/figures')
data <- read.table("data_for_dotplot.tsv", sep = '\t', row.names = 1, header = TRUE)
metadata <- read.table("metadata.tsv", sep = '\t', row.names = 1, header = TRUE)

# Filter and preprocess data
row_means <- rowMeans(data)
data_10up <- data[row_means >= 10, ]
col_sums <- colSums(data_10up)
data_10up_cpm <- (data_10up / col_sums) * 1000000  

# Add pseudocount to avoid zeros before log transformation
data_10up_cpm <- data_10up_cpm + 1

process_group <- function (data, metadata, prefix, color1, color2){

  # Split data by group
  meta <- paste0("meta_", prefix)
  meta <- metadata[metadata$group == prefix, ]
  cpm  <- data_10up_cpm[, rownames(meta)]
  data <- log(cpm)
  data$Taxa <- rownames(data)
  
  # Reshape data to long format
  long_data <- melt(data, id.vars = "Taxa", variable.name = "Sample", value.name = "Abundance")
  
  # Merge 'housing' into long_data_litter
  long_data <- merge(long_data, metadata[, c("Sample", "housing")], by = "Sample", all.x = TRUE)
  
  # Extract the original sample order from data_litter_cpm
  sample_order <- colnames(cpm)
  
  # Ensure the sample order is consistent across both plots
  long_data$Sample <- factor(long_data$Sample, levels = sample_order)
  meta$Sample <- factor(meta$Sample, levels = sample_order)
  
  # Perform Wilcoxon tests for each Taxa on log-transformed data
  housing_info <- meta$housing
  
  # Remove the Taxa column before adding the pseudocount
  data_numeric <- data[, sample_order]
  
  # Add a small pseudocount to avoid zeros before performing the Wilcoxon test
  pseudocount <- 1e-5
  data_adj <- data_numeric + pseudocount
  
  # Perform Wilcoxon tests for each Taxa on log-transformed data
  p_values <- apply(data_adj, 1, function(x) {
    group1 <- x[housing_info == "floor"]
    group2 <- x[housing_info == "cage"]
    if (length(group1) > 0 & length(group2) > 0) {
      if (all(group1 == group2)) {
        return(1)  # No difference between groups
      } else {
        test_result <- wilcox.test(group1, group2, exact = FALSE)
        return(test_result$p.value)
      }
    } else {
      return(NA)
    }
  })
  
  
  
  # Adjust p-values for multiple testing
  p_values_adj <- p.adjust(p_values, method = "BH")
  
  # Create data frame with adjusted p-values
  taxa_pvalues <- data.frame(Taxa = names(p_values_adj), adj_pvalue = p_values_adj)
  
  # Merge adjusted p-values into long_data_litter
  long_data_litter <- merge(long_data_litter, taxa_pvalues, by = "Taxa", all.x = TRUE)
  
  # Calculate -log10 of adjusted p-values
  long_data_litter$neg_log10_adj_pvalue <- -log10(long_data_litter$adj_pvalue)
  
  # Handle infinite values (if any p-values are zero)
  max_finite <- max(long_data_litter$neg_log10_adj_pvalue[is.finite(long_data_litter$neg_log10_adj_pvalue)], na.rm = TRUE)
  long_data_litter$neg_log10_adj_pvalue[is.infinite(long_data_litter$neg_log10_adj_pvalue)] <- max_finite + 1
  
  # Calculate the midpoint for the color scale
  midpoint_value <- median(long_data_litter$neg_log10_adj_pvalue, na.rm = TRUE)
  
  # Adjust the factor levels of 'Taxa' (without 'Housing')
  long_data_litter$Taxa <- factor(long_data_litter$Taxa)
  
  # Create a data frame for the housing tile bar
  top_tiles <- unique(long_data_litter[, c("Sample", "housing")])
  
  # Ensure top_tiles Sample factor levels are set
  top_tiles$Sample <- factor(top_tiles$Sample, levels = sample_order)
  
  # Create the housing tile bar plot
  
  # Set font size for legends
  legend_text_size <- 12  # Adjust this value as needed
  
  # Modify the housing_plot with increased font size for legend
  housing_plot <- ggplot(top_tiles, aes(x = Sample, y = 1, fill = housing)) +
    geom_tile() + 
    scale_fill_manual(values = c("cage" = "salmon", "floor" = "turquoise", "NA" = "gray"), name = "Housing") +
    theme_bw() +
    labs(fill = "Housing") +
    theme(legend.position = "top",
          legend.text = element_text(size = legend_text_size),
          plot.margin = unit(c(0.0, 0, -0.3, 0), "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  # Modify the main_plot with increased font size for legend
  main_plot <- ggplot(long_data_litter, aes(x = Sample, y = Taxa)) +
    geom_point(aes(size = Abundance, color = neg_log10_adj_pvalue)) +
    scale_size_continuous(range = c(1, 10)) +
    scale_color_gradient2(low = "black", mid = "green", high = "red", midpoint = midpoint_value, na.value = "gray", name = "-log10(adj p-value)") +
    theme_bw() +
    labs(title = "", 
         x = NULL, y = "Taxa", size = "Abundance", color = "-log10(adj p-value)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = y_axis_text_size),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  # Combine the plots with adjusted relative height for the tile bar
  combined_plot_litter <- plot_grid(housing_plot, main_plot, ncol = 1, align = "v", rel_heights = c(0.06, 1))
  
  # Save the combined plot
  ggsave(filename = "litter_selected_taxa_dotplot.png",
         plot = combined_plot_litter,
         width = 10, 
         height = 15,
         units = "in",
         dpi = 300)
}
