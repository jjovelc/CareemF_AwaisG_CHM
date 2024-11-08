library(RColorBrewer)
library(plotly)
library(reshape2)
library(dplyr)
library(tibble)

# Set working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/git_repos/EcoliGenomics/SQLite_database/scripts')

# Read the TSV file
data <- read.table("mge_results_report.tsv", row.names = 1, header = TRUE, sep = '	')
plot_outfile <- "mge_combined_plot.html"

# Read metadata
metadata <- read.table("metadata.tsv", header = TRUE, sep = '	')

# Ensure 'sample_id' column is in the correct format
metadata$sample_id <- gsub('-', '.', metadata$sample_id)

# Transpose the data
data_t <- t(data)

# Melt the data for plotly
melted_data <- melt(as.matrix(data_t))
colnames(melted_data) <- c('Sample', 'AMR', 'Count')

# Merge metadata with melted data using 'sample_id'
melted_data <- merge(melted_data, metadata, by.x = "Sample", by.y = "sample_id", all.x = TRUE)

# Check the merged data
head(melted_data, n=100)

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
histogram <- plot_ly(abundance_data, x = ~Sample, y = ~Abundance, type = 'bar', marker = list(color = 'dodgerblue4')) %>%
  layout(
    xaxis = list(title = ''),
    yaxis = list(title = 'Abundance'),
    barmode = 'stack',
    margin = list(b = 150),
    title = 'AMR Abundance per Sample'
  )

# Create the heatmap for AMR counts
heatmap <- plot_ly(
  data = melted_data, 
  x = ~Sample, 
  y = ~AMR, 
  z = ~Count, 
  type = 'heatmap', 
  colors = c("white", RColorBrewer::brewer.pal(8, "Set1")),
  showscale = TRUE
) %>%
  layout(
    xaxis = list(title = 'Sample', tickangle = 45),
    yaxis = list(title = 'Mobile Genetic Elements'),
    margin = list(b = 150),
    title = 'AMR Heatmap'
  )

# Create bottom bars for 'mash_group' and 'source'
mash_group_bar <- plot_ly(
  data = melted_data, 
  x = ~Sample, 
  y = ~mash_group, 
  type = 'heatmap', 
  colors = mash_group_colors,
  showscale = FALSE
) %>%
  layout(
    xaxis = list(title = ''),
    yaxis = list(title = 'Mash Group', showticklabels = FALSE),
    margin = list(b = 150),
    height = 100
  )

source_bar <- plot_ly(
  data = melted_data, 
  x = ~Sample, 
  y = ~source, 
  type = 'heatmap', 
  colors = source_colors,
  showscale = FALSE
) %>%
  layout(
    xaxis = list(title = 'Sample', tickangle = 45),
    yaxis = list(title = 'Source', showticklabels = FALSE),
    margin = list(b = 150),
    height = 100
  )

# Combine the plots into one using subplot from plotly
combined_plot <- subplot(histogram, heatmap, mash_group_bar, source_bar, nrows = 4, shareX = TRUE, titleY = TRUE)

# Display the combined plot
combined_plot

# Optionally, save the plot
htmlwidgets::saveWidget(combined_plot, file = plot_outfile)
