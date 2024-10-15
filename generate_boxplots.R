# Load necessary libraries
library(dplyr)
library(ggplot2)
library(wesanderson)
library(stats)

# Set up file paths and import data
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/awais_ghaffar/phyloseq')
file <- file.path(getwd(), 'all_samples_shortNames.tsv')
metaFile <- file.path(getwd(), 'metadata_kraken2.txt')
counts_data <- read.table(file, header = TRUE, sep = '\t', row.names = 1)
metadata <- read.table(metaFile, header = TRUE, sep = '\t', row.names = 1)

# Filter rows where the average abundance across all columns is >= 10
counts_data_filtered <- counts_data[rowMeans(counts_data) >= 10, ]

#write.table(counts_data_filtered, "all_samples_chicken_241014_10up.tsv", sep = '\t', quote = F)

# Calculate CPM
counts_data_cpm <- t(t(counts_data_filtered) / lib_sizes) * 1e6

# Extract transcript IDs
transc_ids <- rownames(counts_data_filtered)

# Initialize counter and result storage
counter <- 0
result_list <- data.frame(File = character(), Status = character(), stringsAsFactors = FALSE)

# Define the plotting function for different sample types (litter/dust)
plot_group <- function(group_type, transc_ids, counts_data_df, metadata, color1, color2) {
  for (i in 1:length(transc_ids)) {
    counter <<- counter + 1
    trID <- transc_ids[i]
    sel_row <- counts_data_df[trID, ]
    row_names <- colnames(counts_data_df)
    
    # Filter metadata for the current group type (litter/dust)
    metadata_subset <- metadata %>% filter(group == group_type)
    
    # Match the condition and CPM counts for the current group type
    group <- metadata_subset$housing
    counts <- as.numeric(sel_row[colnames(counts_data_df) %in% rownames(metadata_subset)])
    
    df <- data.frame(row.names = row_names[row_names %in% rownames(metadata_subset)],
                     condition = group,
                     counts = counts)
    
    # Use trID (transcript ID) for the gene name and file name
    gene.name <- paste(trID, "Expression (CPM)", sep = " ")
    fileName <- paste0(group_type, "_", gsub("/", "_", trID), ".pdf")
    
    message <- paste("Plot#", counter, "Current plot:", fileName, sep = ' ')
    print(message)
    
    # Perform U test (Wilcoxon test)
    p_val <- wilcox.test(counts ~ condition, data = df)$p.value
    p_adj <- p.adjust(p_val, method = "BH")  # Apply Benjamini-Hochberg (BH) correction
    
    # Determine significance
    status <- ifelse(p_adj < 0.05, "significant", "non-significant")
    
    # Store the result in the result_list data frame
    result_list <<- rbind(result_list, data.frame(File = fileName, Status = status))
    
    # Generate the boxplot
    pdf(file = fileName)
    p <- ggplot(df, aes(factor(condition), counts)) +
      geom_boxplot(notch = FALSE, aes(fill = factor(condition))) +
      theme_bw() +
      scale_fill_manual(values = c(color1, color2)) +
      ggtitle(gene.name) +
      xlab("Housing (Cage/Floor)") +
      ylab("Abundance (CPM)") +
      theme(axis.text.x = element_text(colour = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(colour = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            axis.title.x = element_text(colour = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
            axis.title.y = element_text(colour = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold")) +
      labs(subtitle = paste("Adjusted p-value (BH):", signif(p_adj, digits = 3)))  # Display adjusted p-value
    
    print(p)
    dev.off()
  }
}



# Generate plots for 'litter' group
plot_group("litter", transc_ids, counts_data_cpm, metadata, "#3160f7", "#ebbe1c")

# Generate plots for 'dust' group
plot_group("dust", transc_ids, counts_data_cpm, metadata, "#cc1d06", "#04b5a6")

# Save the result list with file names and significance status
write.table(result_list, file = "significance_status_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)

print("Significance status list saved as significance_status_list.txt")