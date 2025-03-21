# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(forcats)


install.packages("ggpubr")  # If not installed
library(ggpubr)  # Load the package

getwd()


abricate_dir <- "C:/Users/Awais/Desktop/litter data/New folder/newresults"
metadata_file <- "C:/Users/Awais/Desktop/litter data/New folder/newresults/metadata.txt"

# Function to read and process Abricate files
read_abricate_files <- function(sample_file) {
  if(grepl("metadata", sample_file)) return(NULL)
  data <- read.csv(sample_file, sep="\t", header=TRUE)
  sample_name <- str_extract(basename(sample_file), "(litter|dust)\\d+")
  data$Sample <- sample_name
  return(data)
}

# Read and combine all abricate files
abricate_files <- list.files(path = abricate_dir, pattern = "*.txt", full.names = TRUE) %>%
  .[!grepl("metadata", .)]
all_data <- lapply(abricate_files, read_abricate_files) %>% bind_rows()

# Read metadata and join with resistance data
metadata <- read.csv(metadata_file, sep="\t")
resistome_data <- all_data %>% left_join(metadata, by = c("Sample" = "samples"))

# Split data into dust and litter samples
dust_data <- resistome_data %>% filter(grepl("dust", Sample))
litter_data <- resistome_data %>% filter(grepl("litter", Sample))

# Function to get top genes with complete housing coverage
get_top_genes_combined <- function(data, n) {
  samples_per_housing <- data %>% group_by(housing) %>% summarise(total_samples = n_distinct(Sample))
  gene_stats <- data %>%
    group_by(Gene.symbol, Class, housing) %>%
    summarise(Frequency = n(), Samples = n_distinct(Sample), .groups = 'drop') %>%
    left_join(samples_per_housing, by = "housing") %>%
    mutate(Prevalence = (Samples / total_samples) * 100) %>%
    select(-total_samples)
  
  top_genes <- gene_stats %>% 
    group_by(Gene.symbol) %>% 
    summarise(total_freq = sum(Frequency)) %>%
    top_n(n, total_freq) %>% 
    pull(Gene.symbol)
  
  complete_gene_stats <- expand.grid(
    Gene.symbol = top_genes, 
    housing = unique(data$housing),
    stringsAsFactors = FALSE
  ) %>%
    left_join(gene_stats, by = c("Gene.symbol", "housing")) %>%
    mutate(
      Frequency = replace_na(Frequency, 0),
      Samples = replace_na(Samples, 0),
      Prevalence = replace_na(Prevalence, 0)
    ) %>%
    left_join(gene_stats %>% select(Gene.symbol, Class) %>% distinct(), by = "Gene.symbol")
  
  return(complete_gene_stats)
}

# Function to create bubble plots with increased font size for legends and axis labels
create_bubble_plot <- function(data, title) {
  p <- ggplot(data, aes(x = reorder(Gene.symbol, -Frequency), y = Prevalence, 
                        size = Frequency, color = housing)) +
    geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
    scale_size_continuous(range = c(1, 12)) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = title,
      x = "Gene Symbol",
      y = "Prevalence (%)",
      size = "Frequency",
      color = "Housing System"
    )
  return(p)
}

# Function to create heatmaps with 300 DPI
create_heatmap <- function(data, sample_type, top_n = 20) {
  top_genes <- data %>%
    count(Gene.symbol, sort = TRUE) %>%
    slice_max(n, n = top_n) %>%
    pull(Gene.symbol)
  
  heatmap_data <- data %>%
    filter(Gene.symbol %in% top_genes) %>%
    group_by(Sample, Gene.symbol, housing) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene.symbol, values_from = Count, values_fill = 0)
  
  row_annotation <- data.frame(
    Housing = heatmap_data$housing, 
    row.names = heatmap_data$Sample
  )
  
  heatmap_matrix <- heatmap_data %>%
    column_to_rownames("Sample") %>%
    select(-housing) %>%
    as.matrix()
  
  ann_colors <- list(Housing = c(cage = "#66C2A5", floor = "#FC8D62"))
  
  pheatmap(
    heatmap_matrix,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    main = paste("Top", top_n, "ARG Patterns in", sample_type, "Samples"),
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    filename = paste0("heatmap_top", top_n, "_", tolower(sample_type), ".png"),
    width = 12,
    height = 8,
    res = 300
  )
}

# Adjusted class analysis function with larger axis and legend text
create_class_analysis <- function(data, sample_type) {
  class_dist <- data %>%
    group_by(Class, housing) %>%
    summarise(
      count = n(),
      unique_genes = n_distinct(Gene.symbol),
      .groups = 'drop'
    )
  
  p1 <- ggplot(class_dist, aes(x = reorder(Class, -count), y = count, fill = housing)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("cage" = "#66C2A5", "floor" = "#FC8D62")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 08),
      legend.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = paste("AMR Class Distribution in", sample_type, "Samples"),
      x = "Resistance Class",
      y = "Count",
      fill = "Housing System"
    )
  
  class_prop <- class_dist %>%
    group_by(housing) %>%
    mutate(proportion = count / sum(count) * 100)
  
  p2 <- ggplot(class_prop, aes(x = reorder(Class, -proportion), y = proportion, fill = housing)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("cage" = "#66C2A5", "floor" = "#FC8D62")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 08),
      legend.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = paste("AMR Class Proportions in", sample_type, "Samples"),
      x = "Resistance Class",
      y = "Percentage (%)",
      fill = "Housing System"
    )
  
  return(list(distribution = p1, proportion = p2))
}

# Generate all visualizations
for(n in c(20, 50)) {
  # Dust bubble plot
  dust_genes <- get_top_genes_combined(dust_data, n)
  p <- create_bubble_plot(dust_genes, paste("Top", n, "ARGs in Dust Samples"))
  ggsave(
    paste0("bubble_plot_top", n, "_dust.png"),
    p,
    width = 10, # or any appropriate width
    height = 6, # or any appropriate height
    dpi = 300
  )
  
  # Litter bubble plot
  litter_genes <- get_top_genes_combined(litter_data, n)
  p <- create_bubble_plot(litter_genes, paste("Top", n, "ARGs in Litter Samples"))
  ggsave(
    paste0("bubble_plot_top", n, "_litter.png"),
    p,
    width = 10, # or any appropriate width
    height = 6, # or any appropriate height
    dpi = 300
  )
  
  # Heatmaps
  create_heatmap(dust_data, "Dust", n)
  create_heatmap(litter_data, "Litter", n)
}


# Generate AMR class analysis plots for dust samples
dust_class_plots <- create_class_analysis(dust_data, "Dust")
ggsave("class_distribution_dust.png", dust_class_plots$distribution, width = 8, height = 6, dpi = 300)
ggsave("class_proportion_dust.png", dust_class_plots$proportion, width = 8, height = 6, dpi = 300)

# Generate AMR class analysis plots for litter samples
litter_class_plots <- create_class_analysis(litter_data, "Litter")
ggsave("class_distribution_litter.png", litter_class_plots$distribution, width = 10, height = 8, dpi = 300)
ggsave("class_proportion_litter.png", litter_class_plots$proportion, width = 10, height = 8, dpi = 300)




####stat housing system wise
# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(tibble)  # Load the package

# Function to calculate and plot alpha diversity
calculate_alpha_diversity <- function(data, sample_type) {
  # Prepare abundance matrix
  abundance_matrix <- data %>%
    group_by(Sample, Gene.symbol) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene.symbol, values_from = Count, values_fill = 0) %>%
    column_to_rownames("Sample")
  
  # Calculate diversity indices
  shannon <- diversity(abundance_matrix, index = "shannon")
  simpson <- diversity(abundance_matrix, index = "simpson")
  richness <- specnumber(abundance_matrix)
  
  # Combine indices with metadata
  alpha_div <- data.frame(
    Sample = rownames(abundance_matrix),
    Shannon = shannon,
    Simpson = simpson,
    Richness = richness
  ) %>%
    left_join(data %>% select(Sample, housing) %>% distinct(), by = "Sample")
  
  # Statistical tests with Kruskal-Wallis
  shannon_test <- kruskal.test(Shannon ~ housing, data = alpha_div)
  simpson_test <- kruskal.test(Simpson ~ housing, data = alpha_div)
  richness_test <- kruskal.test(Richness ~ housing, data = alpha_div)
  
  # Function to add significance stars
  add_significance_stars <- function(p_value) {
    if(is.na(p_value)) return("ns")
    if(p_value <= 0.001) return("***")
    if(p_value <= 0.01) return("**")
    if(p_value <= 0.05) return("*")
    return("ns")
  }
  
  # Create plots with statistics
  # Modify theme settings within the create_diversity_plot function
  create_diversity_plot <- function(data, y_var, y_lab, title, test_result) {
    p_val <- round(test_result$p.value, 3)
    sig_stars <- add_significance_stars(p_val)
    
    y_max <- max(data[[y_var]], na.rm = TRUE)
    y_pos <- y_max * 1.1
    
    ggplot(data, aes(x = housing, y = !!sym(y_var), fill = housing)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      annotate("text", x = 1.5, y = y_pos,
               label = paste0("p = ", p_val, " ", sig_stars)) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste(title, "-", sample_type),
           x = "Housing System", 
           y = y_lab) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
  }
  
  
  # Generate plots
  shannon_plot <- create_diversity_plot(alpha_div, "Shannon", "Shannon Index", 
                                        "Shannon Diversity", shannon_test)
  simpson_plot <- create_diversity_plot(alpha_div, "Simpson", "Simpson Index", 
                                        "Simpson Diversity", simpson_test)
  richness_plot <- create_diversity_plot(alpha_div, "Richness", "Number of Genes", 
                                         "Species Richness", richness_test)
  
  # Save plots
  ggsave(paste0("alpha_diversity_shannon_", tolower(sample_type), ".png"), 
         shannon_plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0("alpha_diversity_simpson_", tolower(sample_type), ".png"), 
         simpson_plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0("alpha_diversity_richness_", tolower(sample_type), ".png"), 
         richness_plot, width = 8, height = 6, dpi = 300)
  
  # Save statistical results with effect sizes
  stats_results <- data.frame(
    Metric = c("Shannon", "Simpson", "Richness"),
    P_value = c(shannon_test$p.value, simpson_test$p.value, richness_test$p.value),
    Significance = sapply(c(shannon_test$p.value, simpson_test$p.value, richness_test$p.value), 
                          add_significance_stars)
  )
  write.csv(stats_results, 
            paste0("alpha_diversity_stats_", tolower(sample_type), ".csv"), 
            row.names = FALSE)
  
  return(list(
    alpha_div = alpha_div,
    plots = list(shannon = shannon_plot, simpson = simpson_plot, richness = richness_plot),
    stats = stats_results
  ))
}

# Function to calculate and plot beta diversity
calculate_beta_diversity <- function(data, sample_type) {
  tryCatch({
    # Prepare abundance matrix
    abundance_matrix <- data %>%
      group_by(Sample, Gene.symbol) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = Gene.symbol, values_from = Count, values_fill = 0) %>%
      column_to_rownames("Sample")
    
    # Calculate distance matrices
    bray_dist <- vegdist(abundance_matrix, method = "bray")
    jaccard_dist <- vegdist(abundance_matrix, method = "jaccard")
    
    # Perform PERMANOVA
    metadata <- data %>% 
      select(Sample, housing) %>% 
      distinct() %>%
      arrange(match(Sample, rownames(abundance_matrix)))
    
    permanova_bray <- adonis2(bray_dist ~ housing, data = metadata, permutations = 999)
    permanova_jaccard <- adonis2(jaccard_dist ~ housing, data = metadata, permutations = 999)
    
    # Perform NMDS with error handling
    nmds_bray <- metaMDS(bray_dist, trymax = 100)
    nmds_jaccard <- metaMDS(jaccard_dist, trymax = 100)
    
    # Create NMDS plots with statistics
    create_nmds_plot <- function(nmds_obj, dist_type, permanova_result) {
      nmds_data <- data.frame(
        NMDS1 = nmds_obj$points[,1],
        NMDS2 = nmds_obj$points[,2],
        Sample = rownames(abundance_matrix)
      ) %>%
        left_join(metadata, by = "Sample")
      
      p_val <- round(permanova_result$`Pr(>F)`[1], 3)
      r2_val <- round(permanova_result$R2[1], 3)
      stress <- round(nmds_obj$stress, 3)
      
      ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = housing)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95) +
        annotate("text", x = min(nmds_data$NMDS1), y = max(nmds_data$NMDS2),
                 label = paste0("Stress = ", stress, "\n",
                                "R² = ", r2_val, "\n",
                                "p = ", p_val),
                 hjust = 0, vjust = 1) +
        scale_color_brewer(palette = "Set2") +
        labs(title = paste("NMDS Plot (", dist_type, ") -", sample_type),
             color = "Housing System") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)
        )
    }
    
    
    # Generate plots
    bray_plot <- create_nmds_plot(nmds_bray, "Bray-Curtis", permanova_bray)
    jaccard_plot <- create_nmds_plot(nmds_jaccard, "Jaccard", permanova_jaccard)
    
    # Save plots
    ggsave(paste0("beta_diversity_bray_", tolower(sample_type), ".png"), 
           bray_plot, width = 8, height = 6, dpi = 300)
    ggsave(paste0("beta_diversity_jaccard_", tolower(sample_type), ".png"), 
           jaccard_plot, width = 8, height = 6, dpi = 300)
    
    # Save PERMANOVA results
    write.csv(data.frame(
      Test = c("Bray-Curtis", "Jaccard"),
      P_value = c(permanova_bray$`Pr(>F)`[1], permanova_jaccard$`Pr(>F)`[1]),
      R2 = c(permanova_bray$R2[1], permanova_jaccard$R2[1])
    ), paste0("beta_diversity_stats_", tolower(sample_type), ".csv"), row.names = FALSE)
    
    return(list(
      permanova = list(bray = permanova_bray, jaccard = permanova_jaccard),
      plots = list(bray = bray_plot, jaccard = jaccard_plot)
    ))
  }, error = function(e) {
    message("Error in beta diversity calculation: ", e$message)
    return(NULL)
  })
}

# Run analyses for dust samples
dust_alpha <- calculate_alpha_diversity(dust_data, "Dust")
dust_beta <- calculate_beta_diversity(dust_data, "Dust")

# Run analyses for litter samples
litter_alpha <- calculate_alpha_diversity(litter_data, "Litter")
litter_beta <- calculate_beta_diversity(litter_data, "Litter")

# Create combined plot for each sample type
create_combined_plot <- function(alpha_results, beta_results, sample_type) {
  if(!is.null(beta_results)) {
    combined_plot <- ggarrange(
      alpha_results$plots$shannon,
      alpha_results$plots$simpson,
      alpha_results$plots$richness,
      beta_results$plots$bray,
      beta_results$plots$jaccard,
      ncol = 2, nrow = 3,
      common.legend = TRUE, legend = "right"
    )
    
    ggsave(paste0("combined_diversity_analysis_", tolower(sample_type), ".png"),
           combined_plot, width = 16, height = 20, dpi = 300)
  }
}

create_combined_plot(dust_alpha, dust_beta, "Dust")
create_combined_plot(litter_alpha, litter_beta, "Litter")
