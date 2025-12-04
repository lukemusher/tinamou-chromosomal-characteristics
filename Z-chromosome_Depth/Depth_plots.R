library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)

setwd("~/Documents/ANSDU/Tinamous/cladeA_paper/z-chrom-depth/")

# Define your female depth files and corresponding sample names
depth_files <- c(
  "Cstr9577_chrz_depth.txt",
  "Crnoc59478_chrz_depth.txt"
  # Add more files as needed
)

sample_names <- c(
  "C. strigulosus (female)",
  "C. noctivagus (female)"
  # Add corresponding sample names
)

# Read and combine all depth files
read_depth_file <- function(file_path, sample_name) {
  read_tsv(file_path, 
           col_names = c("chr", "pos", "depth"),
           col_types = "cii") %>%
    mutate(sample = sample_name) %>%
    filter(!is.na(depth) & depth >= 0)
}

# Combine all samples into one dataframe
depth_data <- map2_dfr(depth_files, sample_names, read_depth_file)

# Define window size
window_size <- 50000

# Calculate window assignments and averages
window_avg <- depth_data %>%
  mutate(window = floor((pos - 1) / window_size) * window_size + 1,
         window_center = window + window_size/2) %>%
  group_by(sample, window, window_center) %>%
  summarise(
    avg_depth = mean(depth, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    n_sites = n(),
    .groups = 'drop'
  ) %>%
  filter(n_sites >= 5) %>%  # Filter windows with at least 5 sites
  filter(avg_depth < 50) %>%  # Remove windows with depth > 100
  group_by(sample) %>%
  mutate(
    normalized_depth = (avg_depth - min(avg_depth, na.rm = TRUE)) / 
      (max(avg_depth, na.rm = TRUE) - min(avg_depth, na.rm = TRUE))
  ) %>%
  ungroup()

# Print summary statistics
cat("Summary statistics for ChrZ_RagTag:\n")
cat("Samples found:", paste(unique(window_avg$sample), collapse = ", "), "\n")
cat("Total 50kb windows per sample (after filtering depth >100):", nrow(window_avg)/length(unique(window_avg$sample)), "\n")
cat("Scaffold length (approx):", max(depth_data$pos), "bp\n")

# Print per-sample statistics
sample_stats <- window_avg %>%
  group_by(sample) %>%
  summarise(
    mean_depth = round(mean(avg_depth), 2),
    median_depth = round(median(avg_depth), 2),
    min_depth = round(min(avg_depth), 2),
    max_depth = round(max(avg_depth), 2),
    mean_normalized = round(mean(normalized_depth), 3),
    median_normalized = round(median(normalized_depth), 3),
    .groups = 'drop'
  )
print(sample_stats)

# Create stacked plots using facet_wrap
p1_females <- ggplot(window_avg, aes(x = window_center/1e6, y = avg_depth)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(
    title = "Normalized Read Depth Across ChrZ_RagTag Scaffold",
    subtitle = paste("50kb windows,", length(unique(window_avg$sample)), "samples (depth >100 removed)"),
    x = "Position (Mb)",
    y = "Normalized Depth (0-1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p1_females)

# Alternative: Create overlaid plot with different colors per sample
p2_females <- ggplot(window_avg, aes(x = window_center/1e6, y = avg_depth, color = sample)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=25), 
              linewidth = 1, se = FALSE) +
  labs(
    title = "Female sequence coverage across the Z-chromosome",
    subtitle = paste("50kb windows,", length(unique(window_avg$sample)), "samples (depth >75 removed, with smoothed splines)"),
    x = "Position (Mb)",
    y = "Normalized read depth",
    color = "Sample"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d()  # Use a colorblind-friendly palette

print(p2_females)

# Save the processed data
write_csv(window_avg, "females_chrz_50kb_windows_multi_sample.csv")

#######do the same for males
# Define your male depth files and corresponding sample names
depth_files <- c(
  "Ccgol2213_chrz_depth.txt",
  "Ceery21146_chrz_depth.txt",
  "Cratr320360_chrz_depth.txt"
  # Add more files as needed
)

sample_names <- c(
  "C. cinnamomeus (male)",
  "C. erythropus (male)",
  "C. atrocapillus (male)"
  # Add corresponding sample names
)

# Read and combine all depth files
read_depth_file <- function(file_path, sample_name) {
  read_tsv(file_path, 
           col_names = c("chr", "pos", "depth"),
           col_types = "cii") %>%
    mutate(sample = sample_name) %>%
    filter(!is.na(depth) & depth >= 0)
}

# Combine all samples into one dataframe
depth_data <- map2_dfr(depth_files, sample_names, read_depth_file)

# Define window size
window_size <- 50000

# Calculate window assignments and averages
window_avg <- depth_data %>%
  mutate(window = floor((pos - 1) / window_size) * window_size + 1,
         window_center = window + window_size/2) %>%
  group_by(sample, window, window_center) %>%
  summarise(
    avg_depth = mean(depth, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    n_sites = n(),
    .groups = 'drop'
  ) %>%
  filter(n_sites >= 5) %>%  # Filter windows with at least 5 sites
  filter(avg_depth < 100) %>%  # Remove windows with depth > 100
  group_by(sample) %>%
  mutate(
    normalized_depth = (avg_depth - min(avg_depth, na.rm = TRUE)) / 
      (max(avg_depth, na.rm = TRUE) - min(avg_depth, na.rm = TRUE))
  ) %>%
  ungroup()

# Print summary statistics
cat("Summary statistics for ChrZ_RagTag:\n")
cat("Samples found:", paste(unique(window_avg$sample), collapse = ", "), "\n")
cat("Total 50kb windows per sample (after filtering depth >100):", nrow(window_avg)/length(unique(window_avg$sample)), "\n")
cat("Scaffold length (approx):", max(depth_data$pos), "bp\n")

# Print per-sample statistics
sample_stats <- window_avg %>%
  group_by(sample) %>%
  summarise(
    mean_depth = round(mean(avg_depth), 2),
    median_depth = round(median(avg_depth), 2),
    min_depth = round(min(avg_depth), 2),
    max_depth = round(max(avg_depth), 2),
    mean_normalized = round(mean(normalized_depth), 3),
    median_normalized = round(median(normalized_depth), 3),
    .groups = 'drop'
  )
print(sample_stats)

# Create stacked plots using facet_wrap
p1_males <- ggplot(window_avg, aes(x = window_center/1e6, y = avg_depth)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(
    title = "Normalized Read Depth Across ChrZ_RagTag Scaffold",
    subtitle = paste("50kb windows,", length(unique(window_avg$sample)), "samples (depth >100 removed)"),
    x = "Position (Mb)",
    y = "Normalized Depth (0-1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p1_males)

# Alternative: Create overlaid plot with different colors per sample
p2_males <- ggplot(window_avg, aes(x = window_center/1e6, y = avg_depth, color = sample)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=25), 
              linewidth = 1, se = FALSE) +
  labs(
    title = "Male sequence coverage across the Z-chromosome",
    subtitle = paste("50kb windows,", length(unique(window_avg$sample)), "samples (depth >100 removed, with smoothed splines)"),
    x = "Position (Mb)",
    y = "Read Depth",
    color = "Sample"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d()  # Use a colorblind-friendly palette

print(p2_males)

# Save the processed data
write_csv(window_avg, "males_chrz_50kb_windows_multi_sample.csv")


###introgression_plot


tab2<-read.csv("../abba-baba-windows/ABBABABAwindows.w100k.T2-cinP2.csv")
names(tab2)

tab3<-tab2[tab2$scaffold==40,]
p3_smooth <- ggplot(tab3, aes(x = mid/1e6, y = fdM)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=15), 
              color = "red", linewidth = 1, se = F, alpha = 0.3) +
  labs(
    title = "Introgression Across Chromosome Z",
    subtitle = paste("100kb windows (with smoothed spline)"),
    x = "Position (Mb)",
    y = "fdM"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p3_smooth)


library(ggpubr)

pdf(file = "../figs/depth-introgression-plot3.pdf", width = 6.25, height = 8.25, bg = "white")
ggarrange(p2_females, p2_males, p3_smooth,
          ncol = 1, nrow = 2)
dev.off()



#####other plots

# Define your depth files and corresponding sample names
depth_files <- c(
  "Cratr320360_chrz_depth.txt",
  "Cstr9577_chrz_depth.txt",
  "Crnoc59478_chrz_depth.txt",
  "Cvar27099_chrz_depth.tx"
  
  # Add more files as needed
)

sample_names <- c(
  "C. atrocapillus (male)",
  "C. strigulosus (female)",
  "C. noctivagus (female)",
  "C. variegatus (female)"
  
  # Add corresponding sample names
)

# Read and combine all depth files
read_depth_file <- function(file_path, sample_name) {
  read_tsv(file_path, 
           col_names = c("chr", "pos", "depth"),
           col_types = "cii") %>%
    mutate(sample = sample_name) %>%
    filter(!is.na(depth) & depth >= 0)
}

# Combine all samples into one dataframe
depth_data <- map2_dfr(depth_files, sample_names, read_depth_file)

# Define window size
window_size <- 50000

# Calculate window assignments and averages
window_avg <- depth_data %>%
  mutate(window = floor((pos - 1) / window_size) * window_size + 1,
         window_center = window + window_size/2) %>%
  group_by(sample, window, window_center) %>%
  summarise(
    avg_depth = mean(depth, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    n_sites = n(),
    .groups = 'drop'
  ) %>%
  filter(n_sites >= 5) %>%  # Filter windows with at least 5 sites
  filter(avg_depth < 100) %>%  # Remove windows with depth > 100
  group_by(sample) %>%
  mutate(
    normalized_depth = (avg_depth - min(avg_depth, na.rm = TRUE)) / 
      (max(avg_depth, na.rm = TRUE) - min(avg_depth, na.rm = TRUE))
  ) %>%
  ungroup()

p1_smooth <- ggplot(window_avg, aes(x = window_center/1e6, y = avg_depth)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=50), 
              color = "red", linewidth = 1, se = TRUE, alpha = 0.3) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(
    title = "Normalized Read Depth Across ChrZ_RagTag Scaffold",
    subtitle = paste("50kb windows,", length(unique(window_avg$sample)), "samples (depth >100 removed, with smoothed spline)"),
    x = "Position (Mb)",
    y = "Normalized Depth (0-1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p1_smooth)

# Create a comparison plot showing depth ratios (useful for CNV detection)
if(length(unique(window_avg$sample)) >= 2) {
  # Calculate ratios relative to first sample
  reference_sample <- unique(window_avg$sample)[1]
  
  ratio_data <- window_avg %>%
    select(window_center, sample, avg_depth) %>%
    pivot_wider(names_from = sample, values_from = avg_depth) %>%
    pivot_longer(cols = -c(window_center, all_of(reference_sample)), 
                 names_to = "sample", 
                 values_to = "depth") %>%
    rename(ref_depth = all_of(reference_sample)) %>%
    mutate(depth_ratio = depth / ref_depth) %>%
    filter(!is.na(depth_ratio) & is.finite(depth_ratio))
  
  p3 <- ggplot(ratio_data, aes(x = window_center/1e6, y = depth_ratio)) +
    geom_point(color = "steelblue", linewidth = 0.8) +
    geom_smooth(method = "gam", formula = y ~ s(x, k=50), 
                color = "red", linewidth = 1, se = TRUE, alpha = 0.3)+
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = c(0.5, 2), linetype = "dotted", alpha = 0.5) +
    facet_wrap(~sample, ncol = 1) +
    labs(
      title = paste("Depth Ratios Relative to", reference_sample),
      subtitle = "Dashed line = 1:1 ratio, dotted lines = 2-fold change",
      x = "Position (Mb)",
      y = paste("Depth Ratio (relative to", reference_sample, ")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA)
    )
  
  print(p3)
  ggsave("chrz_depth_ratios.png", p3, width = 12, height = 8, dpi = 300)
}

pdf(file = "../figs/depth_ratio-introgression-plot2.pdf", width = 6.5, height = 8.5, bg = "white")
ggarrange(p3, p3_smooth,
          ncol = 1, nrow = 2)
dev.off()
