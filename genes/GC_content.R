rm(list=ls())
library(Biostrings)
library(ggplot2)
library(dplyr)

setwd("~/Documents/ANSDU/Tinamous/cladeA_paper/genes/")

# Function to calculate GC content in sliding windows
calculate_gc_windows <- function(fasta_file, window_size = 10000, step_size = 5000) {
  
  # Read FASTA file
  cat("Reading FASTA file...\n")
  genome <- readDNAStringSet(fasta_file)
  
  # Initialize results list
  results_list <- list()
  
  # Process each chromosome/sequence
  for (i in seq_along(genome)) {
    chrom_name <- names(genome)[i]
    cat("Processing", chrom_name, "...\n")
    
    seq <- genome[[i]]
    seq_length <- length(seq)
    
    # Skip sequences shorter than window size
    if (seq_length < window_size) {
      cat("  Skipping", chrom_name, "- too short (", seq_length, "bp)\n")
      next
    }
    
    # Calculate number of windows
    starts <- seq(1, seq_length - window_size + 1, by = step_size)
    
    # Calculate GC content for each window
    gc_values <- numeric(length(starts))
    positions <- numeric(length(starts))
    st <- numeric(length(starts))
    en <- numeric(length(starts))
    
    for (j in seq_along(starts)) {
      start <- starts[j]
      end <- min(start + window_size - 1, seq_length)
      
      # Extract window sequence
      window_seq <- subseq(seq, start, end)
      
      # Calculate GC content
      freq <- alphabetFrequency(window_seq, baseOnly = TRUE)
      gc_count <- freq["G"] + freq["C"]
      total_bases <- freq["A"] + freq["T"] + freq["G"] + freq["C"]
      
      # Store results
      gc_values[j] <- ifelse(total_bases > 0, (gc_count / total_bases) * 100, NA)
      positions[j] <- (start + end) / 2  # Middle of window
      st[j] <- start
      en[j] <-end
    }
    
    # Store chromosome results
    results_list[[i]] <- data.frame(
      chromosome = chrom_name,
      position = positions,
      gc_content = gc_values,
      start_pos = st,
      end_pos = en,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all chromosomes
  results <- do.call(rbind, results_list)
  
  return(results)
}

# Main execution
# Set parameters
fasta_file <- "Cuund34614.ragtag.fasta"  # Replace with your FASTA file
window_size <- 100000  # 100kb windows
step_size <- 100000     # 100kb steps 

# Calculate GC content
gc_data <- calculate_gc_windows(fasta_file, window_size, step_size)

write.csv(gc_data,file = "GC_100kb.csv") #manually remove non-chromosome-level scaffolds for downstream analysis

window_size <- 10000  # 10kb windows
step_size <- 10000     # 10kb steps 

# Calculate GC content
gc_data <- calculate_gc_windows(fasta_file, window_size, step_size)

write.csv(gc_data,file = "GC_10kb.csv") #manually remove non-chromosome-level scaffolds for downstream analysis

# Summarize GC results for 10kb windows

summary(gc_data$gc_content)
gc_data_fil<-gc_data[gc_data$chromosome %in% unique(gc_data$chromosome)[1:49],]
unique(gc_data_fil$chromosome)

# Summary by chromosome
chrom_summary <- gc_data_fil %>%
  group_by(chromosome) %>%
  summarize(
    mean_gc = mean(gc_content, na.rm = TRUE),
    median_gc = median(gc_content, na.rm = TRUE),
    sd_gc = sd(gc_content, na.rm = TRUE),
    min_gc = min(gc_content, na.rm = TRUE),
    max_gc = max(gc_content, na.rm = TRUE),
    n_windows = n()
  )

print(chrom_summary)

# MANHATTAN PLOT

# Filter data for specific chromosomes
chrom_pattern <- paste0("Chr", c(1:39, "Z"), "_RagTag")
gc_filtered <- gc_data %>%
  filter(chromosome %in% chrom_pattern)

# Order chromosomes 
gc_filtered <- gc_filtered %>%
  mutate(
    chrom_order = case_when(
      chromosome == "ChrZ_RagTag" ~ 40,
      TRUE ~ as.numeric(gsub("Chr([0-9]+)_RagTag", "\\1", chromosome))
    )
  ) %>%
  arrange(chrom_order, position)

# Calculate cumulative positions for Manhattan plot
gc_filtered_manhattan <- gc_filtered %>%
  group_by(chromosome) %>%
  mutate(max_pos = max(position)) %>%
  ungroup() %>%
  arrange(chrom_order, position) %>%
  group_by(chromosome) %>%
  mutate(
    cum_pos = position + sum(c(0, head(max_pos, -1)))
  ) %>%
  ungroup() %>%
  group_by(chromosome) %>%
  mutate(cum_pos = cum_pos - min(cum_pos) + 
           ifelse(chrom_order == 1, 0, 
                  lag(cumsum(max_pos), default = 0)[1])) %>%
  ungroup()

# Recalculate cumulative positions properly
gc_filtered_manhattan <- gc_filtered %>%
  group_by(chromosome, chrom_order) %>%
  summarise(max_pos = max(position), .groups = "drop") %>%
  arrange(chrom_order) %>%
  mutate(add_pos = lag(cumsum(max_pos), default = 0)) %>%
  select(chromosome, add_pos) %>%
  left_join(gc_filtered, by = "chromosome") %>%
  mutate(cum_pos = position + add_pos)

# Get chromosome centers for x-axis labels
axis_df_filtered <- gc_filtered_manhattan %>%
  group_by(chromosome, chrom_order) %>%
  summarize(center = mean(cum_pos), .groups = "drop") %>%
  arrange(chrom_order) %>%
  mutate(chrom_label = gsub("_RagTag", "", chromosome))

# Create filtered Manhattan plot
manhattan_filtered <- ggplot(gc_filtered_manhattan, 
                             aes(x = cum_pos, y = gc_content)) +
  geom_point(aes(color = as.factor(chrom_order %% 2)), 
             alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("#276FBF", "#183059")) +
  scale_x_continuous(label = axis_df_filtered$chrom_label, 
                     breaks = axis_df_filtered$center) +
  labs(
    title = "GC Content: Chr1-Chr39 and ChrZ",
    x = "Chromosome",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

print(manhattan_filtered)
ggsave("gc_content_manhattan_filtered.png", manhattan_filtered, 
       width = 16, height = 6, dpi = 300)

###quadratic models:
library(gridExtra)

# Read GC content data (from previous script output)
gc_data <- read.csv("GC_10kb.csv")

# Function to fit and plot quadratic model for a chromosome
fit_quadratic_model <- function(data, chrom_name) {
  
  # Filter data for specific chromosome
  chrom_data <- data %>%
    filter(chromosome == chrom_name) %>%
    filter(!is.na(gc_content))
  
  if (nrow(chrom_data) < 10) {
    cat("Not enough data points for", chrom_name, "\n")
    return(NULL)
  }
  
  # Normalize position for better numerical stability
  chrom_data <- chrom_data %>%
    mutate(position_scaled = (position - min(position)) / 1e6)  # Scale to Mb
  
  # Fit quadratic model: GC ~ a + b*pos + c*pos^2
  model <- lm(gc_content ~ position_scaled + I(position_scaled^2), 
              data = chrom_data)
  
  # Get model summary
  model_summary <- summary(model)
  
  # Extract coefficients
  coeffs <- coef(model)
  r_squared <- model_summary$r.squared
  adj_r_squared <- model_summary$adj.r.squared
  p_value <- pf(model_summary$fstatistic[1], 
                model_summary$fstatistic[2], 
                model_summary$fstatistic[3], 
                lower.tail = FALSE)
  
  # Print model statistics
  cat("\n=================================\n")
  cat("Chromosome:", chrom_name, "\n")
  cat("=================================\n")
  cat("Model: GC = a + b*pos + c*pos²\n")
  cat(sprintf("Coefficients:\n"))
  cat(sprintf("  Intercept (a): %.4f\n", coeffs[1]))
  cat(sprintf("  Linear (b):    %.6f\n", coeffs[2]))
  cat(sprintf("  Quadratic (c): %.8f\n", coeffs[3]))
  cat(sprintf("\nR-squared:     %.4f\n", r_squared))
  cat(sprintf("Adj R-squared: %.4f\n", adj_r_squared))
  cat(sprintf("P-value:       %.2e\n", p_value))
  
  # Generate predictions
  chrom_data$predicted <- predict(model, chrom_data)
  
  # Create plot
  plot <- ggplot(chrom_data, aes(x = position / 1e6)) +
    geom_point(aes(y = gc_content), alpha = 0.3, size = 1, color = "#276FBF") +
    geom_line(aes(y = predicted), color = "#D32F2F", size = 1.2) +
    labs(
      title = paste("Quadratic Model Fit:", gsub("_RagTag", "", chrom_name)),
      subtitle = sprintf("GC = %.2f + %.4f*pos + %.6f*pos²  (R² = %.4f)", 
                         coeffs[1], coeffs[2], coeffs[3], r_squared),
      x = "Position (Mb)",
      y = "GC Content (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  
  return(list(model = model, plot = plot, data = chrom_data, 
              summary = model_summary))
}

# Fit models for specific chromosomes
chromosomes_to_analyze <- all_chroms #c("Chr1_RagTag", "Chr2_RagTag", "Chr5_RagTag", "Chr10_RagTag", "ChrZ_RagTag")

results <- list()
plots <- list()

for (chrom in chromosomes_to_analyze) {
  if (chrom %in% gc_data$chromosome) {
    result <- fit_quadratic_model(gc_data, chrom)
    if (!is.null(result)) {
      results[[chrom]] <- result
      plots[[chrom]] <- result$plot
    }
  }
}


# ===== Compare all chromosomes =====
# Fit models for all Chr1-39 + ChrZ
chrom_pattern <- paste0("Chr", c(1:39, "Z"), "_RagTag")
all_chroms <- intersect(chrom_pattern, unique(gc_data$chromosome))

model_comparison <- data.frame()

for (chrom in all_chroms) {
  chrom_data <- gc_data %>%
    filter(chromosome == chrom, !is.na(gc_content)) %>%
    mutate(position_scaled = (position - min(position)) / 1e6)
  
  if (nrow(chrom_data) >= 10) {
    model <- lm(gc_content ~ position_scaled + I(position_scaled^2), 
                data = chrom_data)
    summary_stats <- summary(model)
    coeffs <- coef(model)
    
    model_comparison <- rbind(model_comparison, data.frame(
      chromosome = chrom,
      intercept = coeffs[1],
      linear_coef = coeffs[2],
      quadratic_coef = coeffs[3],
      r_squared = summary_stats$r.squared,
      adj_r_squared = summary_stats$adj.r.squared,
      mean_gc = mean(chrom_data$gc_content, na.rm = TRUE)
    ))
  }
}

# Display comparison table
print(model_comparison)

# Save comparison table
write.csv(model_comparison, "gc_quadratic_model_comparison.csv", row.names = FALSE)


####normalized window position on macrochromosomes

# Select chromosomes Chr1-8 and ChrZ
chromosomes_to_plot <- paste0("Chr", c(1:8, "Z"), "_RagTag")

# Filter and normalize positions within each chromosome
gc_normalized <- gc_data %>%
  filter(chromosome %in% chromosomes_to_plot) %>%
  filter(!is.na(gc_content)) %>%
  group_by(chromosome) %>%
  mutate(
    # Normalize position from 0 to 1
    position_norm = (position - min(position)) / (max(position) - min(position))
  ) %>%
  ungroup()

# Fit a single quadratic model to all combined data
model <- lm(gc_content ~ position_norm + I(position_norm^2), 
            data = gc_normalized)

# Print model summary
summary(model)
coeffs <- coef(model)
r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

# Generate smooth predictions
pred_positions <- seq(0, 1, length.out = 300)
predictions <- predict(model, 
                       newdata = data.frame(position_norm = pred_positions),
                       interval = "confidence", level = 0.95)

pred_df <- data.frame(
  position_norm = pred_positions,
  gc_predicted = predictions[, "fit"],
  lower = predictions[, "lwr"],
  upper = predictions[, "upr"]
)

# Create the plot
quadratic_plot <- ggplot() +
  # Plot raw data points
  geom_point(data = gc_normalized, 
             aes(x = position_norm, y = gc_content),
             alpha = 0.2, size = 0.8, color = "#276FBF") +
  # Add confidence interval
  geom_ribbon(data = pred_df,
              aes(x = position_norm, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "#D32F2F") +
  # Plot quadratic model fit
  geom_line(data = pred_df,
            aes(x = position_norm, y = gc_predicted),
            size = 1.5, color = "#D32F2F") +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    title = "Quadratic Model: GC Content by Normalized Position",
    subtitle = sprintf("Chr1-8 and ChrZ combined | GC = %.2f + %.2f*pos + %.2f*pos²  (R² = %.4f)",
                       coeffs[1], coeffs[2], coeffs[3], r_squared),
    x = "Normalized Position (% along chromosome)",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    panel.grid.minor = element_blank()
  )

print(quadratic_plot)

# Save plot
ggsave("gc_normalized_quadratic_combined.png", quadratic_plot, 
       width = 10, height = 7, dpi = 300)


#Final plot

# Select chromosomes Chr1-29 and ChrZ
chromosomes_to_plot <- paste0("Chr", c(1:29, "Z"), "_RagTag")

# Filter and normalize positions within each chromosome
gc_normalized <- gc_data %>%
  filter(chromosome %in% chromosomes_to_plot) %>%
  filter(!is.na(gc_content)) %>%
  group_by(chromosome) %>%
  mutate(
    # Normalize position from 0 to 1
    position_norm = (position - min(position)) / (max(position) - min(position)),
    # Create clean chromosome labels and ordering
    chrom_label = gsub("_RagTag", "", chromosome),
    chrom_order = ifelse(chromosome == "ChrZ_RagTag", 30,
                         as.numeric(gsub("Chr([0-9]+)_RagTag", "\\1", chromosome)))
  ) %>%
  ungroup() %>%
  arrange(chrom_order)

# Fit quadratic models for each chromosome and generate predictions
predictions_list <- list()
model_stats <- data.frame()

for (chrom in chromosomes_to_plot) {
  chrom_data <- gc_normalized %>%
    filter(chromosome == chrom)
  
  if (nrow(chrom_data) >= 10) {
    # Fit quadratic model
    model <- lm(gc_content ~ position_norm + I(position_norm^2), 
                data = chrom_data)
    
    # Generate smooth predictions
    pred_positions <- seq(0, 1, length.out = 100)
    predictions <- predict(model, 
                           newdata = data.frame(position_norm = pred_positions))
    
    # Store predictions
    chrom_label <- gsub("_RagTag", "", chrom)
    chrom_order <- ifelse(chrom == "ChrZ_RagTag", 30,
                          as.numeric(gsub("Chr([0-9]+)_RagTag", "\\1", chrom)))
    
    pred_df <- data.frame(
      chromosome = chrom,
      chrom_label = chrom_label,
      chrom_order = chrom_order,
      position_norm = pred_positions,
      gc_predicted = predictions
    )
    predictions_list[[chrom]] <- pred_df
    
    # Store model statistics
    coeffs <- coef(model)
    r_squared <- summary(model)$r.squared
    
    # Calculate p-value
    model_summary <- summary(model)
    f_stat <- model_summary$fstatistic
    p_value <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    
    model_stats <- rbind(model_stats, data.frame(
      chrom_label = chrom_label,
      chrom_order = chrom_order,
      intercept = coeffs[1],
      linear = coeffs[2],
      quadratic = coeffs[3],
      r_squared = r_squared,
      p_value = p_value
    ))
  }
}

# Combine all predictions
all_predictions <- do.call(rbind, predictions_list)

# Add R² and p-value labels for facets
label_data <- model_stats %>%
  mutate(facet_label = sprintf("%s\nR² = %.3f, p = %.2e", 
                               chrom_label, r_squared, p_value))

# Add labels to data frames
gc_normalized <- gc_normalized %>%
  left_join(label_data %>% select(chrom_label, facet_label), by = "chrom_label")

all_predictions <- all_predictions %>%
  left_join(label_data %>% select(chrom_label, facet_label), by = "chrom_label")

# Reorder factor levels for proper facet ordering
gc_normalized$facet_label <- factor(gc_normalized$facet_label,
                                    levels = label_data$facet_label[order(label_data$chrom_order)])
all_predictions$facet_label <- factor(all_predictions$facet_label,
                                      levels = label_data$facet_label[order(label_data$chrom_order)])

# Create faceted plot
facet_plot <- ggplot() +
  geom_point(data = gc_normalized, 
             aes(x = position_norm, y = gc_content),
             alpha = 0.3, size = 0.5, color = "#276FBF") +
  geom_line(data = all_predictions,
            aes(x = position_norm, y = gc_predicted),
            size = 0.8, color = "#D32F2F") +
  facet_wrap(~ facet_label, ncol = 5) +
  scale_x_continuous(labels = scales::percent_format(),
                     breaks = c(0, 0.5, 1)) +
  labs(
    title = "Quadratic Models of GC Content by Normalized Position",
    subtitle = "Chromosomes 1-29 and Z",
    x = "Normalized Position (% along chromosome)",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 8, lineheight = 1.1),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )

print(facet_plot)

# Save plot (adjusted dimensions for 5 columns)
ggsave("gc_quadratic_faceted_chr1-29_chrz.png", facet_plot, 
       width = 15, height = 24, dpi = 300)

# Print model statistics
cat("\n=================================\n")
cat("Quadratic Model Statistics\n")
cat("=================================\n")
model_stats <- model_stats %>% arrange(chrom_order)
print(model_stats)

# Save model statistics
write.csv(model_stats, "gc_quadratic_model_stats.csv", row.names = FALSE)


