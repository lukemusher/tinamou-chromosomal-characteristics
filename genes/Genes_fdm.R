setwd("~/Documents/ANSDU/Tinamous/cladeA_paper/github/genes/")

library(dplyr)
library(readr)
library(ggplot2)
genes<-read_tsv("full_table.tsv")

intr<-na.omit(read.csv("../abba-baba-windows/ABBABABAwindows.w100k.T2-cinP2_error_rm.csv"))
names(intr)
gc.tab<-read.csv("GC_100kb.csv")


# Add Z-scores and P-values to existing intr table
# Assumes your table is called 'intr' and has a column 'fdM'

add_fdM_statistics <- function(intr_table) {
  
  # Extract fdM values
  fdM <- intr_table$fdM
  
  # Calculate z-scores
  # Z = (value - mean) / standard deviation
  intr_table$fdM_zscore <- (fdM - mean(fdM, na.rm = TRUE)) / sd(fdM, na.rm = TRUE)
  
  # Calculate two-tailed p-values from z-scores
  # P = 2 * (1 - pnorm(|z|))
  intr_table$fdM_pvalue <- 2 * (1 - pnorm(abs(intr_table$fdM_zscore)))
  
  # Add absolute fdM for easier identification of extreme values
  intr_table$fdM_abs <- abs(fdM)
  
  # Apply multiple testing correction (FDR by default)
  intr_table$fdM_pvalue_adj <- p.adjust(intr_table$fdM_pvalue, method = "fdr")
  
  # Add significance flags
  intr_table$fdM_sig_005 <- intr_table$fdM_pvalue_adj < 0.05
  intr_table$fdM_sig_001 <- intr_table$fdM_pvalue_adj < 0.01
  
  return(intr_table)
}

# Simple one-liner version if you prefer:
add_fdM_stats_simple <- function(intr_table) {
  intr_table$fdM_zscore <- (intr_table$fdM - mean(intr_table$fdM, na.rm = TRUE)) / sd(intr_table$fdM, na.rm = TRUE)
  intr_table$fdM_pvalue <- 2 * (1 - pnorm(abs(intr_table$fdM_zscore)))
  intr_table$fdM_pvalue_adj <- p.adjust(intr_table$fdM_pvalue, method = "fdr")
  return(intr_table)
}

intr 

# Summary and exploration functions:
summarize_fdM_results <- function(intr_table) {
  cat("fdM Analysis Summary\n")
  cat("===================\n")
  cat("Total windows:", nrow(intr_table), "\n")
  cat("Mean fdM:", round(mean(intr_table$fdM, na.rm = TRUE), 4), "\n")
  cat("SD fdM:", round(sd(intr_table$fdM, na.rm = TRUE), 4), "\n")
  cat("Range fdM:", round(range(intr_table$fdM, na.rm = TRUE), 4), "\n\n")
  
  if("fdM_pvalue_adj" %in% names(intr_table)) {
    n_sig_005 <- sum(intr_table$fdM_pvalue_adj < 0.05, na.rm = TRUE)
    n_sig_001 <- sum(intr_table$fdM_pvalue_adj < 0.01, na.rm = TRUE)
    
    cat("Significant results (FDR adjusted):\n")
    cat("α = 0.05:", n_sig_005, "windows (", round(100*n_sig_005/nrow(intr_table), 1), "%)\n")
    cat("α = 0.01:", n_sig_001, "windows (", round(100*n_sig_001/nrow(intr_table), 1), "%)\n\n")
  }
  
  # Show most extreme values
  extreme_idx <- order(abs(intr_table$fdM), decreasing = TRUE)[1:min(10, nrow(intr_table))]
  cat("Top 10 most extreme fdM values:\n")
  cols_to_show <- intersect(c("fdM", "fdM_zscore", "fdM_pvalue", "fdM_pvalue_adj"), names(intr_table))
  print(intr_table[extreme_idx, cols_to_show])
}

# Get significant windows
get_significant_windows <- function(intr_table, alpha = 0.05, use_adjusted = TRUE) {
  p_col <- if(use_adjusted) "fdM_pvalue_adj" else "fdM_pvalue"
  
  if(!p_col %in% names(intr_table)) {
    stop("P-value column not found. Run add_fdM_statistics() first.")
  }
  
  sig_windows <- intr_table[intr_table[[p_col]] < alpha, ]
  sig_windows <- sig_windows[order(abs(sig_windows$fdM), decreasing = TRUE), ]
  
  return(sig_windows)
}

# Alternative multiple testing correction methods
apply_different_correction <- function(intr_table, method = "bonferroni") {
  # Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
  if(!"fdM_pvalue" %in% names(intr_table)) {
    stop("Raw p-values not found. Calculate p-values first.")
  }
  
  correction_col <- paste0("fdM_pvalue_", method)
  intr_table[[correction_col]] <- p.adjust(intr_table$fdM_pvalue, method = method)
  
  return(intr_table)
}


# Add statistics to your table
intr <- add_fdM_statistics(intr)

# Summarize results  
summarize_fdM_results(intr)

# Get significant windows
significant_windows <- get_significant_windows(intr, alpha = 0.05)

# Look at the most significant
head(significant_windows)

# Alternative corrections if needed
intr <- apply_different_correction(intr, method = "bonferroni")

# Count significant windows:
table(intr$fdM_sig_005)  # TRUE = significant at 0.05
table(intr$fdM_sig_001)  # TRUE = significant at 0.01

# New columns
intr$busco_ids <- vector("list", nrow(intr))
intr$busco_descriptions <- vector("list", nrow(intr))
intr$gene_count <- c()

# Loop through each row of intr
for(i in 1:nrow(intr)) {
  intr$busco_ids[[i]] <- genes %>%
    filter(scaffold == intr$scaffold[i] & 
             genes$start <= intr$end[i] & 
             genes$end >= intr$start[i]) %>%
    pull(busco.id)
  intr$gene_count[i]<-length(intr$busco_ids[[i]])
  intr$busco_descriptions[[i]] <- genes %>%
    filter(scaffold == intr$scaffold[i] & 
             genes$start <= intr$end[i] & 
             genes$end >= intr$start[i]) %>%
    pull(Description)
}

intr$busco_ids.500kb <- vector("list", nrow(intr))
intr$gene_count.500kb <- c()

for(i in 1:nrow(intr)) {
  intr$busco_ids.500kb[[i]] <- genes %>%
    filter(scaffold == intr$scaffold[i] & 
             genes$start <= intr$end[i]+500000 & 
             genes$end >= intr$start[i]-500000) %>%
    pull(busco.id)
  intr$gene_count.500kb[i]<-length(intr$busco_ids.500kb[[i]])
}

intr.sig<- intr[intr$fdM_pvalue_bonferroni < 0.0001,] %>% arrange(fdM_pvalue)
head(intr.sig)
nrow(intr.sig)
intr.sig

busco.list<-c()
counts=0

for (i in 1:nrow(intr.sig)){
  for (j in 1:length(intr.sig$busco_ids.500kb[[i]])){
    counts=counts+1
    print(intr.sig$busco_ids.500kb[[i]][j])
    if(length(intr.sig$busco_ids.500kb[[i]])>0){
      busco.list[counts]<-intr.sig$busco_ids.500kb[[i]][j]
    }
  }
}

busco.list<-na.omit(busco.list)

busco.tab.sig<-genes[genes$busco.id %in% busco.list,]
write.csv(busco.tab.sig, "busco.p0.0001.csv")

intr.sort<-intr %>% arrange(fdM_pvalue_bonferroni)

intr.sort$busco_ids <- vapply(intr.sort$busco_ids, paste, collapse = ",", character(1L))
intr.sort$busco_descriptions <- vapply(intr.sort$busco_descriptions, paste, collapse = ",", character(1L))
intr.sort$busco_ids.500kb <- vapply(intr.sort$busco_ids.500kb, paste, collapse = ",", character(1L))
write_tsv(intr.sort,file="./fdm_windows_buscos.tsv")

hist(intr.sort$gene_count.500kb)
hist(intr.sort$fdM)
plot(intr.sort$gene_count.500kb,intr.sort$fdM)
abline((lm(intr.sort$fdM~intr.sort$gene_count.500kb)))
summary(lm(intr.sort$fdM~intr.sort$gene_count.500kb))


# fdM plotting with stricter z-score thresholds

# Function to calculate z-score thresholds from Bonferroni correction
calculate_correct_thresholds <- function(intr_table) {
  
  # Calculate the number of tests (total windows with valid data)
  valid_data <- intr_table[complete.cases(intr_table[, c("fdM_pvalue_bonferroni", "fdM_zscore")]), ]
  n_tests <- nrow(valid_data)
  
  # Calculate raw p-value thresholds that correspond to Bonferroni thresholds
  # Bonferroni correction: raw_p = bonferroni_p / n_tests
  raw_p_0001 <- 0.0001 / n_tests
  raw_p_00001 <- 0.00001 / n_tests
  
  # Calculate corresponding z-score thresholds (two-tailed)
  z_0001 <- qnorm(1 - raw_p_0001/2)
  z_00001 <- qnorm(1 - raw_p_00001/2)
  
  # Handle cases where thresholds are extreme
  z_0001 <- ifelse(is.finite(z_0001), z_0001, 10)  # Cap at reasonable value
  z_00001 <- ifelse(is.finite(z_00001), z_00001, 10)
  
  return(list(
    n_tests = n_tests,
    z_0001 = z_0001,
    z_00001 = z_00001,
    raw_p_0001 = raw_p_0001,
    raw_p_00001 = raw_p_00001
  ))
}

# Manhattan plots function with threshold calculation
create_manhattan_plots <- function(intr_table, plot_width = 12, plot_height = 10) {
  
  # Ensure we have the required columns
  required_cols <- c("scaffold", "mid", "fdM_zscore", "fdM_pvalue_bonferroni", "gene_count")
  missing_cols <- required_cols[!required_cols %in% names(intr_table)]
  if(length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate correct z-score thresholds
  thresholds <- calculate_correct_thresholds(intr_table)
  
  # Remove rows with missing values
  intr_clean <- intr_table[complete.cases(intr_table[, required_cols]), ]
  
  # Create a factor for scaffolds to ensure consistent ordering
  intr_clean$scaffold <- factor(intr_clean$scaffold)
  scaffolds <- levels(intr_clean$scaffold)
  
  # Assign alternating colors to scaffolds
  scaffold_colors <- rep(c("#2E86AB", "#A23B72"), length.out = length(scaffolds))
  names(scaffold_colors) <- scaffolds
  
  # Add color column to data
  intr_clean$color <- scaffold_colors[as.character(intr_clean$scaffold)]
  
  # Create cumulative positions for continuous x-axis
  intr_clean <- intr_clean[order(intr_clean$scaffold, intr_clean$mid), ]
  
  # Calculate cumulative positions
  scaffold_lengths <- aggregate(mid ~ scaffold, intr_clean, function(x) max(x) - min(x) + 1)
  scaffold_lengths$cum_start <- c(0, cumsum(scaffold_lengths$mid[-nrow(scaffold_lengths)]))
  
  # Merge cumulative starts back to main data
  intr_clean <- merge(intr_clean, scaffold_lengths[, c("scaffold", "cum_start")], by = "scaffold")
  intr_clean$pos_cum <- intr_clean$mid + intr_clean$cum_start
  
  # Calculate scaffold midpoints for labeling
  scaffold_mids <- aggregate(pos_cum ~ scaffold, intr_clean, function(x) (max(x) + min(x))/2)
  
  # Set up the plotting area for three plots
  par(mfrow = c(3, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Plot 1: Z-scores
  plot(intr_clean$pos_cum, intr_clean$fdM_zscore, 
       col = intr_clean$color, 
       pch = 20, 
       cex = 1.2,
       xlab = "Genomic Position", 
       ylab = "fdM Z-score",
       main = "fdM Z-scores Across Scaffolds",
       xaxt = "n")
  
  # Add horizontal lines using calculated thresholds
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  abline(h = c(-thresholds$z_0001, thresholds$z_0001), col = "pink", lty = 2, lwd = 1)
  abline(h = c(-thresholds$z_00001, thresholds$z_00001), col = "red", lty = 2, lwd = 1)
  
  # Add scaffold labels
  axis(1, at = scaffold_mids$pos_cum, labels = scaffold_mids$scaffold, 
       las = 2, cex.axis = 1.2)
  
  # Add legend for z-score plot
  legend("topleft", 
         legend = c("p < 0.0001", "p < 0.00001"),
         col = c("pink", "red"),
         lty = 2,
         cex = 1.2)
  
  # Plot 2: -log10(Bonferroni p-values)
  neg_log_p <- -log10(intr_clean$fdM_pvalue_bonferroni)
  
  plot(intr_clean$pos_cum, neg_log_p,
       col = intr_clean$color, 
       pch = 20, 
       cex = 1.2,
       xlab = "Genomic Position", 
       ylab = "-log10(Bonferroni P-value)",
       main = "fdM Significance Across Scaffolds (Bonferroni Corrected)",
       xaxt = "n")
  
  # Add Bonferroni-adjusted significance threshold lines
  abline(h = -log10(0.0001), col = "pink", lty = 2, lwd = 2)
  abline(h = -log10(0.00001), col = "red", lty = 2, lwd = 2)
  
  # Add scaffold labels
  axis(1, at = scaffold_mids$pos_cum, labels = scaffold_mids$scaffold, 
       las = 2, cex.axis = 1.2)
  
  # Add legend for p-value plot
  legend("topleft",
         legend = c("p < 0.0001", "p < 0.00001"),
         col = c("pink", "red"),
         lty = 2,
         lwd = 2,
         cex = 1.2)
  
  # Plot 3: Gene count
  plot(intr_clean$pos_cum, intr_clean$gene_count,
       type = "n",
       xlab = "Genomic Position", 
       ylab = "Gene Count",
       main = "Gene Count Across Scaffolds",
       xaxt = "n")
  
  # Draw lines for each scaffold separately to maintain colors
  for(scaffold in scaffolds) {
    scaffold_data <- intr_clean[intr_clean$scaffold == scaffold, ]
    scaffold_data <- scaffold_data[order(scaffold_data$pos_cum), ]
    lines(scaffold_data$pos_cum, scaffold_data$gene_count, 
          col = scaffold_colors[scaffold], lwd = 1.5)
  }
  
  # Add scaffold labels
  axis(1, at = scaffold_mids$pos_cum, labels = scaffold_mids$scaffold, 
       las = 2, cex.axis = 1.2)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  # Return summary statistics
  n_scaffolds <- length(scaffolds)
  n_sig_0001 <- sum(intr_clean$fdM_pvalue_bonferroni < 0.0001)
  n_sig_00001 <- sum(intr_clean$fdM_pvalue_bonferroni < 0.00001)
  
  cat("Manhattan Plot Summary:\n")
  cat("======================\n")
  cat("Number of scaffolds:", n_scaffolds, "\n")
  cat("Total windows plotted:", nrow(intr_clean), "\n")
  cat("Number of tests (for Bonferroni):", thresholds$n_tests, "\n")
  cat("Z-score thresholds: ±", round(thresholds$z_0001, 2), " and ±", round(thresholds$z_00001, 2), "\n")
  cat("Significant at Bonferroni p < 0.0001:", n_sig_0001, "\n")
  cat("Significant at Bonferroni p < 0.00001:", n_sig_00001, "\n")
}

# Zoom plots function with automatic threshold calculation
create_zoom_plots <- function(intr_table, scaffold_name, start_pos = NULL, end_pos = NULL, 
                              plot_width = 12, plot_height = 10) {
  
  # Ensure we have the required columns
  required_cols <- c("scaffold", "mid", "fdM_zscore", "fdM_pvalue_bonferroni", "gene_count")
  missing_cols <- required_cols[!required_cols %in% names(intr_table)]
  if(length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate correct z-score thresholds from full dataset
  thresholds <- calculate_correct_thresholds(intr_table)
  
  # Filter to specific scaffold
  scaffold_data <- intr_table[intr_table$scaffold == scaffold_name, ]
  
  if(nrow(scaffold_data) == 0) {
    stop("No data found for scaffold: ", scaffold_name)
  }
  
  # Filter to window range if specified
  if(!is.null(start_pos)) {
    scaffold_data <- scaffold_data[scaffold_data$mid >= start_pos, ]
  }
  if(!is.null(end_pos)) {
    scaffold_data <- scaffold_data[scaffold_data$mid <= end_pos, ]
  }
  
  if(nrow(scaffold_data) == 0) {
    stop("No data found in the specified window range")
  }
  
  # Remove rows with missing values
  scaffold_data <- scaffold_data[complete.cases(scaffold_data[, required_cols]), ]
  
  # Order by position
  scaffold_data <- scaffold_data[order(scaffold_data$mid), ]
  
  # Determine color (use first color from alternating scheme)
  plot_color <- "#2E86AB"
  
  # Set up the plotting area for three plots
  par(mfrow = c(3, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Create title with range info
  range_text <- paste0("Scaffold: ", scaffold_name)
  if(!is.null(start_pos) || !is.null(end_pos)) {
    range_start <- if(is.null(start_pos)) min(scaffold_data$mid) else start_pos
    range_end <- if(is.null(end_pos)) max(scaffold_data$mid) else end_pos
    range_text <- paste0(range_text, " (", format(range_start, big.mark = ","), 
                         " - ", format(range_end, big.mark = ","), ")")
  }
  
  # Plot 1: Z-scores
  plot(scaffold_data$mid, scaffold_data$fdM_zscore, 
       col = plot_color, 
       type='l',
       pch = 20, 
       cex = 1.2,
       xlab = "Genomic Position (bp)", 
       ylab = "fdM Z-score",
       main = paste("fdM Z-scores -", range_text))
  
  # Add horizontal lines using calculated thresholds
  abline(h = 0, col = "black", lty = 2, lwd = 1)
  abline(h = c(-thresholds$z_0001, thresholds$z_0001), col = "pink", lty = 2, lwd = 1)
  abline(h = c(-thresholds$z_00001, thresholds$z_00001), col = "red", lty = 2, lwd = 1)
  
  # # Add legend for z-score plot
  # legend("bottomright", 
  #        legend = c("Bonferroni p < 0.0001", "Bonferroni p < 0.00001"),
  #        col = c("pink", "red"),
  #        lty = 2,
  #        cex = 1.2)
  
  # Plot 2: -log10(Bonferroni p-values)
  neg_log_p <- -log10(scaffold_data$fdM_pvalue_bonferroni)
  
  plot(scaffold_data$mid, neg_log_p,
       col = plot_color,
       type='l',
       pch = 20, 
       cex = 1.2,
       xlab = "Genomic Position (bp)", 
       ylab = "-log10(Bonferroni P-value)",
       main = paste("fdM Significance -", range_text))
  
  # Add Bonferroni-adjusted significance threshold lines
  abline(h = -log10(0.0001), col = "pink", lty = 2, lwd = 2)
  abline(h = -log10(0.00001), col = "red", lty = 2, lwd = 2)
  
  # # Add legend for p-value plot
  # legend("topright", 
  #        legend = c("Bonferroni p < 0.0001", "Bonferroni p < 0.00001"),
  #        col = c("pink", "red"),
  #        lty = 2,
  #        lwd = 2,
  #        cex = 0.8)
  
  # Plot 3: Gene count
  plot(scaffold_data$mid, scaffold_data$gene_count,
       type = "l",
       col = plot_color,
       #lwd = 1,
       xlab = "Genomic Position (bp)", 
       ylab = "Gene Count",
       main = paste("Gene Count -", range_text))
  
  # Add points to show actual data points
  # points(scaffold_data$mid, scaffold_data$gene_count,
  #        col = plot_color, pch = 20, cex = 0.6)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  # Return summary statistics for the region
  n_windows <- nrow(scaffold_data)
  n_sig_0001 <- sum(scaffold_data$fdM_pvalue_bonferroni < 0.0001)
  n_sig_00001 <- sum(scaffold_data$fdM_pvalue_bonferroni < 0.00001)
  
  cat("Zoom Plot Summary:\n")
  cat("==================\n")
  cat("Scaffold:", scaffold_name, "\n")
  cat("Position range:", format(min(scaffold_data$mid), big.mark = ","), 
      "to", format(max(scaffold_data$mid), big.mark = ","), "\n")
  cat("Windows plotted:", n_windows, "\n")
  cat("Mean fdM Z-score:", round(mean(scaffold_data$fdM_zscore), 3), "\n")
  cat("Empirical z-score thresholds: ±", round(thresholds$z_0001, 2), " and ±", round(thresholds$z_00001, 2), "\n")
  cat("Significant at Bonferroni p < 0.0001:", n_sig_0001, "\n")
  cat("Significant at Bonferroni p < 0.00001:", n_sig_00001, "\n")
  cat("Gene count range:", min(scaffold_data$gene_count), "to", max(scaffold_data$gene_count), "\n")
  
  # Return the filtered data for further analysis if needed
  invisible(scaffold_data)
}

# Plot across all scaffolds
create_manhattan_plots(intr)

# "zoom" into scaffold
create_zoom_plots(intr, scaffold_name = "10")

# "zoom" into more
create_zoom_plots(intr, scaffold_name = "10", start_pos = 5e06, end_pos = 1.5e07)
create_zoom_plots(intr, scaffold_name = "7", start_pos = 0, end_pos = 5e06)


# # Save plots:
# pdf(file = "../figs/fdM_manhattan_plots.pdf", width = 13, height = 9)
# create_manhattan_plots(intr)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr10.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "10", start_pos = 5e06, end_pos = 1.5e07)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr7.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "7", start_pos = 0, end_pos = 5e06)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr2a.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "2", start_pos = 1e08, end_pos = 1.5e08)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr2b.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "2", start_pos = 0, end_pos = 5e07)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr21.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "21", start_pos = 0, end_pos=810000)
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr15.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "15")
# dev.off()
# 
# pdf(file = "../figs/fdM_manhattan_plots_chr28.pdf", width = 6.5, height = 9)
# create_zoom_plots(intr, scaffold_name = "28")
# dev.off()


###add gc content to intr table

gc.cont <- c()
for (i in 1:length(intr[,1])){
  match_value <- gc.tab$gc_content[gc.tab$chromosome==intr$scaffold[i] & gc.tab$start_pos == intr$start[i]]
  
  if(length(match_value) > 0){
    gc.cont[i] <- match_value[1]  # Take first match if multiple exist
  } else {
    gc.cont[i] <- NA  # Assign NA when no match found
  }
}

intr.gc<-data.frame(cbind(intr,gc.cont))

plot(intr.gc$gc.cont,intr.gc$fdM)
abline(lm(intr.gc$fdM~intr.gc$gc.cont), col="red")
summary(lm(intr.gc$fdM~intr.gc$gc.cont))


summary(glm(fdM~gc.cont+gene_count.500kb, data=intr.gc))

