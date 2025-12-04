# R script to extract chromosome and position information from tree file names and generate chromosome bar plot
# Expected filename format: "Scaffold_name.windowstart-windowend-aln.fasta.treefile"; e.g., "Chr1_RagTag.2600001-2700000-aln.fasta.treefile"

library(stringr)
library(ape)

setwd("~/Documents/ANSDU/Tinamous/cladeA_paper/github/windowed-gene-trees/")

# Function to parse filename, extract genomic coordinates, and read file contents
parse_treefile_name <- function(filename, directory_path) {
  # Remove the file extension (.treefile)
  base_name <- str_remove(filename, "\\.treefile$")
  
  # Extract chromosome (everything before the first underscore after "Chr")
  chr_match <- str_extract(base_name, "^Chr[^_]+")
  
  # Extract the position range (numbers separated by hyphen)
  pos_match <- str_extract(base_name, "\\d+-\\d+")
  
  if (is.na(chr_match) || is.na(pos_match)) {
    return(data.frame(filename = filename, chromosome = NA, start_pos = NA, end_pos = NA, tree = NA))
  }
  
  # Split the position range
  positions <- str_split(pos_match, "-")[[1]]
  start_pos <- as.numeric(positions[1])
  end_pos <- as.numeric(positions[2])
  
  # Read the file contents
  file_path <- file.path(directory_path, filename)
  tree_content <- NA
  
  tryCatch({
    tree_content <- readLines(file_path, warn = FALSE)
    tree_content <- paste(tree_content, collapse = "\n")
  }, error = function(e) {
    warning(paste("Could not read file:", filename, "-", e$message))
  })
  
  return(data.frame(
    filename = filename,
    chromosome = chr_match,
    start_pos = start_pos,
    end_pos = end_pos,
    tree = tree_content
  ))
}

# Main function to process all files in directory
process_treefile_directory <- function(directory_path = ".") {
  # Get all .treefile files in the directory
  treefile_pattern <- "\\.treefile$"
  files <- list.files(directory_path, pattern = treefile_pattern, full.names = FALSE)
  
  if (length(files) == 0) {
    warning("No .treefile files found in the specified directory.")
    return(data.frame(filename = character(0), chromosome = character(0), 
                      start_pos = numeric(0), end_pos = numeric(0), tree = character(0)))
  }
  
  cat("Found", length(files), "treefile(s) to process.\n")
  
  # Process each file and combine results
  results_list <- lapply(files, function(f) parse_treefile_name(f, directory_path))
  results_df <- do.call(rbind, results_list)
  
  # Sort by chromosome and start position
  results_df <- results_df[order(results_df$chromosome, results_df$start_pos), ]
  
  # Reset row names
  rownames(results_df) <- NULL
  
  return(results_df)
}

#path to directory with sliding window trees
directory_path <- "./all_gts"  # Modify this path to your directory

# Process the files and create the table
genomic_coords_table <- process_treefile_directory(directory_path)

head(genomic_coords_table)

# Save to CSV file
write.csv(genomic_coords_table, "tree_genomic_coordinates.csv", row.names = FALSE)

#Define T1 through T5
library(ape)
par(mfrow=c(2,2))
topo<-c() 
for (i in 1:length(genomic_coords_table$filename)){
  t<-drop.tip(root(read.tree(text=genomic_coords_table$tree[i]),outgroup = "Cvar27099"), "Crnoc59478")
  if(is.monophyletic(phy = t, tips = c("Cratr320360","Cstr9577","Ceery21146")) && is.monophyletic(phy = t, tips = c("Cratr320360","Cstr9577","Ceery21146","Ccgol2213")) && is.monophyletic(phy = t, tips = c("Cstr9577","Ceery21146"))){
    topo[i]<-"T1"
  }
  if(is.monophyletic(phy = t, tips = c("Ccgol2213","Cstr9577","Ceery21146")) && is.monophyletic(phy = t, tips = c("Cratr320360","Ccgol2213","Cstr9577","Ceery21146")) && is.monophyletic(phy = t, tips = c("Cstr9577","Ceery21146"))){
    topo[i]<-"T2"
  }
  # if(is.monophyletic(phy = t, tips = c("Ccgol2213","Cratr320360")) && is.monophyletic(phy = t, tips = c("Ceery21146","Cstr9577"))){
  #   topo[i]<-"T3"
  # }
  if(is.monophyletic(phy = t, tips = c("Cstr9577","Cratr320360")) && is.monophyletic(phy = t, tips = c("Ceery21146","Cstr9577","Cratr320360"))){
    topo[i]<-"T4"
  }
  if(is.monophyletic(phy = t, tips = c("Ceery21146","Cratr320360","Ccgol2213"))&& is.monophyletic(phy = t, tips = c("Crnoc59478","Cstr9577"))){
    topo[i]<-"T5"
  }
}

#All other topologies<-"Tx"
for (i in 1:length(topo)){
  if(is.na(topo[i])){
    topo[i]<-"Tx"
  }
}

#generate table of topology and window positions
tab1<-data.frame(cbind(genomic_coords_table$chromosome,as.numeric(genomic_coords_table$start_pos),as.numeric(genomic_coords_table$end_pos),topo))
names(tab1)<-c("chromosome","start_pos","end_pos","topo")
head(tab1)

tab1$start_pos<-as.numeric(tab1$start_pos)
tab1$end_pos<-as.numeric(tab1$end_pos)

tab1$chromosome<-factor(tab1$chromosome,levels = 
                         c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7",
                           "Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14",
                           "Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21",
                           "Chr22","Chr23","Chr24","Chr25","Chr26","Chr27","Chr28",
                           "Chr29","Chr30","Chr31","Chr32","Chr33","Chr34","Chr35",
                           "Chr36","Chr37","Chr38","Chr39","ChrZ"))


write.csv(tab1,file="topo.position.csv")
# R script to create chromosome visualization with topology coloring
# Assumes you have a table called "tab1" with columns: chromosome, start_pos, end_pos, topo

library(ggplot2)
library(dplyr)

# Function to create chromosome topology plot
create_chromosome_plot <- function(data) {
  
  # Define custom chromosome order (1-39, then Z)
  chr_order <- c(paste0("Chr", 1:39), "ChrZ")
  
  # Calculate chromosome lengths (max end position for each chromosome)
  chr_lengths <- data %>%
    group_by(chromosome) %>%
    summarise(max_end = max(end_pos, na.rm = TRUE), .groups = 'drop') %>%
    mutate(chromosome = factor(chromosome, levels = chr_order)) %>%
    arrange(chromosome)
  
  # Create segments for plotting - each row becomes a colored segment
  plot_data <- data %>%
    left_join(chr_lengths, by = "chromosome") %>%
    mutate(
      # Create y-axis positions for each chromosome with custom order
      chromosome = factor(chromosome, levels = chr_order),
      chr_num = as.numeric(chromosome)
    ) %>%
    arrange(chromosome, start_pos)
  
  # Get unique topologies for color mapping
  unique_topos <- sort(unique(plot_data$topo))
  
  # Create the plot
  p <- ggplot() +
    # Add chromosome background bars (optional - shows full chromosome length)
    geom_rect(data = chr_lengths %>% 
                mutate(chr_num = as.numeric(chromosome)),
              aes(xmin = 0, xmax = max_end, 
                  ymin = chr_num - 0.4, ymax = chr_num + 0.4),
              fill = "lightgray", color = "black", alpha = 0.3) +
    
    # Add colored segments for each topology region
    geom_rect(data = plot_data,
              aes(xmin = start_pos, xmax = end_pos, 
                  ymin = chr_num - 0.35, ymax = chr_num + 0.35,
                  fill = factor(topo)),, size = 0.1) +
    
    # Customize the plot
    scale_y_continuous(
      breaks = 1:nrow(chr_lengths),
      labels = levels(chr_lengths$chromosome),
      expand = c(0.02, 0)
    ) +
    scale_x_continuous(
      labels = scales::comma_format(scale = 1e-6, suffix = "M"),
      expand = c(0.01, 0)
    ) +
    scale_fill_manual(
      values = c("#1F78B4", "#E31A1C", "#FFD700", "#000000", "#FB9A99", 
                 "#FF1493", "#33A02C", "#B2DF8A", "#FF7F00", "#FDBF6F", 
                 "#CAB2D6", "#B15928", "#6A3D9A", "#808080", "#A6CEE3"),
      name = "Topology"
    ) +
    
    # Labels and theme
    labs(
      title = "Chromosome Topology Distribution",
      x = "Position (Mb)",
      y = "Chromosome",
      caption = "Each bar represents a chromosome with regions colored by topology"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  
  return(p)
}

# Create the plot
chromosome_plot <- create_chromosome_plot(tab1)

# Display the plot
print(chromosome_plot)

# save the plot
pdf(file = "200kb_trees_chromosomes.pdf", width = 9, height = 6.5)
print(chromosome_plot)
dev.off()
