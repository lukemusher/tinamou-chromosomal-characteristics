rm(list = ls())
###distance from center
library(ggplot2)
library(ggpubr)
library(ape)
library(phytools)
library(phangorn)
library(FSA)
library(dplyr)
library(readr)
library(stringr)

#laptop
setwd("~/Documents/ANSDU/Tinamous/cladeA_paper/github/")

# Read abba-baba-output files
intr.100kb <- read_csv("abba-baba-windows/ABBABABAwindows.w100k.T2-cinP2_error_rm.csv") #100kb windows assuming T2 after removing erroneous T4 windows
intr.200kb <- read_csv("abba-baba-windows/ABBABABAwindows.w200k.T2-cinP2.csv") #200kb windows assumint T2


#standard error function
se<-function(x){
  std.err<-sd(x)/sqrt(length(x))
  print(std.err)
}


# The first section of this R script is to extract chromosome and 
# position information from tree file names and generate chromosome bar plot
# Expected filename format: "Scaffold_name.windowstart-windowend-aln.fasta.treefile"; 
# e.g., "Chr1_RagTag.2600001-2700000-aln.fasta.treefile"

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

# Function to process all files in directory
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

#Execute the script
# Change the directory path as needed - use "." for current directory
directory_path <- "windowed-gene-trees/all_gts"  # Modify this path to your directory

# Process the files and create the table
# genomic_coords_table <- process_treefile_directory(directory_path)


# Save csv file
# write.csv(genomic_coords_table, "windowed-gene-trees/window_trees_genomic_coordinates.csv", row.names = FALSE)

#read csv file
genomic_coords_table<-read.csv("windowed-gene-trees/window_trees_genomic_coordinates.csv")

#define topologies T1 through T5
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
  if(is.monophyletic(phy = t, tips = c("Cstr9577","Cratr320360")) && is.monophyletic(phy = t, tips = c("Ceery21146","Cstr9577","Cratr320360"))){
    topo[i]<-"T3"
  }
  if(is.monophyletic(phy = t, tips = c("Ceery21146","Cratr320360","Ccgol2213"))&& is.monophyletic(phy = t, tips = c("Crnoc59478","Cstr9577"))){
    topo[i]<-"T4"
  }
  if(is.monophyletic(phy = t, tips = c("Cratr320360","Cstr9577"))&& is.monophyletic(phy = t, tips = c("Ceery21146","Ccgol2213"))){
    topo[i]<-"T5"
  }
}

#alt topos <- "Tx"
for (i in 1:length(topo)){
  if(is.na(topo[i])){
    topo[i]<-"Tx"
  }
}

tree.table<-data.frame(cbind(genomic_coords_table$chromosome,as.numeric(genomic_coords_table$start_pos),as.numeric(genomic_coords_table$end_pos),topo))
names(tree.table)<-c("chromosome","start_pos","end_pos","topo")
head(tree.table)

tree.table$start_pos<-as.numeric(tree.table$start_pos)
tree.table$end_pos<-as.numeric(tree.table$end_pos)

tree.table$chromosome<-factor(tree.table$chromosome,levels = 
                          c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7",
                            "Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14",
                            "Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21",
                            "Chr22","Chr23","Chr24","Chr25","Chr26","Chr27","Chr28",
                            "Chr29","Chr30","Chr31","Chr32","Chr33","Chr34","Chr35",
                            "Chr36","Chr37","Chr38","Chr39","ChrZ"))
chr<-c()
fdM.val<-c()
st<-c()
en<-c()
micro.macro<-c()
mid<-c()
ABBA<-c()

for (i in 1:length(tree.table$chromosome)){
  chr[i]<-str_split(tree.table$chromosome[i],"Chr")[[1]][2]
  if(chr[i]=="Z"){
    chr[i]<-40
  }
  st[i]<-tree.table$start_pos[i]
  en[i]<-tree.table$end_pos[i]
  fdM.val[i]<-intr.200kb$fdM[intr.200kb$start==st[i] & intr.200kb$end==en[i] & intr.200kb$scaffold==chr[i]]
  ABBA[i]<-intr.200kb$ABBA[intr.200kb$start==st[i] & intr.200kb$end==en[i] & intr.200kb$scaffold==chr[i]]
  mid[i]<-intr.200kb$mid[intr.200kb$start==st[i] & intr.200kb$end==en[i] & intr.200kb$scaffold==chr[i]]
  micro.macro[i]<-"micro"
  if(as.numeric(chr[i])<=8){
    micro.macro[i]<-"macro"
  }
  if(chr[i]=="40"){
    micro.macro[i]<-"Z-chromosome"
  }
}

tree.fdm.tab<-data.frame(cbind(tree.table,chr,mid,ABBA,fdM.val,micro.macro))

############################################################
##################### Figure1 barplots #####################
############################################################

# Define custom colors for topo values
topo_colors <- c("T1" = "#1F78B4", 
                 "T2" = "#E31A1C", 
                 "T3" = "#FFD700", 
                 "T4" = "#000000", 
                 "T5" = "#FB9A99",
                 "Tx" = "lightgray")

# Create stacked bar plot with custom colors
ggplot(tree.fdm.tab, aes(x = chromosome, fill = topo)) +
  geom_bar(position = "stack") +
  labs(
    title = "Topologies by Chromosome",
    x = "Chromosome",
    y = "Count",
    fill = "Topology"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors)

# Create stacked bar plot with custom colors
ggplot(tree.fdm.tab, aes(x = micro.macro, fill = topo)) +
  geom_bar(position = "stack") +
  labs(
    title = "Topologies by Chromosome",
    x = "Chromosome",
    y = "Count",
    fill = "Topology"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors)

a<-ggplot(tree.fdm.tab, aes(x = micro.macro, fill = topo)) +
  geom_bar(position = "fill") +
  labs(
    title = "Topologies by Chromosome",
    x = "Chromosome",
    y = "Percentage",
    fill = "Topology"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors) +
  scale_y_continuous(labels = scales::percent)


#get percentages of each topology on autosomes
sum(tree.fdm.tab$topo=="T1" & tree.fdm.tab$chromosome!="ChrZ")/sum(tree.fdm.tab$chromosome!="ChrZ")
sum(tree.fdm.tab$topo=="T2" & tree.fdm.tab$chromosome!="ChrZ")/sum(tree.fdm.tab$chromosome!="ChrZ")
sum(tree.fdm.tab$topo=="T3" & tree.fdm.tab$chromosome!="ChrZ")/sum(tree.fdm.tab$chromosome!="ChrZ")
sum(tree.fdm.tab$topo=="T4" & tree.fdm.tab$chromosome!="ChrZ")/sum(tree.fdm.tab$chromosome!="ChrZ")
sum(tree.fdm.tab$topo=="T5" & tree.fdm.tab$chromosome!="ChrZ")/sum(tree.fdm.tab$chromosome!="ChrZ")

#get percentages of each topology on Z
sum(tree.fdm.tab$topo=="T1" & tree.fdm.tab$chromosome=="ChrZ")/sum(tree.fdm.tab$chromosome=="ChrZ")
sum(tree.fdm.tab$topo=="T2" & tree.fdm.tab$chromosome=="ChrZ")/sum(tree.fdm.tab$chromosome=="ChrZ")
sum(tree.fdm.tab$topo=="T3" & tree.fdm.tab$chromosome=="ChrZ")/sum(tree.fdm.tab$chromosome=="ChrZ")
sum(tree.fdm.tab$topo=="T4" & tree.fdm.tab$chromosome=="ChrZ")/sum(tree.fdm.tab$chromosome=="ChrZ")
sum(tree.fdm.tab$topo=="T5" & tree.fdm.tab$chromosome=="ChrZ")/sum(tree.fdm.tab$chromosome=="ChrZ")

sum(tree.fdm.tab$topo=="T1" )/length(tree.fdm.tab$chromosome)
sum(tree.fdm.tab$topo=="T2" )/length(tree.fdm.tab$chromosome)
sum(tree.fdm.tab$topo=="T3" )/length(tree.fdm.tab$chromosome)
sum(tree.fdm.tab$topo=="T4" )/length(tree.fdm.tab$chromosome)
sum(tree.fdm.tab$topo=="T5" )/length(tree.fdm.tab$chromosome)


# Create bar plot with topo on x-axis showing percentages for whole genome
b<-ggplot(tree.fdm.tab, aes(x = topo, fill = topo)) +
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count)) * 100)) +
  labs(
    title = "Percentage Distribution of Topo Values",
    x = "Topo",
    y = "Percentage (%)",
    fill = "Topo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

# Figure 2

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
                 "#808080", "#33A02C", "#B2DF8A", "#FF7F00", "#FDBF6F", 
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
chromosome_plot <- create_chromosome_plot(tree.fdm.tab)
print(chromosome_plot)

# Save Figure 2
# pdf(file = "windowed-gene-trees/200kb_trees_chromosomes.pdf", width = 9, height = 6.5)
# print(chromosome_plot)
# dev.off()

###this next section focuses on Figure3 and relevent statistical analyses, first with 
###results assuming T2, then with results assuming T1
###analyses looking including fdm and topology use abba-baba-results from 200kb windows
###all other fdm results are based on 100kb windows

#######plot topo~fdm ASSUMING T2
tree.fdm.tab.T4.rm<-tree.fdm.tab[tree.fdm.tab$topo!="T4",]

ggplot(barrier.loci, aes(x = topo, fill = topo)) +
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count)) * 100)) +
  labs(
    title = "Percentage Distribution of Topo Values",
    x = "Topo",
    y = "Percentage (%)",
    fill = "Topo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

# Create binary variable for topo=="T3"
tree.fdm.tab.T4.rm$is_T3 <- as.numeric(tree.fdm.tab.T4.rm$topo == "T3")

# Fit logistic regression model
model <- glm(is_T3 ~ fdM.val, data = tree.fdm.tab.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T4.rm)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT3<-ggplot(tree.fdm.tab.T4.rm, aes(x = fdM.val, y = is_T3)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T3",
       title = "Logistic Regression: fdM vs T3",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT3

# Create binary variable for topo=="T1"
tree.fdm.tab.T4.rm$is_T1 <- as.numeric(tree.fdm.tab.T4.rm$topo == "T1")

# Fit logistic regression model
model <- glm(is_T1 ~ fdM.val, data = tree.fdm.tab.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T4.rm)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT1<-ggplot(tree.fdm.tab.T4.rm, aes(x = fdM.val, y = is_T1)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T1",
       title = "Logistic Regression: fdM vs T1",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT1

# Create binary variable for topo=="T2"
tree.fdm.tab.T4.rm$is_T2 <- as.numeric(tree.fdm.tab.T4.rm$topo == "T2")

# Fit logistic regression model
model <- glm(is_T2 ~ fdM.val, data = tree.fdm.tab.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T4.rm)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT2<-ggplot(tree.fdm.tab.T4.rm, aes(x = fdM.val, y = is_T2)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T2",
       title = "Logistic Regression: fdM vs T2",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT2

###boxplots examining differences in fdm~topo
kwt<-kruskal.test(fdM.val~topo,data=tree.fdm.tab.T4.rm)
dunnTest(as.numeric(fdM.val)~topo,data=tree.fdm.tab.T4.rm)

m.topo<-ggplot(tree.fdm.tab.T4.rm[tree.fdm.tab.T4.rm$topo!="Tx",],aes(x=topo, y=fdM.val, group=topo))+
  ggtitle("Introgression ~ topology", subtitle = paste("K-W X² = ",round(kwt$statistic[[1]], digits = 2), ", df = ", kwt$parameter[[1]],", p ", ifelse(kwt$p.value < 0.0001, "< 0.0001", round(kwt$p.value, 3)), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Gene tree topology") +
  ylab("fdM") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

m.topo

#what is distribution of topos after removing introgressed regions?
fdm.tab2<-tree.fdm.tab.T4.rm[tree.fdm.tab.T4.rm$fdM<0.02,]

# Create bar plot with topo on x-axis showing percentages
ggplot(fdm.tab2, aes(x = topo, fill = topo)) +
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count)) * 100)) +
  labs(
    title = "Distribution of Topologies (introgressed regions removed)",
    x = "Topo",
    y = "Percentage (%)",
    fill = "Topo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))


#########supplemental figure##########
#######showing logistic regressions###
########assuming T2 as the st#########

#pdf(file = "logistic_regression.T2.pdf", width = 12, height = 4)
ggarrange(logisticT1,logisticT2,logisticT3, nrow = 1, ncol = 3)
#dev.off()

####100kb windows

intr.100kb$scaffold<-factor(intr.100kb$scaffold,levels = 
                         c("1","2","3","4","5","6","7",
                           "8","9","10","11","12","13","14",
                           "15","16","17","18","19","20","21",
                           "22","23","24","25","26","27","28",
                           "29","30","31","32","33","34","35",
                           "36","37","38","39","40"))

data_clean <- intr.100kb %>%
  filter(!is.na(fdM) & !is.na(mid)) %>%
  mutate(scaffold = as.character(scaffold))

# Function to extract chromosome number/letter for proper sorting
extract_chr <- function(scaffold) {
  # Extract the chromosome part (number or letter after "Chr")
  chr_match <- regexpr("Chr(\\d+|[A-Z])", scaffold, perl = TRUE)
  if (chr_match > 0) {
    chr_part <- regmatches(scaffold, chr_match)
    chr_value <- gsub("Chr", "", chr_part)
    return(chr_value)
  }
  return(scaffold)
}

# Add chromosome extraction and create sorting key
data_clean <- intr.100kb %>%
  mutate(
    chr_extract = sapply(scaffold, extract_chr),
    # Create numeric version for sorting (letters get high numbers)
    chr_numeric = ifelse(
      grepl("^\\d+$", chr_extract),
      as.numeric(chr_extract),
      1000 + utf8ToInt(chr_extract)  # Letters get values > 1000
    )
  ) %>%
  # Sort by chromosome then by mid position
  arrange(chr_numeric, mid) %>%
  # Create sequential x-position for plotting
  mutate(x_position = row_number())

# Get scaffold labels for x-axis (sample for readability)
n_scaffolds <- length(unique(data_clean$scaffold))
label_step <- max(1, floor(n_scaffolds / 20))  # Show ~20 labels max

scaffold_labels <- data_clean %>%
  group_by(scaffold) %>%
  summarise(x_pos = first(x_position), .groups = 'drop') %>%
  slice(seq(1, n(), by = label_step))

# Create alternating colors for chromosomes
data_clean <- data_clean %>%
  mutate(
    # Create unique chromosome identifiers
    chr_id = dense_rank(chr_numeric),
    # Alternate colors for every other chromosome
    chr_color = ifelse(chr_id %% 2 == 1, "Color1", "Color2")
  )

min(data_clean$fdM[data_clean$scaffold=="40"])

# Create the scatter plot with alternating colors (Top of Fig 4)
# introgression across the genome
p.T2 <- ggplot(data_clean, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T2",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
p.T2

#Show chromosome Z with a smoothed spline: the misalinged region shows a dubious pattern of high positive fdm
data_clean_z<-data_clean[data_clean$scaffold==40,]
data_clean_z$scaffold<-factor(data_clean_z$scaffold, levels = "40")
p2.T2 <- ggplot(data_clean_z, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method="gam", formula = y ~ s(x, k=15),
              linewidth=1, se=F, color="darkred")+
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T1",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(p2.T2)


#Kruskal Wallis and Dunn tests for significance
kwt<-kruskal.test(fdM.val~chromosome,data=tree.fdm.tab.T4.rm)
kwt1<-kruskal.test(fdM.val~micro.macro,data=tree.fdm.tab.T4.rm)
dunnTest(as.numeric(fdM.val)~micro.macro,data=tree.fdm.tab.T4.rm)

#More boxplots
m<-ggplot(tree.fdm.tab.T4.rm,aes(x=chromosome, y=fdM.val, group=chromosome))+
  ggtitle("Introgression ~ chromosome", subtitle = paste("X²=",round(kwt$statistic[[1]], digits = 2), ", df=", kwt$parameter[[1]],", p=", round(kwt$p.value, digits = 2), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Chromosome") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

m

m1<-ggplot(tree.fdm.tab.T4.rm,aes(x=micro.macro, y=fdM.val, group=micro.macro))+
  ggtitle("Introgression ~ chromosome type", subtitle = paste("X²=",round(kwt1$statistic[[1]], digits = 2), ", df=", kwt1$parameter[[1]],", p=", round(kwt1$p.value, digits = 2), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Chromosome type") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))
m1

####Get new table to look at fdM X chrom. length and position on chrom
names(tree.fdm.tab.T4.rm)
counts=0
fdm.val<-c()
chrom.log.length<-c()
chrom.length<-c()
chrom.number<-c()
abba.count<-c()
abba.prop<-c()
position<-c()
chr.ctr<-c()
per.dist.fr.ctr<-c()
micro.macro<-c()
tree.fdm.tab.T4.rm<-na.omit(tree.fdm.tab.T4.rm)
tree.fdm.tab.T4.rm$chr<-as.numeric(tree.fdm.tab.T4.rm$chr)
tree.fdm.tab.T4.rm$end<-as.numeric(tree.fdm.tab.T4.rm$end)

for (i in 1:length(tree.fdm.tab.T4.rm[,1])){
  counts=counts+1
  abba.count[counts]=tree.fdm.tab.T4.rm$ABBA[i]
  abba.prop[counts]=tree.fdm.tab.T4.rm$ABBA[i]/100000
  fdm.val[counts]=tree.fdm.tab.T4.rm$fdM[i]
  chrom.log.length[counts]<-log10(max(tree.fdm.tab.T4.rm$end[tree.fdm.tab.T4.rm$chr==tree.fdm.tab.T4.rm$chr[i]])) #log length
  chrom.length[counts]<-max(tree.fdm.tab.T4.rm$end[tree.fdm.tab.T4.rm$chr==tree.fdm.tab.T4.rm$chr[i]]) #length
  chrom.number[counts]<-tree.fdm.tab.T4.rm$chr[i]
  position[counts]<-tree.fdm.tab.T4.rm$mid[i]
  chr.ctr[counts]<-chrom.length[counts]/2
  per.dist.fr.ctr[counts]<-abs(position[counts]-chr.ctr[counts])/chr.ctr[counts]
  micro.macro[counts]<-tree.fdm.tab.T4.rm$micro.macro[i]
}

fdm.chr.tab<-na.omit(data.frame(cbind(chrom.length,chrom.log.length,chrom.number,position,abba.count,abba.prop,fdm.val,per.dist.fr.ctr,micro.macro)))
fdm.chr.tab$fdm.val<-as.numeric(fdm.chr.tab$fdm.val)
fdm.chr.tab$per.dist.fr.ctr<-as.numeric(fdm.chr.tab$per.dist.fr.ctr)
fdm.chr.tab$abba.prop<-as.numeric(fdm.chr.tab$abba.prop)
fdm.chr.tab$chrom.log.length<-as.numeric(fdm.chr.tab$chrom.log.length)
fdm.chr.tab$position<-as.numeric(fdm.chr.tab$position)
fdm.chr.tab$chrom.length<-as.numeric(fdm.chr.tab$chrom.length)

unique(fdm.chr.tab$chrom.number[fdm.chr.tab$chrom.length>25000000])

fdm.chr.tab <- fdm.chr.tab %>%
  group_by(chrom.number) %>%
  mutate(position_normalized = (position - min(position, na.rm = TRUE)) / (max(position, na.rm = TRUE) - min(position, na.rm = TRUE))) %>%
  ungroup()

quadratic_model <- lm(fdm.val ~ position_normalized + I(position_normalized^2), data=fdm.chr.tab[fdm.chr.tab$chrom.length>25000000,])

# Get model statistics
model_summary <- summary(quadratic_model)
r_squared <- model_summary$r.squared
p_val_linear <- coef(model_summary)[2,4]  # p-value for linear term
p_val_quad <- coef(model_summary)[3,4]    # p-value for quadratic term
linear_coef <- coef(model_summary)[1,1]   # intercept
quad_coef <- coef(model_summary)[3,1]     # quadratic coefficient

# Create the plot (Fig 4E)
quadratic_plot <- ggplot(fdm.chr.tab[fdm.chr.tab$chrom.length>25000000,], aes(x=position_normalized, y=fdm.val)) +
  ggtitle("Introgression ~ position on macrochromosome", 
          subtitle = paste("R² =", round(r_squared, digits = 3), " coef=", round(quad_coef, digits = 2)," p=", round(p_val_quad, digits = 2), sep = " ")) +
  geom_point(alpha=0.3, size=0.8) +
  stat_smooth(method = "lm", 
              formula = y ~ x + I(x^2),
              colour="red", se=TRUE, size=1.2) +
  theme_bw() +
  xlab("Normalized Window Position (0-1)") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

# Plot fdM~normalized window position macrochromosomes
print(quadratic_plot)

# Print model coefficients and equation
print("Model Coefficients:")
print(coef(model_summary))
print(paste("Equation: fdm.val =", round(coef(quadratic_model)[1], digits = 4), 
            "+", round(coef(quadratic_model)[2], digits = 4), "* mid_normalized +", 
            round(coef(quadratic_model)[3], digits = 4), "* mid_normalized^2"))





mean.fdm<-c()
chrom.length<-c()
chrom.log.length<-c()
abba.prop<-c()

for(i in 1:40){
  mean.fdm[i]<-mean(fdm.chr.tab$fdm.val[fdm.chr.tab$chrom.number==i])
  chrom.length[i]<-max(fdm.chr.tab$chrom.length[fdm.chr.tab$chrom.number==i])
  chrom.log.length[i]<-max(fdm.chr.tab$chrom.log.length[fdm.chr.tab$chrom.number==i])
  abba.prop[i]<-max(fdm.chr.tab$abba.prop[fdm.chr.tab$chrom.number==i])
}

tab3<-na.omit(data.frame(cbind(mean.fdm,chrom.length,chrom.log.length, abba.prop)))
tab3$mean.fdm<-as.numeric(tab3$mean.fdm)
tab3$chrom.log.length<-as.numeric(tab3$chrom.log.length)
tab3$abba<-as.numeric(tab3$abba.prop)

model <- lm(mean.fdm ~ chrom.log.length, data=tab3)
summary(model)
pval<-coef(summary(model))[2,4]
coef<-coef(summary(model))[2,1]
r_squared <- summary(model)$r.squared

#fdm~chromosome length
length.plot<-ggplot(tab3,aes(x=chrom.log.length,y=mean.fdm)) +
  geom_point(alpha=0.5) +
  ggtitle("Introgression ~ chromosome length",subtitle = paste("R²=",round(r_squared, digits = 2)," coef=", round(coef, digits = 2), ", p=", round(pval,digits = 2), sep = ""))+
  stat_smooth(colour="red", method="glm",se=TRUE)+
  theme_bw() +
  xlab("log10(chromosome length)") +
  ylab("Mean FDM") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

length.plot

#####################################
#####Quadratic Chromosome Plots######
#####################################

tab1<-intr.100kb
names(tab1)

counts=0
fdm.val<-c()
chrom.log.length<-c()
chrom.length<-c()
chrom.number<-c()
abba.count<-c()
abba.prop<-c()
position<-c()
chr.ctr<-c()
per.dist.fr.ctr<-c()
mid<-c()

for (i in 1:length(tab1$scaffold)){
  counts=counts+1
  abba.count[counts]=tab1$ABBA[i]
  abba.prop[counts]=tab1$ABBA[i]/100000
  fdm.val[counts]=tab1$fdM[i]
  chrom.log.length[counts]<-log10(max(tab1$end[tab1$scaffold==tab1$scaffold[i]])) #log length
  chrom.length[counts]<-max(tab1$end[tab1$scaffold==tab1$scaffold[i]]) #length
  chrom.number[counts]<-tab1$scaffold[i]
  position[counts]<-tab1$mid[i]
  chr.ctr[counts]<-chrom.length[counts]/2
  per.dist.fr.ctr[counts]<-abs(position[counts]-chr.ctr[counts])/chr.ctr[counts]
  mid[counts]<-tab1$mid[i]
}

tab<-na.omit(data.frame(cbind(chrom.length,chrom.log.length,mid,chrom.number,abba.count,abba.prop,fdm.val,per.dist.fr.ctr,micro.macro)))
tab$fdm.val<-as.numeric(tab$fdm.val)
tab$per.dist.fr.ctr<-as.numeric(tab$per.dist.fr.ctr)
tab$abba.prop<-as.numeric(tab$abba.prop)
tab$chrom.log.length<-as.numeric(tab$chrom.log.length)
tab$mid<-as.numeric(tab$mid)
tab.pos<-tab#[tab$fdm.val<=0,]

tab.pos <- tab.pos %>%
  group_by(chrom.number) %>%
  mutate(mid_normalized = (mid - min(mid, na.rm = TRUE)) / (max(mid, na.rm = TRUE) - min(mid, na.rm = TRUE))) %>%
  ungroup()

# Print model coefficients and equation
print("Model Coefficients:")
print(coef(model_summary))
print(paste("Equation: fdm.val =", round(coef(quadratic_model)[1], digits = 4), 
            "+", round(coef(quadratic_model)[2], digits = 4), "* mid_normalized +", 
            round(coef(quadratic_model)[3], digits = 4), "* mid_normalized^2"))



####Linear model
for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$baba.prop<-as.numeric(tab2$abba.prop)
  model <- lm(fdm.val ~ per.dist.fr.ctr, data=tab2)
  #summary(model)
  pval<-coef(summary(model))[2,4]
  slope<-coef(summary(model))[2,1]
  nam<-paste("p",i,sep="")
  assign(nam,ggplot(tab2,aes(x=per.dist.fr.ctr,y=fdm.val)) + ggtitle(paste("Chr",i,sep=""),subtitle = paste("p=", round(pval,digits = 2),sep = ""))+
           geom_point(alpha=0.5) +
           stat_smooth(colour="red", method="glm",se=TRUE)+
           theme_bw() + xlab("% dist. from center") + ylab("fdm"))
}



for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$abba.prop<-as.numeric(tab2$abba.prop)
  
  
  # Skip if not enough data points
  if(nrow(tab2) < 10) {
    cat("Skipping chromosome", i, "- insufficient data\n")
    next
  }
  
  # Try multiple exponential approaches
  model_fitted <- FALSE
  
  # Quadratic
  if(!model_fitted) {
    tryCatch({
      model <- lm(fdm.val ~ mid + I(mid^2), data=tab2)
      pval <- coef(summary(model))[3,4]  # p-value for quadratic term
      r_squared <- summary(model)$r.squared
      quad_coef <- coef(model)[3]
      
      nam<-paste("r",i,sep="")
      assign(nam, ggplot(tab2, aes(x=mid, y=fdm.val)) + 
               ggtitle(paste("Chromosome ",i,sep=""), 
                       subtitle = paste("R² = ", round(r_squared, digits = 3), 
                                        ", P= ", round(pval, digits = 3), sep = "")) +
               geom_point(alpha=0.5) +
               stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                           colour="red", se=TRUE) +
               theme_bw() + xlab("window position") + ylab("fdm"))
      
      model_fitted <- TRUE
      cat("Quadratic fitted for chromosome", i, "\n")
      
    }, error = function(e) {})
  }
  
}


#linear models of fdm~%distance from chromosome center (not in paper)
# pdf(file = "../figs/fdmXdist.from.ctr.T2.pdf", width = 13, height = 15)
ggarrange(p1,p2, p3,p4, 
          p5, p6, p7, p8,
          p9,p10,p11,p12,
          p13,p14,p15,p16,p17,p18,p19,p20,p21,
          p22,p23,p24,p25,p26,p27,p28,p29,p40,
          ncol = 5, nrow = 6)

# dev.off()

#quadratic models of fdm~position on chromosome
# pdf(file = "../figs/quadratic_fdmXmid.T2.pdf", width = 13, height = 9)
ggarrange(r1,r2, r3,r4, 
          r5, r6, r7, r8,
          r9,r10,r11,r12,
          r13,r14,r15,r16,r17,r18,r19,r20,
          r22,r23,r24,r25,r26,r27,r28,r29,r40,
          ncol = 5, nrow = 6)

# dev.off()

##########################################################
######################Repeat for T1#######################
##########################################################

intr.100kb.T1 <- read_csv("abba-baba-windows/ABBABABAwindows.w100k.T1-cinp3.csv")
intr.200kb.T1 <- read_csv("abba-baba-windows/ABBABABAwindows.w200k.T1-cinp3.csv")

intr.100kb.T1$scaffold<-factor(intr.100kb.T1$scaffold,levels = 
                              c("1","2","3","4","5","6","7",
                                "8","9","10","11","12","13","14",
                                "15","16","17","18","19","20","21",
                                "22","23","24","25","26","27","28",
                                "29","30","31","32","33","34","35",
                                "36","37","38","39","40"))


intr.100kb.T1.err.rm <- intr.100kb.T1 |> semi_join(intr.100kb, by = c("scaffold", "end"))
length(intr.100kb.T1$scaffold) #get length of original table
length(intr.100kb$scaffold) #new table is shorter
length(intr.100kb.T1.err.rm$scaffold) #new table matches expected number of rows after T4 windows removed
chr<-c()
fdM.val<-c()
st<-c()
en<-c()
micro.macro<-c()
mid<-c()
ABBA<-c()

for (i in 1:length(tree.table$chromosome)){
  chr[i]<-str_split(tree.table$chromosome[i],"Chr")[[1]][2]
  if(chr[i]=="Z"){
    chr[i]<-40
  }
  st[i]<-tree.table$start_pos[i]
  en[i]<-tree.table$end_pos[i]
  fdM.val[i]<-intr.200kb.T1$fdM[intr.200kb.T1$start==st[i] & intr.200kb.T1$end==en[i] & intr.200kb.T1$scaffold==chr[i]]
  ABBA[i]<-intr.200kb.T1$ABBA[intr.200kb.T1$start==st[i] & intr.200kb.T1$end==en[i] & intr.200kb.T1$scaffold==chr[i]]
  mid[i]<-intr.200kb.T1$mid[intr.200kb.T1$start==st[i] & intr.200kb.T1$end==en[i] & intr.200kb.T1$scaffold==chr[i]]
  micro.macro[i]<-"micro"
  if(as.numeric(chr[i])<=8){
    micro.macro[i]<-"macro"
  }
  if(chr[i]=="40"){
    micro.macro[i]<-"Z-chromosome"
  }
}

tree.fdm.tab.T1<-data.frame(cbind(tree.table,chr,mid,ABBA,fdM.val,micro.macro))

#######plot top~fdm
tree.fdm.tab.T1.T4.rm<-tree.fdm.tab.T1[tree.fdm.tab.T1$topo!="T4",]

# Create binary variable for topo=="T3"
tree.fdm.tab.T1.T4.rm$is_T3 <- as.numeric(tree.fdm.tab.T1.T4.rm$topo == "T3")

# Fit logistic regression model
model <- glm(is_T3 ~ fdM.val, data = tree.fdm.tab.T1.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T1.T4.rm)
aic_value <- AIC(model)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT3.T1<-ggplot(tree.fdm.tab.T1.T4.rm, aes(x = fdM.val, y = is_T3)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T3",
       title = "Logistic Regression: fdM vs T3",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT3.T1

# Create binary variable for topo=="T1"
tree.fdm.tab.T1.T4.rm$is_T1 <- as.numeric(tree.fdm.tab.T1.T4.rm$topo == "T1")

# Fit logistic regression model
model <- glm(is_T1 ~ fdM.val, data = tree.fdm.tab.T1.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T1.T4.rm)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT1.T1<-ggplot(tree.fdm.tab.T1.T4.rm, aes(x = fdM.val, y = is_T1)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T1",
       title = "Logistic Regression: fdM vs T1",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT1.T1

# Create binary variable for topo=="T2"
tree.fdm.tab.T1.T4.rm$is_T2 <- as.numeric(tree.fdm.tab.T1.T4.rm$topo == "T2")

# Fit logistic regression model
model <- glm(is_T2 ~ fdM.val, data = tree.fdm.tab.T1.T4.rm, family = binomial)

# Print model summary
summary(model)

# Calculate pseudo R-squared (McFadden's)
null_deviance <- model$null.deviance
residual_deviance <- model$deviance
pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

# Get model statistics
n_obs <- nrow(tree.fdm.tab.T4.rm)
p_value <- summary(model)$coefficients[2, 4]  # p-value for fdM coefficient
odds_ratio<-exp(coef(model)[2])
# Create subtitle with key statistics
subtitle_text <- paste0("Pseudo R² = ", round(pseudo_r_squared, 3),
                        ", n = ", n_obs,
                        #", AIC = ", round(aic_value, 1),
                        ", p ", ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 3)))

# Create predictions for smooth curve
fdM_range <- seq(min(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 max(tree.fdm.tab.T1.T4.rm$fdM.val, na.rm = TRUE), 
                 length.out = 100)
predictions <- predict(model, 
                       newdata = data.frame(fdM.val = fdM_range), 
                       type = "response")

# Create the plot
logisticT2.T1<-ggplot(tree.fdm.tab.T1.T4.rm, aes(x = fdM.val, y = is_T2)) +
  # Add raw data points with some jitter to avoid overplotting
  geom_point(alpha = 0.6, 
             position = position_jitter(height = 0.02, width = 0)) +
  # Add logistic regression curve
  geom_line(data = data.frame(fdM = fdM_range, predicted = predictions),
            aes(x = fdM, y = predicted), 
            color = "red", size = 1.2) +
  # Labels and theme
  labs(x = "fdM",
       y = "Probability of topo = T2",
       title = "Logistic Regression: fdM vs T2",
       subtitle = subtitle_text) +
  theme_bw() +
  # Set y-axis limits to show full probability range
  ylim(0, 1) +
  # Add horizontal reference lines
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dashed", alpha = 0.3)

logisticT2.T1

#supplemental figure showing logistic regression assuming T1
pdf(file = "logistic_regression.T1.pdf", width = 12, height = 4)
ggarrange(logisticT1.T1,logisticT2.T1,logisticT3.T1, nrow = 1, ncol = 3)
dev.off()

kwt<-kruskal.test(fdM.val~topo,data=tree.fdm.tab.T1.T4.rm)
dunnTest(as.numeric(fdM.val)~topo,data=tree.fdm.tab.T1.T4.rm)

m.topo.T1<-ggplot(tree.fdm.tab.T1.T4.rm[tree.fdm.tab.T1.T4.rm$topo!="Tx",],aes(x=topo, y=fdM.val, group=topo))+
  ggtitle("Introgression ~ topology", subtitle = paste("K-W X² = ",round(kwt$statistic[[1]], digits = 2), ", df = ", kwt$parameter[[1]],", p ", ifelse(kwt$p.value < 0.0001, "< 0.0001", round(kwt$p.value, 3)), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Gene tree topology") +
  ylab("fdM") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

m.topo.T1


fdm.tab2.T1<-tree.fdm.tab.T1.T4.rm[tree.fdm.tab.T1.T4.rm$fdM<0.02,]

# Create bar plot with topo on x-axis showing percentages
ggplot(fdm.tab2.T1, aes(x = topo, fill = topo)) +
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count)) * 100)) +
  labs(
    title = "Distribution of Topologies (introgressed regions removed)",
    x = "Topo",
    y = "Percentage (%)",
    fill = "Topo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = topo_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

####100kb windows

intr.100kb.T1.err.rm$scaffold<-factor(intr.100kb.T1.err.rm$scaffold,levels = 
                              c("1","2","3","4","5","6","7",
                                "8","9","10","11","12","13","14",
                                "15","16","17","18","19","20","21",
                                "22","23","24","25","26","27","28",
                                "29","30","31","32","33","34","35",
                                "36","37","38","39","40"))

data_clean <- intr.100kb.T1.err.rm %>%
  filter(!is.na(fdM) & !is.na(mid)) %>%
  mutate(scaffold = as.character(scaffold))

# Function to extract chromosome number/letter for proper sorting
extract_chr <- function(scaffold) {
  # Extract the chromosome part (number or letter after "Chr")
  chr_match <- regexpr("Chr(\\d+|[A-Z])", scaffold, perl = TRUE)
  if (chr_match > 0) {
    chr_part <- regmatches(scaffold, chr_match)
    chr_value <- gsub("Chr", "", chr_part)
    return(chr_value)
  }
  return(scaffold)
}

# Add chromosome extraction and create sorting key
data_clean <- intr.100kb.T1.err.rm %>%
  mutate(
    chr_extract = sapply(scaffold, extract_chr),
    # Create numeric version for sorting (letters get high numbers)
    chr_numeric = ifelse(
      grepl("^\\d+$", chr_extract),
      as.numeric(chr_extract),
      1000 + utf8ToInt(chr_extract)  # Letters get values > 1000
    )
  ) %>%
  # Sort by chromosome then by mid position
  arrange(chr_numeric, mid) %>%
  # Create sequential x-position for plotting
  mutate(x_position = row_number())

# Get scaffold labels for x-axis (sample for readability)
n_scaffolds <- length(unique(data_clean$scaffold))
label_step <- max(1, floor(n_scaffolds / 20))  # Show ~20 labels max

scaffold_labels <- data_clean %>%
  group_by(scaffold) %>%
  summarise(x_pos = first(x_position), .groups = 'drop') %>%
  slice(seq(1, n(), by = label_step))

# Create alternating colors for chromosomes
data_clean <- data_clean %>%
  mutate(
    # Create unique chromosome identifiers
    chr_id = dense_rank(chr_numeric),
    # Alternate colors for every other chromosome
    chr_color = ifelse(chr_id %% 2 == 1, "Color1", "Color2")
  )

min(data_clean$fdM[data_clean$scaffold=="40"])

# Create the scatter plot with alternating colors (Top of Fig 4)
p.T1 <- ggplot(data_clean, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T2",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
p.T1

#Show chromosome Z with a smoothed spline: the misalinged region shows a dubious pattern of high positive fdm
data_clean_z<-data_clean[data_clean$scaffold==40,]
data_clean_z$scaffold<-factor(data_clean_z$scaffold, levels = "40")
p2.T1 <- ggplot(data_clean_z, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method="gam", formula = y ~ s(x, k=15),
              linewidth=1, se=F, color="darkred")+
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T1",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(p2.T1)


#Kruskal Wallis and Dunn tests for significance
kwt<-kruskal.test(fdM.val~chromosome,data=tree.fdm.tab.T1.T4.rm)
kwt1<-kruskal.test(fdM.val~micro.macro,data=tree.fdm.tab.T1.T4.rm)
dunnTest(as.numeric(fdM.val)~micro.macro,data=tree.fdm.tab.T1.T4.rm)

#boxplots
m2<-ggplot(tree.fdm.tab.T1.T4.rm,aes(x=chromosome, y=fdM.val, group=chromosome))+
  ggtitle("Introgression ~ chromosome", subtitle = paste("X²=",round(kwt$statistic[[1]], digits = 2), ", df=", kwt$parameter[[1]],", p=", round(kwt$p.value, digits = 2), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Chromosome") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

m2

m3<-ggplot(tree.fdm.tab.T1.T4.rm,aes(x=micro.macro, y=fdM.val, group=micro.macro))+
  ggtitle("Introgression ~ chromosome type", subtitle = paste("X²=",round(kwt1$statistic[[1]], digits = 2), ", df=", kwt1$parameter[[1]],", p=", round(kwt1$p.value, digits = 2), sep = ""))+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +
  xlab("Chromosome type") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))
m3

####Get new table to look at fdM X chrom. length and position on chrom
names(tree.fdm.tab.T1.T4.rm)
counts=0
fdm.val<-c()
chrom.log.length<-c()
chrom.length<-c()
chrom.number<-c()
abba.count<-c()
abba.prop<-c()
position<-c()
chr.ctr<-c()
per.dist.fr.ctr<-c()
micro.macro<-c()
tree.fdm.tab.T1.T4.rm<-na.omit(tree.fdm.tab.T1.T4.rm)
tree.fdm.tab.T1.T4.rm$chr<-as.numeric(tree.fdm.tab.T1.T4.rm$chr)
tree.fdm.tab.T1.T4.rm$end<-as.numeric(tree.fdm.tab.T1.T4.rm$end)

for (i in 1:length(tree.fdm.tab.T1.T4.rm[,1])){
  counts=counts+1
  abba.count[counts]=tree.fdm.tab.T1.T4.rm$ABBA[i]
  abba.prop[counts]=tree.fdm.tab.T1.T4.rm$ABBA[i]/100000
  fdm.val[counts]=tree.fdm.tab.T1.T4.rm$fdM[i]
  chrom.log.length[counts]<-log10(max(tree.fdm.tab.T1.T4.rm$end[tree.fdm.tab.T1.T4.rm$chr==tree.fdm.tab.T1.T4.rm$chr[i]])) #log length
  chrom.length[counts]<-max(tree.fdm.tab.T1.T4.rm$end[tree.fdm.tab.T1.T4.rm$chr==tree.fdm.tab.T1.T4.rm$chr[i]]) #length
  chrom.number[counts]<-tree.fdm.tab.T1.T4.rm$chr[i]
  position[counts]<-tree.fdm.tab.T1.T4.rm$mid[i]
  chr.ctr[counts]<-chrom.length[counts]/2
  per.dist.fr.ctr[counts]<-abs(position[counts]-chr.ctr[counts])/chr.ctr[counts]
  micro.macro[counts]<-tree.fdm.tab.T1.T4.rm$micro.macro[i]
}

fdm.chr.tab.T1<-na.omit(data.frame(cbind(chrom.length,chrom.log.length,chrom.number,position,abba.count,abba.prop,fdm.val,per.dist.fr.ctr,micro.macro)))
fdm.chr.tab.T1$fdm.val<-as.numeric(fdm.chr.tab.T1$fdm.val)
fdm.chr.tab.T1$per.dist.fr.ctr<-as.numeric(fdm.chr.tab.T1$per.dist.fr.ctr)
fdm.chr.tab.T1$abba.prop<-as.numeric(fdm.chr.tab.T1$abba.prop)
fdm.chr.tab.T1$chrom.log.length<-as.numeric(fdm.chr.tab.T1$chrom.log.length)
fdm.chr.tab.T1$position<-as.numeric(fdm.chr.tab.T1$position)
fdm.chr.tab.T1$chrom.length<-as.numeric(fdm.chr.tab.T1$chrom.length)

unique(fdm.chr.tab.T1$chrom.number[fdm.chr.tab.T1$chrom.length>25000000])

fdm.chr.tab.T1 <- fdm.chr.tab.T1 %>%
  group_by(chrom.number) %>%
  mutate(position_normalized = (position - min(position, na.rm = TRUE)) / (max(position, na.rm = TRUE) - min(position, na.rm = TRUE))) %>%
  ungroup()

quadratic_model <- lm(fdm.val ~ position_normalized + I(position_normalized^2), data=fdm.chr.tab.T1[fdm.chr.tab.T1$chrom.length>25000000,])

# Get model statistics
model_summary <- summary(quadratic_model)
r_squared <- model_summary$r.squared
p_val_linear <- coef(model_summary)[2,4]  # p-value for linear term
p_val_quad <- coef(model_summary)[3,4]    # p-value for quadratic term
linear_coef <- coef(model_summary)[1,1]   # intercept
quad_coef <- coef(model_summary)[3,1]     # quadratic coefficient

# Create the plot (Fig 4E)
quadratic_plot.T1 <- ggplot(fdm.chr.tab.T1[fdm.chr.tab.T1$chrom.length>25000000,], aes(x=position_normalized, y=fdm.val)) +
  ggtitle("Introgression ~ position on macrochromosome", 
          subtitle = paste("R² =", round(r_squared, digits = 3), " coef=", round(quad_coef, digits = 2)," p=", round(p_val_quad, digits = 2), sep = " ")) +
  geom_point(alpha=0.3, size=0.8) +
  stat_smooth(method = "lm", 
              formula = y ~ x + I(x^2),
              colour="red", se=TRUE, size=1.2) +
  theme_bw() +
  xlab("Normalized Window Position (0-1)") +
  ylab("FDM Value") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

# plot
print(quadratic_plot.T1)

# Print model coefficients and equation
print("Model Coefficients:")
print(coef(model_summary))
print(paste("Equation: fdm.val =", round(coef(quadratic_model)[1], digits = 4), 
            "+", round(coef(quadratic_model)[2], digits = 4), "* mid_normalized +", 
            round(coef(quadratic_model)[3], digits = 4), "* mid_normalized^2"))





mean.fdm<-c()
chrom.length<-c()
chrom.log.length<-c()
abba.prop<-c()

for(i in 1:40){
  mean.fdm[i]<-mean(fdm.chr.tab.T1$fdm.val[fdm.chr.tab.T1$chrom.number==i])
  chrom.length[i]<-max(fdm.chr.tab.T1$chrom.length[fdm.chr.tab.T1$chrom.number==i])
  chrom.log.length[i]<-max(fdm.chr.tab.T1$chrom.log.length[fdm.chr.tab.T1$chrom.number==i])
  abba.prop[i]<-max(fdm.chr.tab.T1$abba.prop[fdm.chr.tab.T1$chrom.number==i])
}

tab4<-na.omit(data.frame(cbind(mean.fdm,chrom.length,chrom.log.length, abba.prop)))
tab4$mean.fdm<-as.numeric(tab4$mean.fdm)
tab4$chrom.log.length<-as.numeric(tab4$chrom.log.length)
tab4$abba<-as.numeric(tab4$abba.prop)

model <- lm(mean.fdm ~ chrom.log.length, data=tab4)
summary(model)
pval<-coef(summary(model))[2,4]
coef<-coef(summary(model))[2,1]
r_squared <- summary(model)$r.squared

length.plot.T1<-ggplot(tab4,aes(x=chrom.log.length,y=mean.fdm)) +
  geom_point(alpha=0.5) +
  ggtitle("Introgression ~ chromosome length",subtitle = paste("R²=",round(r_squared, digits = 2)," coef=", round(coef, digits = 2), ", p=", round(pval,digits = 2), sep = ""))+
  stat_smooth(colour="red", method="glm",se=TRUE)+
  theme_bw() +
  xlab("log10(chromosome length)") +
  ylab("Mean FDM") +
  theme(plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12))

length.plot.T1

#####################################
#####Quadratic Chromosome Plots######
#####################################

tab1<-intr.100kb.T1.err.rm
names(tab1)

counts=0
fdm.val<-c()
chrom.log.length<-c()
chrom.length<-c()
chrom.number<-c()
abba.count<-c()
abba.prop<-c()
position<-c()
chr.ctr<-c()
per.dist.fr.ctr<-c()
mid<-c()

for (i in 1:length(tab1$scaffold)){
  counts=counts+1
  abba.count[counts]=tab1$ABBA[i]
  abba.prop[counts]=tab1$ABBA[i]/100000
  fdm.val[counts]=tab1$fdM[i]
  chrom.log.length[counts]<-log10(max(tab1$end[tab1$scaffold==tab1$scaffold[i]])) #log length
  chrom.length[counts]<-max(tab1$end[tab1$scaffold==tab1$scaffold[i]]) #length
  chrom.number[counts]<-tab1$scaffold[i]
  position[counts]<-tab1$mid[i]
  chr.ctr[counts]<-chrom.length[counts]/2
  per.dist.fr.ctr[counts]<-abs(position[counts]-chr.ctr[counts])/chr.ctr[counts]
  mid[counts]<-tab1$mid[i]
}

tab<-na.omit(data.frame(cbind(chrom.length,chrom.log.length,mid,chrom.number,abba.count,abba.prop,fdm.val,per.dist.fr.ctr,micro.macro)))
tab$fdm.val<-as.numeric(tab$fdm.val)
tab$per.dist.fr.ctr<-as.numeric(tab$per.dist.fr.ctr)
tab$abba.prop<-as.numeric(tab$abba.prop)
tab$chrom.log.length<-as.numeric(tab$chrom.log.length)
tab$mid<-as.numeric(tab$mid)
tab.pos<-tab#[tab$fdm.val<=0,]

tab.pos <- tab.pos %>%
  group_by(chrom.number) %>%
  mutate(mid_normalized = (mid - min(mid, na.rm = TRUE)) / (max(mid, na.rm = TRUE) - min(mid, na.rm = TRUE))) %>%
  ungroup()

# Print model coefficients and equation
print("Model Coefficients:")
print(coef(model_summary))
print(paste("Equation: fdm.val =", round(coef(quadratic_model)[1], digits = 4), 
            "+", round(coef(quadratic_model)[2], digits = 4), "* mid_normalized +", 
            round(coef(quadratic_model)[3], digits = 4), "* mid_normalized^2"))



####Linear model
for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$baba.prop<-as.numeric(tab2$abba.prop)
  model <- lm(fdm.val ~ per.dist.fr.ctr, data=tab2)
  #summary(model)
  pval<-coef(summary(model))[2,4]
  slope<-coef(summary(model))[2,1]
  nam<-paste("T1.p",i,sep="")
  assign(nam,ggplot(tab2,aes(x=per.dist.fr.ctr,y=fdm.val)) + ggtitle(paste("Chr",i,sep=""),subtitle = paste("p=", round(pval,digits = 2),sep = ""))+
           geom_point(alpha=0.5) +
           stat_smooth(colour="red", method="glm",se=TRUE)+
           theme_bw() + xlab("% dist. from center") + ylab("fdm"))
}



for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$abba.prop<-as.numeric(tab2$abba.prop)
  
  
  # Skip if not enough data points
  if(nrow(tab2) < 10) {
    cat("Skipping chromosome", i, "- insufficient data\n")
    next
  }
  
  # Try multiple exponential approaches
  model_fitted <- FALSE
  
  # Quadratic
  if(!model_fitted) {
    tryCatch({
      model <- lm(fdm.val ~ mid + I(mid^2), data=tab2)
      pval <- coef(summary(model))[3,4]  # p-value for quadratic term
      r_squared <- summary(model)$r.squared
      quad_coef <- coef(model)[3]
      
      nam<-paste("T1.r",i,sep="")
      assign(nam, ggplot(tab2, aes(x=mid, y=fdm.val)) + 
               ggtitle(paste("Chromosome ",i,sep=""), 
                       subtitle = paste("R² = ", round(r_squared, digits = 3), 
                                        ", P= ", round(pval, digits = 3), sep = "")) +
               geom_point(alpha=0.5) +
               stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                           colour="red", se=TRUE) +
               theme_bw() + xlab("window position") + ylab("fdm"))
      
      model_fitted <- TRUE
      cat("Quadratic fitted for chromosome", i, "\n")
      
    }, error = function(e) {})
  }
  
}


# Print model coefficients
print("Model Coefficients:")
print(coef(model_summary))

#linear model fdm~%distance from chromosome center (not in paper)
#pdf(file = "fdmXdist.from.ctr.T1.pdf", width = 13, height = 15)
ggarrange(T1.p1,T1.p2, T1.p3,T1.p4, 
          T1.p5, T1.p6, T1.p7, T1.p8,
          T1.p9,T1.p10,T1.p11,T1.p12,
          T1.p13,T1.p14,T1.p15,T1.p16,T1.p17,T1.p18,T1.p19,T1.p20,
          T1.p22,T1.p23,T1.p24,T1.p25,T1.p26,T1.p27,T1.p28,T1.p29,T1.p40,
          ncol = 5, nrow = 6)

#dev.off()

#quadratic models
pdf(file = "quadratic_fdmXmid.T1.pdf", width = 13, height = 9)
ggarrange(T1.r1,T1.r2, T1.r3,T1.r4, 
          T1.r5, T1.r6, T1.r7, T1.r8,
          T1.r9,T1.r10,T1.r11,T1.r12,
          T1.r13,T1.r14,T1.r15,T1.r16,T1.r17,T1.r18,T1.r19,T1.r20,
          T1.r22,T1.r23,T1.r24,T1.r25,T1.r26,T1.r27,T1.r28,T1.r29,T1.r40,
          ncol = 5, nrow = 6)

dev.off()



#simplified Figure 3 showing subset of comparisons
pdf(file = "200kb_broad_patterns_T1+T2.pdf", width = 13, height = 9)
ggarrange(m3,m.topo.T1,length.plot.T1,quadratic_plot.T1,logisticT3,
          m1,m.topo,length.plot,quadratic_plot, logisticT3.T1,
          ncol = 5, nrow = 2)
dev.off()



###Figure 3 (analyses assuming T2 as the species tree)
pdf(file = "../figs/fdm.patterns.pdf", width = 13, height = 6)
ggarrange(m1,m.topo,length.plot,logisticT3,
          quadratic_plot,r3, r7,r15,
          ncol = 4, nrow = 2)
dev.off()

###Figure 3A
pdf(file = "fdm.patterns.T1.pdf", width = 13, height = 6)
ggarrange(m3,m.topo.T1,length.plot.T1,logisticT2.T1,
          quadratic_plot.T1,T1.r3, T1.r7,T1.r15,
          ncol = 4, nrow = 2)
dev.off()




