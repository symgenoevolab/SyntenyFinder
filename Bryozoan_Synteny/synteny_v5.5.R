# R script for metazoan synteny analysis
# Isabel Jiah-Yih Liao, Tom Lewin, Yi-Jyun Luo
# Symbiosis Genomics & Evolution Lab 2024
# version 5.5 | 2024-01-16
#
# This script performs an extensive pairwise comparison of macrosynteny among selected species,
# and calculates both the genome-wide and specific chromosome pair synteny mixing rates using 
# Spearman's rank correlation coefficient (ρ).
# Key outputs include:
# - 'source_target_synteny.pdf': A PDF file showcasing in-depth synteny information.
# - 'source_target_oxford_plot.pdf': A PDF file featuring an Oxford dot plot, useful for
#                                    visualizing inversion and translocation events (scaled in Mb).
#
# Additional functionalities of the script:
# - (a) Generating a heatmap to visualize the genome-wide synteny mixing rate.
# - (b) Analyzing specific chromosome subsets within the Oxford dot plot.
# - (c) Providing detailed synteny mixing rate analyses for orthologous chromosome pairs.


# Setting the working directory
# The './input' directory should contain karyotype and coordinate files for each species.
# To set the working directory to the parent folder of 'input', use setwd('path_to_input').
# For example, if the full path to 'input' is 'path_to_input/input', use setwd('path_to_input').
setwd('/Users/isabel/Desktop/Bioinformatics/Synteny/bryozoans')

# Load required packages
library(RIdeogram)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ape)

# List of species
# mme, Membranipora membranacea (Darwin Tree of Life) (n=11)
# bst, Bugulina stolonifera (Darwin Tree of Life) (n=11)
# wsu, Watersipora subatra (Darwin Tree of Life) (n=11)
# cpa, Cryptosula pallasiana (Darwin Tree of Life) (n=12)
# cmu, Cristatella mucedo (NCBI GEO) (n=8)
species_list <- c("Mme", "Bst", "Wsu", "Cpa", "Cmu")

# Assign source and target species 
source <- "Mme"
target <- "Cmu"

# Initialize a matrix to store mixing rate results
mixing_rate_matrix <- matrix(NA, nrow = length(species_list), ncol = length(species_list),
                             dimnames = list(species_list, species_list))

# Initialize a datafram to store mixing rate results
chr_mixing_rate_data <- data.frame(Species1 = character(), Species2 = character(), OrthoChromosome = character(),
                                   MixingRate = character(), stringsAsFactors = FALSE)

# Loop over each pair of species
for (source in species_list) {
  for (target in species_list) {
    if (source != target) {  # Skip if source and target are the same
      
      # Setup filter for jumping genes (default = 5)
      jumping_genes = 5
      
      # If you would like to reorder the chromosomes, put the chromosome order here.
      # Otherwise, set reorder to FALSE
      reorder <- FALSE
      Mme_order <- c(1,4,2,3,7,8,6,9,10,5,11)
      Bst_order <- c(1,2,3,4,11,8,5,7,9,10,6)
      Wsu_order <- c(1,3,4,2,6,5,7,9,11,8,10)
      Cpa_order <- c(10,12,4,6,8,5,2,1,11,7,3,9)
      Cmu_order <- c(5,1,2,3,7,4,8,6)
      source_order <- get(paste0(source, "_order"))
      target_order <- get(paste0(target, "_order"))
      
      # List any chromosomes in either the source or the target genome which should be flipped
      Mme_flip <- c()
      Bst_flip <- c()
      Wsu_flip <- c(1,4,2,6,7,11,10)
      Cpa_flip <- c(10,12,2,11,7,3,9)
      Cmu_flip <- c()
      source_flip <- get(paste0(source, "_flip"))
      target_flip <- get(paste0(target, "_flip"))
      
      ########################################
      ### 1. Set karyotypes
      ########################################
      
      # Load karyotype files
      source_karyotype <- read.table(paste0("input/", source, "_karyotype.txt"), 
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      target_karyotype <- read.table(paste0("input/", target, "_karyotype.txt"), 
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      
      # Sort karyotypes such that chromosomes are ordered from largest to smallest
      source_karyotype <- source_karyotype[order(source_karyotype$End, 
                                                 decreasing = TRUE), ]
      target_karyotype <- target_karyotype[order(target_karyotype$End,
                                                 decreasing = TRUE), ]
      
      # Create dictionaries such that each karyotype name can be mapped to an integer
      source_mapping <- setNames(1:nrow(source_karyotype), source_karyotype$Chr) 
      target_mapping <- setNames(1:nrow(target_karyotype), target_karyotype$Chr)
      
      # Save the original chromosome names to a new column
      source_karyotype$orig_ids <- source_karyotype$Chr
      target_karyotype$orig_ids <- target_karyotype$Chr
      
      # Remap the karyotypes such that the names are integers
      source_karyotype$Chr <- source_mapping[source_karyotype$Chr]
      target_karyotype$Chr <- target_mapping[target_karyotype$Chr]
      
      # Reorder karyotypes as specified. 
      # Create a mapping to link karyotype to order of appearance 
      if (reorder == TRUE) {
        source_karyotype <- source_karyotype[order(match(source_karyotype$Chr, source_order)),]
        target_karyotype <- target_karyotype[order(match(target_karyotype$Chr, target_order)),]
        source_remap <- setNames(1:nrow(source_karyotype), source_order)
        target_remap <- setNames(1:nrow(target_karyotype), target_order)
      }
      
      # Combine the karyotypes into one matrix 
      karyotype_all <- rbind(source_karyotype, target_karyotype, 
                             make.row.names=FALSE)
      karyotype_all$orig_ids <- NULL
      # Set fill colours for nodes 
      karyotype_all$fill <- "cccccc"
      
      # Reorder the columns 
      karyotype_all <- karyotype_all[, c(1,2,3,7,4,5,6)]
      
      # Change the datatype of of the variables into a data frame 
      karyotype_all <- as.data.frame(karyotype_all, stringsAsFactors = F )
      head(karyotype_all)
      
      ########################################
      ### 2. Set orthologous gene coordinates
      ########################################
      
      # Load coordinates files 
      source_coords <- read.table(paste0("input/", source, "_alg_coordinates.tsv"), 
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      target_coords <- read.table(paste0("input/", target, "_alg_coordinates.tsv"), 
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      
      # Check file length
      nrow(source_coords)
      nrow(target_coords)
      
      # Get rid of coordinates that aren't on the selected chromosomes
      source_coords <- source_coords[source_coords$V3 %in% source_karyotype$orig_ids,]
      target_coords <- target_coords[target_coords$V3 %in% target_karyotype$orig_ids,]
      
      # Recheck file length
      nrow(source_coords)
      nrow(target_coords)
      
      # Check coordinates file
      head(source_coords)
      head(target_coords)
      
      # Map the chromosome names to appropriate karyotype numbers for both coordinate sets 
      source_coords$V3 <- source_mapping[source_coords$V3]
      target_coords$V3 <- target_mapping[target_coords$V3]
      
      # Flip any specified chromosomes 
      for (chromosome in source_flip){
        chromosome_length = source_karyotype[source_karyotype$Chr == chromosome, 'End']
        new_coords <- (chromosome_length - source_coords[source_coords$V3 == chromosome, c("V4", "V5")]) 
        source_coords[source_coords$V3 == chromosome, c("V4", "V5")] <- new_coords
      }
      for (chromosome in target_flip){
        chromosome_length = target_karyotype[target_karyotype$Chr == chromosome, 'End']
        new_coords <- (chromosome_length - target_coords[target_coords$V3 == chromosome, c("V4", "V5")]) 
        target_coords[target_coords$V3 == chromosome, c("V4", "V5")] <- new_coords
      }
      
      # Remap the chromosomes to the reshuffled chromosome order if necessary
      if (reorder == TRUE) {
        source_coords$V3 <- source_remap[as.character(source_coords$V3)]
        target_coords$V3 <- target_remap[as.character(target_coords$V3)]
      }
      
      # Create connections
      coords_all <- (merge(source_coords, target_coords, by = "V1")
                     [, c(-1, -2, -6, -7)])
      head(coords_all)
      
      # Relabel columns 
      colnames(coords_all) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2", "End_2", "ALG")
      
      # Create ALG colour map 
      # The 'XXXX' designation indicates orthology assignments made independently of ALG orthologs.
      alg_hexmap <- c("A1a" = "A6CDE3", "A1b" = "A6CDE3", "A2" = "197AB4", 
                      "B1" = "089E78", "B2" = "31A149", "B3" = "B2D287", 
                      "C1" = "F19799", "C2" = "E42027", 
                      "D"  = "F7F09B", 
                      "Ea" = "AF5B2A", "Eb" = "AF5B2A", 
                      "F"  = "DA6227", 
                      "G"  = "7671B0", 
                      "H"  = "E42C8A", 
                      "I"  = "65A744", 
                      "J1" = "C9B2D3", "J2" = "664292", 
                      "K"  = "E6AB22", 
                      "L"  = "A4782B", 
                      "M"  = "90CDC3", 
                      "N"  = "686868", 
                      "O1" = "EF7E25", "O2" = "F9BD6F", 
                      "P"  = "BDB9DA", 
                      "Qa" = "FED832", "Qb" = "FED832", "Qc" = "FED832", "Qd" = "FED832", 
                      "R"  = "00008B",
                      "XXXX" = "27A9E0")

      # Reassign line colours based on the mapping above
      coords_all$fill <- alg_hexmap[coords_all$ALG]
      
      # Check entries
      head(coords_all)
      
      # Filter out the jumping genes to reduce noise
      coords_all <- coords_all %>%
        group_by(Species_1, Species_2) %>%
        filter(n() > jumping_genes) %>%
        ungroup()
      
      # Check file
      head(coords_all)
      
      # Re-order connections by colour 
      coords_all <- coords_all[with(coords_all, order(fill)),]
      
      # Reformat into data frame to plot
      synteny_coord <- as.data.frame(coords_all[,-7])
      head(synteny_coord)
      head(karyotype_all)
      
      ########################################
      ### 3. Create synteny plot
      ########################################
      
      # Create synteny plot 
      ideogram(karyotype = karyotype_all, synteny = synteny_coord,  
               output = paste0(source, "_", target, "_synteny.svg"))
      
      # Convert svg file to pdf
      convertSVG(paste0(source, "_", target, "_synteny.svg"), 
                 file = paste0(source, "_", target, "_synteny.pdf"), 
                 device = "pdf")
      
      ########################################
      ### 4. Generate Oxford dot plot
      ########################################
      
      # Calculate cumulative sums
      # Initialize the Cumulative_End column with NA values
      source_karyotype$Cumulative_End <- 0
      target_karyotype$Cumulative_End <- 0
      
      # Iterate over the rows, starting from the second row
      for (i in 2:nrow(source_karyotype)) {
        source_karyotype$Cumulative_End[i] <- sum(source_karyotype$End[1:(i-1)])
      }
      for (i in 2:nrow(target_karyotype)) {
        target_karyotype$Cumulative_End[i] <- sum(target_karyotype$End[1:(i-1)])
      }
      
      # Separate cumulative sums by species (unit in Mb)
      source_cumsum <- source_karyotype[source_karyotype$species == source, "Cumulative_End"]/1e6
      target_cumsum <- target_karyotype[target_karyotype$species == target, "Cumulative_End"]/1e6
      
      # Adjust start sites
      synteny_coord$adjusted_start_1 <- synteny_coord$Start_1 + source_karyotype$Cumulative_End[synteny_coord$Species_1]
      synteny_coord$adjusted_start_2 <- synteny_coord$Start_2 + target_karyotype$Cumulative_End[synteny_coord$Species_2]
      
      # Create the plot
      # Check if all values in the ALG column are 'XXXX'
      if (all(coords_all$ALG == 'XXXX')) {
        # If all values are 'XXXX', use a fixed color for points
        dp <- ggplot(synteny_coord, aes(x=adjusted_start_1/1e6, y=adjusted_start_2/1e6)) + 
          geom_point(size=0.5, color="#27A9E0") +
          theme_classic()
      } else {
        # If there are other values, use the 'fill' column for coloring
        # Add '#' to the color codes in the 'fill' column
        synteny_coord$fill <- paste0("#", synteny_coord$fill)
        dp <- ggplot(synteny_coord, aes(x=adjusted_start_1/1e6, y=adjusted_start_2/1e6, color=fill)) + 
          scale_color_identity() +
          geom_point(size=0.5) + 
          theme_classic()
      }
      for (x in source_cumsum) {
        dp <- dp + geom_vline(xintercept=x, color="gray", size=0.3)
      }
      for (y in target_cumsum) {
        dp <- dp + geom_hline(yintercept=y, color="gray", size=0.3)
      }
      dp <- dp + labs(x=source, y=target) + theme(legend.position = "none")
      
      # Remove '#' from the fill column for downstream analysis
      synteny_coord$fill <- gsub("#", "", synteny_coord$fill)
      
      # Save Oxford plot as PDF
      ggsave(paste0(source, "_", target, "_oxford_plot.pdf"), dp, device = "pdf")
    
      ########################################
      ### 5. Calculate genome-wide synteny 
      ###    mixing rate 
      ########################################
      
      # Rank the 'adjusted_start' column and add as 'rank'
      synteny_coord$rank_1 <- rank(synteny_coord$adjusted_start_1)
      synteny_coord$rank_2 <- rank(synteny_coord$adjusted_start_2)

      # Calculate Spearman's coefficient for the adjusted start sites
      spearman_coefficient <- cor(synteny_coord$adjusted_start_1, synteny_coord$adjusted_start_2, method = "spearman")
      
      # Calculate synteny mixing rate
      mixing_rate <- 1 - abs(spearman_coefficient)
      
      # Store the mixing rate result in the matrix
      mixing_rate_matrix[source, target] <- mixing_rate
      
      ########################################
      ### 6. Calculate synteny mixing rate for 
      ###    orthologous chromosome pairs
      ########################################
      
      # Select orthologous chromosomes based on source species
      chromosomes <- sort(unique(synteny_coord$Species_1))

      # Loop through each chromosome and create a plot
      for (chr in chromosomes) {
        # Filter the data for the current chromosome
        chr_data <- subset(synteny_coord, Species_1 == chr)
        
        # Rank the 'adjusted_start' column and add as 'rank'
        chr_data$rank_1 <- rank(chr_data$adjusted_start_1)
        chr_data$rank_2 <- rank(chr_data$adjusted_start_2)
        
        # Calculate Spearman's coefficient for the adjusted start sites
        spearman_coefficient <- cor(chr_data$adjusted_start_1, chr_data$adjusted_start_2, method = "spearman")
        
        # Calculate synteny mixing rate
        chr_mixing_rate <- 1 - abs(spearman_coefficient)
        
        # Store the correlation result in the table
        new_row <- data.frame(
          Species1 = source,
          Species2 = target,
          OrthoChromosome = chr,
          MixingRate = chr_mixing_rate,
          stringsAsFactors = FALSE
        )
        chr_mixing_rate_data <- rbind(chr_mixing_rate_data, new_row)
      }
    }
  }
}

# End of loop. Please continue with the subsequent steps below for additional analyses.

########################################################################
### (a) Create a heatmap visualizing the genome-wide synteny mixing rate
########################################################################

# Check the matrix
head(mixing_rate_matrix)

# Convert the matrix to a long format for ggplot
mixing_rate_data <- melt(mixing_rate_matrix, varnames = c("Source", "Target"), value.name = "Mixing_rate")

# Create a heatmap
heatmap1 <- ggplot(mixing_rate_data, aes(x = Source, y = Target, fill = Mixing_rate)) +
  geom_tile() +
  scale_fill_gradient2(low = "royalblue", mid = "white",high = "darkorange", 
                       midpoint = 0.3, limit = c(0, 0.6)) +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  xlab("Source Species") +
  ylab("Target Species") +
  ggtitle("Heatmap of syteny mixing rate")

# Print the heatmap
print(heatmap1)

# Create a heatmap with hierarchical clustering
my_palette <- colorRampPalette(c("royalblue", "white", "darkorange"))(n = 64)
heatmap2 <- pheatmap(mixing_rate_matrix, 
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         cluster_rows=TRUE, col=my_palette, scale="none", 
                         key=TRUE, symkey=FALSE, density.info="none", trace="none", 
                         cexCol=1, cexRow=1, cellwidth=30, cellheight=30)

# Print the heatmap
print(heatmap2)

# Save the heatmap
ggsave("mixing_rate_heatmap_1.pdf", heatmap1, device = "pdf")
ggsave("mixing_rate_heatmap_2.pdf", heatmap2, device = "pdf")

########################################################################
### (b) Examine specific chromosome subsets within the Oxford dot plot
########################################################################

# Assign source and target species 
source <- "wsu"
target <- "cpa"

# (Rerun the script from line 59 to line 294)

# Replace chromosome number with the target chromosome number you are interested in
specific_chromosome <- 1

# Plot highlighting points from a specific target chromosome
dp1 <- ggplot(synteny_coord, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(size = 0.3, color = 'royalblue') +  # Plot all points in blue
  geom_point(data = subset(synteny_coord, Species_1 == specific_chromosome), 
             aes(color = "specific_chromosome"), size = 0.3) +  # Highlight specific chromosome in red
  geom_vline(xintercept = source_cumsum, color = "gray", size = 0.3) +
  geom_hline(yintercept = target_cumsum, color = "gray", size = 0.3) +
  theme_classic() +
  labs(x = source, y = target) +
  scale_color_manual(values = c("specific_chromosome" = "red")) + theme(legend.position = "none")

# Print the plot
print(dp1)

# Subset the data for chromosomes 1 and 2, replace them with target chromosomes
subset <- c(1,2)
subset_chrs <- synteny_coord[synteny_coord$Species_1 %in% subset, ]

# Create color mapping for the subset
color_mapping <- setNames(c("royalblue", "red"), subset)

# Plotting the chromosome subset
dp2 <- ggplot(subset_chrs, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(aes(color = as.factor(Species_1)), size = 0.3) +
  theme_classic() +
  labs(x = source, y = target, color = "Chromosome") +
  scale_color_manual(values = color_mapping)

# Print the plot
print(dp2)

# Plot highlighting points from a specific target ALG (B2=31A149)
dp3 <- ggplot(synteny_coord, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(size = 0.3, color = "#27A9E0") +  # Plot all points in blue
  geom_point(data = subset(synteny_coord, fill == "31A149"), 
             aes(color = "specific_ALG"), size = 0.3) +  # Highlight specific chromosome in red
  geom_vline(xintercept = source_cumsum, color = "gray", size = 0.3) +
  geom_hline(yintercept = target_cumsum, color = "gray", size = 0.3) +
  theme_classic() +
  labs(x = source, y = target) +
  scale_color_manual(values = c("specific_ALG" = "red")) + theme(legend.position = "none")

# Print the diagnostic plot
print(dp3)

# Loop through each chromosome and create a plot
chromosomes <- sort(unique(synteny_coord$Species_1))

for (chr in chromosomes) {
  # Filter the data for the current chromosome
  chr_data <- subset(synteny_coord, Species_1 == chr)
  # Create the plot for the current chromosome
  # The 'Number of OrthoChromosome' corresponds to the order in the synteny plot
  plot_name <- paste0(source, "_", "chr_", chr, "_", target, "_plot.pdf")
  dp4 <- ggplot(chr_data, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
    geom_point(size = 0.3, color = 'royalblue') +
    theme_classic() +
    labs(x = source, y = target, color = "Chromosome") +
    scale_color_manual(values = color_mapping) +
    ggtitle(paste("OrthoChromosome", chr))
  # Save the plot to a file
  ggsave(plot_name, dp4)
}

########################################################################
### (c) Calculate mixing rates for orthologous chromosome pairs
########################################################################

# In this section, we merge the evolutionary distance data with synteny 
# mixing rates for orthologous chromosome pairs.

# Check the mixing rate data
head(chr_mixing_rate_data)

# Parse the Newick format phylogenetic tree. 
# This tree is rooted, with 'Pecten maximus' (pma) as the outgroup.
tree_newick <- "(pma:0.502254533,(cmu:0.8598396118,(mme:0.2989849558,(bst:0.2983355281,(wsu:0.2429769398,cpa:0.1794422116):0.1702015822[100]):0.0824568567[100]):0.6346684761[100]):0.502254533);"
tree <- read.tree(text = tree_newick)

# Calibrate the tree based on the split of pma and cmu around the Cambrian Period (540-480 mya)
cal <- makeChronosCalib(tree, interactive = T)

# Check the calibration table
head(cal)
# node age.min age.max soft.bounds
# 1    7     480     540       FALSE

# Create chronogram with estimated divergence times
chronogram <- chronos(tree, calibration = cal)

# Calculate pairwise divergence times
divergence_times <- cophenetic(chronogram)

# Convert the distance matrix to a table with unique pairs
divergence_time_data <- melt(as.matrix(divergence_times))
divergence_time_data <- divergence_time_data[divergence_time_data$Var1 != divergence_time_data$Var2, ]
names(divergence_time_data) <- c("Species1", "Species2", "DivergenceTime")

# Check the phylogenetic distance data
head(divergence_time_data)

# Merge mixing rate and phylogenetic distance data
merged_data <- merge(chr_mixing_rate_data, divergence_time_data, by = c("Species1", "Species2"), all.x = TRUE)

# Choose the species located deepest within the tree for further analysis (e.g., wsu or cpa)
# Exclude cmu since it doesn't have one-to-one orthologous chromosome pairs with other species.
merged_data <- subset(merged_data, Species1 == "cpa" & Species2 != "cmu")

# Fit the linear model
linear_model <- lm(MixingRate ~ DivergenceTime, data = merged_data)
equation <- paste("y =",
                  format(coef(linear_model)[1], digits = 4),
                  "+",
                  format(coef(linear_model)[2], digits = 4),
                  "x")
model_summary <- summary(linear_model)
r_squared_linear <- model_summary$r.squared

# Fit the logarithmic model
log_model <- lm(MixingRate ~ log(DivergenceTime), data = merged_data)
intercept <- coef(log_model)[1]
slope <- coef(log_model)[2]
equation <- paste("y =", round(intercept, 3), "+", round(slope, 3), "* log(x)")
model_summary <- summary(log_model)
r_squared_log <- model_summary$r.squared

# Create a dot plot with a linear model
dp6 <- ggplot(merged_data, aes(x = DivergenceTime, y = MixingRate)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  theme_classic() +
  labs(title = "Relationship between phylogenetic distance and mixing rate",
       x = "Divergence Time (MYA)",
       y = "Mixing rate")

# Create a dot plot with a logarithmic model
dp7 <- ggplot(merged_data, aes(x = DivergenceTime, y = MixingRate)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ log(x), color = "black", se = FALSE) +
  theme_classic() +
  labs(title = "Relationship between phylogenetic distance and mixing rate",
       x = "Divergence Time (MYA)",
       y = "Mixing rate")

# Print the dot plot
print(dp6)
print(dp7)

# Save the data
write.csv(merged_data, file="mixing_rate.csv")