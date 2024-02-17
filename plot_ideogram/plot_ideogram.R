# R script for metazoan synteny analysis
# Isabel Jiah-Yih Liao, Tom Lewin, Yi-Jyun Luo
# Symbiosis Genomics & Evolution Lab July 2023

# This code works with two subdirectories within the Synteny working directory: 
# Synteny/input contains karyotype and coordinate files for each organism
# Synteny/ideograms is the output folder containing the ideogram svg and pdf files. 

# Setting the working directory
# The './input' directory should contain karyotype and coordinate files for each species.
# To set the working directory to the parent folder of 'input', use setwd('path_to_input').
# For example, if the full path to 'input' is 'path_to_input/input', use setwd('path_to_input').

# Set working directory
#setwd("Desktop/Bioinformatics/Synteny/")

# Load required packages 

require(RIdeogram)
require(dplyr)

# Assign source and target species 
source <- "lan"
target <- "pma"

# If you would like to reorder the chromosomes, put the chromosome order here. 
# Otherwise, set reorder to FALSE
reorder <- TRUE
bfl_order <- c(4,11,2,15,19,9,3,13,18,1,16,8,7,6,10,5,17,14,12)
pma_order <- c(10,12,8,9,4,5,17,1,19,3,13,2,6,11,15,7,18,16,14)
llo_order <- c(8,6,2,4,5,14,1,17,3,18,13,16,9,11,15,7,19,12,10)
lan_order <- c(1,2,7,3,4,5,6,8,9,10)
source_order <- lan_order
target_order <- pma_order
# List any chromosomes in either the source or the target genome which should be flipped 
source_flip <- c()
target_flip <- c()

############################
### Set karyotypes
############################

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

############################
### Set genes 
############################

# Load coordinates files 
source_coords <- read.table(paste0("input/", source, "_coordinates.tsv"), 
                            sep = "\t", header = FALSE, stringsAsFactors = F)
target_coords <- read.table(paste0("input/", target, "_coordinates.tsv"), 
                            sep = "\t", header = FALSE, stringsAsFactors = F)

# Check file length
nrow(source_coords)
nrow(target_coords)

# Get rid of coordinates that aren't on the selected chromosomes
source_coords <- source_coords[source_coords$V3 %in% source_karyotype$orig_ids,]
target_coords <-target_coords[target_coords$V3 %in% target_karyotype$orig_ids,]

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
# Relabel columns 
colnames(coords_all) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2", "End_2", "ALG")

# Create ALG colour map 
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
                "R"  = "00008B")

# Reassign line colours based on the mapping above
coords_all$fill <- alg_hexmap[coords_all$ALG]

# Check entries
head(coords_all)

# Filter out the jumping genes to reduce noise
coords_all <- coords_all %>%
    group_by(Species_1, Species_2) %>%
    filter(n() > 5) %>%
    ungroup()

# Check file
head(coords_all)

# Re-order connections by colour 
coords_all <- coords_all[with(coords_all, order(fill)),]

# Remove the 7th column and reformat into data frame to plot
synteny_coord <- as.data.frame(coords_all[,-7])
head(synteny_coord)
head(karyotype_all)

# Create synteny plot 
ideogram(karyotype = karyotype_all, synteny = synteny_coord,  
         output = paste0("ideograms/svg/", source, "_", target, "_synteny.svg"))

# Convert svg file to pdf
convertSVG(paste0("ideograms/svg/", source, "_", target, "_synteny.svg"), 
           file = paste0("ideograms/pdf/", source, "_", target, "_synteny.pdf"), 
           device = "pdf")



