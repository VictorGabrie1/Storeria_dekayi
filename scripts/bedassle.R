##############################################################
# R script for preparing allele frequency and distance data
# for BEDASSLE analyses
##############################################################

# Load required packages
library(BEDASSLE)    # main analysis
library(fossil)      # distance utilities
library(geosphere)   # geographic distance calculations
library(writexl)     # optional: write Excel output
library(dplyr)       # data manipulation
library(readxl)      # read Excel files


##############################################################
# 1. Convert STRUCTURE format file to allele/sample matrices
##############################################################

# Input files
structure_file <- "dekayistructure"   # STRUCTURE file (edited: pop column removed, spaces → tabs)
popmap_file    <- "popmap.txt"        # popmap file: sampleID   population

# Read STRUCTURE-like file
x    <- read.table(structure_file, sep = "\t")
pops <- read.table(popmap_file, sep = "\t")

# Extract individual IDs
inds   <- as.character(x[2:nrow(x), 1])
n.pops <- length(unique(pops[,2]))

# Remove population ID column if present
if (is.na(x[1, 2]) == TRUE) {
  x <- cbind(x[,1], x[,3:ncol(x)])
}

# SNP names
snp.names <- as.vector(as.character(x[1, 2:ncol(x)]))


##############################################################
# 2. Filter SNPs: retain only loci present in all populations
##############################################################

for (a in 1:ncol(x)) {
  if (a == 1) {
    x2 <- as.character(x[,1])
  } else {
    x.test <- x[2:nrow(x), a]
    skip <- 0
    for (b in 1:n.pops) {
      test <- match(inds, as.character(pops[pops[,2] == b, 1]))
      test[is.na(test)] <- 0
      test <- x.test[test > 0]
      test <- test[test > 0]
      if (length(test) == 0) {
        skip <- skip + 1
      }
    }
    if (skip == 0) {
      x.test <- c(snp.names[a-1], x.test)
      x2 <- cbind(x2, x.test)
    }
  }
}

# Update dataset
x <- x2
snp.names <- x[1, 2:ncol(x)] |> as.vector()


##############################################################
# 3. Filter SNPs: retain only biallelic loci
##############################################################

for (a in 1:ncol(x)) {
  if (a == 1) {
    x2 <- x[,1]
  } else {
    x.test <- x[2:nrow(x), a]
    x.test <- x.test[x.test > 0]
    if (length(unique(x.test)) == 2) {
      x.test <- c(snp.names[a-1], x[2:nrow(x), a])
      x2 <- cbind(x2, x.test)
    }
  }
}

# Update dataset
x <- x2
snp.names <- x[1, 2:ncol(x)] |> as.vector()


##############################################################
# 4. Build allele frequency and sample size matrices
##############################################################

for (a in 1:ncol(x)) {
  if (a == 1) {
    next
  } else {
    x.rep   <- as.vector(x[2:nrow(x), a])
    allele.a <- x.rep[x.rep > 0][1]  # reference allele
    alleles <- c()
    samples <- c()
    for (b in 1:n.pops) {
      test <- match(inds, as.character(pops[pops[,2] == b, 1]))
      test[is.na(test)] <- 0
      test <- x.rep[test > 0]
      alleles <- c(alleles, length(na.omit(match(test, allele.a))))
      samples <- c(samples, length(test[test > 0]))
    }
    if (a == 2) {
      allele.matrix      <- alleles
      sample.size.matrix <- samples
    } else {
      allele.matrix      <- cbind(allele.matrix, alleles)
      sample.size.matrix <- cbind(sample.size.matrix, samples)
    }
  }
}

# Assign row/column names
colnames(allele.matrix)      <- snp.names
rownames(allele.matrix)      <- 1:n.pops
colnames(sample.size.matrix) <- snp.names
rownames(sample.size.matrix) <- 1:n.pops

# Save matrices
write.table(allele.matrix, file = "bedassle.allele.counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(sample.size.matrix, file = "bedassle.sample.sizes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##############################################################
# 5. Build geographic distance matrix (Haversine distance)
##############################################################

# Load coordinates file with columns: Sample | LAT | LON
coords <- read.table("sample_info.txt", sep = "\t", header = TRUE)

# Extract latitude/longitude
latlong <- as.matrix(coords[, c("LON", "LAT")])

# Calculate pairwise geographic distances in meters
geo_dist_matrix <- distm(latlong, fun = distHaversine)

# Convert to kilometers
geo_dist_matrix_km <- geo_dist_matrix / 1000

# Save as CSV
write.csv(geo_dist_matrix_km, file = "matriz_distancias.csv", row.names = FALSE)

# Reload distance matrix as BEDASSLE-compatible object
D <- read.csv("matriz_distancias.csv", header = TRUE) |> as.matrix()


##############################################################
# 6. Import allele/sample matrices (edited manually in Excel)
##############################################################
# NOTE: before importing, edit allele.counts and sample.sizes files:
# - Add row/column headers (population IDs, locus names)
# - Save as tab-delimited text.

counts <- read.table("bedassle.allele.counts.txt", sep = "\t", header = TRUE, row.names = 1) |> as.matrix()
sample_size <- read.table("bedassle.sample.sizes.txt", sep = "\t", header = TRUE, row.names = 1) |> as.matrix()


##############################################################
# 7. Save BEDASSLE input objects for downstream analysis
##############################################################

save(counts, sample_size, D, file = "Storeria.Robj")

##############################################################
# R script to calculate genetic (Fst) and geographic distances
# and visualize Isolation by Distance (IBD)
##############################################################

# Load required packages
library(BEDASSLE)   # Fst and allele frequency analyses
library(reshape2)    # for melting matrices
library(ggplot2)     # plotting

##############################################################
# 1. Load prepared data
##############################################################

# Load previously prepared BEDASSLE objects
# counts       = allele counts per population/locus
# sample_size  = number of sampled alleles per population/locus
# D            = geographic distance matrix (km)
load("Storeria.Robj")


##############################################################
# 2. Function to calculate pairwise distances and organize data
##############################################################

get_dist <- function() {
  
  # Assign populations to major groups
  # NOTE: ensure the sample ordering matches the counts matrix
  sample_info <- data.frame(SAMPLE_ID = rownames(counts))
  sample_info$POP_GROUP <- c(rep("Mexico", 19),
                             rep("USAEast", 2),
                             rep("USAWest", 12),
                             rep("VCC", 2))
  
  # Geographic distances
  geo_dist <- D
  colnames(geo_dist) <- sample_info$SAMPLE_ID
  rownames(geo_dist) <- sample_info$SAMPLE_ID
  diag(geo_dist) <- NA            # Remove self-comparisons
  geo_dist[geo_dist == 0] <- NA  # Remove zero distances (same location)
  
  # Genetic distances: pairwise Fst using BEDASSLE
  gen_dist <- calculate.all.pairwise.Fst(allele.counts = counts,
                                         sample.sizes = sample_size)
  gen_dist <- as.matrix(gen_dist)
  
  # Ensure matrix is symmetric
  if (!isSymmetric(gen_dist)) stop("Genetic distance matrix is not symmetric.")
  
  # Assign row/column names and remove diagonal
  colnames(gen_dist) <- rownames(counts)
  rownames(gen_dist) <- rownames(counts)
  diag(gen_dist) <- NA
  
  # Convert matrices to long format for plotting/analysis
  geo_vector <- melt(geo_dist)
  gen_vector <- melt(gen_dist)
  
  # Combine into a single dataframe
  dist_df <- data.frame(
    SAMPLE1 = gen_vector$Var1,
    SAMPLE2 = gen_vector$Var2,
    Genetic_distance   = gen_vector$value,
    Geographic_distance = geo_vector$value
  )
  
  # Add population group info
  dist_df$POP1 <- sample_info$POP_GROUP[match(dist_df$SAMPLE1, sample_info$SAMPLE_ID)]
  dist_df$POP2 <- sample_info$POP_GROUP[match(dist_df$SAMPLE2, sample_info$SAMPLE_ID)]
  
  # Label intra- vs inter-group comparisons
  dist_df$comparison <- ifelse(dist_df$POP1 == dist_df$POP2, "intragroup", "intergroup")
  
  # Save results for downstream use
  write.csv(dist_df, "distances_results.csv", row.names = FALSE)
  
  # Return dataframe for plotting or further analyses
  return(dist_df)
}


##############################################################
# 3. Function to plot Isolation by Distance (IBD)
##############################################################

plot_IBD <- function(dist_df) {
  
  ggplot(data = dist_df) +
    geom_point(aes(x = Geographic_distance, y = Genetic_distance, fill = comparison),
               color = "black", size = 2, shape = 21, alpha = 0.5, stroke = 0.2) +
    xlab("Geographic distance (km)") +
    ylab("Genetic distance (Fst)") +
    ylim(0, 1) +
    guides(fill = "none") +
    scale_fill_manual(values = c("black", "gray")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 14),
      panel.border = element_rect(linewidth = 0.75, colour = "gray20"),
      panel.grid = element_blank()
    )
}


##############################################################
# 4. Run analysis and plot
##############################################################

dist_df <- get_dist()       # compute pairwise distances
plot_IBD(dist_df)           # visualize Isolation by Distance

