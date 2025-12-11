##############################################################
# R script to plot admixture proportions on geographic map
##############################################################

# Load required packages
library(mapmixture)  # plotting admixture proportions as pie charts
library(ggplot2)     # saving high-quality figures

##############################################################
# 1. Load data
##############################################################

# Admixture proportions per individual or population
admixture1 <- read.csv("Storeria_clusters.csv")

# Geographic coordinates of each individual or population
coordinates <- read.csv("Coordenadas.csv")

##############################################################
# 2. Generate admixture map
##############################################################

# Create map with admixture pie charts
# pie_size       : size of pies
# boundary       : map limits (latitude & longitude)
# scalebar       : add scale bar to map
# cluster_cols   : colors for each genetic cluster
map1 <- mapmixture(
  admixture1,
  coordinates,
  pie_size = 1.2,
  boundary = c(ymin = 10, ymax = 45, xmin = -105, xmax = -70),
  scalebar = TRUE,
  cluster_cols = c("purple", "blue", "red", "green")
)

# Display the map
plot(map1)

##############################################################
# 3. Save high-quality SVG figure
##############################################################

ggsave(
  filename = "map1.svg",
  plot = map1,
  device = "svg",
  scale = 1,
  units = "cm",
  dpi = 300,
  width = 23.5,
  height = 12.5
)
