#############################################
# R script for divergence time estimation 
# and visualization using bppr and BPP output
#############################################

# Load required libraries
library(ape)        # phylogenetic trees
library(bppr)       # BPP post-processing
library(coda)       # posterior summarization
library(bpptools)   # visualization utilities


#############################################
# 1. Import MCMC posterior samples
#############################################

# Load MCMC chain produced by BPP
# The file should be a tab-delimited text file (e.g., "mcmc.txt")
chain <- read.delim("mcmc.txt", header = TRUE, sep = "\t")

# Convert coalescent units to absolute time (millions of years)
# u.m  = mean mutation rate per site per generation
# u.sd = standard deviation of mutation rate (here set to ~10% of u.m)
# g.m  = mean generation time in years
# g.sd = standard deviation of generation time
time_scaled <- bppr::msc2time.r(
  chain,
  u.m  = 8.1e-10,
  u.sd = 0.81e-10,
  g.m  = 2,
  g.sd = 0.2
)


#############################################
# 2. Posterior summary statistics
#############################################

# Posterior mean values for each parameter
posterior_means <- apply(time_scaled, 2, mean)

# 95% Highest Posterior Density (HPD) intervals
posterior_HPD <- coda::HPDinterval(coda::as.mcmc(time_scaled))

# Combine into a summary table
summary_table <- cbind(posterior_means, posterior_HPD)


#############################################
# 3. Phylogenetic tree input
#############################################

# Read the species tree inferred with BPP
# (exported in FigTree/Nexus format)
tree_raw <- read.tree("FigTree.tre")

# Some BPP trees in Nexus format store the actual tree in the field "UTREE1="
# Extract the tree for visualization
tree <- tree_raw$`UTREE1=`


#############################################
# 4. Divergence time visualization
#############################################

# Plot density tree of divergence times
# Parameters:
#  - thin    : thinning of posterior samples for visualization
#  - alpha   : transparency of density
#  - pfract  : proportion of posterior samples plotted
#  - y.offset: vertical spacing of lineages
mcmc2densitree(
  tree,
  time_scaled,
  "t_",
  thin    = 0.09,
  alpha   = 0.01,
  pfract  = 0.5,
  y.offset= 0.5
)


#############################################
# 5. Export figure
#############################################

# Save divergence time density tree as SVG for editing/publication
svg(filename = "storeria/arbol_divergencia.svg", width = 15, height = 6)

mcmc2densitree(
  tree,
  time_scaled,
  "t_",
  thin    = 0.09,
  alpha   = 0.01,
  pfract  = 0.1,
  y.offset= 0.5
)

title(xlab = "Divergence time (Ma)")
dev.off()
