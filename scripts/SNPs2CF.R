##############################################################
# R script for converting SNP data to character frequencies
# using SNPs2CF functions
##############################################################

# Load required packages
library(pegas)       # handling genetic data
library(foreach)     # parallel processing
library(doMC)        # multicore backend for foreach

##############################################################
# 1. Load SNPs2CF functions
##############################################################

# Source the main functions from the SNPs2CF package
source("SNPs2CF-master/functions_v1.6.R")

##############################################################
# 2. Run SNPs2CF
##############################################################

# Arguments:
# seqMatrix      : input sequence matrix (e.g., .phy file)
# ImapName       : species/population mapping file
# rm.outgroup    : remove outgroup sequences? TRUE/FALSE
# outgroupSp     : name of outgroup species (if any)
# between.sp.only: only include SNPs between species?
# max.SNPs       : maximum number of SNPs to process
# bootstrap      : perform bootstrapping? TRUE/FALSE
# outputName     : name of the output CSV file
# save.progress  : save progress while running? TRUE/FALSE
# cores          : number of cores for parallel processing

SNPs2CF(
  seqMatrix       = "dekayi92.phy",
  ImapName        = "Imap4.txt",
  rm.outgroup     = FALSE,
  outgroupSp      = "S_occipitomaculata",
  between.sp.only = TRUE,
  max.SNPs        = 100000,
  bootstrap       = TRUE,
  outputName      = "Storeria.csv",
  save.progress   = FALSE,
  cores           = 6
)
