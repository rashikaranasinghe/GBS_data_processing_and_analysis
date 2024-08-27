###################################################################
## Flameback_GBS_analysis_script.R
# started by Rashika W. Ranasinghe 23 July 2024
 
# This file includes R scripts for generating SNP-based PCA, genomic cline analysis, and phenotypic cline analysis, as used in the manuscript titled: "Cryptic Hybridization Dynamics in a Three-Way Hybrid Zone of Dinopium Flamebacks on a Tropical Island".


############################################################################
################## 1.  Generate SNPs-based PCA ############################
############################################################################

# Make this PCA

# Requred data files :
# 1. .indv 
# 2. .pos
# 3. .012
# Please refer to the Bioinformatics file ******* for instructions on how to produce these files on Linux.

# set the working drectory 
setwd("/Users/rashikaranasinghe/Library/CloudStorage/OneDrive-UBC/Dinopium_GBS_on_laptop/012NA_150_samples/")

# Load functions
# source file can be downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4j2662g
source("/Users/rashikaranasinghe/Library/CloudStorage/OneDrive-UBC/Dinopium_GBS_on_laptop/012NA_150_samples/Dinopium_GBS_R_project/genomics_R_functions_V2.R")

# loada the packages
library(tidyverse)

# Set up requred variables
groups_and_colors <- rbind(c("Dino_R_Sym", "red"),
                           c("Dino_R_Allo", "red"),
                           c("Dino_Y_J", "green"),
                           c("Dino_Y_M", "blue"),
                           c("Dino_Y_Sym", "yellow"),
                           c("Dino_H", "orange"))
groups.to.plot.PCA <- groups_and_colors[,1]
group.colors.PCA <- groups_and_colors[,2]
groups <- c("Dino_R_Sym", "Dino_R_Allo", "Dino_Y_J", "Dino_Y_M", "Dino_Y_Sym", "Dino_H") # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)

base.file.name <- "DW.chrom_150.n_Oct.2022.genotypes.variant_only.whole_genome.noIndels_biallelic.GATKfiltered.maxmiss60.mac3.70pcntIndMiss.SNP30miss.tab"
pos <- read.table(paste0(base.file.name, ".012.pos"), col.names = c("chrom", "position"))
# Change chromosome numbers
pos$chrom <- str_replace_all(pos$chrom, c("CM025994.1"="1", "CM025995.1"="2", "CM025996.1"="3", "CM025997.1"="4", "CM025998.1"="5", "CM025999.1"="6", "CM026000.1"="7", "CM026001.1"="8", "CM026002.1"="9", "CM026003.1"="10", "CM026004.1"="11", "CM026005.1"="12", "CM026006.1"="13", "CM026007.1"="14", "CM026008.1"="15", "CM026009.1"="16", "CM026010.1"="17", "CM026011.1"="18", "CM026012.1"="19", "CM026013.1"="20", "CM026014.1"="21", "CM026015.1"="22", "CM026016.1"="23", "CM026017.1"="24", "CM026018.1"="25", "CM026019.1"="26", "CM026020.1"="27", "CM026021.1"="28", "CM026022.1"="29", "CM026023.1"="30", "CM026024.1"="31", "CM026025.1"="32", "CM026026.1"="33", "CM026027.1"="34", "CM026028.1"="35", "CM026029.1"="36", "CM026030.1"="37", "CM026031.1"="38", "CM026032.1"="39", "CM026033.1"="40", "CM026034.1"="41", "CM026035.1"="42", "CM026036.1"="43", "CM026037.1"="44", "CM026038.1"="W", "CM026039.1"="Z"))
column_names <- c("null", paste("c", pos$chrom, pos$position, sep="_"))
#geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
#save(geno, file = "geno")
load(file = "geno")
dim(geno) # [1]      114 1251423
SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
SNPnum
##SNPs = 1,251,422
ind <- read.table(paste0(base.file.name, ".012.indv"))
locations <- read.table("Fst_groups_III.txt", header=TRUE, sep = "\t", quote = "")
num_loc_cols <- ncol(locations)
ind_with_locations <- cbind(ind,locations)
combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
X <- 70   # this is the percentage threshold
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo[(num_loc_cols+1):ncol(combo)]))
numNAs_by_ID <- data.frame(combo$ID, numNAs)  # useful to see numNAs per sample: numNAs_by_ID
selection <- (numNAs < threshold_NA)
if(any(is.na(selection))) cat("selection contains NA values\n")  # this is a check for noticing errors / bugs
combo.NApass.all <- combo[selection,]
dim(combo.NApass.all)   # [1]    114 1251427
combo$ID[which(selection==F)]  # Non-left out for 70%

# filter out SNPs with too many missing genotypes:
SNP_NAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):ncol(combo.NApass.all)]))
X <- 10   # this is the percentage threshold
threshold_SNP_NAs <- length(combo.NApass.all[,1]) * X/100
selection <- (SNP_NAs <= threshold_SNP_NAs)
if(any(is.na(selection))) cat("selection contains NA values\n")  # this is a check for noticing errors / bugs
combo.NApass.subset <- combo.NApass.all[, c(rep(TRUE, times=num_loc_cols),selection)]
pos.subset <- pos[selection,]
dim(pos.subset)

# option to filter out all except selected chromosome (or set of them):
choose.chrom <- F
if (choose.chrom == TRUE) {
  chrom <- "1"
  # selection <- (pos.subset$chrom == chrom)
  selection <- (pos.subset$chrom == chrom)
  if(any(is.na(selection))) cat("selection contains NA values\n")  # this is a check for noticing errors / bugs
  #pos.subset.one.chr <- pos.subset[selection,]
  #loci.selection <- c(rep(TRUE, times=num_loc_cols), selection)  # add placeholders for info columns
  # which(loci.selection == T)    # to check which entries are TRUE
  combo.NApass <- combo.NApass.subset[, c(rep(TRUE, times=num_loc_cols), selection)]
  pos.NApass <- pos.subset[selection,]
  region.text <- paste0("chr", chrom)
}	else {
  region.text <- "whole_genome"
  combo.NApass <- combo.NApass.subset
  pos.NApass <- pos.subset
}


# Calculate allele freqs and sample sizes (use column Fst_group)
temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)
# calculate WC84_Fst 
temp.list <- getWC84Fst(freqs, sample_size, groups, among=TRUE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
WC84_Fst <- temp.list$WC84_Fst
rm(temp.list)
#dim(WC84_Fst)
# make the figure:
Fst.filter <- F   # option to filter to high-Fst markers only, using cutoff below
Fst.cutoff <- 0.5  # has no effect if Fst.filter is FALSE
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
axes <- 3

PCA_DW.chrom_n108_734806SNP <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text, groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=T)
## The PCA has some individuals that have genotype and phenotype missmatch. seems like genotype doesn't really match the phenotype.

####### Variations explained by each PCs for n108_734806SNP PCA
PCA.table_PCA_DW.chrom_n108_734806SNP <- data.frame(PCA_DW.chrom_n108_734806SNP$scores, PCA_DW.chrom_n108_734806SNP$data[,1:3])

############################################################################
############################ END SNPs-based PCA ############################
############################################################################



############################################################################
######################## 2. Genomic cline analysis ########################
############################################################################

###### Notes:
# The distances were calculated from Point Calimere (10.2845° N, 79.8241° E) in India to each sampling location using ruler tool in google earth Pro.
# Admixture values (for K=3) calculated from ADMIXTURE(v.) are used as the genomic data. 
# When multiple birds were captured at a single sampling location, their ancestry proportions were averaged to represent the genomic data for that site.

#set working directory
setwd("~/Documents/Dinopium_GBS_on_laptop/012NA_150_samples/Cline_analysis_n.150_2022/")

# load the requred libraries
library(ggplot2)
library(sf)
library(stringr)
library(raster)
library(sp)
library(splines)
library(dplyr)
library(hzar)
library(viridis)
library(forcats) 


## Read the data file
admix_data <- read.csv("Geno cline/admix_data.k3_locality_data_added_with_mean_ancestry.csv")

### make data object of allele frequency data
all.ind.k3_Q <-
  hzar.doMolecularData1DPops(admix_data$distance_km,
                             admix_data$ancestry_D.psarodes,
                             admix_data$nSamples)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(all.ind.k3_Q)

## Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
all.ind.k3_Q_model_free_both <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="free", tails="both")
all.ind.k3_Q_model_free_none <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="free", tails="none")
all.ind.k3_Q_model_free_right <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="free", tails="right")
all.ind.k3_Q_model_free_left <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="free", tails="left")
all.ind.k3_Q_model_free_mirror <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="free", tails="mirror")

all.ind.k3_Q_model_fixed_both <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="fixed", tails="both")
all.ind.k3_Q_model_fixed_none <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="fixed", tails="none")
all.ind.k3_Q_model_fixed_right <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="fixed", tails="right")
all.ind.k3_Q_model_fixed_left <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="fixed", tails="left")
all.ind.k3_Q_model_fixed_mirror <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="fixed", tails="mirror")

all.ind.k3_Q_model_none_both <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="none", tails="both")
all.ind.k3_Q_model_none_none <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="none", tails="none")
all.ind.k3_Q_model_none_right <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="none", tails="right")
all.ind.k3_Q_model_none_left <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="none", tails="left")
all.ind.k3_Q_model_none_mirror <- hzar.makeCline1DFreq(all.ind.k3_Q, scaling="none", tails="mirror")


## The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
all.ind.k3_Q_model_free_both <- hzar.model.addBoxReq(all.ind.k3_Q_model_free_both,-100,750)
all.ind.k3_Q_model_free_left <- hzar.model.addBoxReq(all.ind.k3_Q_model_free_left,-100,750)
all.ind.k3_Q_model_free_right <- hzar.model.addBoxReq(all.ind.k3_Q_model_free_right,-100,750)
all.ind.k3_Q_model_free_mirror <- hzar.model.addBoxReq(all.ind.k3_Q_model_free_mirror,-100,750)
all.ind.k3_Q_model_free_none <- hzar.model.addBoxReq(all.ind.k3_Q_model_free_none,-100,750)

all.ind.k3_Q_model_fixed_both <- hzar.model.addBoxReq(all.ind.k3_Q_model_fixed_both,-100,750)
all.ind.k3_Q_model_fixed_left <- hzar.model.addBoxReq(all.ind.k3_Q_model_fixed_left,-100,750)
all.ind.k3_Q_model_fixed_right <- hzar.model.addBoxReq(all.ind.k3_Q_model_fixed_right,-100,750)
all.ind.k3_Q_model_fixed_mirror <- hzar.model.addBoxReq(all.ind.k3_Q_model_fixed_mirror,-100,750)
all.ind.k3_Q_model_fixed_none <- hzar.model.addBoxReq(all.ind.k3_Q_model_fixed_none,-100,750)

all.ind.k3_Q_model_none_both <- hzar.model.addBoxReq(all.ind.k3_Q_model_none_both,-100,750)
all.ind.k3_Q_model_none_left <- hzar.model.addBoxReq(all.ind.k3_Q_model_none_left,-100,750)
all.ind.k3_Q_model_none_right <- hzar.model.addBoxReq(all.ind.k3_Q_model_none_right,-100,750)
all.ind.k3_Q_model_none_mirror <- hzar.model.addBoxReq(all.ind.k3_Q_model_none_mirror,-100,750)
all.ind.k3_Q_model_none_none <- hzar.model.addBoxReq(all.ind.k3_Q_model_none_none,-100,750)


#### cline model fitting
## generate an hzar.fitRequest object suitable for hzar.doFit
all.ind.k3_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_free_both, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_free_left, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_free_right, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_free_mirror, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_free_none, all.ind.k3_Q, verbose=FALSE)

all.ind.k3_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_fixed_both, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_fixed_left, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_fixed_right, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_fixed_mirror, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_fixed_none, all.ind.k3_Q, verbose=FALSE)

all.ind.k3_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_none_both, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_none_left, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_none_right, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_none_mirror, all.ind.k3_Q, verbose=FALSE)
all.ind.k3_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=all.ind.k3_Q_model_none_none, all.ind.k3_Q, verbose=FALSE)

## set mcmc chain length and buRAllo in
all.ind.k3_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_free_bothFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_free_leftFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_free_rightFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_free_mirrorFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_free_noneFitR$mcmcParam$buRAlloin <- 5e5

all.ind.k3_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_fixed_bothFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_fixed_leftFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_fixed_rightFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_fixed_mirrorFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_fixed_noneFitR$mcmcParam$buRAlloin <- 5e5

all.ind.k3_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_none_bothFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_none_leftFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_none_rightFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_none_mirrorFitR$mcmcParam$buRAlloin <- 5e5
all.ind.k3_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
all.ind.k3_Q_model_none_noneFitR$mcmcParam$buRAlloin <- 5e5

# ### Run the optimizer using the parameters listed in the hzar.fitRequest given.
all.ind.k3_Q_model_free_bothFit <- hzar.doFit(all.ind.k3_Q_model_free_bothFitR)
all.ind.k3_Q_model_free_leftFit <- hzar.doFit(all.ind.k3_Q_model_free_leftFitR)
all.ind.k3_Q_model_free_rightFit <- hzar.doFit(all.ind.k3_Q_model_free_rightFitR)
all.ind.k3_Q_model_free_mirrorFit <- hzar.doFit(all.ind.k3_Q_model_free_mirrorFitR)
all.ind.k3_Q_model_free_noneFit <- hzar.doFit(all.ind.k3_Q_model_free_noneFitR)

all.ind.k3_Q_model_fixed_bothFit <- hzar.doFit(all.ind.k3_Q_model_fixed_bothFitR)
all.ind.k3_Q_model_fixed_leftFit <- hzar.doFit(all.ind.k3_Q_model_fixed_leftFitR)
all.ind.k3_Q_model_fixed_rightFit <- hzar.doFit(all.ind.k3_Q_model_fixed_rightFitR)
all.ind.k3_Q_model_fixed_mirrorFit <- hzar.doFit(all.ind.k3_Q_model_fixed_mirrorFitR)
all.ind.k3_Q_model_fixed_noneFit <- hzar.doFit(all.ind.k3_Q_model_fixed_noneFitR)

all.ind.k3_Q_model_none_bothFit <- hzar.doFit(all.ind.k3_Q_model_none_bothFitR)
all.ind.k3_Q_model_none_leftFit <- hzar.doFit(all.ind.k3_Q_model_none_leftFitR)
all.ind.k3_Q_model_none_rightFit <- hzar.doFit(all.ind.k3_Q_model_none_rightFitR)
all.ind.k3_Q_model_none_mirrorFit <- hzar.doFit(all.ind.k3_Q_model_none_mirrorFitR)
all.ind.k3_Q_model_none_noneFit <- hzar.doFit(all.ind.k3_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and retuRAllo the mcmc data with an added a log likelihood column.
par(mar=c(1,1,1,1))
plot(hzar.mcmc.bindLL(all.ind.k3_Q_model_free_bothFit))
plot(hzar.mcmc.bindLL(all.ind.k3_Q_model_free_leftFit))
plot(hzar.mcmc.bindLL(all.ind.k3_Q_model_free_rightFit))
plot(hzar.mcmc.bindLL(all.ind.k3_Q_model_free_mirrorFit))
plot(hzar.mcmc.bindLL(all.ind.k3_Q_model_free_noneFit))
# 
##group multiple fits of the same model and the same observation data into a single object
all.ind.k3_Q_model_free_bothData <- hzar.dataGroup.add(all.ind.k3_Q_model_free_bothFit)
all.ind.k3_Q_model_free_leftData <- hzar.dataGroup.add(all.ind.k3_Q_model_free_leftFit)
all.ind.k3_Q_model_free_rightData <- hzar.dataGroup.add(all.ind.k3_Q_model_free_rightFit)
all.ind.k3_Q_model_free_mirrorData <- hzar.dataGroup.add(all.ind.k3_Q_model_free_mirrorFit)
all.ind.k3_Q_model_free_noneData <- hzar.dataGroup.add(all.ind.k3_Q_model_free_noneFit)

all.ind.k3_Q_model_fixed_bothData <- hzar.dataGroup.add(all.ind.k3_Q_model_fixed_bothFit)
all.ind.k3_Q_model_fixed_leftData <- hzar.dataGroup.add(all.ind.k3_Q_model_fixed_leftFit)
all.ind.k3_Q_model_fixed_rightData <- hzar.dataGroup.add(all.ind.k3_Q_model_fixed_rightFit)
all.ind.k3_Q_model_fixed_mirrorData <- hzar.dataGroup.add(all.ind.k3_Q_model_fixed_mirrorFit)
all.ind.k3_Q_model_fixed_noneData <- hzar.dataGroup.add(all.ind.k3_Q_model_fixed_noneFit)

all.ind.k3_Q_model_none_bothData <- hzar.dataGroup.add(all.ind.k3_Q_model_none_bothFit)
all.ind.k3_Q_model_none_leftData <- hzar.dataGroup.add(all.ind.k3_Q_model_none_leftFit)
all.ind.k3_Q_model_none_rightData <- hzar.dataGroup.add(all.ind.k3_Q_model_none_rightFit)
all.ind.k3_Q_model_none_mirrorData <- hzar.dataGroup.add(all.ind.k3_Q_model_none_mirrorFit)
all.ind.k3_Q_model_none_noneData <- hzar.dataGroup.add(all.ind.k3_Q_model_none_noneFit)

### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
all.ind.k3_Q_modelNull <- hzar.dataGroup.null(all.ind.k3_Q)

##make list of cline models and null models
all.ind.k3_Q_dGs <- list(cline_free_bothModel = all.ind.k3_Q_model_free_bothData,
                         cline_free_leftModel = all.ind.k3_Q_model_free_leftData,
                         cline_free_rightModel = all.ind.k3_Q_model_free_rightData,
                         cline_free_mirrorModel = all.ind.k3_Q_model_free_mirrorData,
                         cline_free_noneModel = all.ind.k3_Q_model_free_noneData,
                         cline_fixed_bothModel = all.ind.k3_Q_model_fixed_bothData,
                         cline_fixed_leftModel = all.ind.k3_Q_model_fixed_leftData,
                         cline_fixed_rightModel = all.ind.k3_Q_model_fixed_rightData,
                         cline_fixed_mirrorModel = all.ind.k3_Q_model_fixed_mirrorData,
                         cline_fixed_noneModel = all.ind.k3_Q_model_fixed_noneData,
                         cline_none_bothModel = all.ind.k3_Q_model_none_bothData,
                         cline_none_leftModel = all.ind.k3_Q_model_none_leftData,
                         cline_none_rightModel = all.ind.k3_Q_model_none_rightData,
                         cline_none_mirrorModel = all.ind.k3_Q_model_none_mirrorData,
                         cline_none_noneModel = all.ind.k3_Q_model_none_noneData,
                         nullModel = all.ind.k3_Q_modelNull)

##Collect optimizer output based on the same hzar.obsData object
all.ind.k3_Q_oDG <- hzar.make.obsDataGroup(all.ind.k3_Q_dGs) #

##Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
all.ind.k3_Q_oDG <- hzar.copyModelLabels(all.ind.k3_Q_dGs,all.ind.k3_Q_oDG)

##Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(all.ind.k3_Q_oDG)

##Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(all.ind.k3_Q_oDG))
#                            AICc
# cline_free_bothModel    36.37253
# cline_free_leftModel    32.56917
# cline_free_rightModel   32.62327
# cline_free_mirrorModel  32.69089
# cline_free_noneModel    28.04608
# cline_fixed_bothModel   31.87585
# cline_fixed_leftModel   27.90802
# cline_fixed_rightModel  27.31332
# cline_fixed_mirrorModel 27.67267
# cline_fixed_noneModel   23.63733
# cline_none_bothModel    31.77375
# cline_none_leftModel    27.90845
# cline_none_rightModel   27.47133
# cline_none_mirrorModel  27.58143
# cline_none_noneModel    23.63724
# nullModel               89.88868
##### The best model has the lowers AIC value

### make em purty, get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
min(hzar.AICc.hzar.obsDataGroup(all.ind.k3_Q_oDG))
# the best fit model is cline_free_noneModel    25.24257
#cline_fixed_noneModel  = all.ind.k3_Q_model_fixed_noneData

cline.center <- all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$center
print(all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$center) # Centre of the cline = 172.8894
cline.width <- all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$width
print(all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$width) # Width of the cline = 130.7505
print(hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))) # get the boundaries
#   center2LLLow center2LLHigh width2LLLow width2LLHigh
#1     153.8569      191.0664    92.69923     197.2396

### Caclculate likelyhoods of clide widths
Northern.border.Upper.bound.Likelihood <- cline.center - (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[4]/2) # 74.26966
Northern.border.of.cline.width <- (cline.center-(cline.width/2)) # 107.5142
Northern.border.Lower.bound.Likelihood <- cline.center - (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[3]/2) # 126.5398
Southern.border.Upper.bound.Likelihood <- cline.center + (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[3]/2) # 219.2391
Southern.border.of.cline.width <- (cline.center+(cline.width/2)) # 238.2647
Sourthern.border.Lower.bound.Likelihood <-  cline.center + (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[4]/2) # 271.5092

# Plot the final cline
cline.center <- all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$center
cline.width <- all.ind.k3_Q_model_fixed_noneData$ML.cline$param.free$width
quartz(width = 6, height = 4)
hzar.plot.fzCline(all.ind.k3_Q_model_fixed_noneData, fzCol = "transparent", pch = 20, col="saddlebrown", main = "Genomic cline of all lindividuals (n=108) - Admixture K=3", xlab = "Distance (km)", ylab = "Ancestry proportion of D. psarodes", cex.main=1, cex.axis=0.8, cex.lab=0.8, xlim=c(30, 500))
abline(v=cline.center, lty=2, col="red") # to reprecent center of the cline
abline(v=(cline.center-(cline.width/2)), lty=3, col="saddlebrown") # lower boundary of cline width = starting point of the cline
abline(v=(cline.center+(cline.width/2)), lty=3, col="saddlebrown") # upper boundary of cline width = end point of the cline
rect((cline.center - (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[4])/2), -3.5, (cline.center + (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[4])/2), 4, border = "transparent", col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.1))
rect((hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[1]), -0.5, (hzar.getLLCutParam(all.ind.k3_Q_model_fixed_noneData,c("center","width"))[2]), 1.5, col=rgb(red = 0.95, green = 0.8, blue = 0.3, alpha = 0.3), border = "transparent") ## shows the likelyhood of cline center (CI of the cline center)

############################################################################
###################### END Genomic cline analysis ########################
############################################################################


############################################################################
######################## 3. Phenotypic cline analysis ########################
############################################################################

###### Notes:
# Using the HZAR example R script (men12209-2up-002-datas2.pdf) to conduct this analysis.
# The distances were calculated from Point Calimere (10.2845° N, 79.8241° E) in India to each sampling location using ruler tool in google earth Pro.
# PC1 of the Lab-PCA was considered as the phenotype data for the analysis. (L, a and b values were taken from photoshop --> make a PCA on L,a and b values -> get the PC1 values as phenotype scores)
# samples that were caught from the same location were lumped together into the "Locality" column.

# set working directory
setwd("~/Documents/Dinopium_GBS_on_laptop/012NA_150_samples/Cline_analysis_n.150_2022/Pheno cline/")

# load the requred libraries
library(ggplot2)
library(sf)
library(stringr)
library(raster)
library(sp)
library(splines)
library(dplyr)
library(hzar)
library(viridis)
library(forcats) 


####### Visualize some data before generating the cline #######

# ####### Plot the PCA for LAB values
# lab_PCA <- read.csv("LAB.PCA.sample.data.for.HZAR.n109.csv",stringsAsFactors = T)
# 
# # plot the Lab-PCA
# lab_PCA$Group_phenotype <- factor(lab_PCA$Group_phenotype, levels = c("Dino_R_Allo", "Dino_R_Sym", "Dino_H", "Dino_Y_Sym", "Dino_Y_J", "Dino_Y_M"))
# 
# quartz(width = 6, height = 4)
# ggplot(lab_PCA, aes(x=pc.values.PC1, y=pc.values.PC2, fill=Group_phenotype)) + geom_point(shape=21, size=2, colour="gray25") +
#   scale_fill_manual(values = c("red", "red", "orange", "yellow", "yellow","yellow"))+ theme_test() +
#   xlab("PC1") + ylab("PC2") +
#   labs(title = "PCA based on L, a and b values (n=108)") +
#   theme(plot.title=element_text(lineheight = 5, size = 9))
# 
# 
# ### Plot phenotype scores agianst the geographic location
# ######### load the phenotype data file
# phenotype_data <- read.csv("PhenoCline_sample_data_n=109_new.csv",stringsAsFactors = T) 
# 
# # plot the data
# phenotype_data$Group_phenotype <- factor(phenotype_data$Group_phenotype, levels = c("Dino_R_Allo", "Dino_R_Sym", "Dino_H", "Dino_Y_Sym", "Dino_Y_J", "Dino_Y_M"))
# 
# quartz(width = 6, height = 4)
# ggplot(phenotype_data, aes(x=googleearth.dist_km, y=pc.values.PC1, fill=Group_phenotype)) + geom_point(shape=21, size=2, colour="gray25") + 
#   scale_fill_manual(values = c("red", "red", "orange", "yellow", "yellow","yellow"))+ theme_test() +
#   xlab("Distance (km)") + #ylab("LAB PC1 values") +
#   labs(title = "LAB PC1 values alogn the geographical distance for all individuals (n=108) - (K=3)") + 
#   theme(plot.title=element_text(lineheight = 5, size = 9))
# 
# # Plot the distribution again with a non - linear model filt 
# fit <- loess(pc.values.PC1 ~ googleearth.dist_km, phenotype_data)
# quartz(width = 6, height = 4)
# ggplot(phenotype_data, aes(x=googleearth.dist_km, y=pc.values.PC1, fill=Group_phenotype)) + geom_point(shape=21, size=2, colour="gray25") + 
#   geom_line(aes(googleearth.dist_km, predict(fit))) +
#   scale_fill_manual(values = c("red", "red", "orange", "yellow", "yellow","yellow"))+ theme_test() +
#   xlab("Distance (km)") + ylab("LAB PC1 values") +
#   labs(title = "LAB PC1 values alogn the geographical distance for all individuals (n=108) - (K=3)") + 
#   theme(plot.title=element_text(lineheight = 5, size = 9))


##### load the phenotype data file
phenotype_data <- read.csv("PhenoCline_sample_data_n=109_new.csv",stringsAsFactors = T)

##### Load the Locality data - This data set inclde locations and distance for each location
phenoLocations <- read.csv("PhenoCline_sample_data_n=109_new_LocalityData.csv", header = T) 

## Make the cline
n_108_LAB.PC1 <- hzar.doNormalData1DRaw(hzar.mapSiteDist(phenoLocations$Locality_new,
                                                         phenoLocations$googleearth.dist_km),
                                        phenotype_data$Locality_new,
                                        phenotype_data$pc.values.PC1)


hzar.plot.obsData(n_108_LAB.PC1)

# make the models
n_108_LAB.PC1_model_none <- hzar.makeCline1DNormal(n_108_LAB.PC1, tails="none")
n_108_LAB.PC1_model_right <- hzar.makeCline1DNormal(n_108_LAB.PC1, tails="right")
n_108_LAB.PC1_model_left <- hzar.makeCline1DNormal(n_108_LAB.PC1, tails="left")
n_108_LAB.PC1_model_mirror <- hzar.makeCline1DNormal(n_108_LAB.PC1, tails="mirror")
n_108_LAB.PC1_model_both <- hzar.makeCline1DNormal(n_108_LAB.PC1, tails="both")

## The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.

n_108_LAB.PC1_model_none <- hzar.model.addBoxReq(n_108_LAB.PC1_model_none, -100,750)
n_108_LAB.PC1_model_right <- hzar.model.addBoxReq(n_108_LAB.PC1_model_right, -100,750)
n_108_LAB.PC1_model_left <- hzar.model.addBoxReq(n_108_LAB.PC1_model_left, -100,750)
n_108_LAB.PC1_model_mirror <- hzar.model.addBoxReq(n_108_LAB.PC1_model_mirror, -100,750)
n_108_LAB.PC1_model_both <- hzar.model.addBoxReq(n_108_LAB.PC1_model_both, -100,750)

##### I did not set the model setting but in the example it has been set up

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines
n_108_LAB.PC1_model_noneFitR <- hzar.first.fitRequest.gC(gModel=n_108_LAB.PC1_model_none, obsData=n_108_LAB.PC1, verbose = F)
n_108_LAB.PC1_model_rightFitR <- hzar.first.fitRequest.gC(gModel=n_108_LAB.PC1_model_right, obsData=n_108_LAB.PC1, verbose = F)
n_108_LAB.PC1_model_leftFitR <- hzar.first.fitRequest.gC(gModel=n_108_LAB.PC1_model_left, obsData=n_108_LAB.PC1, verbose = F)
n_108_LAB.PC1_model_mirrorFitR <- hzar.first.fitRequest.gC(gModel=n_108_LAB.PC1_model_mirror, obsData=n_108_LAB.PC1, verbose = F)
n_108_LAB.PC1_model_bothFitR <- hzar.first.fitRequest.gC(gModel=n_108_LAB.PC1_model_both, obsData=n_108_LAB.PC1, verbose = F)

## set mcmc chain length and buRAllo for each model fit
n_108_LAB.PC1_model_noneFitR$mcmcParam$chainLength <- 1e5 # this is the default value
n_108_LAB.PC1_model_noneFitR$mcmcParam$burnin <- 5e5
n_108_LAB.PC1_model_rightFitR$mcmcParam$chainLength <- 1e5 
n_108_LAB.PC1_model_rightFitR$mcmcParam$burnin <- 5e5
n_108_LAB.PC1_model_leftFitR$mcmcParam$chainLength <- 1e5 
n_108_LAB.PC1_model_leftFitR$mcmcParam$burnin <- 5e5
n_108_LAB.PC1_model_mirrorFitR$mcmcParam$chainLength <- 1e5 
n_108_LAB.PC1_model_mirrorFitR$mcmcParam$burnin <- 5e5
n_108_LAB.PC1_model_bothFitR$mcmcParam$chainLength <- 1e5 
n_108_LAB.PC1_model_bothFitR$mcmcParam$burnin <- 5e5

# ### Run the optimizer using the parameters listed in the hzar.fitRequest given.
n_108_LAB.PC1_model_noneFit <- hzar.doFit(n_108_LAB.PC1_model_noneFitR)
n_108_LAB.PC1_model_rightFit <- hzar.doFit(n_108_LAB.PC1_model_rightFitR)
n_108_LAB.PC1_model_leftFit <- hzar.doFit(n_108_LAB.PC1_model_leftFitR)
n_108_LAB.PC1_model_mirrorFit <- hzar.doFit(n_108_LAB.PC1_model_mirrorFitR)
n_108_LAB.PC1_model_bothFit <- hzar.doFit(n_108_LAB.PC1_model_bothFitR) # Fitting failed here

# ### plot model to look for run stability and convergence and retuRAllo the mcmc data with an added a log likelihood column.
quartz()
par(mar=c(1,1,1,1))
plot(hzar.mcmc.bindLL(n_108_LAB.PC1_model_noneFit))
plot(hzar.mcmc.bindLL(n_108_LAB.PC1_model_rightFit))
plot(hzar.mcmc.bindLL(n_108_LAB.PC1_model_leftFit))
plot(hzar.mcmc.bindLL(n_108_LAB.PC1_model_mirrorFit))
plot(hzar.mcmc.bindLL(n_108_LAB.PC1_model_bothFit))

# 
##group multiple fits of the same model and the same observation data into a single object
n_108_LAB.PC1_model_noneData <- hzar.dataGroup.add(n_108_LAB.PC1_model_noneFit)
n_108_LAB.PC1_model_rightData <- hzar.dataGroup.add(n_108_LAB.PC1_model_rightFit)
n_108_LAB.PC1_model_leftData <- hzar.dataGroup.add(n_108_LAB.PC1_model_leftFit)
n_108_LAB.PC1_model_mirrorData <- hzar.dataGroup.add(n_108_LAB.PC1_model_mirrorFit)
n_108_LAB.PC1_model_bothData <- hzar.dataGroup.add(n_108_LAB.PC1_model_bothFit)

##make list of cline models
n_108_LAB.PC1_dGs <- list(cline_nonModel=n_108_LAB.PC1_model_noneData,
                          cline_rightModel=n_108_LAB.PC1_model_rightData,
                          cline_leftModel=n_108_LAB.PC1_model_leftData,
                          cline_mirrorModel=n_108_LAB.PC1_model_mirrorData,
                          cline_bothModel=n_108_LAB.PC1_model_bothData)

##Collect optimizer output based on the same hzar.obsData object
n_108_LAB.PC1_oDG <- hzar.make.obsDataGroup(n_108_LAB.PC1_dGs)

##Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
n_108_LAB.PC1_oDG <- hzar.copyModelLabels(n_108_LAB.PC1_dGs,n_108_LAB.PC1_oDG)

##Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(n_108_LAB.PC1_oDG)

##Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(n_108_LAB.PC1_oDG))
#                     AICc
# cline_nonModel    182.1636
# cline_rightModel  233.2423
# cline_leftModel   234.8826
# cline_mirrorModel 237.1733
# cline_bothModel   187.1353

### make em purty, get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
## Lower the AIC value is better the model fit is
min(hzar.AICc.hzar.obsDataGroup(n_108_LAB.PC1_oDG)) # [1] 182.1636 - cline_nonModel =n_108_LAB.PC1_model_noneData

cline.center <- n_108_LAB.PC1_model_noneData$ML.cline$param.free$center
print(n_108_LAB.PC1_model_noneData$ML.cline$param.free$center) # Centre of the cline = 184.7431
cline.width <- n_108_LAB.PC1_model_noneData$ML.cline$param.free$width
print(n_108_LAB.PC1_model_noneData$ML.cline$param.free$width) # Width of the cline = 167.6345
print(hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))) # get the boundaries
#      center2LLLow center2LLHigh width2LLLow width2LLHigh
# 1     142.1823      209.2158    76.93385     256.8465

### Caclculate likelyhoods of clide widths
Northern.border.Upper.bound.Likelihood <- cline.center - (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[4]/2) # 56.31983
Northern.border.of.cline.width <- (cline.center-(cline.width/2)) # 100.9258
Northern.border.Lower.bound.Likelihood <- cline.center - (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[3]/2) # 146.2762
Southern.border.Upper.bound.Likelihood <- cline.center + (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[3]/2) # 223.21
Southern.border.of.cline.width <- (cline.center+(cline.width/2)) # 268.5603
Sourthern.border.Lower.bound.Likelihood <-  cline.center + (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[4]/2) # 313.1663


# Plot the final cline
cline.center <- n_108_LAB.PC1_model_noneData$ML.cline$param.free$center
cline.width <- n_108_LAB.PC1_model_noneData$ML.cline$param.free$width
quartz(width = 6, height = 4)
hzar.plot.fzCline(n_108_LAB.PC1_model_noneData, fzCol = "transparent", pch = 20, col="maroon4",main = "Phenotypic cline of all lindividuals (n=108)-googleeath distances",  xlab = "Distance (km)", ylab = "Plumage score", cex.main=1, cex.axis=0.8, cex.lab=0.8, xlim=c(30, 500))
abline(v=cline.center, lty=2, col="red") # to reprecent center of the cline
abline(v=(cline.center-(cline.width/2)), lty=3, col="maroon4") # lower boundary of cline width = starting point of the cline
abline(v=(cline.center+(cline.width/2)), lty=3, col="maroon4") # upper boundary of cline width = end point of the cline
rect((cline.center - (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[4])/2), -3.5, (cline.center + (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[4])/2), 4, border = "transparent", col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.1))
rect((hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[1]), -3.5, (hzar.getLLCutParam(n_108_LAB.PC1_model_noneData,c("center","width"))[2]), 4, border = "transparent", col=rgb(red = 0.95, green = 0.6, blue = 0.95, alpha = 0.2))

############################################################################
######################## END Phenotypic cline analysis ########################
############################################################################

