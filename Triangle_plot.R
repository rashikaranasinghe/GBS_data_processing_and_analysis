## Install the packages
#devtools::install_github("omys-omics/triangulaR")

# load the package
library(triangulaR)
library(vcfR)

##### Read the data
#data <- read.csv(file = 'Interspecific.heterozygosity_II.csv', stringsAsFactors = T)
data <- read.vcfR("DW.chrom_150.n_Oct.2022.genotypes.variant_only.whole_genome.noIndels_biallelic.GATKfiltered.maxmiss60.mac3.70pcntIndMiss.SNP30miss.recode.vcf", verbose = F)

# make a popmap
info <- read.table("Populaiton_info.txt", header = F)
popmap <- data.frame(colnames(data@gt) ,info[,2])
popmap <- data.frame(colnames(data@gt)[-1] ,info[-1,2])
colnames(popmap) <- c( "id", "pop")
rownames(popmap) <- NULL # Resetting the row names of popmap
head(popmap)

# Choose sites above an allele frequency difference threshold
# choose fixed differences between species
fixed.snps <- alleleFreqDiff(vcfR = data, pm=popmap, p1 = "D.psarodes", p2 = "D.b.jaffnense", difference = 0.7)
# for 0.8 - [1] "4 sites passed allele frequency difference threshold"
# For 0.7 - [1] "8 sites passed allele frequency difference threshold"
# for 0.6 - [1] "35 sites passed allele frequency difference threshold"
# for 0.5 - [1] "263 sites passed allele frequency difference threshold"

#  Calculate hybrid index and heterozygosity for each sample
hi.het <- hybridIndex(vcfR = fixed.snps, pm = popmap, p1 = "D.psarodes", p2 = "D.b.jaffnense")

# make the triangle plt
# Generate colors manually
cols <- c("white", "white", "yellow", "red3", "orange")

# View triangle plot
quartz()
triangle.plot(hi.het, colors = cols)

info <- read.csv("Pop_info.csv", stringsAsFactors = T)
hi.het$phenotype <- info$Group_phenotype
hi.het$Red_ancestry <- info$q_Red

# make the plot
library(dplyr)
quartz()
ggplot(filter(hi.het, pop != "Chrys_R", pop != "Chrys_Y"), aes(y = heterozygosity , x = hybrid.index, fill = phenotype)) +
  geom_polygon(data = data.frame(x = c(0, 1, 0.5), y = c(0, 0, 1)), 
               aes(x = x, y = y), 
               fill = NA, 
               color = "grey60", 
               size = 0.4) +
  stat_function(fun = function(x) -2 * (x - 0.5)^2 + 0.5, color = "grey60", size = 0.4,linetype = "dashed") +
  geom_point(shape=23, size=2.5) +
 #scale_fill_manual(values = c("yellow", "red3", "orange")) + for pop
scale_fill_manual(values = c("orange", "red3", "red3", "yellow", "yellow", "yellow")) + # for phenotype
# scale_fill_manual(values = c("orange", "red3", "red3", "green3", "blue", "yellow")) + # for phenotype
   ylim(c(0, 1)) +
 labs(x = "Hybrid index",  y = "Interspecific heterozygosity (0.7)") + 
   theme_test() + theme(legend.position = "left")

# how to add a curved line starting form 0,0 and end at 1,0 and the highest plot is at 0.5,0.5 on this plot


# Color triangle plot by missing data
quartz()
missing.plot(hi.het)
