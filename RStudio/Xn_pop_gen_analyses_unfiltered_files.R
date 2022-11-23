---
title: "Population genetic analyses of **Xylaria necrophora**"
author: "Teddy Garcia-Aroca"
date: 'Last update: October 31st, 2022'
output:
  html_document: rmdformats::downcute
pdf_document: default
---
  
# Description.

#Complete analyses for X. necrophora 153 genomes
#Using code from https://grunwaldlab.github.io/Population_Genetics_in_R/
#And some phylogeography analyses based on Shakya et. al. 2021 (Phylogeography of P. cinnanomi)

#Install packages needed
# Package names
packages <- c("adegenet", "ape", "dartR", "cowplot", "devtools", "dplyr", "ggplot2", "hierfstat", "igraph", "knitr", "lattice", "magrittr", "mmod", "pegas", "pinfsc50", "poppr", "RColorBrewer", "reshape2", "treemap", "vcfR", "vegan", "viridisLite")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Install SNPRelate dependency with BiocManager

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.15")

#BiocManager::install("SNPRelate")

#library("devtools")
#install_github("zhengxwen/SNPRelate",force = TRUE)

#SNPRelate can be tricky to install on a mac. If none of that works, then try these steps:

#1. Install homebrew
#Open the Terminal app and type:
 # > ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
#2. Install Command-Line tools
#In the Terminal app type:
 # > xcode-select --install
#3. Install compiler
#In the Terminal app type:
  #> brew install gcc
#4. Make a symlink from your gcc version to gcc
#In the Terminal app type:
  #> ln -s gcc-11 gcc
#5. Create a makevars file and put the path to libgfortran into FLIBS
#In the R console type 
#> dir.create('~/.R')
#> write.table("FLIBS=`gfortran -print-search-dirs | grep '^libraries:' | sed 's|libraries: =||' | sed 's|:| -L|g' | sed 's|^|-L|'`",file='~/.R/Makevars',col.names = FALSE,row.names = FALSE,quote = FALSE)
#6. Install SNPRelate
#In the R console type:
#  > install.packages("SNPRelate")
#7. Install dartR
#In the R console type:
 # > install.packages("dartR")
#BiocManager::install("SNPRelate")

#Load specific packages needed for analyses (If preferred)
library("poppr")
library(dartR)
library(mmod)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)


#vcfR
#Using the R package vcfR, we can read VCF format files into memory using the function read.vcfR(). Once in memory we can use the head() method to summarize the information in the three VCF regions.

#Personal laptop
#setwd("~/Documents/SUNY_Albany/TGA_Postdoc/TGA_PhD/Xylaria_WGS/R_analyses_X_157_genomes/")

#FEL laptop
setwd("~/Dropbox/work/TGA_PhD/Xylaria_WGS_2023/RStudio/")

vcf <- read.vcfR("~/Desktop/work_fungeco/phd/X.necrophora.popgen/necrophora_153_filtering/homozygous_only/Xn_153_genomes.MinQG40.MinDP5.MaxMiss0.7.recode.Homozygous.recode.vcf")
head(vcf)

#Load population data (175 genomes)
pop.data <- read.table("../../Xylaria_WGS/metadata/Xylaria_175_genomes.popdata.Updated.05.29.2022.txt", sep = "\t", header = TRUE)

#Check that all the data is contained in our pop.data matches the data contained in our vcf
all(colnames(vcf@gt)[-1] == pop.data$No)

#Convert to genind
my_genind <- vcfR2genind(vcf)
#Individuals with no scored loci have been removed, when converting to genind

#Check new number of individuals
nInd(my_genind)

#Print new names for individuals that remained after conversion
adegenet::indNames(my_genind)

#Extract labels as vector
newNames <- as.vector(adegenet::indNames(my_genind))

#Subset those from original pop data
pop.data2 <- pop.data[pop.data$ID %in% newNames,]

#Add ploidy to our genind
ploidy(my_genind) <- 1

#Add new strata to genind
strata(my_genind) <- pop.data2
setPop(my_genind) <- ~Species_Type

#Convert to genclone object
my_genclone <- poppr::as.genclone(my_genind)

#Using hypothesized pop2 in our pop.data
#pop(my_genclone) <- pop.data2$Hyp_Pop2
ploidy(my_genclone) <- 1

#Subsetting data (removing historical)
setPop(my_genclone) <- ~Species_Type
my_genclone.subset <- popsub(my_genclone, sublist = "X.necrophora_Contemp", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

#Subset to remove populations with less than 10 ind (TN is the threshold)
setPop(my_genclone.subset) <- ~State
my_genclone.subset.2 <- popsub(my_genclone.subset, sublist = c("AL","AR","LA", "MS","TN"), exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

#Removing 6 X. necrophora isolates because they showed a lot of missing data
#setPop(my_genind.subset) <- ~ID
setPop(my_genclone) <- ~ID
my_genclone.subset.3 <- popsub(my_genclone, sublist = c("DMCC1467","DMCC1468","DMCC1469","DMCC1485","DMCC1486","DMCC1487","DMCC1509","DMCC2089","DMCC2091","DMCC2097","DMCC2098","DMCC2099","DMCC2100","DMCC2101","DMCC2102","DMCC2103","DMCC2104","DMCC2105","DMCC2106","DMCC2108","DMCC2109","DMCC2110","DMCC2111","DMCC2113","DMCC2115","DMCC2116","DMCC2117","DMCC2118","DMCC2119","DMCC2120","DMCC2121","DMCC2122","DMCC2123","DMCC2124","DMCC2125","DMCC2126","DMCC2127","DMCC2165","DMCC2166","DMCC2456","DMCC2458","DMCC2463","DMCC2470","DMCC2471","DMCC2477","DMCC2478","DMCC2514","DMCC2515","DMCC2517","DMCC2524","DMCC2525","DMCC2526","DMCC2528","DMCC2529","DMCC2530","DMCC2531","DMCC2532","DMCC2533","DMCC2534","DMCC2615","DMCC2616","DMCC2617","DMCC2618","DMCC2619","DMCC2620","DMCC2621","DMCC2622","DMCC2624","DMCC2625","DMCC3164","DMCC3165","DMCC3166","DMCC3168","DMCC3171","DMCC3173","DMCC3175","DMCC3176","DMCC3177","DMCC3180","DMCC3182","DMCC3183","DMCC3186","DMCC3205","DMCC3222","DMCC3229","DMCC3230","DMCC3247","DMCC3255","DMCC3260","DMCC3267","DMCC3270","DMCC3273","DMCC3274","DMCC3275","DMCC3276","DMCC3282","DMCC3285","DMCC3289","DMCC3291","DMCC3299","DMCC3300","DMCC3404","DMCC3405","DMCC3412","DMCC3413","DMCC3414","DMCC3415","DMCC3416","DMCC3417","DMCC3418","DMCC3419","DMCC3420","DMCC3423","DMCC3424","DMCC3425","DMCC3426","DMCC3427","DMCC3429","DMCC3430","DMCC3431","DMCC3432","DMCC3433","DMCC3434","DMCC3435","DMCC3436","DMCC3437","DMCC3438","DMCC3439","DMCC3440","DMCC3441","DMCC3442","DMCC3443","DMCC3444","DMCC3445","DMCC3446","DMCC3448","DMCC3449","DMCC3450","DMCC3451","DMCC3452","DMCC3453","DMCC3455","DMCC3456","DMCC3457","DMCC3458","DMCC3459","DMCC3460","DMCC3461","DMCC3463","DMCC3464","DMCC3829","DMCC3852","DMCC3861","DMCC3866"), exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

#Filtering datasets for missing data
my_genclone.filtered <- missingno(my_genclone.subset, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE)
my_genclone.subset.filtered <- missingno(my_genclone.subset, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE)
my_genclone.subset.2.filtered <- missingno(my_genclone.subset.2, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE)
my_genclone.subset.3.filtered <- missingno(my_genclone.subset.3, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE)

#Convert to genlight and genind
my_genind <- genclone2genind(my_genclone)
my_genind.filtered <- genclone2genind(my_genclone.filtered)
my_genind.subset <- genclone2genind(my_genclone.subset)
my_genind.subset.filtered <- genclone2genind(my_genclone.subset.filtered)
my_genind.subset.2 <- genclone2genind(my_genclone.subset.2)
my_genind.subset.2.filtered <- genclone2genind(my_genclone.subset.2.filtered)
my_genind.subset.3 <- genclone2genind(my_genclone.subset.3)
my_genind.subset.3.filtered <- genclone2genind(my_genclone.subset.3.filtered)

gl.necrophora <- gi2gl(my_genind)
gl.necrophora.filtered <- gi2gl(my_genind.filtered)
gl.necrophora.subset.filtered <- gi2gl(my_genind.subset.filtered)
gl.necrophora.subset <- gi2gl(my_genind.subset)
gl.necrophora.subset.2 <- gi2gl(my_genind.subset.2)
gl.necrophora.subset.2.filtered <- gi2gl(my_genind.subset.2.filtered)
gl.necrophora.subset.3 <- gi2gl(my_genind.subset.3)
gl.necrophora.subset.3.filtered <- gi2gl(my_genind.subset.3.filtered)

#Add pop data to genlight objects.
names.filtered <- gl.necrophora$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.unfiltered <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora) <- pop.data.unfiltered
setPop(gl.necrophora) <- ~State

#Set pop for subset
names.filtered <- gl.necrophora.filtered$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.filtered <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.filtered) <- pop.data.filtered
setPop(gl.necrophora.filtered) <- ~State

#Set pop for subset
names.filtered <- gl.necrophora.subset$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset) <- pop.data.subset
setPop(gl.necrophora.subset) <- ~State

#Set pop for subset
names.filtered <- gl.necrophora.subset.filtered$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset.filtered <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset.filtered) <- pop.data.subset.filtered
setPop(gl.necrophora.subset.filtered) <- ~State

#Set pop for subset 2
names.filtered <- gl.necrophora.subset.2$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset.2 <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset.2) <- pop.data.subset.2
setPop(gl.necrophora.subset.2) <- ~State

#Set pop for subset 2 filtered
names.filtered <- gl.necrophora.subset.2.filtered$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset.2.filtered <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset.2.filtered) <- pop.data.subset.2.filtered
setPop(gl.necrophora.subset.2.filtered) <- ~State

#Export as csv
#write.csv(pop.data.subset.2.filtered, "~/Dropbox/work/TGA_PhD/Xylaria_WGS/raxml_Xn_156_genomes/metadata.pop.data.subset.2.filtered.csv")


#Set pop for subset 3
names.filtered <- gl.necrophora.subset.3$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset.3 <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset.3) <- pop.data.subset.3
setPop(gl.necrophora.subset.3) <- ~State

#Set pop for subset 3 filtered
names.filtered <- gl.necrophora.subset.3.filtered$ind.names
newNames <- as.vector(names.filtered)
#Subset those from original pop data
pop.data.subset.3.filtered <- pop.data[pop.data$ID %in% newNames,]
strata(gl.necrophora.subset.3.filtered) <- pop.data.subset.3.filtered
setPop(gl.necrophora.subset.3.filtered) <- ~State


#Distance tree
#Let’s start our analysis by building a genetic distance tree that represents the genetic relatedness of the samples. The similarity between samples and groups of samples is represented by the branch length. In most trees, the branch length is represented by the number of substitutions per site for a cluster or a sample. When samples are very similar, they are grouped by short branches. The longer the branch, the higher the number of substitutions and the higher the genetic distance is between samples or clusters.
gl.list <- list(gl.necrophora, gl.necrophora.filtered, gl.necrophora.subset, gl.necrophora.subset.2, gl.necrophora.subset.2.filtered, gl.necrophora.subset.3, gl.necrophora.subset.3.filtered)

#######Run function for unrooted trees
lapply(gl.list, function (x){
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})

############Run function for rooted trees
lapply(gl.list, function (x){
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T, root = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})


###########Run same analyses by lineages determined by ADM and Structure
#######Run function for unrooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hy_PopAdmixture1
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})

############Run function for rooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hy_PopAdmixture1
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T, root = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})


###########Run same analyses by lineages determined by ADM and Structure and post assessment after PCA
#######Run function for unrooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Admix_postPCA
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})

############Run function for rooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Admix_postPCA
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T, root = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})


###########Run same analyses by lineages determined by ADM and Structure and post assessment after PCA
########### Adding both lineages and admixed genotypes
#######Run function for unrooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Admix_postPCA_Final
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})


############Run function for rooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Admix_postPCA_Final
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T, root = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})



###########Run same analyses by lineages determined by ADM and Structure and post assessment after PCA
####FINAL pops after removing outliers
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Final
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})


############Run function for rooted trees
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Final
  tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T, root = T)
  
  #Next, we will color the tips of the tree based on the population of origin of the samples, and draw conclusions from what we observe in the tree:
  #Plotting by Year by state
  cols <- c(brewer.pal(n = nPop(x), name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 6))
  #cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
  
  plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color =  cols[pop(x)])
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 0.2, xpd = TRUE)
  #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
  legend('topleft', legend = as.vector(unique(x$pop)), fill = cols, border = FALSE, bty = "n", cex = 0.5)
  axis(side = 1)
  title(xlab = "Genetic distance (proportion of loci that are different)")
})



###########Plot happlotype network by Lineage
lapply(gl.list, function (x){
  setPop(x) <- ~Hy_PopAdmixture1
  necrophora.dist <- bitwise.dist(x)
  necrophora.msn <- poppr.msn(x, necrophora.dist, showplot = FALSE, include.ties = T)
  
  node.size <- rep(2, times = nInd(x))
  #names(node.size) <- indNames(x)
  #vertex.attributes(necrophora.msn$graph)$size <- node.size
  
  set.seed(9)
  plot_poppr_msn(x, necrophora.msn , palette = brewer.pal(n = 8, name = "Dark2"), gadj = 70)
})



###########Plot happlotype network by Lineage after PCA assessment
lapply(gl.list, function (x){
  setPop(x) <- ~~Hyp_Admix_postPCA
  necrophora.dist <- bitwise.dist(x)
  necrophora.msn <- poppr.msn(x, necrophora.dist, showplot = FALSE, include.ties = T)
  
  node.size <- rep(2, times = nInd(x))
  #names(node.size) <- indNames(x)
  #vertex.attributes(necrophora.msn$graph)$size <- node.size
  
  set.seed(9)
  plot_poppr_msn(x, necrophora.msn , palette = brewer.pal(n = 8, name = "Dark2"), gadj = 70)
})


###########Plot happlotype network by Lineage after PCA assessment and including admixed
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Admix_postPCA_Final
  necrophora.dist <- bitwise.dist(x)
  necrophora.msn <- poppr.msn(x, necrophora.dist, showplot = FALSE, include.ties = T)
  
  node.size <- rep(2, times = nInd(x))
  #names(node.size) <- indNames(x)
  #vertex.attributes(necrophora.msn$graph)$size <- node.size
  
  set.seed(9)
  plot_poppr_msn(x, necrophora.msn , palette = brewer.pal(n = 8, name = "Dark2"), gadj = 70)
})


###########Plot happlotype network by Lineage after PCA assessment and including admixed (FINAL)
lapply(gl.list, function (x){
  setPop(x) <- ~Hyp_Final
  necrophora.dist <- bitwise.dist(x)
  necrophora.msn <- poppr.msn(x, necrophora.dist, showplot = FALSE, include.ties = T)
  
  node.size <- rep(2, times = nInd(x))
  #names(node.size) <- indNames(x)
  #vertex.attributes(necrophora.msn$graph)$size <- node.size
  
  set.seed(9)
  plot_poppr_msn(x, necrophora.msn , palette = brewer.pal(n = 8, name = "Dark2"), gadj = 70)
})


###Run and plot PCA and DAPC analyses (Unfiltered and filtered subset 2 only)

#######################Visualizing K- means clustering
####If you peak at the bottom of this document you’ll see that our goal is a multi-panel ggplot. Each panel will be a different ggplot object, so we’ll have to give them unique names as we make them.
#setPop(gl.necrophora.subset.2.filtered) <- ~Hy_PopAdmixture1

setPop(my_genclone) <- ~Hy_PopAdmixture1

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(my_genclone, n.pca = 20, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}


#Plot the results from above
library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


##DAPC
#In general, it is recommended to explore several values of K. The find.clusters() function includes some stochasticity. When we’re at the figure creation step we’ll need consistency, so we’ll set a seed. If you’re at an earlier stage in your analysis you should comment the set.seed() call out to explore how sensitive your results are to the seed.
my_k <- 2:3

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(my_genclone, n.pca = 20, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(my_genclone, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
  #dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

#A nice perspective is to create a scatterplot based on the discriminant functions. This helps us see how diffferent the resulting clusters are and if we may have chosen too high of a value for K.

#LD1 and LD2 are discriminant functions
my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)


my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#The above palette only takes 8 colors maximum
#create a new palette:
#my_pal <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))

#Plot LD 1 and LD2
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#Plot LD2 and LD3
p2 <- ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2


#Run PCA for the above
setPop(gl.necrophora) <-~Hy_PopAdmixture1

necrophora.pca <- glPca(gl.necrophora, nf = 3)
barplot(100*necrophora.pca$eig/sum(necrophora.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

necrophora.pca.scores <- as.data.frame(necrophora.pca$scores)
necrophora.pca.scores$pop <- pop(gl.necrophora)

cols <- brewer.pal(n = nPop(gl.necrophora), name = "Dark2")

library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC1, y=PC3, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC2, y=PC3, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p

s.class(necrophora.pca$scores, pop(gl.necrophora), col=colors()[c(131,134)])
add.scatter.eig(necrophora.pca$eig,4,1,2)


#Extraction values from PCA
#pca<-gl.pcoa(gl.necrophora.subset,nfactors=5)
gl.pcoa.plot(necrophora.pca, gl.necrophora, ellipse=TRUE, p=0.99, labels="pop",hadjust=1.5,
             vadjust=1)
if (requireNamespace("plotly", quietly = TRUE)) {
  #interactive plot to examine labels
  gl.pcoa.plot(necrophora.pca, gl.necrophora, labels="interactive")
}



#######################Visualizing K- means clustering
####If you peak at the bottom of this document you’ll see that our goal is a multi-panel ggplot. Each panel will be a different ggplot object, so we’ll have to give them unique names as we make them.
setPop(gl.necrophora.subset.2.filtered) <- ~Hy_PopAdmixture1

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

#Plot the results from above
library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


##DAPC
#In general, it is recommended to explore several values of K. The find.clusters() function includes some stochasticity. When we’re at the figure creation step we’ll need consistency, so we’ll set a seed. If you’re at an earlier stage in your analysis you should comment the set.seed() call out to explore how sensitive your results are to the seed.
my_k <- 2:3

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
  #dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

#A nice perspective is to create a scatterplot based on the discriminant functions. This helps us see how diffferent the resulting clusters are and if we may have chosen too high of a value for K.

#LD1 and LD2 are discriminant functions
my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)


my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#The above palette only takes 8 colors maximum
#create a new palette:
#my_pal <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))

#Plot LD 1 and LD2
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#Plot LD2 and LD3
p2 <- ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#A principal components analysis (PCA) converts the observed SNP data into a set of values of linearly uncorrelated variables called principal components that summarize the variation between samples. We can perform a PCA on our genlight object by using the glPCA function.

setPop(gl.necrophora.subset.2.filtered) <-~Hy_PopAdmixture1
necrophora.pca <- glPca(gl.necrophora.subset.2.filtered, nf = 3)
barplot(100*necrophora.pca$eig/sum(necrophora.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#The barplot indicates that we will need to only retain the first 3 PCAs, which cumulatively explain explain 63.183 percent of the variance of the data.
#To view the results of the PCA we can use the package ggplot2. We need to convert the data frame that contains the principal components (necrophora.pca$scores) into the new object necrophora.pca.scores. In addition, we will add the population values as a new column in our necrophora.pca.scores object, in order to be able to color samples by population.

#ggplot2 will plot the PCA, color the samples by population, and create ellipses that include 95% of the data for each the population:

necrophora.pca.scores <- as.data.frame(necrophora.pca$scores)
necrophora.pca.scores$pop <- pop(gl.necrophora.subset.2.filtered)

cols <- brewer.pal(n = nPop(gl.necrophora.subset.2.filtered), name = "Dark2")

library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


s.class(necrophora.pca$scores, pop(gl.necrophora.subset.2.filtered), col=colors()[c(131,134)])
add.scatter.eig(necrophora.pca$eig,2,1,2)

#Extraction values from PCA
#pca<-gl.pcoa(gl.necrophora.subset,nfactors=5)
gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, ellipse=TRUE, p=0.99, labels="pop",hadjust=1.5,
             vadjust=1)
if (requireNamespace("plotly", quietly = TRUE)) {
  #interactive plot to examine labels
  gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, labels="interactive")
}

#DAPC
pnw.dapc <- dapc(gl.necrophora.subset.2.filtered, n.pca = 3, n.da = 2)

#To confirm that the DAPC is similar to the PCA we can plot the data in a scatter plot.

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)


#We see that the results of the PCA and DAPC are very similar. The DAPC object we created includes the population membership probability for each sample to each of the predetermined populations. To visualize the posterior assignments of each sample, we use a composite stacked bar plot (compoplot). A compoplot illustrates the probability of population membership on the y-axis. Each sample is a bin on the x-axis, and the assigned probability of population membership is shown as a stacked bar chart with clusters or populations shown in color.
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'bottom')


#These plots are hard to interpret and we will thus separate the samples by population.
#ggplot2 can be used to reconstruct these plots, but we need to convert the data into a ggplot2 friendly object. We will extract the DAPC calculated population membership assignments (pnw.dapc$posterior) into a new data frame (dapc.results), include the original population assignment as a new column in the data frame (dapc.results$pop), and add a column that includes the sample names (dapc.results$indNames).

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.necrophora.subset.2.filtered)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 has specific requirements for the structure of the data frame format, as it requires each observation in rows, and all different values of these observations in columns (i.e., a long format data frame). To transform the data frame we use pivot_longer from the package tidyr.

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#We have now reorganized the data frame into the required format, where each membership probability observation for a given population is a row with the sample name, original population, and assigned population as columns.

head(dapc.results, n = 6)

#Then, we rename the columns into more familiar terms:
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#library(devtools)
#install_github("ronammar/randomcoloR")
#library(randomcoloR)
#n <- 33
#palette <- distinctColorPalette(n)

#cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
cols <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#ggplot2 will plot the dapc.results data frame we reorganized using pivot_longer, using the samples on the X-axis and membership probabilities on the Y-axis. The fill color will indicate the original population assignments. Each facet represents the original population assignment for each sample:
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))
p


###Do it again with 3 populations
setPop(gl.necrophora.subset.2.filtered) <- ~Hyp_Admix_postPCA_Final

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

#Plot the results from above
library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


##DAPC
#In general, it is recommended to explore several values of K. The find.clusters() function includes some stochasticity. When we’re at the figure creation step we’ll need consistency, so we’ll set a seed. If you’re at an earlier stage in your analysis you should comment the set.seed() call out to explore how sensitive your results are to the seed.
my_k <- 2:3

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
  #dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

#A nice perspective is to create a scatterplot based on the discriminant functions. This helps us see how diffferent the resulting clusters are and if we may have chosen too high of a value for K.

#LD1 and LD2 are discriminant functions
my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)


my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#The above palette only takes 8 colors maximum
#create a new palette:
#my_pal <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))

#Plot LD 1 and LD2
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#Plot LD2 and LD3
p2 <- ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#A principal components analysis (PCA) converts the observed SNP data into a set of values of linearly uncorrelated variables called principal components that summarize the variation between samples. We can perform a PCA on our genlight object by using the glPCA function.

setPop(gl.necrophora.subset.2.filtered) <-~Hyp_Admix_postPCA_Final
necrophora.pca <- glPca(gl.necrophora.subset.2.filtered, nf = 3)
barplot(100*necrophora.pca$eig/sum(necrophora.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#The barplot indicates that we will need to only retain the first 3 PCAs, which cumulatively explain explain 63.183 percent of the variance of the data.
#To view the results of the PCA we can use the package ggplot2. We need to convert the data frame that contains the principal components (necrophora.pca$scores) into the new object necrophora.pca.scores. In addition, we will add the population values as a new column in our necrophora.pca.scores object, in order to be able to color samples by population.

#ggplot2 will plot the PCA, color the samples by population, and create ellipses that include 95% of the data for each the population:

necrophora.pca.scores <- as.data.frame(necrophora.pca$scores)
necrophora.pca.scores$pop <- pop(gl.necrophora.subset.2.filtered)

cols <- brewer.pal(n = nPop(gl.necrophora.subset.2.filtered), name = "Dark2")

library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.99, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


s.class(necrophora.pca$scores, pop(gl.necrophora.subset.2.filtered), col=colors()[c(131,134)])
add.scatter.eig(necrophora.pca$eig,3,1,2)

#Extraction values from PCA
#pca<-gl.pcoa(gl.necrophora.subset,nfactors=5)
gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, ellipse=TRUE, plevel=0.95,hadjust=1.5,
             vadjust=1)
if (requireNamespace("plotly", quietly = TRUE)) {
  #interactive plot to examine labels
  gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, labels="interactive")
}

#DAPC
pnw.dapc <- dapc(gl.necrophora.subset.2.filtered, n.pca = 3, n.da = 2)

#To confirm that the DAPC is similar to the PCA we can plot the data in a scatter plot.

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)


#We see that the results of the PCA and DAPC are very similar. The DAPC object we created includes the population membership probability for each sample to each of the predetermined populations. To visualize the posterior assignments of each sample, we use a composite stacked bar plot (compoplot). A compoplot illustrates the probability of population membership on the y-axis. Each sample is a bin on the x-axis, and the assigned probability of population membership is shown as a stacked bar chart with clusters or populations shown in color.
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'bottom')


#These plots are hard to interpret and we will thus separate the samples by population.
#ggplot2 can be used to reconstruct these plots, but we need to convert the data into a ggplot2 friendly object. We will extract the DAPC calculated population membership assignments (pnw.dapc$posterior) into a new data frame (dapc.results), include the original population assignment as a new column in the data frame (dapc.results$pop), and add a column that includes the sample names (dapc.results$indNames).

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.necrophora.subset.2.filtered)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 has specific requirements for the structure of the data frame format, as it requires each observation in rows, and all different values of these observations in columns (i.e., a long format data frame). To transform the data frame we use pivot_longer from the package tidyr.

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#We have now reorganized the data frame into the required format, where each membership probability observation for a given population is a row with the sample name, original population, and assigned population as columns.

head(dapc.results, n = 6)

#Then, we rename the columns into more familiar terms:
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#library(devtools)
#install_github("ronammar/randomcoloR")
#library(randomcoloR)
#n <- 33
#palette <- distinctColorPalette(n)

#cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
cols <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#ggplot2 will plot the dapc.results data frame we reorganized using pivot_longer, using the samples on the X-axis and membership probabilities on the Y-axis. The fill color will indicate the original population assignments. Each facet represents the original population assignment for each sample:
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))
p




###############Separate in 3 lineages#######################

###Do it again with 3 populations
setPop(gl.necrophora.subset.2.filtered) <- ~Hyp_Final

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

#Plot the results from above
library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


##DAPC
#In general, it is recommended to explore several values of K. The find.clusters() function includes some stochasticity. When we’re at the figure creation step we’ll need consistency, so we’ll set a seed. If you’re at an earlier stage in your analysis you should comment the set.seed() call out to explore how sensitive your results are to the seed.
my_k <- 2:3

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 20, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
  #dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

#A nice perspective is to create a scatterplot based on the discriminant functions. This helps us see how diffferent the resulting clusters are and if we may have chosen too high of a value for K.

#LD1 and LD2 are discriminant functions
my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)


my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#The above palette only takes 8 colors maximum
#create a new palette:
#my_pal <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))

#Plot LD 1 and LD2
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#Plot LD2 and LD3
p2 <- ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

#A principal components analysis (PCA) converts the observed SNP data into a set of values of linearly uncorrelated variables called principal components that summarize the variation between samples. We can perform a PCA on our genlight object by using the glPCA function.

setPop(gl.necrophora.subset.2.filtered) <-~Hyp_Final
necrophora.pca <- glPca(gl.necrophora.subset.2.filtered, nf = 3)
barplot(100*necrophora.pca$eig/sum(necrophora.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#The barplot indicates that we will need to only retain the first 3 PCAs, which cumulatively explain explain 63.183 percent of the variance of the data.
#To view the results of the PCA we can use the package ggplot2. We need to convert the data frame that contains the principal components (necrophora.pca$scores) into the new object necrophora.pca.scores. In addition, we will add the population values as a new column in our necrophora.pca.scores object, in order to be able to color samples by population.

#ggplot2 will plot the PCA, color the samples by population, and create ellipses that include 95% of the data for each the population:

necrophora.pca.scores <- as.data.frame(necrophora.pca$scores)
necrophora.pca.scores$pop <- pop(gl.necrophora.subset.2.filtered)

cols <- brewer.pal(n = nPop(gl.necrophora.subset.2.filtered), name = "Dark2")

library(ggplot2)
set.seed(9)
p <- ggplot(necrophora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.99, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


s.class(necrophora.pca$scores, pop(gl.necrophora.subset.2.filtered), col=colors()[c(120,131,134)])
add.scatter.eig(necrophora.pca$eig,3,1,3)

#Extraction values from PCA
#pca<-gl.pcoa(gl.necrophora.subset,nfactors=5)
gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, ellipse=TRUE, plevel=0.95,hadjust=1.5,
             vadjust=1)
if (requireNamespace("plotly", quietly = TRUE)) {
  #interactive plot to examine labels
  gl.pcoa.plot(necrophora.pca, gl.necrophora.subset.2.filtered, labels="interactive")
}

#DAPC
pnw.dapc <- dapc(gl.necrophora.subset.2.filtered, n.pca = 3, n.da = 2)

#To confirm that the DAPC is similar to the PCA we can plot the data in a scatter plot.

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)


#We see that the results of the PCA and DAPC are very similar. The DAPC object we created includes the population membership probability for each sample to each of the predetermined populations. To visualize the posterior assignments of each sample, we use a composite stacked bar plot (compoplot). A compoplot illustrates the probability of population membership on the y-axis. Each sample is a bin on the x-axis, and the assigned probability of population membership is shown as a stacked bar chart with clusters or populations shown in color.
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'bottom')


#These plots are hard to interpret and we will thus separate the samples by population.
#ggplot2 can be used to reconstruct these plots, but we need to convert the data into a ggplot2 friendly object. We will extract the DAPC calculated population membership assignments (pnw.dapc$posterior) into a new data frame (dapc.results), include the original population assignment as a new column in the data frame (dapc.results$pop), and add a column that includes the sample names (dapc.results$indNames).

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.necrophora.subset.2.filtered)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 has specific requirements for the structure of the data frame format, as it requires each observation in rows, and all different values of these observations in columns (i.e., a long format data frame). To transform the data frame we use pivot_longer from the package tidyr.

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#We have now reorganized the data frame into the required format, where each membership probability observation for a given population is a row with the sample name, original population, and assigned population as columns.

head(dapc.results, n = 6)

#Then, we rename the columns into more familiar terms:
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#library(devtools)
#install_github("ronammar/randomcoloR")
#library(randomcoloR)
#n <- 33
#palette <- distinctColorPalette(n)

#cols <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(name="Paired", n = 12), brewer.pal(name="Accent", n = 8), brewer.pal(name="Pastel1", n = 9))
cols <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

#ggplot2 will plot the dapc.results data frame we reorganized using pivot_longer, using the samples on the X-axis and membership probabilities on the Y-axis. The fill color will indicate the original population assignments. Each facet represents the original population assignment for each sample:
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))
p


#######Run this analysis with several values of K

setPop(gl.necrophora.subset.2.filtered) <- ~Hyp_Final


maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}


my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)


p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


my_k <- 2:8

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.necrophora.subset.2.filtered, n.pca = 8, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 8, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl.necrophora.subset.2.filtered, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2


p2 <- ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2



tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Lin <- pop.data.subset.2.filtered$Hyp_Final
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Lin <- pop.data.subset.2.filtered$Hyp_Final
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Lin, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3


library("ggpubr")
#tiff('dapc__k3_6_dapc.tiff', width=6.5, height=6.5, units='in', compression='lzw', res=300)
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)
#dev.off()


####Export datasets for iTol

#Extract ID and State

pop.data.subset.2.filtered.State <- dplyr::select(pop.data.subset.2.filtered, "ID", "State") 

#pop.data.subset.2.filtered.State <- select(pop.data.subset.2.filtered, "ID", "State") 

pop.data.subset.2.filtered.State$color <- recode(pop.data.subset.2.filtered.State$State,
        AL="#f1c232", AR="#0000ff", Martinique="#741b47", LA="#ff0000", MO="#f44336", MS="#00ff00", TN="#6a329f")


pop.data.subset.2.filtered.State <- select(pop.data.subset.2.filtered.State, "ID", "color", "State")

#write.csv(pop.data.subset.2.filtered.State, "pop.data.subset.2.filtered.State.iTol.csv", sep = "\t")


#Extract ID and Lineage

pop.data.subset.2.filtered.Lin <- select(pop.data.subset.2.filtered, "ID", "Hyp_Final") 

pop.data.subset.2.filtered.Lin$color <- revalue(pop.data.subset.2.filtered.Lin$Hyp_Final,
                                                  c("Lin1"="#1B9E77", "Lin2"="#D95F02", "Lin3"="#7570B3"))


pop.data.subset.2.filtered.Lin <- select(pop.data.subset.2.filtered.Lin, "ID", "color", "Hyp_Final" )

#write.table(pop.data.subset.2.filtered.Lin, "pop.data.subset.2.filtered.Lin.iTol.txt", sep = "\t")



#Export vcf for analyses outside of R

#remotes::install_github("thierrygosselin/radiator")
#library(radiator)

#radiator::detect_genomic_format(data = "gl.necrophora.subset.2.filtered")

#genomic_converter(gl.necrophora.subset.2.filtered, 
                 # strata = NULL, output = "vcf", 
                 # filename = "X.necrophora.153.genomes.filtered.2.dataset.vcf", 
                #  parallel.core = parallel::detectCores() - 1,
                #  verbose = TRUE)


#write.vcf(gl.necrophora.subset.2.filtered, "X.necrophora.153.genomes.filtered.2.dataset.vcf")

#gl2fasta(gl.necrophora.subset.2.filtered, method=1, outfile='X.necrophora.153.genomes.filtered.2.dataset.fasta',verbose=3)


###Standardized index of association by Lineages
setPop(my_genclone.subset.2.filtered) <- ~Hyp_Final

#Lin 1 Standardized index of association
Lin1.subpop <- popsub(my_genclone.subset.2.filtered, "Lin1")
gg1 <- ia(Lin1.subpop, sample = 999)

#Lin 1 Standardized index of association (Clone corrected)
Lin1.subpop.cloneCorr <- popsub(my_genclone.subset.2.filtered, "Lin1")
ia(Lin1.subpop.cloneCorr, sample = 999)


#Lin 2 Standardized index of association
Lin2.subpop <- popsub(my_genclone.subset.2.filtered, "Lin2")
gg3 <- ia(Lin2.subpop, sample = 999)

#Lin 2 Standardized index of association (Clone corrected)
Lin2.subpop.cloneCorr <- popsub(my_genclone.subset.2.filtered, "Lin2")
ia(Lin2.subpop.cloneCorr, sample = 999)


#Lin 3 Standardized index of association
Lin3.subpop <- popsub(my_genclone.subset.2.filtered, "Lin3")
gg5 <- ia(Lin3.subpop, sample = 999)

#Lin 3 Standardized index of association (Clone corrected)
Lin3.subpop.cloneCorr <- popsub(my_genclone.subset.2.filtered, "Lin3")
ia(Lin3.subpop.cloneCorr, sample = 999)




######Simulate and calculate index of association per lineage

setPop(gl.necrophora.subset.2.filtered) <-  ~Hyp_Final


#Lineage 1
t <- seppop(gl.necrophora.subset.2.filtered)
t1 <- clonecorrect(as.snpclone(t$Lin1))
#t2 <- clonecorrect(as.snpclone(t$Lin2))

set.seed(100)
sex <- glSim(n.ind = nInd(t1), n.snp.nonstruc = ceiling(0.9*nLoc(t1)), n.snp.struc = floor(0.1*nLoc(t1)), ploidy=1, LD=TRUE)
### Structure (clonal pops)
clone <- glSim(nInd(t1), n.snp.nonstruc = floor(0.1*nLoc(t1)), n.snp.struc=ceiling(0.9*nLoc(t1)), ploidy=1, LD = T)
### Semi-clonal 
semi_clone <- glSim(nInd(t1),n.snp.nonstruc = 0.5*nLoc(t1), n.snp.struc= 0.5*nLoc(t1), ploidy=1, LD=T)
### Most-clonal 
most_clone <- glSim(nInd(t1), n.snp.nonstruc = ceiling(nLoc(t1)/3), n.snp.struc=2*nLoc(t1)/3, ploidy=1, LD=T)

## IA sex
ia.sex <- samp.ia(sex,quiet = T, reps = 100, n.snp = 500)
## IA clone
ia.clone <- samp.ia(clone, quiet = T, reps = 100, n.snp = 500)
## IA.semiclone
ia.semi <- samp.ia(semi_clone, quiet = T,reps = 100, n.snp = 500)
## IA.mostclone
ia.most <- samp.ia(most_clone, quiet = T, reps = 100, n.snp = 500)

ia.lin1 <- samp.ia(t1,  reps = 100, quiet = T, n.snp = 500)
#ia.lin2 <- samp.ia(t2,  reps = 100, quiet = T, n.snp = 300)


# Summarizing data frames
d1 <- data.frame(ia.lin1, rep("Lineage1", length(ia.lin1)))
d3 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d4 <- data.frame(ia.clone, rep("clonal", length(ia.clone)))
d5 <- data.frame(ia.semi, rep("most-clonal", length(ia.semi)))
d6 <- data.frame(ia.most, rep("semi-clonal", length(ia.semi)))
colnames(d1) <- c("ia","dset")
#colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")
colnames(d6) <- c("ia","dset")
ia.total <- rbind(d4, d6, d5, d3, d1)
#ia.total <- rbind(d1, d2, d3, d4, d5)

# Normality tests
frames <- list(as.data.frame(d1), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5), as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

# Analysis of variance
anova.ia <- aov(lm(ia ~ dset, ia.total))
library(agricolae)
tukey <- HSD.test(anova.ia, "dset", alpha = 0.001)
tukey
# Kluskal wallis test
#kruskal.test(ia ~ dset, ia.total), trt="dset")
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))

# Plot
gg2 <- ggplot(ia.total,aes(x=reorder(dset,-ia), y=ia,fill=dset)) + geom_boxplot() + 
  xlab("Dataset") + ylab("Index of association") + 
  scale_fill_manual(values=c("#FFFFFF", "#1B9E77", "#C1C1C1", "#6D6D6D", "#1E1E1E"))
gg2

#Lineage 2

#setPop(gl.necrophora.subset.2.filtered) <-  ~Hyp_Final
t <- seppop(gl.necrophora.subset.2.filtered)
t1 <- clonecorrect(as.snpclone(t$Lin2))
#t2 <- clonecorrect(as.snpclone(t$Lin2))

set.seed(100)
sex <- glSim(n.ind = nInd(t1), n.snp.nonstruc = ceiling(0.9*nLoc(t1)), n.snp.struc = floor(0.1*nLoc(t1)), ploidy=1, LD=TRUE)
### Structure (clonal pops)
clone <- glSim(nInd(t1), n.snp.nonstruc = floor(0.1*nLoc(t1)), n.snp.struc=ceiling(0.9*nLoc(t1)), ploidy=1, LD = T)
### Semi-clonal 
semi_clone <- glSim(nInd(t1),n.snp.nonstruc = 0.5*nLoc(t1), n.snp.struc= 0.5*nLoc(t1), ploidy=1, LD=T)
### Most-clonal 
most_clone <- glSim(nInd(t1), n.snp.nonstruc = ceiling(nLoc(t1)/3), n.snp.struc=2*nLoc(t1)/3, ploidy=1, LD=T)

## IA sex
ia.sex <- samp.ia(sex,quiet = T, reps = 100, n.snp = 500)
## IA clone
ia.clone <- samp.ia(clone, quiet = T, reps = 100, n.snp = 500)
## IA.semiclone
ia.semi <- samp.ia(semi_clone, quiet = T,reps = 100, n.snp = 500)
## IA.mostclone
ia.most <- samp.ia(most_clone, quiet = T, reps = 100, n.snp = 500)

ia.lin2 <- samp.ia(t1,  reps = 100, quiet = T, n.snp = 500)
#ia.lin2 <- samp.ia(t2,  reps = 100, quiet = T, n.snp = 300)


# Summarizing data frames
d1 <- data.frame(ia.lin2, rep("Lineage2", length(ia.lin2)))
d3 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d4 <- data.frame(ia.clone, rep("clonal", length(ia.clone)))
d5 <- data.frame(ia.semi, rep("semi-clonal", length(ia.semi)))
d6 <- data.frame(ia.most, rep("most-clonal", length(ia.semi)))
colnames(d1) <- c("ia","dset")
#colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")
colnames(d6) <- c("ia","dset")
ia.total <- rbind(d4, d6, d5, d3, d1)
#ia.total <- rbind(d1, d2, d3, d4, d5)

# Normality tests
frames <- list(as.data.frame(d1), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5), as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

# Analysis of variance
anova.ia <- aov(lm(ia ~ dset, ia.total))
library(agricolae)
tukey <- HSD.test(anova.ia, "dset", alpha = 0.001)
tukey
# Kluskal wallis test
#kruskal.test(ia ~ dset, ia.total), trt="dset")
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))

# Plot
gg4 <- ggplot(ia.total,aes(x=reorder(dset,-ia), y=ia,fill=dset)) + 
  geom_boxplot() + xlab("Dataset") + ylab("Index of association") + 
  scale_fill_manual(values=c("#FFFFFF", "#D95F02", "#C1C1C1", "#6D6D6D", "#1E1E1E")) + 
  ylim(NA, 0.15)
gg4


#Lineage 3

#setPop(gl.necrophora.subset.2.filtered) <-  ~Hyp_Final
t <- seppop(gl.necrophora.subset.2.filtered)
t1 <- clonecorrect(as.snpclone(t$Lin3))
#t2 <- clonecorrect(as.snpclone(t$Lin2))

set.seed(100)
sex <- glSim(n.ind = nInd(t1), n.snp.nonstruc = ceiling(0.9*nLoc(t1)), n.snp.struc = floor(0.1*nLoc(t1)), ploidy=1, LD=TRUE)
### Structure (clonal pops)
clone <- glSim(nInd(t1), n.snp.nonstruc = floor(0.1*nLoc(t1)), n.snp.struc=ceiling(0.9*nLoc(t1)), ploidy=1, LD = T)
### Semi-clonal 
semi_clone <- glSim(nInd(t1),n.snp.nonstruc = 0.5*nLoc(t1), n.snp.struc= 0.5*nLoc(t1), ploidy=1, LD=T)
### Most-clonal 
most_clone <- glSim(nInd(t1), n.snp.nonstruc = ceiling(nLoc(t1)/3), n.snp.struc=2*nLoc(t1)/3, ploidy=1, LD=T)

## IA sex
ia.sex <- samp.ia(sex,quiet = T, reps = 100, n.snp = 500)
## IA clone
ia.clone <- samp.ia(clone, quiet = T, reps = 100, n.snp = 500)
## IA.semiclone
ia.semi <- samp.ia(semi_clone, quiet = T,reps = 100, n.snp = 500)
## IA.mostclone
ia.most <- samp.ia(most_clone, quiet = T, reps = 100, n.snp = 500)

ia.lin3 <- samp.ia(t1,  reps = 100, quiet = T, n.snp = 500)
#ia.lin2 <- samp.ia(t2,  reps = 100, quiet = T, n.snp = 300)


# Summarizing data frames
d1 <- data.frame(ia.lin3, rep("Lineage3", length(ia.lin3)))
d3 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d4 <- data.frame(ia.clone, rep("clonal", length(ia.clone)))
d5 <- data.frame(ia.semi, rep("semi-clonal", length(ia.semi)))
d6 <- data.frame(ia.most, rep("most-clonal", length(ia.semi)))
colnames(d1) <- c("ia","dset")
#colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")
colnames(d6) <- c("ia","dset")
ia.total <- rbind(d4, d6, d5, d3, d1)
#ia.total <- rbind(d1, d2, d3, d4, d5)

# Normality tests
frames <- list(as.data.frame(d1), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5), as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

# Analysis of variance
anova.ia <- aov(lm(ia ~ dset, ia.total))
library(agricolae)
tukey <- HSD.test(anova.ia, "dset", alpha = 0.001)
tukey
# Kluskal wallis test
#kruskal.test(ia ~ dset, ia.total), trt="dset")
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))

# Plot
gg6 <- ggplot(ia.total,aes(x=reorder(dset,-ia), y=ia,fill=dset)) + geom_boxplot() + 
  xlab("Dataset") + ylab("Index of association") + 
  scale_fill_manual(values=c("#FFFFFF", "#7570B3", "#C1C1C1", "#6D6D6D", "#1E1E1E"))
gg6



#gridExtra::grid.arrange(gg1, gg3, gg5,
                   # gg2,
                   # gg4,
                   # gg6,
                   # ncol = 2 
                    #labels = c("A", "B", "C", "D", "E", "F")
                   # )
#dev.off()



