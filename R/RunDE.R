################################################################################
### R script based on SARTools by Hugo Varet
### Heath O'Brien
### 12 April, 2017
### designed to be executed with a modified version of SARTools 1.3.0 
### (https://github.com/hobrien/SARTools)
################################################################################
#setwd("../")

rm(list=ls())                                        # remove all the objects from the R session
library("optparse")

option_list <- list(
  make_option(c("-v", "--varInt"), type="character", default="condition", 
              help="Variable of interest (either Sex or PCW)"),
  make_option(c("-r", "--ref"), type="character", default="Uninfected", 
              help="Reference level of varInt"),
  make_option(c("-b", "--batch"), type="character", default=NULL, 
              help="batch correction"),
  make_option(c("-i", "--interaction"), type="character", default=NULL, 
              help="Cofactors to interact with varInt", metavar = 'interact'),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.1, 
              help="corrected pvalue cutoff", metavar="pvalue"),
  make_option(c("-c", "--conditions"), type="character", default="conditions.txt", 
              help="Cofactors (either Sex or PCW)"),
  make_option(c("-e", "--exclude"), type="character", default="Reference/Candida_genes.txt", 
              help="File with list of genes to exclude"),
  make_option(c("-f", "--feature"), type="character", default="genes", 
              help="Type of feature to Analyse (genes, junctions, transcripts)"),
  make_option(c("-k", "--kallisto"), action='store_true', type="logical", default=FALSE, 
              help="Use counts derived from Kallisto"),
  make_option(c("-s", "--sva"), type="integer", default=0, 
              help="Number of Surrogate Variables to estimate"),
  make_option(c("-n", "--name"), type='character', default="ManducaInfection", 
              help="Project name"),
  make_option(c("-o", "--output"), type='character', default="Results", 
              help="Folder for output files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=TRUE)
opt$args<-paste0("Mappings/", c("A4", "A5", "A6", "B5", "B6", "C5", "C6", "D4", "E4"), ".counts.txt")
if (! opt$options$feature %in% c('genes', 'junctions', 'transcripts')) {
  stop("feature type not recognised. Must be one of (genes, junctions, transcripts)")
}

projectName <- opt$options$name                         # name of the project
print(paste("Saving output to", projectName))

author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- opt$options$output     # working directory for the R session                                     # path to the directory containing raw counts files

targetFile <- opt$options$conditions


condRef <- opt$options$ref                                     # reference biological condition

pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

# DESeq parameters
fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
#if numeric, features with maxCooks values above this number are removed 
cooksCutoff <-  TRUE #0.75                          # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
# p-value adjustment method: "BH" (default) or "BY"
testMethod <- 'Wald'

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

################################################################################
###                             running script                               ###
################################################################################

library(devtools)
load_all(pkg = "R/SARTools")
library(tidyverse)
library(stringr)
library(RColorBrewer)

featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove
if (! is.null(opt$options$exclude)) {
  featuresToRemove <- c(featuresToRemove, read_tsv(opt$options$exclude, col_names = FALSE)$X1)
}


# loading target file
LibraryInfo <- read_tsv(opt$options$conditions, 
                        col_types = cols(.default = col_character())
)
LibraryInfo <- tibble(files = opt$args) %>% 
  mutate(SampleID=str_extract(basename(files), '^[^.]+')) %>% full_join(LibraryInfo)

LibraryInfo <- as.data.frame(LibraryInfo)

# loading counts
if (opt$options$kallisto) {
    library(tximport)
    tx2gene <- read_tsv("Data/tx2gene.txt")
    files <- file.path("Kallisto", LibraryInfo$Sample, "abundance.tsv")
    names(files) <- LibraryInfo$Sample
    if ( opt$options$feature == 'genes' ) {
      counts <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader=read_tsv)
    } else if ( opt$options$feature == 'transcripts' ) {
      counts <- tximport(files, type = "kallisto", txOut = TRUE, reader=read_tsv)
    } else {
      stop("Kallisto output cannot be used to analyse junctions")
    }
} else {
    counts <- loadCountData(labels=LibraryInfo$SampleID, files=LibraryInfo$files, featuresToRemove=featuresToRemove)
}


dir.create(workDir)
setwd(workDir)

# description plots
if (opt$options$kallisto) {
  majSequences <- descriptionPlots(counts=counts$counts, group=LibraryInfo[,opt$options$varInt], col=colors)
} else {
  majSequences <- descriptionPlots(counts=counts, group=LibraryInfo[,opt$options$varInt], col=colors)
}

# DEseq analysis
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                                              featuresToRemove=featuresToRemove,varInt=opt$options$varInt,
                                              condRef=opt$options$ref,batch=opt$options$batch,fitType=fitType,cooksCutoff=cooksCutoff,
                                              independentFiltering=independentFiltering,alpha=opt$options$pvalue,pAdjustMethod=pAdjustMethod,
                                              typeTrans=typeTrans,locfunc=locfunc,colors=colors)


out.DESeq2 <- run.DESeq2(counts=counts, target=LibraryInfo, varInt=opt$options$varInt, batch=opt$options$batch, 
                           interact=opt$options$interact, num_sva=opt$options$sva, kallisto=opt$options$kallisto,
                           locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                           cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=opt$options$pvalue)
# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=LibraryInfo[,opt$options$varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=LibraryInfo[,opt$options$varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, kallisto=opt$options$kallisto, alpha=opt$options$pvalue)

################################################################################

write_tsv(as.data.frame(colData(out.DESeq2$dds)), "tables/col_data.txt")
# save image of the R session
save.image(file=paste0(projectName, ".RData"))


# generating HTML report
writeReport.DESeq2(target=LibraryInfo, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, featuresToRemove=featuresToRemove, varInt=opt$options$varInt,
                   condRef=opt$options$ref, batch=opt$options$batch, interact=opt$options$interact, 
                   fitType=fitType, cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, 
                   alpha=opt$options$pvalue, pAdjustMethod=pAdjustMethod, typeTrans=typeTrans, 
                   locfunc=locfunc, colors=colors)
