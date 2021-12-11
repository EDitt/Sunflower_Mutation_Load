#!/usr/bin/env Rscript

### A script that outputs a table with the number of derived variants for 3 frequency classes (Deleterious, Tolerated, Synonymous) for all genotypes

options(warn=1)

### Command line arguments:
# 1.) Directory output from "grep PSC" bcftools stats command (tab delimited table with numbers of variants in each class)
# 2.) File + directory info for "grep PSC" from bcftools stats on the whole VCF file (to get total number of genotype calls for each sample)

#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Directory <- args[1]
AllStats <- args[2]
OutputFile <- args[3]

#########################
####### FUNCTIONS #######
#########################

# import tables with positions for various variant classes and edit variant class names for clarity
ImportTxts <- function (DirPath) {
  my_files <- list.files(path = DirPath, pattern = "*.txt", full.names = TRUE)
  my_data <- lapply(my_files, read.table)
  names(my_data) <- gsub("\\.txt$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  names(my_data) <- gsub("_perSampleCounts", "", names(my_data))
  my_data_RefDerived <- my_data[grep("RefDerived", names(my_data))]
  for (i in seq_along(my_data_RefDerived)) {
    name <- names(my_data_RefDerived[i])
    colnames(my_data_RefDerived[[i]]) <- c("PSC", "id", "sample", "nDerivedHom", "nAncestralHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")
  }
  my_data_AltDerived <- my_data[grep("AltDerived", names(my_data))]
  for (i in seq_along(my_data_AltDerived)) {
    name <- names(my_data_AltDerived[i])
    colnames(my_data_AltDerived[[i]]) <- c("PSC", "id", "sample", "nAncestralHom", "nDerivedHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")
  }
  return(c(my_data_RefDerived, my_data_AltDerived))
}

CombineRefAlt <- function(DataList, Category, AllData) {
	my_data1 <- DataList[grep(Category, names(DataList))]
	Combined_data <- merge(my_data1[[1]][,c("sample", "nAncestralHom", "nDerivedHom", "nHets", "nMissing")],
	my_data1[[2]][,c("sample", "nAncestralHom", "nDerivedHom", "nHets", "nMissing")],
	by="sample", suffixes=c("_Ref", "_Alt"))
	Combined_data$NumAncestralHom <- Combined_data$nAncestralHom_Ref + Combined_data$nAncestralHom_Alt
	Combined_data$NumDerivedHom <- Combined_data$nDerivedHom_Ref + Combined_data$nDerivedHom_Alt
	Combined_data$NumHet <- Combined_data$nHets_Ref + Combined_data$nHets_Alt
	Combined_data$NumMissing <- Combined_data$nMissing_Ref + Combined_data$nMissing_Alt
	Combined_data$Consequence <- Category
	Combined_data2 <- merge(Combined_data, 
		AllData[,c("sample", "CalledGenotypes", "nMissingTotal")],
		by="sample")
	return(Combined_data2[,c("sample", "NumAncestralHom", "NumDerivedHom", "NumHet", "NumMissing", "Consequence", "CalledGenotypes", "nMissingTotal")])
}

#########################
# APPLY FUNCTIONS TO DATA #
#########################

SampleLists <- ImportTxts(Directory)

Full_stats <- read.table(AllStats)
colnames(Full_stats) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissingTotal")
Full_stats$CalledGenotypes <- Full_stats$nRefHom + Full_stats$nNonRefHom + Full_stats$nHets

Categories <- list("Deleterious", "Tolerated", "Synonymous")

AnnotationList <- lapply(Categories, function(x) {CombineRefAlt(SampleLists, x, Full_stats)})

All_Derived <- do.call("rbind", AnnotationList)

write.table(All_Derived, OutputFile, sep = "\t", quote=FALSE, row.names=FALSE)




