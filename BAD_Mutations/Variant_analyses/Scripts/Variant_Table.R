#!/usr/bin/env Rscript

### A script that outputs a large table combining allele frequencies, reference and alternate allele calls, ancestral state, annotation class

options(warn=1)

### Command line arguments:
# 1.) Filepath with SNP Position and annotation
# 2.) File output from bcftools that has the number of alternate alleles and total count of alleles in called genotypes
# 3.) Text file with chromosome, position, and ancestral allele
# 4.) Output filename for Table


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Annotations <- args[1]
Ancestral_file <- args[2]
Output_Dir <- args[3]
Group_name <- args[4]

Input_file <- paste0(Output_Dir, "/intermediates/", Group_name, "_SNP_freq.txt")
Output_file <- paste0(Output_Dir, "/", Group_name, "_SNP_info.txt")

#########################
####### FUNCTIONS #######
#########################


# calculate allele frequencies from bcftools output & merge with variant type information & ancestral allele information
# outputs a list of dataframes split by variant type
SNP_Info <- function (bcftools_file, ancestral_allele_file, annotation_file) {
  annotation_table <- read.table(annotation_file, sep = "\t", header=TRUE)
  data <- read.table(bcftools_file, sep = "\t", header=FALSE)
  colnames(data) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele",
    "Num_Alt_alleles", "Num_alleles", "Alt_Freq")
  data$Altfreq <- data$Num_Alt_alleles / data$Num_alleles
  data$Reffreq <- (data$Num_alleles - data$Num_Alt_alleles) / data$Num_alleles
  data$MAF <- ifelse(data$Reffreq < data$Altfreq, data$Reffreq, data$Altfreq)
  Ancestral_table <- read.table(ancestral_allele_file, sep="\t", header=FALSE)
  colnames(Ancestral_table) <- c("Chromosome", "Position", "Ancestral_Allele")
  data_Anc <- merge(data, Ancestral_table, by=c("Chromosome", "Position"), all.x=TRUE)
  data_Anc$Derived_Freq <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
  data_Anc$Altfreq, ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
    data_Anc$Reffreq, NA))
  data_Anc <- with(data_Anc, data_Anc[order(Chromosome, Position),])
  All_data <-  merge(data_Anc, annotation_table, by=c("Chromosome", "Position"))
  return(All_data)
}

#########################
# APPLY FUNCTIONS TO DATA #
#########################


# Combine allele frequency data, ancestral state calls, and variant annotation class and output a table
AllInfo <- SNP_Info(Input_file, Ancestral_file, Annotations)

write.table(AllInfo, Output_File, sep = "\t", quote=FALSE, row.names=FALSE)
