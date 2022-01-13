#!/usr/bin/env Rscript

### A script that outputs a table of frequencies for two groups

options(warn=1)

### Command line arguments:
# 1.) bcftools output for group 1
# 2.) name of group 1
# 3.) bcftools output for group 2
# 4.) name of group2
# 5.) SNP info table
# 6.) Filepath and output filename for combined table


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

bcftools1 <- args[1]
group1 <- args[2]
bcftools2 <- args[3]
group2 <- args[4]
SNP_table <- args[5]
Output_File <- args[6]

#########################
####### FUNCTIONS #######
#########################


SNP_freq <- function (bcftools_file, SNP_info, group_name) {
	data <- read.table(bcftools_file, sep = "\t", header=FALSE,
		col.names=c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles", "Alt_Freq"))
	data$Altfreq <- data$Num_Alt_alleles / data$Num_alleles
	data$Reffreq <- (data$Num_alleles - data$Num_Alt_alleles) / data$Num_alleles
	data$MAF <- ifelse(data$Reffreq < data$Altfreq, data$Reffreq, data$Altfreq)
	SNP_table <- read.table(SNP_info, sep="\t", header=TRUE)
	data_Anc <- merge(data, SNP_table[,c("Chromosome", "Position", "MAF", "Ancestral_Allele", "Derived_Freq", "Variant_type")], 
		by=c("Chromosome", "Position"))
	Colname1 <- paste0(group_name, "_Derived_Freq")
	data_Anc[,Colname1] <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
	data_Anc$Altfreq, ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
		data_Anc$Reffreq, NA))
	Colname2 <- paste0(group_name, "_Derived_Num")
	data_Anc[,Colname2] <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
	data_Anc$Num_Alt_alleles, ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
		data_Anc$Num_alleles - data_Anc$Num_Alt_alleles, NA))
	data_Anc <- with(data_Anc, data_Anc[order(Chromosome, Position),])
  	return(data_Anc[which(!is.na(data_Anc$Derived_Freq)),c("Chromosome", "Position", "Derived_Freq", "Variant_type", Colname1, Colname2)])
}

#########################
# APPLY FUNCTIONS TO DATA #
#########################

Group1_table <- SNP_freq(bcftools1, SNP_table, group1)
Group2_table <- SNP_freq(bcftools2, SNP_table, group2)

Combined_table <- merge(Group1_table, Group2_table, by=c("Chromosome", "Position", "Derived_Freq", "Variant_type"))

write.table(Combined_table, Output_File, sep = "\t", quote=FALSE, row.names=FALSE)

