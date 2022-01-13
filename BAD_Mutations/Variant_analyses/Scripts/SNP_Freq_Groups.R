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
	Altnum <- paste0(group_name, "_Num_Alt_alleles")
	data <- read.table(bcftools_file, sep = "\t", header=FALSE,
		col.names=c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		Altnum, "Num_alleles", "Alt_Freq"))
	Altfreq <- paste0(group_name, "_AltFreq")
	data[,Altfreq] <- data[,Altnum] / data$Num_alleles
	Refnum <- paste0(group_name, "_Num_Ref_alleles")
	data[,Refnum] <- (data$Num_alleles - data[,Altnum])
	Reffreq <- paste0(group_name, "_RefFreq")
	data[,Reffreq] <- data[,Refnum] / data$Num_alleles
	SNP_table <- read.table(SNP_info, sep="\t", header=TRUE)
	data_Anc <- merge(data, SNP_table[,c("Chromosome", "Position", "Altfreq", "Reffreq", "MAF", "Ancestral_Allele", "Derived_Freq", "Variant_type")], 
		by=c("Chromosome", "Position"))
	Colname1 <- paste0(group_name, "_Derived_Freq")
	data_Anc[,Colname1] <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
	data[,Altfreq], ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
		data[,Reffreq], NA))
	Colname2 <- paste0(group_name, "_Derived_Num")
	data_Anc[,Colname2] <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
	data_Anc[,Altnum], ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
		data[,Refnum], NA))
	data_Anc <- with(data_Anc, data_Anc[order(Chromosome, Position),])
  	return(data_Anc[,c("Chromosome", "Position", "MAF", "Derived_Freq", 
  		"Variant_type", Altnum, Refnum, Altfreq, Reffreq, Colname1, Colname2)])
}

#########################
# APPLY FUNCTIONS TO DATA #
#########################

Group1_table <- SNP_freq(bcftools1, SNP_table, group1)
Group2_table <- SNP_freq(bcftools2, SNP_table, group2)

Combined_table <- merge(Group1_table, Group2_table, by=c("Chromosome", "Position", "MAF", "Derived_Freq", "Variant_type"))

write.table(Combined_table, Output_File, sep = "\t", quote=FALSE, row.names=FALSE)

