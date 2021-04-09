#!/usr/bin/env Rscript

### A script that saves a table with deleterious vs. tolerated for each SNP, including info merged from the Predict output + VeP summary

options(warn=1)


#########################
####### FUNCTION ########
#########################

TolvDel_sites <- function (Predict_file, VeP_file, P_cutoff, minseq, max_constraint, Name_Pvalues_column) {
  dsnp <- read.table(Predict_file, sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
  vep <- read.table(VeP_file, sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", "-"))
  colnames(vep) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
  dsnp_data <- merge(dsnp, vep, by="VariantID")
  dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/") # split amino acids column
  dsnp_data$RefAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 1))
  dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))
  dsnp_data$Alignment_list <- strsplit(dsnp_data$Alignment, "") # make alignment characters into a list
  # determine whether reference or alternate AA's are private:
  dsnp_data$Refderived <- ifelse(apply(dsnp_data, 1, function(row) {
    row["RefAA"] %in% unlist(row["Alignment_list"])
  }) == FALSE,
  "derived_state", "not_derived")
  dsnp_data$Altderived <- ifelse(apply(dsnp_data, 1, function(row) {
    row["AltAA"] %in% unlist(row["Alignment_list"])
  }) == FALSE,
  "derived_state", "not_derived")
  # p-value using BH/FDR correction:
  dsnp_data$pAdjusted <- p.adjust(dsnp_data[,Name_Pvalues_column], method = "BH", n = length(dsnp_data[,Name_Pvalues_column]))
  dsnp_data$Result <- ifelse(dsnp_data$pAdjusted < P_cutoff & 
                               dsnp_data$SeqCount >= minseq & 
                               dsnp_data$MaskedConstraint < max_constraint & 
                               (dsnp_data$Refderived == "derived_state" | 
                                  dsnp_data$Altderived == "derived_state"),
                             "Deleterious", "Tolerated")
  Resultdf <- subset(dsnp_data, select = -c(Alignment_list, Amino_acidsSPLIT))
  return(Resultdf)
}


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Predict <- args[1]
VEP <- args[2]
Pcutoff <- args[3]
MinSeq <- args[4]
MaxConstraint <- args[5]
MaskedvUnMasked <- args[6]
OutFile <- args[7]

PvalsColumn <- ifelse(MaskedvUnMasked=="Unmasked",
				"LogisticP_Unmasked",
                 "LogisticP_Masked")

Mydf <- TolvDel_sites(Predict, VEP, Pcutoff, MinSeq, MaxConstraint, PvalsColumn)

write.table(Mydf, OutFile, sep = "\t", quote=FALSE, row.names=FALSE)
