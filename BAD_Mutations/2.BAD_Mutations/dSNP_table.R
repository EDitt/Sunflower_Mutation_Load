#!/usr/bin/env Rscript

### A script that saves a table with deleterious vs. tolerated for each SNP, including info merged from the Predict output + VeP summary

options(warn=1)


### Command line arguments:
# 1.) Predict_file - this is the compiled report output from the BAD_Mutations compile subcommand
# 2.) VeP file - the text file report output from running VeP
# 3.) The p-value cutoff you wish to use after correction (e.g. 0.05)
# 4.) The minimum number of sequences that must be represented in the alignment
# 5.) The maximum constraint value
# 6.) "Unmasked" if you want to use the un-masked p-values from the predictions. This defaults to masked for any other string

### Example command line:
#     Rscript ${Sunflower_Mutation_Load}/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
#     ${BAD_Mut_Files_Results_DIR}/Compiled_report.txt \
#     ${VEP_OUTPUTDIR}/SAM_SNP_Final_BiallelicNorm \
#     0.05 \
#     10 \
#     1 \
#     Masked \
#     ${OUTPUTDIR}/dsnp_data.table

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
