
# Filter VeP Report

```bash
module load VEP/101.0-foss-2019b-Perl-5.30.0
INPUT=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
filter_vep -i ${INPUT} -o fullsam_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT} -o fullsam_synon.txt -filter "Consequence is synonymous_variant"

# remove header lines to use with R
#### actually R automatically ignores lines that start with "#"
awk 'NR > 29 {print}' fullsam_missense.txt > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/fullsam_missense_noHEADER.txt
awk 'NR > 29 {print}' fullsam_synon.txt > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/fullsam_synon_noHEADER.txt
```

# Significant Variants from Compiled Predict report

I tested 50,838 codons for Sunflower, so significance threshold= 0.05/50838 <- way too stringent for this dataset

Will also filter out alignments with fewer than 10 species

I ended up writing an R script to do what I was trying to do

```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results

srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# job killed at mem=22g
module load R/4.0.0-foss-2019b
module load R_ML/3.3.3 # for MSI

# update: calling the R script from the command line (as below) gives me a different number of dSNPs than running the function in R. I cannot figure out why!!! (see xArchive/Test_RscriptPredictIssues.md for troubleshooting). I will instead run the function inside R which gives the correct number.
# I also edited the R function due to realizing the alignment column *included* the reference genome (ie inmasked alignment) and therefore this had to be taken into account when 

##### previous command line code: (Put script in xArchive since there is something wrong with it)

# Rscript /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
# /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt \
# /scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm \
# 0.05 \
# 10 \
# 1 \
# Masked \
# /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table

# instead will run inside R
R
```

```R
dsnp <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt", sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
length(dsnp$VariantID) # 645,215
dsnp[which(is.na(dsnp$Alignment)),] # 5 show NA in Alignment column (SeqCount=2 for all?)
dsnp[which(is.na(dsnp$ReferenceAA)),] # 2 show NA in VariantID, ReferenceAA, P-value columns
# remove those 7 variants:
dsnp <- dsnp[-(which(is.na(dsnp$Alignment) |
	is.na(dsnp$ReferenceAA)
	)),]
length(dsnp$VariantID) # 645,208

# make alignment characters into a list
dsnp$Alignment_list <- strsplit(dsnp$Alignment, "")

# check whether reference amino acid is derived-

# first, find number of reference alleles represented in alignment across sites
dsnp$NumRefInAlignment <- apply(dsnp, 1, function(row) {
    length(which(unlist(row["Alignment_list"]) %in% row["ReferenceAA"]))
  }
)

# if not represented by any other species, reference allele will only be observed once (alignment is unmasked)
min(dsnp$NumRefInAlignment) # 1
length(dsnp[which(dsnp$NumRefInAlignment==1), "VariantID"]) # 154,759

# reference allele is derived (relative to other spp. in alignment if it's only represented once)
dsnp$Refderived <- ifelse(dsnp$NumRefInAlignment==1, "derived_state",
	ifelse(dsnp$NumRefInAlignment > 1, "not_derived",
		"NA"))
aggregate(dsnp$VariantID, by=list(dsnp$Refderived), length) # derived: 154,759; not_derived: 490,449


# merge with Vep results to get alternate AA:
vep <- read.table("/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", sep = "\t", header=FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "-"))
colnames(vep) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
library(pryr)
mem_used() # 11.9 GB
vep_missense <- subset(vep, Consequence=="missense_variant")
rm("vep")
dsnp_data <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
length(dsnp_data$VariantID) # 645,208

# split amino acids column
dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/")
# alternate AA is second
dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))
# determine whether alternate AA's are private:
dsnp_data$Altderived <- ifelse(apply(dsnp_data, 1, function(row) {
	row["AltAA"] %in% unlist(row["Alignment_list"])
	}) == FALSE,
"derived_state", "not_derived")

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Refderived, dsnp_data$Altderived), length)
#        Group.1       Group.2      x
#1 derived_state derived_state  87194
#2   not_derived derived_state 197740
#3 derived_state   not_derived  67565
#4   not_derived   not_derived 292709

# p-value using BH/FDR correction:
dsnp_data$pAdjusted <- p.adjust(dsnp_data[,"LogisticP_Masked"], method = "BH", n = length(dsnp_data[,"LogisticP_Masked"]))

dsnp_data$Result <- ifelse(dsnp_data$pAdjusted < 0.05 & 
                               dsnp_data$SeqCount >= 10 & 
                               dsnp_data$MaskedConstraint < 1 & 
                               (dsnp_data$Refderived == "derived_state" | 
                                  dsnp_data$Altderived == "derived_state"),
                             "Deleterious", "Tolerated")
aggregate(dsnp_data$VariantID, by=list(dsnp_data$Result,
	dsnp_data$Refderived, dsnp_data$Altderived), length)


Resultdf <- subset(dsnp_data, select = -c(Alignment_list, Feature_type, Consequence, Existing_variation, Extra, Amino_acidsSPLIT))

# define which allele is deleterious
Resultdf$Result_Allele <- ifelse(Resultdf$Result=="Deleterious" & Resultdf$Refderived=="derived_state",
	"Reference_deleterious", ifelse(Resultdf$Result=="Deleterious" & Resultdf$Altderived=="derived_state",
		"Alternate_deleterious", "Tolerated"))
# check
aggregate(Resultdf$VariantID, by=list(Resultdf$Result,
	Resultdf$Refderived, Resultdf$Altderived, Resultdf$Result_Allele), length)

write.table(Resultdf, "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", quote=FALSE, row.names=FALSE)

```

      Group.1       Group.2       Group.3      x
1   Tolerated derived_state derived_state  87194
2 Deleterious   not_derived derived_state  76095
3   Tolerated   not_derived derived_state 121645
4 Deleterious derived_state   not_derived  11796
5   Tolerated derived_state   not_derived  55769
6   Tolerated   not_derived   not_derived 292709

Write function + Test
```R
TolvDel_sites <- function (Predict_file, VeP_file, P_cutoff, minseq, max_constraint, Name_Pvalues_column) {
  dsnp <- read.table(Predict_file, sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
  vep <- read.table(VeP_file, sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", "-"))
  colnames(vep) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
  vep_missense <- subset(vep, Consequence=="missense_variant")
  rm("vep")
  dsnp_data <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
  # determine whether reference AA's are private:
  dsnp_data$Alignment_list <- strsplit(dsnp_data$Alignment, "") # make alignment characters into a list
  dsnp_data$NumRefInAlignment <- apply(dsnp_data, 1, function(row) {
    length(which(unlist(row["Alignment_list"]) %in% row["ReferenceAA"]))
    }
  )
  dsnp_data$Refderived <- ifelse(dsnp_data$NumRefInAlignment==1, "derived_state",
    ifelse(dsnp_data$NumRefInAlignment > 1, "not_derived",
     "NA"))
  # determine whether alternate AA's are private:
  dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/") # split amino acids column
  dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))
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
  Resultdf <- subset(dsnp_data, select = -c(Alignment_list, Feature_type, Consequence, Existing_variation, Extra, Amino_acidsSPLIT))
  # define which allele is deleterious
  Resultdf$Result_Allele <- ifelse(Resultdf$Result=="Deleterious" & Resultdf$Refderived=="derived_state",
    "Reference_deleterious", ifelse(Resultdf$Result=="Deleterious" & Resultdf$Altderived=="derived_state",
      "Alternate_deleterious", "Tolerated"))
  return(Resultdf)
}

Result_test <- TolvDel_sites("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt", "/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", 0.05, 10, 1, "LogisticP_Masked")

aggregate(Result_test$VariantID, by=list(Result_test$Result,
	Result_test$Refderived, Result_test$Altderived, Result_test$Result_Allele), length) # only difference is NAs not removed
```


