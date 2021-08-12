# ancestral state for VCF

Bedtools to get variant at all positions in VCF
````bash
module load BEDTools/2.30.0-GCC-8.3.0

vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
ancestral=/scratch/eld72413/SAM_seq/ANGSD/Ancestral/SRS2413741_0.03_realigned.fa

#bedtools getfasta -bedOut -name -fi ${ancestral} -bed ${vcf} | head # why did this return a vcf?
# using "-tab" instead of "-bedOut" works

# the '-name' flag reports the reference/alternate allele
# 	- I will keep this for now even though the input for Peter's script doesn't use it
bedtools getfasta -tab -name -fi ${ancestral} -bed ${vcf} > /scratch/eld72413/SAM_seq/Polarized/Ancestral_StatesGenotypes.bed

# Use R to parse
cd /scratch/eld72413/SAM_seq/Polarized/
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load R/4.0.0-foss-2019b
R
````

```R
Ancestral <- read.table("Ancestral_StatesGenotypes.bed", sep = c("\t", "::", " "), header=FALSE)

colnames(Ancestral) <- c("Locus", "Ancestral_Allele")
length(Ancestral$Locus) # 37,129,758
levels(as.factor(Ancestral$Ancestral_Allele)) # upper and lowercase N

aggregate(Ancestral$Locus, by=list(Ancestral$Ancestral_Allele), length)
# upper and lower case letters look similar except for N's (much fewer lower case n's)

# convert to upper-case:
Ancestral$Ancestral_Allele2 <- as.factor(toupper(Ancestral$Ancestral_Allele))

# check
levels(Ancestral$Ancestral_Allele2)
````

checking numbers
```R
#Ancestral_noN <- subset(Ancestral, Ancestral_Allele!="N" & Ancestral_Allele!="n")
#length(Ancestral_noN$Locus) # 13,524,630
length(Ancestral$Locus) # 37129758
length(Ancestral[which(Ancestral$Ancestral_Allele2!="N"),"Locus"]) # 13,524,630

#write.table(Ancestral_noN, "Ancestral_StatesGenotypes_noNs.bed",
#	sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#save test
#head -1000 Ancestral_StatesGenotypes_noNs.bed >Ancestral_StatesGenotypes_noNsTEST.bed

AncestralFileParse <- function (dataframe, LocusCol, AncestralAlleleColName) {
  dataframe$List <- strsplit(as.character(LocusCol), "::")
  dataframe$Alleles <- as.factor(sapply(dataframe$List, "[", 1))
  dataframe$Region <- as.factor(sapply(dataframe$List, "[", 2))
  dataframe$List2 <- strsplit(as.character(dataframe$Alleles), "/")
  dataframe$Reference <- as.factor(sapply(dataframe$List2, "[", 1))
  dataframe$Alternate <- as.factor(sapply(dataframe$List2, "[", 2))
  return(dataframe[,c("Region", "Alleles", "Reference", "Alternate", AncestralAlleleColName)])
}

ancestralDF <- AncestralFileParse(Ancestral, Ancestral$Locus, "Ancestral_Allele2")

str(ancestralDF)
#ancestralDF_noN$Ancestral_Allele2 <- droplevels(ancestralDF_noN$Ancestral_Allele2)

aggregate(ancestralDF$Region, by=list(ancestralDF$Reference,
                                      ancestralDF$Alternate,
                                      ancestralDF$Ancestral_Allele2), length)

ancestralDF$Category <- ifelse(ancestralDF$Ancestral_Allele2=="N",
                               "Ancestral_N",
                               ifelse(as.character(ancestralDF$Ancestral_Allele2)==as.character(ancestralDF$Reference),
                                   "Alt_derived", 
                                   ifelse(as.character(ancestralDF$Ancestral_Allele2)==as.character(ancestralDF$Alternate),
                                                         "Ref_derived",
                                                         "Both_derived")))

aggregate(ancestralDF$Region, by=list(ancestralDF$Category), length)
```
Result:
       Group.1        x
1  Alt_derived 10587974
2  Ancestral_N 23605128
3 Both_derived   269140
4  Ref_derived  2667516

Format to use for Peter's script
```R
AncestralFileFormat <- function (dataframe, RegionCol) {
  dataframe$List <- strsplit(as.character(RegionCol), ":")
  dataframe$Chromosome <- as.factor(sapply(dataframe$List, "[", 1))
  dataframe$BED_pos <- as.factor(sapply(dataframe$List, "[", 2))
  dataframe$List2 <- strsplit(as.character(dataframe$BED_pos), "-")
  dataframe$Position <- as.factor(sapply(dataframe$List2, "[", 2))
  return(
    subset(dataframe, select=-c(List, List2))
          )
}
  #return(dataframe[,c("Chromosome", "Position", AncestralAlleleColName)])

#Ancestral$Locus_list <- strsplit(as.character(Ancestral$Locus), "::")
#Ancestral$Region <- sapply(Ancestral$Locus_list, "[", 2)

AncestralDF <- AncestralFileFormat(ancestralDF, ancestralDF$Region)

# to use in Peter's script
write.table(AncestralDF[,c("Chromosome", "Position", "Ancestral_Allele2")], "AncestralStateTable",
	sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
```

Use Peter Morrell's ancestral_state.py code to pull out ancestral state info for all genotypes at all positions
```bash
# gzip ancestral state list
gzip -c AncestralStateTable > AncestralStateTable.gz

vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
AncTable=/scratch/eld72413/SAM_seq/Polarized/AncestralStateTable.gz

module load Python/3.8.6-GCCcore-10.2.0
python /home/eld72413/Utilities/PeterMorrell/Utilities/ancestral_state.py \
${AncTable} ${vcf} > /scratch/eld72413/SAM_seq/Polarized/GenotypeAncestralStates

# gzip table
gzip -c /scratch/eld72413/SAM_seq/Polarized/GenotypeAncestralStates > /scratch/eld72413/SAM_seq/Polarized/GenotypeAncestralStates.gz

# script to pull out derived allele frequency
python /home/eld72413/Utilities/PeterMorrell/Utilities/count_derived.py \
```

Merge with compiled predictions
```R
# make column for VariantID to merge with compiled predictions

AncestralDF$VariantID <- paste0(AncestralDF$Chromosome, "_", AncestralDF$Position,
                                 "_", AncestralDF$Alleles)
library(pryr)
mem_used() # 22.3 GB
save(AncestralDF, file = "AncestralDF.RData")

dsnp_table <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table",
  sep = "\t", stringsAsFactors = FALSE)
mem_used() # 23.2 GB
colnames(dsnp_table) <- dsnp_table[1,]
dsnp_table <- dsnp_table[-1,]
length(dsnp_table$VariantID) # 645213
length(dsnp_table[which(dsnp_table$ReferenceAA == dsnp_table$RefAA), "VariantID"]) # 645213

dsnp_data <- merge(AncestralDF[,c(1:7,9:10)], dsnp_table, by="VariantID")
length(dsnp_data$VariantID) # 645,213

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Category, dsnp_data$Result), length)

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Category, dsnp_data$Result,
  dsnp_data$Refderived, dsnp_data$Altderived), length)

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Category, dsnp_data$Refderived), length)
```
       Group.1     Group.2      x
1  Alt_derived Deleterious  36413
2  Ancestral_N Deleterious  15115
3 Both_derived Deleterious    293
4  Ref_derived Deleterious   2624

5  Alt_derived   Tolerated 311060
6  Ancestral_N   Tolerated 185013
7 Both_derived   Tolerated   4966
8  Ref_derived   Tolerated  89729

        Group.1     Group.2       Group.3       Group.4      x
1   Ancestral_N   Tolerated derived_state derived_state      4
2   Ref_derived   Tolerated derived_state derived_state      1
3   Alt_derived Deleterious   not_derived derived_state  36413
4   Ancestral_N Deleterious   not_derived derived_state  15115
5  Both_derived Deleterious   not_derived derived_state    293
6   Ref_derived Deleterious   not_derived derived_state   2624
7   Alt_derived   Tolerated   not_derived derived_state 125017
8   Ancestral_N   Tolerated   not_derived derived_state  78213
9  Both_derived   Tolerated   not_derived derived_state   2190
10  Ref_derived   Tolerated   not_derived derived_state  25069
11  Alt_derived   Tolerated   not_derived   not_derived 186043
12  Ancestral_N   Tolerated   not_derived   not_derived 106796
13 Both_derived   Tolerated   not_derived   not_derived   2776
14  Ref_derived   Tolerated   not_derived   not_derived  64659

1  Ancestral_N derived_state      4
2  Ref_derived derived_state      1
3  Alt_derived   not_derived 347473
4  Ancestral_N   not_derived 200124
5 Both_derived   not_derived   5259
6  Ref_derived   not_derived  92352

alternative order of operations
```R
dsnp_data$Result2 <- ifelse(dsnp_data$pAdjusted < 0.05 & 
                               dsnp_data$SeqCount >= 10 & 
                               dsnp_data$MaskedConstraint < 1,
                             "Deleterious", "Tolerated")

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Category, dsnp_data$Result2,
  dsnp_data$Refderived, dsnp_data$Altderived), length)

dsnp_data[which(dsnp_data$Refderived=="derived_state"),]
```
        Group.1     Group.2       Group.3       Group.4      x
1   Ancestral_N   Tolerated derived_state derived_state      4
2   Ref_derived   Tolerated derived_state derived_state      1
3   Alt_derived Deleterious   not_derived derived_state  36413
4   Ancestral_N Deleterious   not_derived derived_state  15115
5  Both_derived Deleterious   not_derived derived_state    293
6   Ref_derived Deleterious   not_derived derived_state   2624
7   Alt_derived   Tolerated   not_derived derived_state 125017
8   Ancestral_N   Tolerated   not_derived derived_state  78213
9  Both_derived   Tolerated   not_derived derived_state   2190
10  Ref_derived   Tolerated   not_derived derived_state  25069
11  Alt_derived Deleterious   not_derived   not_derived  22302
12  Ancestral_N Deleterious   not_derived   not_derived  12247
13 Both_derived Deleterious   not_derived   not_derived    218
14  Ref_derived Deleterious   not_derived   not_derived   8386
15  Alt_derived   Tolerated   not_derived   not_derived 163741
16  Ancestral_N   Tolerated   not_derived   not_derived  94549
17 Both_derived   Tolerated   not_derived   not_derived   2558
18  Ref_derived   Tolerated   not_derived   not_derived  56273

the four for which the reference allele is in "derived state" have NA's for alignment
head(dsnp_data[,c(1,26)])
