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
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
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
Ancestral_noN <- subset(Ancestral, Ancestral_Allele!="N" & Ancestral_Allele!="n")
length(Ancestral_noN$Locus) # 13,524,630
write.table(Ancestral_noN, "Ancestral_StatesGenotypes_noNs.bed",
	sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#save test
#head -1000 Ancestral_StatesGenotypes_noNs.bed >Ancestral_StatesGenotypes_noNsTEST.bed

AncestralFileParse <- function (dataframe, LocusCol, AncestralAlleleColName) {
  dataframe$List <- strsplit(as.character(LocusCol), "::")
  dataframe$Alleles <- as.factor(sapply(dataframe$List, "[", 1))
  dataframe$Region <- as.factor(sapply(dataframe$List, "[", 2))
  dataframe$List2 <- strsplit(as.character(dataframe$Alleles), "/")
  dataframe$Reference <- as.factor(sapply(dataframe$List2, "[", 1))
  dataframe$Alternate <- as.factor(sapply(dataframe$List2, "[", 2))
  return(dataframe[,c("Region", "Reference", "Alternate", AncestralAlleleColName)])
}

ancestralDF_noN <- AncestralFileParse(Ancestral_noN, Ancestral_noN$Locus, "Ancestral_Allele2")

str(ancestralDF_noN)
ancestralDF_noN$Ancestral_Allele2 <- droplevels(ancestralDF_noN$Ancestral_Allele2)

aggregate(ancestralDF_noN$Region, by=list(ancestralDF_noN$Reference,
                                                 ancestralDF_noN$Alternate,
                                                 ancestralDF_noN$Ancestral_Allele2), length)

ancestralDF_noN$Category <- ifelse(ancestralDF_noN$Ancestral_Allele2==ancestralDF_noN$Reference,
                               "Alt_derived", ifelse(ancestralDF_noN$Ancestral_Allele2==ancestralDF_noN$Alternate,
                                                     "Ref_derived",
                                                     "Both_derived"))

aggregate(ancestralDF_noN$Region, by=list(ancestralDF_noN$Category), length)
```
Result:
       Group.1        x
1  Alt_derived 10587974
2 Both_derived   269140
3  Ref_derived  2667516

Format to use for Peter's script
```R
AncestralFileFormat <- function (dataframe, RegionCol, AncestralAlleleColName) {
  dataframe$List <- strsplit(as.character(RegionCol), ":")
  dataframe$Chromosome <- as.factor(sapply(dataframe$List, "[", 1))
  dataframe$BED_pos <- as.factor(sapply(dataframe$List, "[", 2))
  dataframe$List2 <- strsplit(as.character(dataframe$BED_pos), "-")
  dataframe$Position <- as.factor(sapply(dataframe$List2, "[", 2))
  return(dataframe[,c("Chromosome", "Position", AncestralAlleleColName)])
}

Ancestral$Locus_list <- strsplit(as.character(Ancestral$Locus), "::")
Ancestral$Region <- sapply(Ancestral$Locus_list, "[", 2)

AncestralDF <- AncestralFileFormat(Ancestral, Ancestral$Region, "Ancestral_Allele2")

write.table(AncestralDF, "AncestralStateTable",
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
