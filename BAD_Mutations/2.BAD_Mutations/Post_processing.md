
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

Rscript /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt \
/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm \
0.05 \
10 \
1 \
Masked \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table

```

Check numbers of deleterious vs. tolerated
```bash

wc -l SAM_SNP_BadMut_Summary # 704,075
awk '{print $4}' SAM_SNP_BadMut_Summary | sort -u | wc -l #699,805 # unique variants
awk '{print $1}' SAM_SNP_BadMut_Summary | sort -u | wc -l #50,838 # unique regions

wc -l dsnp_data.table # 1,120,517  ## why the difference?
awk 'NR > 1 {print $1}' dsnp_data.table | sort -u | wc -l # 641,519

wc -l Sunflower_SAM_Combined_Report.txt # 645,216 (including header line)

grep -v "#" SAM_SNP_Final_BiallelicNorm | wc -l # 43,019,114
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{print $1}' | sort -u | wc -l # 37,120,112

grep "missense_variant" SAM_SNP_Final_BiallelicNorm | wc -l # 712,631
grep "missense_variant" SAM_SNP_Final_BiallelicNorm | awk '{print $1}' | sort -u | wc -l # 708,311 more than expected..
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="missense_variant") {print $1}}' | wc -l # 704,075
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="missense_variant") {print $1}}' | sort -u | wc -l # 699,805

grep -v "#" fullsam_missense.txt | wc -l # 704,075
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort -u | wc -l # 699,805
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort | uniq -cd | head
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort | uniq -cd | wc -l # 4270 are duplicated

```
Summary: there should be 699,805 unique missense variant positions in the VeP output (not counting some that have multiple consequences including missense).

Need to figure out where there are so many more positions in the dsnp_data.table
```R
missense <- read.table("/scratch/eld72413/SAM_seq/VeP/fullsam_missense.txt", sep = "\t", header=FALSE)
synon <- read.table("/scratch/eld72413/SAM_seq/VeP/fullsam_synon.txt", sep = "\t", header=FALSE)
colnames(missense) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
colnames(synon) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

length(synon$VariantID) # 835,622
length(missense$VariantID) # 704,075

dsnp_table <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", header=TRUE)
length(dsnp_table$VariantID) # 1,120,516

length(which(dsnp_table$VariantID %in% missense$VariantID)) # 1,120,516
length(which(dsnp_table$VariantID %in% synon$VariantID)) # 12,962
length(which(missense$VariantID %in% synon$VariantID)) # 5346

dsnp <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt", sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
length(dsnp$VariantID) # 645,215

length(which(dsnp$VariantID %in% missense$VariantID)) # 645213
dsnp[which(!dsnp$VariantID %in% missense$VariantID),] # 2 variant ID's are "NA"
length(which(dsnp$VariantID %in% synon$VariantID))

length(which(is.na(dsnp$VariantID))) # 2
length(which(is.na(dsnp$LogisticP_Masked)))
dsnp[which(is.na(dsnp$LogisticP_Masked)),]
dsnp[which(is.na(dsnp$VariantID)),]

length(which(is.na(dsnp_table$VariantID))) #0

# full vep output
vep <- read.table("/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", "-"))
colnames(vep) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
length(vep$VariantID) # 43,019,114
length(which(is.na(vep$VariantID))) #0
length(unique(vep$VariantID)) # 37,120,112

#vep$fConsequence <- as.factor(vep$Consequence)
#levels(vep$fConsequence) # there is one category that has 2 consequences (including missense)

vep_missense <- subset(vep, Consequence=="missense_variant")
length(vep_missense$VariantID) # 704,075
length(unique(vep_missense$VariantID) # 699,805

# try merging again
dsnp_data1 <- merge(dsnp, vep, by="VariantID")
length(dsnp_data1$VariantID) # 1,120,516

#dsnp_data2 <- merge(dsnp, vep[unique(vep$VariantID),], by="VariantID")
#length(dsnp_data2$VariantID)

# merge using this code?
dsnp_data2 <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
length(dsnp_data2$VariantID)

```
In my compiled prediction report, there are 2 variant ID's that are listed as "NA".
The total number of compiled variants in this report is: 645,215 (including the 2 NA's). This is out of 704,075 missense variants in the VeP output. The only 2 that are missing are the "NA's".

For some reason, my R code gives me a table with 1,120,516 variants. These are all listed as being in the missense vep output (which only has 704,075 rows). Therefore, there must be duplicates.

To add to R script
```R
vep_missense <- subset(vep, Consequence=="missense_variant")
rm("vep")

dsnp_data <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
length(dsnp_data$VariantID) # 645,213 (this probably excludes the NA's)

#vep_missenseU <- vep_missense[unique(vep_missense$VariantID),] # this took a loooooong time. did not work (everything is NA)
#length(vep_missenseU$VariantID) # 699,805

library(pryr)
mem_used() # with vep, vep_missense, and dsnp: 11.5 GB
```

How many are deleterious vs. tolerated?
Get positions to subset VCF
```bash
grep "Tolerated" dsnp_data.table | wc -l # 993,364
```
# Subset VCF file

First subset vcf to find *all* missense positions
```bash

awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | awk '{$1=$1}1' FS=':' OFS='\t' > Missense_positions.txt


sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Missense_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_missense' Subset_vcf.sh # 2164550
```

Subset this vcf file to make separate vcfs of both deleterious and tolerated




---

I want to subset the Vep file by the ones that I've  

```bash
grep -v "#" ${vep} | head # 1st column is variant ID, 2nd column is location in chromosomal coordinates

module load BCFtools/1.10.2-GCC-8.3.0
bcftools view -H $vcf | head # 1st column is chromosome, 2nd column is position
```
