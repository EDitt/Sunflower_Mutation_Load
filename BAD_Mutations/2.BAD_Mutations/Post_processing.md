
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

Cloned Deleterious_Mutations repository: https://github.com/MorrellLAB/Deleterious_Mutations
```bash
module load Biopython/1.75-foss-2019b-Python-2.7.16
Del_Mut_Path="/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts"
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
vep=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
python ${Del_Mut_Path}/VEP_To_Table.py ${vcf} ${vep} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/VEP_Table

```
^ this did not work, possibly only meant for soybean samples?
Traceback (most recent call last):
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 115, in <module>
    main()
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 90, in main
    vdata = parse_vcf(sys.argv[1])
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 19, in parse_vcf
    snpid = tmp[2]
IndexError: list index out of range

# Significant Variants from Compiled Predict report

I tested 50,838 codons for Sunflower, so significance threshold= 0.05/50838 <- way too stringent for this dataset

Will also filter out alignments with fewer than 10 species

```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
module load R/4.0.0-foss-2019b
R

```

```R
dsnp <- read.table("Sunflower_SAM_Combined_Report_preliminary2.txt", sep = "\t", header=TRUE)

missense <- read.table("fullsam_missense_noHEADER.txt", sep = "\t", header=FALSE)
colnames(missense) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
# columns I need:
# 1) variant ID that matches with the BAD_Mutations report
# 2) snp position
# 3) variant allele
# 4) gene
# 11) amino acids

# merge based on position
dsnp_data <- merge(dsnp, missense, by="VariantID")

# split the Amino_acids column:
dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/")
dsnp_data$RefAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 1))
dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))

# I think(?) the "ReferenceAA" in the BAD_Mutations output is the reference of unrelated Angiosperm genomes
length(which(dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA)) # 425,939
length(which(dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA &
	dsnp_data$MaskedConstraint < 1 )) # 350,284

# Tom's criteria
lrt_sig <- 0.05/50838
length(which(dsnp_data$LogisticP_Masked < lrt_sig)) # 0
min(na.omit(dsnp_data$LogisticP_Masked)) # 3.345465e-05
min(na.omit(dsnp_data$LogisticP_Unmasked)) # 3.310914e-07

# too stringent, try FDR
dsnp_data$pAdjusted <- p.adjust(dsnp_data$LogisticP_Masked, method = "BH", n = length(dsnp_data$LogisticP_Masked))
length(which(dsnp_data$pAdjusted < 0.05)) # 83642

minseq <- 10
length(which(dsnp_data$SeqCount >= minseq)) # 385827 (out of 420,768)

length(which(dsnp_data$SeqCount >= minseq & dsnp_data$pAdjusted < 0.05)) # 83642

max_constraint <- 1

lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq & 
	dsnp_data$MaskedConstraint < max_constraint & (dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA)), ] 
length(lrt$VariantID) # 81565

#lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq & 
#	dsnp_data$MaskedConstraint < max_constraint), ] # 81565

#lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq), ] # 83642

write.table(lrt[,-30], "dsnp_data_PRELIM.table", sep = "\t", quote=FALSE, row.names=FALSE)
```

# Subset VCF file

First subset vcf to find *all* missense positions
```bash

awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | awk '{$1=$1}1' FS=':' OFS='\t' > Missense_positions.txt


sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Missense_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_missense' Subset_vcf.sh # 2164550
```

Subset this vcf file to make separate vcfs of both deleterious and tolerated



```bash
####### scratch
awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt > Missense_positions.txt

REGIONS_ARR=($(cat Missense_positions.txt))

REGIONS_ARR=($(awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | cat))

echo ${#REGIONS_ARR[@]}

for pos in `cat Missense_positions.txt`; do
	REGIONS_ARR=("${REGIONS_ARR[@]}", "$pos")
done

REGIONS_string=$(paste -sd Missense_positions.txt)


# or change delimiter
awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | head
awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | awk '{$1=$1}1' FS=':' OFS='\t' | head
```




---

I want to subset the Vep file by the ones that I've  

```bash
grep -v "#" ${vep} | head # 1st column is variant ID, 2nd column is location in chromosomal coordinates

module load BCFtools/1.10.2-GCC-8.3.0
bcftools view -H $vcf | head # 1st column is chromosome, 2nd column is position
```
