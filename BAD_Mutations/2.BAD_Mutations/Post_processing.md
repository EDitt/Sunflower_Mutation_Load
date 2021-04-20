
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


Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report_preliminary2.txt" \
"/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm" \
"0.05" \
"10" \
"1" \
"Masked" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_PRELIM_newTEST2.table"



```

```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R")

Mydf <- TolvDel_sites("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report_preliminary2.txt",
"/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm",
0.05,
10,
1,
"LogisticP_Masked")

write.table(Mydf[-30], "dsnp_data_PRELIM_NEW.table", sep = "\t", quote=FALSE, row.names=FALSE)
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
