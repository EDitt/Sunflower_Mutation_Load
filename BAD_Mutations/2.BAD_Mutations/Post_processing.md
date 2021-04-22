
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

How many are deleterious vs. tolerated?
Get positions to subset VCF
```bash
grep "Tolerated" dsnp_data.table | wc -l # 590,768
grep "Deleterious" dsnp_data.table | wc -l # 54,445

## deleterious positions
grep "Deleterious" dsnp_data.table | awk '{print $17}' | awk '{$1=$1}1' FS=':' OFS='\t' > Deleterious_positions.txt
## tolerated positions
grep "Tolerated" dsnp_data.table | awk '{print $17}' | awk '{$1=$1}1' FS=':' OFS='\t' > Tolerated_positions.txt
```

# Subset VCF file

First subset vcf to find *all* missense positions
```bash

awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | awk '{$1=$1}1' FS=':' OFS='\t' > Missense_positions.txt

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Missense_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_missense' Subset_vcf.sh # 2164550


grep -v "#" SAM_missense.vcf | wc -l # 699,805
wc -l Missense_positions.txt # 704,075 # not sure why they aren't equal
```

Subset vcf file to make separate vcfs of both deleterious and tolerated
```bash
# deleterious
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Deleterious_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_deleterious' Subset_vcf.sh # Submitted batch job 2245897
# tolerated
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Tolerated_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_tolerated' Subset_vcf.sh # Submitted batch job 2245899
```

