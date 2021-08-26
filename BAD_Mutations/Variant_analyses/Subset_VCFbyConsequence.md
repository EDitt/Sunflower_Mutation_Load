## Subset VCF file by category

### Setup
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

```

### Get positions to subset VCF
```bash
dsnp_data=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_Polarized.table

grep "Tolerated" ${dsnp_data} | wc -l # 557,317
grep "Deleterious" ${dsnp_data} | wc -l # 87,891

# how many duplicates?
awk '{print $19}' ${dsnp_data} | wc -l # 645,209
awk '{print $19}' ${dsnp_data} | sort -u | wc -l # 641,515

awk -F'\t' '{print NF; exit}' ${dsnp_data} # 39 columns

## deleterious positions in reference:
awk '{if ($39=="Reference_deleterious") {print $0}}' ${dsnp_data} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > Reference_DelPositions.txt
wc -l Reference_DelPositions.txt # 11,796
awk '{print $0}' Reference_DelPositions.txt | sort -u | wc -l # 11,796 no duplicates

# alternate deleterious:
awk '{if ($39=="Alternate_deleterious") {print $0}}' ${dsnp_data} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > Alternate_DelPositions.txt
wc -l Alternate_DelPositions.txt # 76,095
awk '{print $0}' Alternate_DelPositions.txt | sort -u | wc -l # 76,016
awk '{print $0}' Alternate_DelPositions.txt | sort -u  > Alternate_DelPositionsNoDups.txt
wc -l Alternate_DelPositionsNoDups.txt # 76,016

# Tolerated:
awk '{if ($39=="Tolerated") {print $0}}' ${dsnp_data} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > ToleratedPositions.txt
wc -l ToleratedPositions.txt # 556,577 (740 duplicate positions represented in deleterious set already removed - see 2.BAD_Mutations/Post_processing.md)
awk '{print $0}' ToleratedPositions.txt | sort -u | wc -l # 553,720
awk '{print $0}' ToleratedPositions.txt | sort -u > ToleratedPositionsNoDups.txt
wc -l ToleratedPositionsNoDups.txt # 553,720

```
Note: there are 18 positions that are the same in the alternate deleterious & reference deleterious sets (not removed)

Duplicates removed:
- 79 duplicates within alternate deleterious category
- 740 duplicates that were represented in deleterious & tolerated categories (removed from tolerated)
- 2,857 duplicates within tolerated category

# Subset VCF file

Subset vcf file to make separate vcfs of deleterious (for both reference and alternate alleles) and tolerated
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

# deleterious in reference
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Reference_DelPositions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Refdeleterious' Subset_vcf.sh # Submitted batch job 4232837

grep -v "#" SAM_Refdeleterious.vcf | wc -l # 11796

# deleterious in alternate
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Alternate_DelPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Altdeleterious' Subset_vcf.sh # Submitted batch job 4232853

grep -v "#" SAM_Altdeleterious.vcf | wc -l # 76016

# tolerated
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/ToleratedPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_tolerated' Subset_vcf.sh # Submitted batch job 4232872

grep -v "#" SAM_tolerated.vcf | wc -l # 553,720
```

# VCF of synonymous positions (excluding duplicates)
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

cd /scratch/eld72413/SAM_seq/VeP

grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | wc -l # 835,622
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u | wc -l # 832,090

grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u | awk '{$1=$1}1' FS=':' OFS='\t' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt

# use R to remove duplicates
module load R/4.0.0-foss-2019b
R
```
With the `sort -u` command, I removed the duplicates contained within synonymous positions (3,532). Using R I will remove the synonymous positions that are also non-synonymous positions (go by the most severe consequence for duplicate positions)

```R
dsnp <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
nonsynon_pos <- dsnp$Position

synon <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(synon) <- c("Chromosome", "Pos")

synon$Position <- paste0(synon$Chromosome,":",synon$Pos)

synon$duplicate <- ifelse(synon$Position %in% nonsynon_pos, "yes", "no")

aggregate(synon$Pos, by=list(synon$duplicate), length)
#  Group.1      x
# 1      no 827102
# 2     yes   4988

synon_nodups <- subset(synon, duplicate=="no")

head(synon_nodups[,c(1:2)])

write.table(synon_nodups[,c(1:2)], "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positionsNoDups.txt", sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
Removed 4,988 positions that are also denoted as nonsynonmous

Use this positions file to subset the VCF
```bash
wc -l Synonymous_positions.txt # 832090
wc -l Synonymous_positionsNoDups.txt #827102

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_synonymous' Subset_vcf.sh # Submitted batch job 4335804

```
