# Number of dSNPs per genotype

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
wc -l ToleratedPositions.txt # 556,577 (740 duplicate positions represented in deleterious set already removed)
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

grep -v "#" SAM_tolerated.vcf | wc -l #
```

# Count Number of Alt/Ref Alleles per genotype
(using bcftools)
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_Refdeleterious.vcf > RefDeleterious_StatsperSample.txt
grep "PSC" RefDeleterious_StatsperSample.txt > RefDeleterious_perSampleCounts.txt

bcftools stats -s - SAM_Altdeleterious.vcf > AltDeleterious_StatsperSample.txt
grep "PSC" AltDeleterious_StatsperSample.txt > AltDeleterious_perSampleCounts.txt
```

# Add Deleterious Mutations from both lists for each genotype
```bash
module load R/4.0.0-foss-2019b
R
```
```R
DelRef_data <- read.table("RefDeleterious_perSampleCounts.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(DelRef_data) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")

#DelRef_data <- subset(DelRef_data, select = -c(PSC, id, nIndels, nHapRef, nHapAlt))

DelRef_data$Total <- DelRef_data$nRefHom + DelRef_data$nNonRefHom + DelRef_data$nHets + DelRef_data$nMissing

DelAlt_data <- read.table("AltDeleterious_perSampleCounts.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(DelAlt_data) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")

#DelAlt_data <- subset(DelAlt_data, select = -c(PSC, id, nIndels, nHapRef, nHapAlt))

DelAlt_data$Total <- DelAlt_data$nRefHom + DelAlt_data$nNonRefHom + DelAlt_data$nHets + DelAlt_data$nMissing

# merge datasets
Deldata <- merge(DelRef_data[,c(3:6,14:15)], 
					DelAlt_data[,c(3:6,14:15)],
					by="sample",
					suffixes=c("_DelRef", "_DelAlt"))
length(Deldata$sample) #288

# Number of Deleterious Alleles:
Deldata$TotNum_dSNPs <- 2*Deldata$nRefHom_DelRef +
						Deldata$nHets_DelRef +
						2*Deldata$nNonRefHom_DelAlt +
						Deldata$nHets_DelAlt

# Number of homozygous deleterious sites:
Deldata$nHom_dSNPs <- Deldata$nRefHom_DelRef +
					Deldata$nNonRefHom_DelAlt

# Number of heterozygous deleterious sites:
Deldata$nHet_dSNPs <- Deldata$nHets_DelRef +
					Deldata$nHets_DelAlt

# total number of SNPs (for proportions?)
Deldata$Total <- Deldata$Total_DelRef + Deldata$Total_DelAlt
Deldata$TotalMissing <- Deldata$nMissing_DelRef + Deldata$nMissing_DelAlt

# save file
write.table(Deldata[,c(1,12:16)], "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Genotype_dSNP_counts.txt", sep = "\t", quote=FALSE, row.names=FALSE)
```



Previous code, not updated yet
#########
## tolerated positions
grep "Tolerated" dsnp_data.table | awk '{print $17}' | awk '{$1=$1}1' FS=':' OFS='\t' > Tolerated_positions.txt

cd /scratch/eld72413/SAM_seq/VeP
grep -v "#" fullsam_synon.txt | awk '{print $2}' | awk '{$1=$1}1' FS=':' OFS='\t' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt
```
