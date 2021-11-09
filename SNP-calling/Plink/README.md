# Preparing SNP files for GWAS Analyses

### 1. Convert to Plink file format
plink requires .ped and .map files

Used VCF_convert.sh script
```bash
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',OUT_PREFIX='/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2' VCF_convert.sh # 1905764 # 5507827

```

### Created a IBS matrix
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l

module load PLINK/1.9b_5-x86_64
plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2 \
--cluster \
--matrix \
--allow-extra-chr \
--out /scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2
```

128914 MB RAM detected; reserving 64457 MB for main workspace.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (37129915 variants, 288 people).
--file: /scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2-temporary.bed
+ /scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2-temporary.bim +
/scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2-temporary.fam
written.
37129915 variants loaded from .bim file.
288 people (0 males, 0 females, 288 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
/scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2.nosex .
Using up to 47 threads (change this with --threads).
Before main variant filters, 288 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.894318.
37129915 variants and 288 people pass filters and QC.
Note: No phenotypes present.
Distance matrix calculation complete.
IBS matrix written to
/scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2.mibs , and IDs to
/scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2.mibs.id .
Clustering... done.                        
Cluster solution written to 100%


### Pairwise IBD estimation

```bash
plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2 \
--allow-extra-chr \
--Z-genome \
--out /scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2
```

### MDS Plot

```bash
plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2 \
--allow-extra-chr \
--read-genome /scratch/eld72413/SAM_seq/Plink/IBS/Sunflower_SAM_HA412v2.genome.gz \
--cluster --mds-plot 4 \
--out /scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2
```

.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (37129915 variants, 288 people).
--file: /scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2-temporary.bed
+ /scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2-temporary.bim +
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2-temporary.fam
written.
37129915 variants loaded from .bim file.
288 people (0 males, 0 females, 288 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 288 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.894318.
37129915 variants and 288 people pass filters and QC.
Note: No phenotypes present.
Clustering... done.                        
Cluster solution written to 100%
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2.cluster1 ,
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2.cluster2 , and
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2.cluster3 .
Performing multidimensional scaling analysis (SVD algorithm, 4
dimensions)... done.
MDS solution written to
/scratch/eld72413/SAM_seq/Plink/MDS/Sunflower_SAM_HA412v2.mds .


```R
library(ggplot2)
library(ggthemes)

setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/Results")
mds <- read.table("Sunflower_SAM_HA412v2.mds", header=T)
length(mds$FID) #288
plot(mds$C1 ~ mds$C2)

# grouping info
line_info <- read.csv("/Users/eld72413/Documents/GitHub/Sunflower_Mutation_Load/BAD_Mutations/Line_Info.csv", header=T)

# fix names of weird samples
mds$name <- mds$IID

mds$name[mds$name=="531071"] <- "PI_531071"
mds$name[mds$name=="PPN285"] <- "Hopi_PPN285"
mds$name[mds$name=="PPN136"] <- "NMS373_PPN136"

#merge datasets
Mds_wInfo <- merge(line_info[,c(8:10)], mds, by.x = "VCF_line_name", by.y = "name")
length(Mds_wInfo$VCF_line_name) #286
Mds_wInfo$Class1 <- factor(Mds_wInfo$Class1)
Mds_wInfo$Class2 <- factor(Mds_wInfo$Class2)
str(Mds_wInfo)


### plot

p <- ggplot(Mds_wInfo, aes(x=C1, y=C2))
p + geom_point(aes(color=Class1, shape=Class2)) + theme_minimal()
p + geom_point(aes(color=Class2)) + theme_minimal()

```

### Missing Data

Generate a list of genotyping/missingness rate statistics
```bash
module load PLINK/1.9b_5-x86_64
plink --file Sunflower_SAM_HA412v2 --missing --allow-extra-chr

awk '{if ($3==0) {print $0}}' plink.lmiss | wc -l #819

# filter at 10% missing data
awk '{if ($5 < 0.1) {print $0}}' plink.lmiss | wc -l # 18,119,803
```
37129915 variants loaded from .bim file.
288 people (0 males, 0 females, 288 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 288 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.894318.
--missing: Sample missing data report written to plink.imiss, and variant-based
missing data report written to plink.lmiss.


Could also first do filtering for missing data and MAF (as done for prior SNP set before imputing)

Previously, I had filtered for missing data > 20%, >20% heterogygous individuals, and no minor allele frequency filter
For GWAS, I will filter out sites with >10% missing data, >10% heterogyzous individuals, and <1% MAF
```bash
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/Plink',\
OUT_PREFIX='Sunflower_SAM_HA412v2' \
GWAS_filters.sh
# Submitted batch job 5508504

# Convert filtered VCF to plink format
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_missing_filtered.vcf',OUT_PREFIX='/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_FILTERED' VCF_convert.sh # Submitted batch job 5511088
```
After filtering out sites with MAF <1% and heterozygosity >10%, there are 15096881 sites left
After filtering out sites with >10% missing data, there are 6899183 sites left


Filter for Genotypes used in GWAS?
```R
# what about the subset in the 261?
Geno261 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Plink/Andries_Info/Genotypes in 261 set.txt", header=FALSE)
colnames(Geno261) <- c("SAM", "SAM2", "Col1", "Col2", "Col3", "Col4")

library(stringr)
Geno261$IID <- str_replace(Geno261$SAM, "SAM", "PPN")

length(plink_geno[which(plink_geno$IID %in% Geno261$IID),"IID"]) # N = 260
plink_geno[which(!Geno261$IID %in% plink_geno$IID),] # PPN047

# how many are in the "GoldStandard"?
length(Geno261[which(Geno261$IID %in% GoldStandard),"IID"]) # 58
```
These aren't necessarily biased towards the most highest-coverage genotypes

-----

### Filters
- 5 % MAF
- less than 10% heterozygous sites

Andries' code to calculate haplotype blocks (the following is for chromosome 1)
```bash
./plink --tped XRQv1_412_285_filtered.tped --tfam XRQv1_412_285_filtered.tfam --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 100000 --blocks-strong-lowci 0.7005 --out CHR1_285 --allow-extra-chr --chr Ha412HOChr01 --blocks-inform-frac 0.9
```

`plink --file mydata` is the same as `plink --ped mydata.ped --map mydata.map`

generate a list of MAF for each SNP
`plink --file data --freq`

IBS Similarity matrix
`plink --file mydata --cluster --matrix`

Pairwise IBD estimation
`plink --file mydata --genome` creates a plink.genome file
OR
`plink --file mydata --Z-genome` to create a compressed file
