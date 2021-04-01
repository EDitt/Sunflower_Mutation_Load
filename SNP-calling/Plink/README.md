# Preparing SNP files for GWAS Analyses

### 1. Convert to Plink file format
plink requires .ped and .map files

Used VCF_convert.sh script
```bash
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',OUT_PREFIX='/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2' VCF_convert.sh # 1905764

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
