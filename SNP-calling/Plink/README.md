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
