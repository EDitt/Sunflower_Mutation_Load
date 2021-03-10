# Preparing SNP files for GWAS Analyses

### 1. Convert to Plink file format
plink requires .ped and .map files

Used VCF_convert.sh script
```bash
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',OUT_PREFIX='/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2' VCF_convert.sh # 1905764

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
