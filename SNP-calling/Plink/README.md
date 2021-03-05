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

