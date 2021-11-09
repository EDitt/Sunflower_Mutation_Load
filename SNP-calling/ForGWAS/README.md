# Preparing Variant Files for GWAS

### Filtering
Slightly different filters (MAF) for GWAS set

Used the following commands to filter recalibrated VCF file with the GWAS_filters.sh script
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS

sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz',\
GEN_FASTA='/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/ForGWAS',\
TEMP_DIR='/scratch/eld72413/Tmp',\
OUT_PREFIX='Sunflower_SAM_SNP_GWAS' -o /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS/GWASfilter.%j.out -e /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS/GWASfilter.%j.err GWAS_filters.sh # Submitted batch job 5526751
```

