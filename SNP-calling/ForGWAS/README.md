# Preparing Variant Files for GWAS

### Filtering
Slightly different filters (e.g. MAF) for GWAS set
- because I'm filtering on MAF and more stringent thresholds for missing genotypes and number of heterozygotes (10% or less), I will not filter on genotype quality, QUAL, or depth

Used the following commands to filter recalibrated VCF file with the GWAS_filters.sh script
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS

sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz',\
GEN_FASTA='/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/ForGWAS',\
TEMP_DIR='/scratch/eld72413/Tmp',\
OUT_PREFIX='Sunflower_SAM_SNP_GWAS' -o /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS/GWASfilter.%j.out -e /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/ForGWAS/GWASfilter.%j.err GWAS_filters.sh # Submitted batch job 5526811
```

No sequence dictionary found. Creating one.
Selecting Pass Sites from Variant Recalibrator
After filtering out sites that failed variant recalibrator, there are 81431704 sites left
Removing sites with more than 0.1 missing genotypes
After filtering out sites with too many low quality or missing variants, there are 54412579 sites left
Filtering out sites with more than 0.1 heterozygous genotypes
After filtering out sites with more than 0.1 heterozygous sites, there are 50283441 sites remaining
Filtering out sites with <1% minor allele frequency
After filtering out sites with <1% minor allele frequency, there are 21473569 sites remaining
Selecting biallelic sites
Indexing VCF
Normalizing VCF
Indexing final VCF

Looking at stats file(`bcftools stats`) generated from the GWAS_filters.sh script:
Final SNP set (after selecting for biallelic sites) has 20026451 variants

### Convert to Plink file format
plink requires .ped and .map files

Used VCF_convert.sh script
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/Plink
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/ForGWAS/Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM.vcf.gz',OUT_PREFIX='/scratch/eld72413/SAM_seq/ForGWAS/Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM' VCF_convert.sh # Submitted batch job 5555917
```

### Missing Data
Generate a list of genotyping/missingness rate statistics
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l
module load PLINK/1.9b_5-x86_64

cd /scratch/eld72413/SAM_seq/ForGWAS
plink --file Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM --missing --allow-extra-chr

awk '{if ($3==0) {print $0}}' plink.lmiss | wc -l # 1,177,639


### what about the 261 genotypes in the salt set?
plink --file /scratch/eld72413/SAM_seq/ForGWAS/Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM \
--keep /scratch/eld72413/SAM_seq/Plink/Genotypes_261set.txt \
--missing --allow-extra-chr

awk '{if ($3==0) {print $0}}' plink.lmiss | wc -l # 1,363,909

# how many are under 1% missing data (missing in only 1 or 2 genotypes)?
awk '{if ($5<0.01) {print $0}}' plink.lmiss | wc -l # 6,054,250
# (3,645,537 less than 0.05% (or 1 genotype))

# how many are under 2% missing data (missing in 1-5 genotypes)?
awk '{if ($5<0.02) {print $0}}' plink.lmiss | wc -l # 11,348,014
```

### Subset without missing data
Will first test out the subset without any missing data (while I work on imputation)
- Need to first filter to 261 genotypes from salt experiment, then filter out missing data
```bash
tmux new -s vcf_filter  ## ss-sub4 partition
module load BCFtools/1.13-GCC-8.3.0
cd /scratch/eld72413/SAM_seq/ForGWAS/

# make text file with genotype names
awk '{print $1}' /scratch/eld72413/SAM_seq/Plink/Genotypes_261set.txt > Genotypes_261set_forBCFtools.txt
# make sure names are similar
bcftools query -l /scratch/eld72413/SAM_seq/ForGWAS/Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM.vcf.gz
sed -i 's/RHA415-4/RHA415-4_PPN251/g' Genotypes_261set_forBCFtools.txt
sed -i 's/NMS373/NMS373_PPN136/g' Genotypes_261set_forBCFtools.txt
#Hopi_PPN285
#PI_531071
#SF_33

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l

bcftools view -S /scratch/eld72413/SAM_seq/ForGWAS/Genotypes_261set_forBCFtools.txt \
/scratch/eld72413/SAM_seq/ForGWAS/Sunflower_SAM_SNP_GWAS_BIALLELIC_NORM.vcf.gz \
--force-samples \
-Ou | \
bcftools view -g ^miss -Oz -o /scratch/eld72413/SAM_seq/ForGWAS/Subset261/Sunflower_SAM_SNP_GWAS_261Subset_noMISS.vcf.gz

# Error: subset called for sample that does not exist in header: "PPN053".  Use "--force-samples" to ignore this error.
### after adding flag:
# Warn: subset called for sample that does not exist in header: "PPN053"... skipping
du -h Sunflower_SAM_SNP_GWAS_261Subset_noMISS.vcf.gz # 1.5G
bcftools stats Sunflower_SAM_SNP_GWAS_261Subset_noMISS.vcf.gz # 1363909 SNPs

# convert to Plink format
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/Plink
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/ForGWAS/Subset261/Sunflower_SAM_SNP_GWAS_261Subset_noMISS.vcf.gz',OUT_PREFIX='/scratch/eld72413/SAM_seq/ForGWAS/Subset261/Sunflower_SAM_SNP_GWAS_261Subset_noMISS' VCF_convert.sh # Submitted batch job 5581114

```