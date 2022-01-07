# PCA using SNP data


### 1.) Filter VCF file
- Remove singletons
- Remove sites with large amount of missing data (>10%)- previously I had filtered to remove sites with >20% missing)

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.13-GCC-8.3.0
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
outputdir=/scratch/eld72413/SAM_seq/PCA

bcftools filter -e 'F_MISSING > 0.1' --threads 4 $vcf -Ou | bcftools view --min-ac 2[:minor] > ${outputdir}/Sunflower_SAM_SNP_Calling_PCAfilter.vcf

```

### 2. Convert to Plink and LD prune
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/Plink
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf',OUT_PREFIX='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter' VCF_convert.sh # Submitted batch job 7977411

# calculate r^2 among pairs of SNPs, parallelizing across chromosomes
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses

sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
File_Prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
Kb_Window_Size='1000',\
NumVariant_Windows='1000',\
MinR2_Window='0',\
Output_Dir='/scratch/eld72413/SAM_seq/PCA' LD_stats.sh # submitted batch job 7983974
```


```bash
module load EIGENSOFT/7.2.1-foss-2019b
```