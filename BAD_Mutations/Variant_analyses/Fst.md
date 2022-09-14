# Calculate Fst across genome

### Subset VCF file:

Starting with VCF file from Step 1 of PCA.md with singletons and sites missing more than 10% of data are removed:
  `bcftools filter -e 'F_MISSING > 0.1' --threads 4 ${VCF} -Ou | bcftools view --min-ac 2[:minor] > /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf`

### Gzip/Index VCF file
```bash
VCF_sub="/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf"

bgzip -c ${VCF_sub} > /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz
```

### Make lists of HA and RHA / Oil and NonOil lines
```bash
# HA lines
awk -F',' '{if ($11=="HA" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_INFO | wc -l # 145
awk -F',' '{if ($11=="HA" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_INFO > ${OUT_DIR}/Fst/HA_lines.txt

# RHA lines
awk -F',' '{if ($11=="RHA" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_INFO | wc -l # 112
awk -F',' '{if ($11=="RHA" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_INFO > ${OUT_DIR}/Fst/RHA_lines.txt

# Oil lines
awk -F',' '{if ($12=="Oil" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_INFO | wc -l #193
awk -F',' '{if ($12=="Oil" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_INFO > ${OUT_DIR}/Fst/Oil_lines.txt

# NonOil lines
awk -F',' '{if ($12=="NonOil" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_INFO | wc -l #95
awk -F',' '{if ($12=="NonOil" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_INFO > ${OUT_DIR}/Fst/NonOil_lines.txt
```

### Calculate Fst for 1 Mbp bins across genome (parallelize among chromosomes)

HA/RHA
```bash
sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz',\
LIST1='/scratch/eld72413/SAM_seq/dSNP_results/Fst/HA_lines.txt',\
LIST2='/scratch/eld72413/SAM_seq/dSNP_results/Fst/RHA_lines.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/dSNP_results/Fst/HA_RHA' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Fst_calc.sh # Submitted batch job 12666566

awk '{if ($5 > 0.1) {print $0}}' HA_RHA_Ha412HOChr01.windowed.weir.fst | wc -l # 3 (out of 161)
awk '{if ($5 > 0.1) {print $0}}' HA_RHA_Ha412HOChr10.windowed.weir.fst | wc -l # 184 (out of 192)
```

Oil/NonOil
```bash
sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz',\
LIST1='/scratch/eld72413/SAM_seq/dSNP_results/Fst/Oil_lines.txt',\
LIST2='/scratch/eld72413/SAM_seq/dSNP_results/Fst/NonOil_lines.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/dSNP_results/Fst/Oil_NonOil' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Fst_calc.sh # Submitted batch job 12666622
```
See script `Plots/Fst_Genome.R` for graphing

### Calculate Nucleotide diversity/Tajima's D across genome for 1 Mbp bins for HA vs. RHA
look at diversity for HA and RHA? & Oil/Non-Oil to identify SS?

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

cd /scratch/eld72413/SAM_seq/dSNP_results/Fst
mkdir Pi_heterotic

group=HA
group=RHA
awk -v var="$group" -F',' '{if ($11==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO > ${OUT_DIR}/IntermediateFiles/${group}_genotypes.txt

sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
GenoList='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/HA_genotypes.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/dSNP_results/Fst/Pi_heterotic/HA_Mbp1' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Pi_TajD_calc.sh # 13838554

sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
GenoList='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/RHA_genotypes.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/dSNP_results/Fst/Pi_heterotic/RHA_Mbp1' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Pi_TajD_calc.sh #  13838643 (prev 12680201)- resubmitted due to missing data for RHA on chr 17

```





##### no longer using the code below
```bash

#VCF_sub=/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz actually I don't want to use this

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

vcftools --gzvcf ${VCF_sub} \
--keep ${OUT_DIR}/IntermediateFiles/${group}_genotypes.txt \
--window-pi 1000000 \
--window-pi-step 1000000 \
--out /scratch/eld72413/SAM_seq/dSNP_results/Fst/Pi_heterotic/Pi_1Mbp_${group}

vcftools --gzvcf ${VCF_sub} \
--keep ${OUT_DIR}/IntermediateFiles/${group}_genotypes.txt \
--TajimaD 1000000 \
--out /scratch/eld72413/SAM_seq/dSNP_results/Fst/Pi_heterotic/Pi_1Mbp_${group}


```

```bash
sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf.gz',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/dSNP_results/Fst/Mbp1' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Pi_TajD_calc.sh # Submitted batch job 12669005, 12669096
```

