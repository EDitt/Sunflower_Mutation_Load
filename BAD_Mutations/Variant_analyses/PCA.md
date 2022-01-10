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
# now contains 13,973,358 variants
```

### 2. Convert to Plink and LD prune
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/Plink
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter.vcf',OUT_PREFIX='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter' VCF_convert.sh # Submitted batch job 7997074

# calculate r^2 among pairs of SNPs, parallelizing across chromosomes
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses

sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
File_Prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
Kb_Window_Size='1000',\
NumVariant_Windows='1000',\
MinR2_Window='0',\
Output_Dir='/scratch/eld72413/SAM_seq/PCA' LD_stats.sh # submitted batch job 7983974

sbatch --export=Output_Dir='/scratch/eld72413/SAM_seq/PCA',\
File_Prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
Window_Size='1',Step_Size='1',Rsquared='0.9' LD_prune.sh # Submitted batch job 7991417
# Got warning message: step size should be 1 when window size is in kb
# Pruning complete.  9037810 of 13973358 variants removed (at R^2 of 0.8)
# Pruning complete.  8460435 of 13973358 variants removed (at R^2 of 0.9)


plink --file /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter \
--extract /scratch/eld72413/SAM_seq/PCA/PrunedLists/Plink_0.9.prune.in \
--allow-extra-chr \
--out /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_Pruned_R2_0.9 \
--recode

#--make-bed
```

### 3. Convert from Plink format to Eigensoft format
```bash
module load EIGENSOFT/7.2.1-foss-2019b

# first need to convert

# created a parfile:
parfile="/scratch/eld72413/SAM_seq/PCA/EigenstratFiles/par.PED.EIGENSTRAT"
# chromosome names MUST be integers
# take out leading text and leading zeros (where applicable); change all scaffold chromosome names to 18

awk 'BEGIN{FS=OFS="\t"}; {gsub("Ha412HOChr","",$1)}1' /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_Pruned_R2_0.9.map | sed 's/^0//' | sed 's/^0c...../18/' > /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_Pruned_R2_0.9_ChrEdit.map

convertf -p $parfile
# bad chrom (before I changed the naming structure)

## snp order check fail; snp list not ordered

## all individuals set ignore.  Likely input problem (col 6)
## resetting all individual...
## genotype file processed
## numvalidind:    288  maxmiss: 288001
## eigenstrat output
## ##end of convertf run



```


do I wnat to change column 6 and alter the "outputgroup" parameter?

```bash
awk '{print $1,$2,$3,$4,$5,$6}' /scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_Pruned_R2_0.9.ped | head

```

### 4. Calculate PCA using Eigensoft (SmartPCA)
```bash
#smartpca -p PARAMS_FILE > smartpca.log
Dir=/scratch/eld72413/SAM_seq/PCA/EigenstratFiles
smartpca.perl -i ${Dir}/Sunflower_SAM.eigenstratgeno \
-a ${Dir}/Sunflower_SAM.snp \
-b ${Dir}/Sunflower_SAM.ind \
-k 10 \
-o ${Dir}/Sunflower_SAM.pca \
-p ${Dir}/Sunflower_Sam_Unlabeled.plot \
-e ${Dir}/Sunflower_SAM.eval \
-l ${Dir}/Sunflower_SAM_PCA.log \
-m 5 \
-t 10 \
-s 6.0

##### 31 genotypes were removed?

# smartpca -p /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM.pca.par >/scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM_PCA.log
# ploteig -i /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM.pca.evec -c 1:2  -p ???  -x  -y  -o /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_Sam_Unlabeled.plot.xtxt 
# sh: ploteig: command not found
# evec2pca.perl 10 /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM.pca.evec /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM.ind /scratch/eld72413/SAM_seq/PCA/EigenstratFiles/Sunflower_SAM.pca
```