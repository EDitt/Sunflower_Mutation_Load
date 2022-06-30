# Number of dSNPs per genotype

### Setup
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh
```

### Average dSNPs per individual
First, get total number of called sites
```bash
# Already ran command on entire vcf file: "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllVariantStats.txt"

### note: the bcftools stats "-s -" flag means to include all samples
bcftools stats -s - ${VCF} | grep "PSC" > ${OUT_DIR}/GenotypeInfo/All_SNP_stats.txt

```

#### Parse total number of dSNPs per genotype
Count Number of Alt/Ref Alleles per genotype
(using bcftools)
```bash
mkdir ${OUT_DIR}/IntermediateFiles/dSNPCounts_perGeno
### issues with passing comma separated list as variable"Argument list too long"
#positions=$(awk 'BEGIN{FS="\t"; OFS=":"}; {if ($39 == "Alternate_deleterious" && $7 == "Alt_derived") {print $8,$9}}' ${DSNP_DATA} | paste -sd,)

awk 'BEGIN{FS=OFS="\t"}; {if ($39 == "Alternate_deleterious") {print $8,$9}}' ${DSNP_DATA} > ${OUT_DIR}/IntermediateFiles/AltDeleteriousPositions.txt # 76,095 (51,009 when including only those derived relative to H. debilis)
positions="${OUT_DIR}/IntermediateFiles/AltDeleteriousPositions.txt"

### note: the bcftools stats "-s -" flag means to include all samples
bcftools view -Oz ${VCF} -R ${positions} | \
bcftools stats -s - | \
grep "PSC" > ${OUT_DIR}/IntermediateFiles/dSNPCounts_perGeno/SampleCounts_AltDeleterious.txt

awk 'BEGIN{FS=OFS="\t"}; {if ($39 == "Reference_deleterious") {print $8,$9}}' ${DSNP_DATA} > ${OUT_DIR}/IntermediateFiles/RefDeleteriousPositions.txt # 11,796 (4954 when including only those derived relative to H. debilis)
positions="${OUT_DIR}/IntermediateFiles/RefDeleteriousPositions.txt"

bcftools view -Oz ${VCF} -R ${positions} | \
bcftools stats -s - | \
grep "PSC" > ${OUT_DIR}/IntermediateFiles/dSNPCounts_perGeno/SampleCounts_RefDeleterious.txt

```

# Add Deleterious Mutations from both lists for each genotype

Use R function to output a text file with counts of dSNPs for all genotypesa
```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
# all_stats[c("sample", "CalledGenotypes", "nMissingTotal")],

Full_stats <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/All_SNP_stats.txt")
colnames(Full_stats) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissingTotal")
Full_stats$CalledGenotypes <- Full_stats$nRefHom + Full_stats$nNonRefHom + Full_stats$nHets

dSNP_table <- Genotype_dSNP_count("/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/dSNPCounts_perGeno",
	"RefDeleterious", "AltDeleterious", 
	Full_stats[,c("sample", "CalledGenotypes", "nMissingTotal", "nRefHom", "nNonRefHom", "nHets")])

write.table(dSNP_table, "/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/All_dSNP_stats.txt", sep = "\t", quote=FALSE, row.names=FALSE)
```


### Total homozygous and heterogzygous load