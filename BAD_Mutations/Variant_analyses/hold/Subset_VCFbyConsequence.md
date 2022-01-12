# Subset VCF file by category

### Minor Allele Frequency (un-polarized)





#########################
# Subset VCF file

Subset vcf file to make separate vcfs of deleterious (for both reference and alternate alleles) and tolerated
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

# deleterious in reference
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Reference_DelPositions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Refdeleterious' Subset_vcf.sh # Submitted batch job 4232837

grep -v "#" SAM_Refdeleterious.vcf | wc -l # 11796

# deleterious in alternate
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Alternate_DelPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Altdeleterious' Subset_vcf.sh # Submitted batch job 4232853

grep -v "#" SAM_Altdeleterious.vcf | wc -l # 76016

# tolerated
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/ToleratedPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_tolerated' Subset_vcf.sh # Submitted batch job 4232872

grep -v "#" SAM_tolerated.vcf | wc -l # 553,720
```

# VCF of synonymous positions (excluding duplicates)
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

cd /scratch/eld72413/SAM_seq/VeP

grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | wc -l # 835,622
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u | wc -l # 832,090

grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u | awk '{$1=$1}1' FS=':' OFS='\t' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt

# use R to remove duplicates
module load R/4.0.0-foss-2019b
R
```
With the `sort -u` command, I removed the duplicates contained within synonymous positions (3,532). Using R I will remove the synonymous positions that are also non-synonymous positions (go by the most severe consequence for duplicate positions)

```R
dsnp <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
nonsynon_pos <- dsnp$Position

synon <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(synon) <- c("Chromosome", "Pos")

synon$Position <- paste0(synon$Chromosome,":",synon$Pos)

synon$duplicate <- ifelse(synon$Position %in% nonsynon_pos, "yes", "no")

aggregate(synon$Pos, by=list(synon$duplicate), length)
#  Group.1      x
# 1      no 827102
# 2     yes   4988

synon_nodups <- subset(synon, duplicate=="no")

head(synon_nodups[,c(1:2)])

write.table(synon_nodups[,c(1:2)], "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positionsNoDups.txt", sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
Removed 4,988 positions that are also denoted as nonsynonmous

Use this positions file to subset the VCF
```bash
wc -l Synonymous_positions.txt # 832090
wc -l Synonymous_positionsNoDups.txt #827102

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_synonymous' Subset_vcf.sh # Submitted batch job 4335804

# check
grep -v "#" SAM_synonymous.vcf | wc -l # 827102

```

# Subset VCF files for Derived SNPs

```bash
dsnp_data=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_Polarized.table
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/PositionsFiles

awk '{print $7,$39}' ${dsnp_data} | head
awk '{print $7}' ${dsnp_data} | head -50
# column 39: Whether allele is deleterious- Tolerated, Alternate_deleterious, Reference_deleterious
# column 7: Polarized relative to H. debilis- Ancestral_N, Alt_derived, Both_derived, Ref_derived 

## *derived* positions in reference:
awk '{if ($7=="Ref_derived" || $7=="Both_derived") {print $0}}' ${dsnp_data} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > ${OUTPUTDIR}/Reference_Derived_Positions.txt
wc -l ${OUTPUTDIR}/Reference_Derived_Positions.txt # 97,611

## *derived* positions in alternate:
awk '{if ($7=="Alt_derived" || $7=="Both_derived") {print $0}}' ${dsnp_data} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > ${OUTPUTDIR}/Alternate_Derived_Positions.txt
wc -l ${OUTPUTDIR}/Alternate_Derived_Positions.txt # 352,732
```

Subset vcf files made earlier (ref deleterious, alt deleterious, tolerated) for derived alleles

```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations
# derived deleterious in reference

# need to bgzip vcf files:
module load HTSlib/1.10.2-GCC-8.3.0
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/PositionsFiles/Reference_Derived_Positions.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',\
name='SAM_RefDerivedDeleterious' Subset_vcf.sh # Submitted batch job 4961118

grep -v "#" /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedDeleterious.vcf | wc -l # 5018

# derived deleterious in alternate
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/PositionsFiles/Alternate_Derived_Positions.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',\
name='SAM_AltDerivedDeleterious' Subset_vcf.sh # Submitted batch job #  4961121

grep -v "#" /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedDeleterious.vcf | wc -l # 51,342

```

Do the same for tolerated vcf files
```bash
# gzip vcf files:
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
bgzip -c --threads 4 SAM_tolerated.vcf > SAM_tolerated.vcf.gz
tabix -p vcf SAM_tolerated.vcf.gz

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations
# tolerated, reference derived
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/PositionsFiles/Reference_Derived_Positions.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_tolerated.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',name='SAM_RefDerivedTolerated' Subset_vcf.sh # Submitted batch job 4961496
grep -v "#" SAM_RefDerivedTolerated.vcf | wc -l # 

# tolerated, alternate derived
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/PositionsFiles/Alternate_Derived_Positions.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_tolerated.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',name='SAM_AltDerivedTolerated' Subset_vcf.sh # Submitted batch job 4961500
grep -v "#" SAM_AltDerivedTolerated.vcf | wc -l 
```

Do the same for synonymous vcf files
#### I need a different positions file!
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
bgzip -c --threads 4 SAM_synonymous.vcf > SAM_synonymous.vcf.gz
tabix -p vcf SAM_synonymous.vcf.gz

# positions files
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
cd /scratch/eld72413/SAM_seq/Polarized
module load R/4.0.0-foss-2019b
R
```

```R
load("/scratch/eld72413/SAM_seq/Polarized/AncestralDF.RData")

Ref_derived <- subset(AncestralDF, Category=="Ref_derived" | Category=="Both_derived")
write.table(Ref_derived[,c("Chromosome", "Position")], "Ref_derived_sites.txt", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

Alt_derived <- subset(AncestralDF, Category=="Alt_derived" | Category=="Both_derived")
write.table(Alt_derived[,c("Chromosome", "Position")], "Alt_derived_sites.txt", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

```
```bash
wc -l /scratch/eld72413/SAM_seq/Polarized/Ref_derived_sites.txt # 2,936,656
wc -l /scratch/eld72413/SAM_seq/Polarized/Alt_derived_sites.txt # 10,857,114

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations
# synonymous, reference derived
sbatch --export=positions='/scratch/eld72413/SAM_seq/Polarized/Ref_derived_sites.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_synonymous.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',name='SAM_RefDerivedSynonymous' Subset_vcf.sh # Submitted batch job 4976697
grep -v "#" SAM_RefDerivedSynonymous.vcf | wc -l # 159553

# synonymous, alternate derived
sbatch --export=positions='/scratch/eld72413/SAM_seq/Polarized/Alt_derived_sites.txt',\
vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_synonymous.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs',name='SAM_AltDerivedSynonymous' Subset_vcf.sh # Submitted batch job 4976698
grep -v "#" SAM_AltDerivedSynonymous.vcf | wc -l # 443172

module load BCFtools/1.10.2-GCC-8.3.0
bcftools view -H SAM_synonymous.vcf.gz | wc -l

```