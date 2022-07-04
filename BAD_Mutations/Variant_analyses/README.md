# dSNP Analysis

## Navigation: Jump to Section

- [Polarize Ancestral State](#polarize-ancestral-state)
- [SNP Annotation Classes](#snp-annotation-classes)
- [Create SNP Table](#create-snp-table)
- [Site Frequency Spectra](#site-frequency-spectra)
- [Germplasm Patterns](#germplasm-patterns)



- [Genomic Patterns](#genomic-patterns)
	- [Recombination across genome](#recombination-across-genome)
	- [Binning SNP classes](#binning-snp-classes)
- 

---

The compiled dSNP report was merged with VeP output (missense variants only) and processed using criteria detailed in `${REPO_DIR}/2.BAD_Mutations/Post_processing.md`

Resulting file: `/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table`

## Polarize Ancestral State

Created ancestral fasta file with ANGSD. See: ANGSD directory

See `${REPO_DIR}/Variant_analyses/Polarize_SNPs.md` for process of creating ancestral state table.
Output a table with ancestral state calls at 13,524,630 variant positions: `/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt`

Merged ancestral states with dsnp_data.table file to get: `/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_Polarized.table`

Subset VCF file to sites with ancestral state calls
```bash
awk 'BEGIN{FS=OFS="\t"}; {if ($3 != "N") {print $1,$2}}' ${ANCESTRAL_STATE} > ${OUT_DIR}/IntermediateFiles/AncestralStateCalls.txt # 13,524,630
sbatch --export=positions='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/AncestralStateCalls.txt',\
vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles',\
name='SAM_AncestralStateCalls' /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/Subset_vcf.sh # Submitted batch job 8178667
```

## SNP Annotation Classes
Parse VeP Table and Predictions output to obtain annotation classes for all SNPs

See: `${REPO_DIR}/Variant_analyses/Variant_class_numbers.md` for commands used to parse the prediction output and VeP output to get the variants for each class
Tab-delimited lists of variants positions for each annotation class found in: `/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/FinalPositionFiles`

Save list of SNPs & annotation
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

Rscript "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SNP_annotation.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/FinalPositionFiles" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_Positions.txt"

```

## Create SNP Table
Table with SNP Frequency, Ancestral State, Genotype Calls, Annotation, for all SNPs
** note: I changed R script after using this code so will have to change the input arguments for the Rscript commands to work again

Use bcftools to output a table with Allele Position, Reference and Alternate alleles, alternate allele count, total allele count, alternate allele frequency
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' ${VCF} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt

# remove N's from ancestral state table to save memory
awk '{if ($3 != "N") {print $0}}' ${ANCESTRAL_STATE} > /scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt

# Use R script to combine information
Rscript "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/Variant_Table.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/FinalPositionFiles" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info_new.txt"

awk 'NR>1 {print $1,$2}' $SNP_INFO | wc -l # 36,708,692
awk 'NR>1 {print $1,$2}' $SNP_INFO | sort -u | wc -l # 36,708,692

awk 'NR>1 {print $1,$2}' All_SNP_Info_new.txt | wc -l # 37,120,112
awk 'NR>1 {print $1,$2}' All_SNP_Info_new.txt | sort -u | wc -l # 37,120,112
```

## Site Frequency Spectra
Get frequency bins for graphing (folded and unfolded)

To look at derived dSNPs, I need to remove variants for which the derived variant is not the deleterious variant
```bash
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

# use new SNP_INFO file (after edits)
SNP_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info_new.txt

## Folded SFS:
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${SNP_INFO}" \
"0.5" \
"0.05" \
"MAF" \
"/scratch/eld72413/SAM_seq/dSNP_results/MAF_Bins.txt"

## UnFolded SFS:
#	=> First remove derived frequency information for deleterious SNPs for variants at which the derived variant != the deleterious variant
awk '{if (($39 == "Alternate_deleterious" && $7 == "Ref_derived") ||  ($39 == "Reference_deleterious" && $7 == "Alt_derived")) {print $0}}' ${DSNP_DATA} | wc -l # 6144
awk 'BEGIN{FS=OFS="\t"}; {if (($39 == "Alternate_deleterious" && $7 == "Ref_derived") ||  ($39 == "Reference_deleterious" && $7 == "Alt_derived")) {print $8,$9}}' ${DSNP_DATA} | sort -V > ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt

# remove these lines from dSNP_data from Unfolded SFS calculations:
grep -w -Ff ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt ${SNP_INFO} | wc -l # 6141 (3 sites not present in the SNP_INFO table)
grep -v -w -f ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt ${SNP_INFO} | sort -V > ${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt

## Unfolded SFS:
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt" \
"1.0" \
"0.05" \
"Derived_Freq" \
"/scratch/eld72413/SAM_seq/dSNP_results/DerivedFreq_Bins.txt"

# rm ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt
```

See Plots directory for SFS plot

Do I want to remove singletons?
```bash
# take out the sites with allele count of 1 (heterozygous in only 1 sample due to sample pooling)
Pos_info=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt
wc -l  $Pos_info # 37,129,915
awk '{if ($7<0.05) {print $0}}' $Pos_info | wc -l # 28,667,994 in this frequency class

awk '{if ($5==1) {print $0}}' $Pos_info | wc -l # 7,805,470
awk '{if ($5==2) {print $0}}' $Pos_info | wc -l # 7,625,046
awk '{if ($5>1) {print $0}}' $Pos_info > All_alleleFreqInfo_noSingleton.txt # 29,324,445

# what about cases where reference exhibits the minor allele only in heterozygote form?
```

----- 

## Germplasm patterns
Spreadsheet `scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/LineKeywINFO.csv` contains grouping information for each genotype to match with Variant ID's from VCF

See `${REPO_DIR}/Variant_analyses/PCA.md` for information on the creation of a PCA for genotypes in the SAM pop

### Heterotic Group Differentiaton

Plotted Fst across the genome. See: `${REPO_DIR}/Variant_analyses/Fst.md`

### Number of dSNPs per Genotype

See GenotypeLoadNums.md for getting the number of dSNPs per genotype

### Numbers of Derived Variants for all Genotypes for different variant classes

Use scripts to output a table with number of derived variants across all variant classes
```bash
awk 'NR>1 {print $13}' ${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt | sort -u #  annotation classes

# script to return nRefHom, nNonRefHom, nHets, etc. for all genotypes (for all variant classes) for those that are polarized
sbatch --export=Table='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/SNPINFO_ForUnfolded.txt',\
vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/SampleCounts' /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/GenotypeStats_Class.sh # Submitted batch job 12443241

# I moved just the deleterious, tolerated and nonsynonymous into a sub-directory "ToUse" to combine this information 
# This script converts the numbers to derived and ancestral and combines across the variant classes
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm/Derived_Variant_Numbers.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/SampleCounts/intermediates/ToUse" \
"/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/Annotation_VariantStats.txt" \
"derivedCounts" \
"_Ref_" \
"_Alt_"
```

In `/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/SampleCounts/intermediates`
	- lists of sample, nRefHom, nNonRefHom, nHets, etc. for each variant class
	- (deleterious, synonymous, and tolerated are in subdirectory 'ToUse'

In `/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo`, file: `Annotation_VariantStats.txt` which contains the number of ancestral homozygous, number derived homozygous, number of heterozygotes, number missing for each line for deleterious, synonymous, and tolerated variant classes

See graph of derived dSNPs/sSNPs for different germplasm groups at: Plots/Germplasm_boxplot.R

### Frequency of derived SNPs for different germplasm groups
```bash
sbatch --export=SAM_INFO='/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/GroupFreqs',\
snps_remove='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/ToRemove_NotDerived.txt' /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Group_Freqs.sh # Submitted batch job 12461287

# this returns tables for all groups of variants with the 
```

In `/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/GroupFreqs` 
	- lists of each SNP with num ref/alt alleles; MAF, derived freq, for each variant class, for each germplasm group
	- in subdirectory `intermediates`, `_info_ForUnfolded.txt` is the same as the above but without deleterious alleles that are not derived relative to Debilis
	- `_DerivedFreq_Bins.txt` has the frequency bins for all variant classes for all germplasm groups

### Proportion of heterozygous dSNPs/sSNPs across different MAF classes
To see if the greater heterozygosity in dSNPs is due only to lower frequencies
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

SNP_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info_new.txt

awk 'BEGIN{FS=OFS="\t"}; {if ($13 == "AllDel") {print $1, $2-1, $2}}' $SNP_INFO > $OUT_DIR/IntermediateFiles/AllDel_Positions.bed

awk 'BEGIN{FS=OFS="\t"}; {if ($13 == "SynonymousNodups") {print $1, $2-1, $2}}' $SNP_INFO > $OUT_DIR/IntermediateFiles/AllSynon_Positions.bed

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8

gatk IndexFeatureFile \
     -F ${VCF}

SNPset=AllDel
SNPset=AllSynon

gatk VariantsToTable \
     -V ${VCF} \
     -L ${OUT_DIR}/IntermediateFiles/${SNPset}_Positions.bed \
     -F CHROM -F POS -F HET -F HOM-REF -F HOM-VAR -F NCALLED \
     -O ${OUT_DIR}/GenotypeInfo/Heterozygosity/${SNPset}_Het_table.txt
```
see: heterozygosity.R

#### Private/Shared
I will look at the proportion of shared versus private SNPs for deleterious and synonymous across the different frequency classes as well as plot the frequencies of all SNPs for both groups
1.) Are there more private deleterious alleles than expected (based on private synonymous SNPs)? (controlling for frequency class)

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l #tmux window in ss-sub4
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

mkdir -p ${OUT_DIR}/GenotypeInfo/HeteroticGroup

awk 'NR>1 {print $13}' ${OUT_DIR}/GenotypeInfo/GroupFreqs/intermediates/RHA-Oil_SNP_info_ForUnfolded.txt | sort -u # Variant_type strings

# reduce to save memory
awk '{if ($13=="AllDel" || $13=="SynonymousNodups" || $13=="Variant_type") {print $0}}' ${OUT_DIR}/GenotypeInfo/GroupFreqs/intermediates/RHA-Oil_SNP_info_ForUnfolded.txt > ${OUT_DIR}/GenotypeInfo/HeteroticGroup/RHA-Oil_SNP_info_ForUnfolded_reduced.txt # N= 1,009,331
awk '{if ($13=="AllDel" || $13=="SynonymousNodups" || $13=="Variant_type") {print $0}}' ${OUT_DIR}/GenotypeInfo/GroupFreqs/intermediates/HA-Oil_SNP_info_ForUnfolded.txt > ${OUT_DIR}/GenotypeInfo/HeteroticGroup/HA-Oil_SNP_info_ForUnfolded_reduced.txt
awk '{if ($13=="AllDel" || $13=="SynonymousNodups" || $13=="Variant_type") {print $0}}' ${OUT_DIR}/GenotypeInfo/GroupFreqs/intermediates/RHA-NonOil_SNP_info_ForUnfolded.txt > ${OUT_DIR}/GenotypeInfo/HeteroticGroup/RHA-NonOil_SNP_info_ForUnfolded_reduced.txt
awk '{if ($13=="AllDel" || $13=="SynonymousNodups" || $13=="Variant_type") {print $0}}' ${OUT_DIR}/GenotypeInfo/GroupFreqs/intermediates/HA-NonOil_SNP_info_ForUnfolded.txt > ${OUT_DIR}/GenotypeInfo/HeteroticGroup/HA-NonOil_SNP_info_ForUnfolded_reduced.txt
```

See: R code at Heterotic_AlleleInfo.R for data wrangling

#### Joint SFS

See: Joint_SFS.R

#### IBS patterns
Scatterplot showing relationship between pairwise IBS for synonymous versus deleterious variants.
Difference between inter and intra heterotic group crosses?

```bash
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

awk '{print $13}' $SNP_INFO | sort -u

### using the "PCA filter" Plink file which was converted from a VCF with singletons removed (see `${REPO_DIR}/Variant_analyses/PCA.md`)

# synonymous positions
type="SynonymousNodups"
awk -v var="$type" 'BEGIN{OFS=":"}; {if ($13==var) {print $1,$2}}' $SNP_INFO > /scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Synonymous_Pos4Plink.txt

sbatch --export=Input_prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
SNP_List='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Synonymous_Pos4Plink.txt',\
Window_Size='1',\
Step_Size='1',\
Rsquared='0.9',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups',\
output_prefix='Synonymous' ${REPO_DIR}/BAD_Mutations/Variant_analyses/IBS_SNP_Subset.sh # Submitted batch job 12569029

# deleterious positions
type="AllDel"
awk -v var="$type" 'BEGIN{OFS=":"}; {if ($13==var) {print $1,$2}}' $SNP_INFO > /scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Deleterious_Pos4Plink.txt

sbatch --export=Input_prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
SNP_List='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Deleterious_Pos4Plink.txt',\
Window_Size='1',\
Step_Size='1',\
Rsquared='0.9',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups',\
output_prefix='Deleterious' ${REPO_DIR}/BAD_Mutations/Variant_analyses/IBS_SNP_Subset.sh # Submitted batch job 12569031
```

Synonymous:
Pruning SNPs with R^2 more than 0.9 in 1 kb windows will remove 128706 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_Plink_0.9.prune.out variants, 
keeping 257340 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_Plink_0.9.prune.in variants
Pairwise IBS will be calculated using 257340 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_SNPsToKeep.txt SNPs

Deleterious:
Pruning SNPs with R^2 more than 0.9 in 1 kb windows will remove 4586 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_Plink_0.9.prune.out variants, 
keeping 32210 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_Plink_0.9.prune.in variants
Pairwise IBS will be calculated using 32210 /scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/Intermediates/TEMP_SNPsToKeep.txt SNPs

Combine and analyze with R
```R
# srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")

IBS_table <- ImportFilesAsDf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups", "_IBS.txt", "Prefix", c("Genotype1", "Genotype2", "IBS")) # no prefix # N=82656

SAM_info <- read.csv("/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/LineKeywINFO.csv", header=T)

IBS_GroupInfo <- merge(IBS_table, SAM_info[,c(9,11)], by.x="Genotype1", by.y="SequenceName")
colnames(IBS_GroupInfo)[5] <- c("HeteroticGroup1")
IBS_GroupInfo <- merge(IBS_GroupInfo, SAM_info[,c(9,11)], by.x="Genotype2", by.y="SequenceName")
colnames(IBS_GroupInfo)[6] <- c("HeteroticGroup2")

# define cross type
IBS_GroupInfo$Cross <- ifelse((IBS_GroupInfo$HeteroticGroup1=="HA" & IBS_GroupInfo$HeteroticGroup2=="HA") |
								(IBS_GroupInfo$HeteroticGroup1=="RHA" & IBS_GroupInfo$HeteroticGroup2=="RHA"),
								"Within_Group", ifelse((IBS_GroupInfo$HeteroticGroup1=="HA" & IBS_GroupInfo$HeteroticGroup2=="RHA") |
								(IBS_GroupInfo$HeteroticGroup1=="RHA" & IBS_GroupInfo$HeteroticGroup2=="HA"),
								"Between_Group", "Other"))
aggregate(IBS_GroupInfo$IBS, by=list(IBS_GroupInfo$Cross), length)
aggregate(IBS_GroupInfo$IBS, by=list(IBS_GroupInfo$HeteroticGroup1, IBS_GroupInfo$HeteroticGroup2, IBS_GroupInfo$Cross), length) # check


Mod1 <- lm(IBS ~ Variant_type + Cross + Variant_type:Cross,
	data=IBS_GroupInfo[which(IBS_GroupInfo$Cross!="Other"),]) # significant Variant type : cross interaction
drop1(Mod1, test = "F") # F=33.887, p<0.0001

# make wide
IBS_GroupWide <- reshape(IBS_GroupInfo[which(IBS_GroupInfo$Cross!="Other"),],
	idvar=c("Genotype1", "HeteroticGroup1", "Genotype2", "HeteroticGroup2", "Cross"),
	timevar="Variant_type",
	direction="wide")

Mod2 <- lm(IBS.Deleterious ~ IBS.Synonymous + Cross +
	IBS.Synonymous:Cross,
	data = IBS_GroupWide)
drop1(Mod2, test = "F") # p <0.0001, F=148.83
library(car)
Anova(Mod2)
summary(Mod2)

## write.table(IBS_GroupWide, file = "/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/IBS_syndel.txt", sep = "\t", quote=FALSE, row.names=FALSE)
```

---


## Genomic Patterns

Continued with `dSNP_DerAncBINS.R` <- for binning across genome


### Recombination across genome

Using John Bowers genetic map information, I used the file I manipulated for the truth SNPs: `SNP_Genetic_Map_Unique.txt` for the cM distance of the markers. This file contains 6984 markers that mapped uniquely. 
- The first column is the linkage group (chromosome number), second column is the locus name, third column is the distance in cM.

I also used the vcf file that I created with SNPutils- `MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf` where I mapped the SNPs to the new genome build. (N=6523)
- Here the first column is the chromosome, 2nd is the position in bp, third is Locus name

##### Create Recombination Data Files
```bash
GeneticMap=/scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
RemappedVCF=/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf

cd /scratch/eld72413/SAM_seq/Recombination

grep -v "#" $RemappedVCF | awk '{print $1,"\t",$2,"\t",$3}' > SNParray_BPpositions.txt

awk '{print $1,"\t",$2,"\t",$3}' $GeneticMap > SNParray_cMpositions.txt


# I need to change the chromosome names to match:

awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' # need to remove the four contigs
awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' | cat -n > ChromNames.txt 


while read line; do
	OldChromName=$(echo $line | awk '{print $1}')
	NewChromName=$(echo $line | awk '{print $2}')
	echo "$OldChromName to $NewChromName"
	sed -i 's/'^"${OldChromName}"'\b/'"${NewChromName}"'/g' SNParray_cMpositions.txt
done < ChromNames.txt
```

To graph Recombination patterns across genome, I used R (see Genomic_Patterns/Recombination.R script).

## Binning SNP classes across Genome


### Number of codons per window
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh
module load BEDTools/2.30.0-GCC-8.3.0

GFF3=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.gff3
#OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

awk '{if ($3=="mRNA") {print $1,$4,$5}}' $GFF3 | wc -l # 72995

# use bedtools to find # bp (and # counts) for each mRNA feature in 10Mbp window- divide by 3 to get number of codons
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.genome.fasta.fai -w 10000000 \
| bedtools intersect -a - -b $GFF3 -wo \
| awk '{if ($6=="mRNA") {print $0}}' \
| bedtools groupby -g 1,2,3,6 -c 13 -o sum,count \
| awk '{print $1,$2,$3,$4, $5,$6, $5/3}' > ${OUT_DIR}/GenomicPatterns/Bins/Mbp_10/mRNA_10MbpNumCodonCounts.txt

```

### Number of dSNPs and sSNPs per window
```bash
SNP_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info_new.txt

# save a windows file to use for bedtools intersect commands

# make genome file with only 17 chromosomes (no scaffolds)
awk -v OFS='\t' {'print $1,$2'} "/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.genome.fasta.fai" | head -17 > "/scratch/eld72413/SunflowerGenome/GenomeFile.bed"

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 > $OUT_DIR/IntermediateFiles/Mbp10_windows.bed

# number of dSNPs and sSNPs
VariantClass="AllDel"
VariantClass="SynonymousNodups"

awk -v var="$VariantClass" 'OFS="\t" {if ($13==var){print $1, $2-1, $2}}' $SNP_INFO \
| bedtools intersect -a $OUT_DIR/IntermediateFiles/Mbp10_windows.bed -b - -c > ${OUT_DIR}/GenomicPatterns/Bins/Mbp_10/Num_${VariantClass}_10MbpCounts.bed

```

High frequency alleles
- high derived frequency (> 0.9) (N=2091 deleterious)
- high MAF (0.1-0.50) (N=12,557 deleterious)

- low MAF (<0.5) N= 67,005 (and not high derived)
	- 31,267 are only present in single individuals
```bash
SNP_INFO_2=${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt # with deleterious alleles where the deleterious is not derived relative to the outgroup are removed (N=6141)
awk '{if ($13=="AllDel") {print $0}}' $SNP_INFO_2 | wc -l #81,653 (including header)

# high derived allele frequency
awk '{if ($12 > 0.9 && $12!="NA") {print $0}}' $SNP_INFO_2 | wc -l # 1,667,416
awk '{if ($12 > 0.9 && $12!="NA" && $13=="AllDel") {print $0}}' $SNP_INFO_2 | wc -l # 2091 (an additional 2617 fit this category bc the ancestral allele is deleterious!)

awk 'NR>1 {if ($12 > 0.9 && $12!="NA") {print $0}}' $SNP_INFO_2 > ${OUT_DIR}/GenomicPatterns/Freq/HighDerivedFreqAlleles.txt

# number of dSNPs and sSNPs
VariantClass="AllDel"
VariantClass="SynonymousNodups"

awk -v var="$VariantClass" 'OFS="\t" {if ($13==var){print $1, $2-1, $2}}' ${OUT_DIR}/GenomicPatterns/Freq/HighDerivedFreqAlleles.txt \
| bedtools intersect -a $OUT_DIR/IntermediateFiles/Mbp10_windows.bed -b - -c > ${OUT_DIR}/GenomicPatterns/Freq/Num_HighDerFreq${VariantClass}_10MbpCounts.bed


# high MAF
awk '{if ($10 >= 0.1 && $10 <= 0.5) {print $0}}' $SNP_INFO_2 | wc -l # #6,542,031
awk '{if ($10 >= 0.1 && $10 <= 0.5 && $13=="AllDel" ) {print $0}}' $SNP_INFO_2 | wc -l # 12,557

awk '{if ($10 >= 0.1 && $10 <= 0.5) {print $0}}' $SNP_INFO_2 > ${OUT_DIR}/GenomicPatterns/Freq/HighMAFAlleles.txt

# number of dSNPs and sSNPs
VariantClass="AllDel"
VariantClass="SynonymousNodups"

awk -v var="$VariantClass" 'OFS="\t" {if ($13==var){print $1, $2-1, $2}}' ${OUT_DIR}/GenomicPatterns/Freq/HighMAFAlleles.txt \
| bedtools intersect -a $OUT_DIR/IntermediateFiles/Mbp10_windows.bed -b - -c > ${OUT_DIR}/GenomicPatterns/Freq/Num_HighMAF${VariantClass}_10MbpCounts.bed


# low MAF to make sure numbers add up
awk '{if (($10 < 0.1 && $12 < 0.9) || ($10 < 0.1 && $12 == "NA")) {print $0}}' $SNP_INFO_2 | awk '{if ($13=="AllDel") {print $0}}' | wc -l #67,005
awk '{if (($10 < 0.1 && $12 < 0.9) || ($10 < 0.1 && $12 == "NA")) {print $0}}' $SNP_INFO_2 > ${OUT_DIR}/GenomicPatterns/Freq/LowMAFAlleles.txt

VariantClass="AllDel"
VariantClass="SynonymousNodups"

awk -v var="$VariantClass" 'OFS="\t" {if ($13==var){print $1, $2-1, $2}}' ${OUT_DIR}/GenomicPatterns/Freq/LowMAFAlleles.txt \
| bedtools intersect -a $OUT_DIR/IntermediateFiles/Mbp10_windows.bed -b - -c > ${OUT_DIR}/GenomicPatterns/Freq/Num_LowMAF${VariantClass}_10MbpCounts.bed
```

Singletons or private doubletons with VCFtools
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l # tmux window in ss-sub4
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

vcftools --gzvcf ${VCF} --singletons --out ${OUT_DIR}/GenomicPatterns/Freq/All_snps
wc -l All_snps.singletons # 12,803,824
awk '{if ($3=="D") {print $0}}' All_snps.singletons | wc -l # 4,997,248 are private doubletons
awk '{if ($3=="S") {print $0}}' All_snps.singletons | wc -l # 7,806,575

awk 'OFS="\t" {print $1,$2}' ${OUT_DIR}/GenomicPatterns/Freq/All_snps.singletons > ${OUT_DIR}/GenomicPatterns/Freq/SingletonPositions.txt

grep -w -Ff ${OUT_DIR}/GenomicPatterns/Freq/SingletonPositions.txt ${SNP_INFO_2} | wc -l # 12,798,551
grep -w -Ff ${OUT_DIR}/GenomicPatterns/Freq/SingletonPositions.txt ${SNP_INFO_2} > ${OUT_DIR}/GenomicPatterns/Freq/SingleIndiv_Alleles.txt

# number of dSNPs and sSNPs
VariantClass="AllDel"
VariantClass="SynonymousNodups"

awk -v var="$VariantClass" 'OFS="\t" {if ($13==var){print $1, $2-1, $2}}' ${OUT_DIR}/GenomicPatterns/Freq/SingleIndiv_Alleles.txt \
| bedtools intersect -a $OUT_DIR/IntermediateFiles/Mbp10_windows.bed -b - -c > ${OUT_DIR}/GenomicPatterns/Freq/SingleIndiv_Alleles${VariantClass}_10MbpCounts.bed
# 265202 synonymous
# 31267 deleterious
```

See R script

how many dSNPs are singletons or private doubletons?
```bash
# how many dSNPs are singletons or private doubletons?
awk 'NR>1 {if ($3=="S") {print $1":"$2}}' ${OUT_DIR}/GenomicPatterns/Freq/All_snps.singletons > ${OUT_DIR}/IntermediateFiles/SingletonPositions.txt
awk 'NR>1 {if ($3=="D") {print $1":"$2}}' ${OUT_DIR}/GenomicPatterns/Freq/All_snps.singletons > ${OUT_DIR}/IntermediateFiles/PrivateDoubletonPositions.txt

grep -w -Ff ${OUT_DIR}/IntermediateFiles/SingletonPositions.txt ${OUT_DIR}/IntermediateFiles/AllDel_Positions.txt | wc -l # 22,824
grep -w -Ff ${OUT_DIR}/IntermediateFiles/PrivateDoubletonPositions.txt ${OUT_DIR}/IntermediateFiles/AllDel_Positions.txt | wc -l # 10,120
```
###### Answer: 
- 22,824 dSNPs are Singletons (26%)
- 10,120 Private doubletons (11.5%)
Total: 32,944 (37.5%)


### Haplotypes blocks across genome
Haplotype block coordinates obtained from Plink by Andries Temme (uploaded to cluster "blocks_with_delmut261.csv")

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l # tmux window in ss-sub4
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

# convert to 0-based numbering with Chromosome, start, end, name
# remove quotation marks from characters
awk 'BEGIN { FS=","; OFS="\t" } { gsub("\"", "") } { $1=$1 } 1' ${OUT_DIR}/Haplotypes/blocks_with_delmut.csv | awk 'NR>1 {{OFS="\t"}; print $3,$4-1,$5,$2}' > ${OUT_DIR}/Haplotypes/blocks_with_delmut.bed

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

vcftools --gzvcf ${VCF} --hapcount ${OUT_DIR}/Haplotypes/blocks_with_delmut.bed --out ${OUT_DIR}/Haplotypes/HaploBlocks
## the above command did not work

```

Find coordinates of the longest block on chromosome 10
```bash
SNP_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info_new.txt
awk '{if ($3-$2 > 50000000) {print $1":"$2+1"-"$3}}' ${OUT_DIR}/Haplotypes/blocks_with_delmut.bed
# 20 larger than 10M, 8 larger than 20M (one on block 10 > 70M)


awk '{if ($1=="Ha412HOChr10" && $2 > 40336916 && $2 <= 113115980) {print $0}}' $SNP_INFO > ${OUT_DIR}/Haplotypes/Chrom10/BranchingHaplotypeSNPs
# 631,584 ()
awk '{if ($1=="Ha412HOChr10" && $2 > 40336916 && $2 <= 113115980 && $13=="AllDel" && $12!="NA") {print $0}}' $SNP_INFO | wc -l # 801 (537 polarized)
awk '{if ($1=="Ha412HOChr10" && $2 > 40336916 && $2 <= 113115980 && $13=="SynonymousNodups" && $12!="NA") {print $0}}' $SNP_INFO | wc -l # 5855 (3968 polarized)
# above numbers are without the snps removed where ancestral allele = deleterious allele

#bcftools view -Ou -r - ${VCF} |\
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${OUT_DIR}/Haplotypes/Chrom10/
```
See code in Joint_SFS.R for getting info for HA/RHA derived freqs for this region

```bash
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${OUT_DIR}/Haplotypes/Chrom10/Haplotype_Chr10_Info.txt" \
"1.0" \
"0.05" \
"Derived_Freq_HA" \
"${OUT_DIR}/Haplotypes/Chrom10/DerivedFreq_Chr10_HA_Bins.txt"

Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${OUT_DIR}/Haplotypes/Chrom10/Haplotype_Chr10_Info.txt" \
"1.0" \
"0.05" \
"Derived_Freq_RHA" \
"${OUT_DIR}/Haplotypes/Chrom10/DerivedFreq_Chr10_RHA_Bins.txt"
```



############# May redo the code below:



Combine Frequency info with Fst and pruned list for plotting:
```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
library(ggplot2)

Fst <- Combine_Chromosomes("/scratch/eld72413/SAM_seq/Fst/HA_RHA", ".weir.fst", "HA_RHA_") # N=36,389,432
threshold <- quantile(Fst$WEIR_AND_COCKERHAM_FST, c(0.995), na.rm = T) # 0.376806

PrunedList <- read.table("/scratch/eld72413/SAM_seq/PCA/PrunedLists/Plink_0.9.prune.in", sep=":",
	col.names=c("CHROM", "POS")) # N=5,512,923

#Fst_pruned <- merge(PrunedList, Fst, by=c("CHROM", "POS")) # N=5,413,348 (pruned list minus scaffold regions)
freqs <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_Freqs.txt", 
	header=T, sep="\t")
# Derived allele frequencies
#freqs <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_DerivedFreqs.txt", 
#	header=T, sep="\t")

FreqFst_Prune <- merge(freqs, 
	merge(PrunedList, Fst, by=c("CHROM", "POS")), 
	by.x=c("Chromosome", "Position"), by.y=c("CHROM", "POS"))
#write.table(FreqFst_Prune, "/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/FreqPrunedFst.txt", sep = "\t", quote=FALSE, row.names=FALSE)

plot <- FreqScatterplot(FreqFst_Prune,
                        threshold, "AllDel", "Synonymous", "HA_AltFreq", "RHA_AltFreq")

pdf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_FreqScatterplot.pdf")
print(plot)
dev.off()

# distributed non-randomly?
Del_outlier <- FreqFst_Prune[which(FreqFst_Prune$WEIR_AND_COCKERHAM_FST > threshold & FreqFst_Prune$Variant_type=="AllDel"),] # N=62
aggregate(Del_outlier$Position, by=list(Del_outlier$Chromosome), length) # only found on chromosomes 8 (N=4), 9 (N=2), 10 (N=44), 12 (N=12)

FreqFst_Prune_sub <- subset(FreqFst_Prune, Variant_type=="AllDel" | Variant_type=="Synonymous")
plot <- FreqScatterplot(FreqFst_Prune_sub,
                        threshold, "AllDel", "Synonymous", "HA_Derived_Freq", "RHA_Derived_Freq")
pdf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_FreqDerivedScatterplot.pdf")
print(plot)
dev.off()
```






----- 



### Number of different SNP classes per window
Using R object saved from SFS_Info.R script
```R
# tmux new -s SNP_bins
# srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# module load R/4.0.0-foss-2019b
load("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/FrequencyInfo.RData") # object is FrequencyInfo

# function to bin number of each SNP class across genome
SNP_bins <- function(dataset, BinSize_bp, SNP_class_Name) {
	dataset <- dataset[order(dataset$Position),]
	Num_windows <- ceiling(max(dataset$Position)/BinSize_bp)
	dataset$bin <- cut(dataset$Position,seq(0,Num_windows*BinSize_bp,BinSize_bp))
	counts <- aggregate(dataset$Position,
                      by=list(dataset$bin), length, drop=FALSE)
	colnames(counts) <- c("Bin", paste0("Number_", SNP_class_Name))
	counts$BinList <- strsplit(as.character(counts$Bin), ",")
	counts$StartPos <- as.numeric(gsub('[(]', '', sapply(counts$BinList, "[", 1)))
	counts$EndPos <- as.numeric(gsub('[]]', '', sapply(counts$BinList, "[", 2)))
	return(counts[,c(4,5,2,1)])
}

# function to apply SNP_bins across all chromosomes
apply_SNP_bins <- function(dataset, BinSize_bp, SNP_class_Name) {
	dataset$Chromosome <- factor(dataset$Chromosome)
	dataset_chrom <- split(dataset, dataset$Chromosome)
	binCounts <- lapply(dataset_chrom, function(x) {
		SNP_bins(x, BinSize_bp, SNP_class_Name)
	})
	binCounts <- lapply(names(binCounts), function(x) {
		binCounts[[x]]["Chromosome"] <- x; return(binCounts[[x]])
	})
	binCounts_df <- do.call("rbind", binCounts)
	return(binCounts_df)
}

### first subset to only derived SNPs
FrequencyInfo_Derived <- lapply(FrequencyInfo, function(x) {
	subset(x, !is.na(x$Derived_Freq))
	})
names(FrequencyInfo_Derived)
names(FrequencyInfo_Derived)[2:5] # not using the deleterious SNPs from this list

Derived_Bins <- lapply(names(FrequencyInfo_Derived)[2:5], function(x) {
	apply_SNP_bins(FrequencyInfo_Derived[[x]], 10000000, x)
	})
names(Derived_Bins) <- names(FrequencyInfo_Derived)[2:5]

## for dSNPs:
setwd("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/AlleleFreqs")
# use dSNPNums function above
dSNP_summary <- dSNPNums("SAM_Altdeleterious_AC_AN.txt", "SAM_Refdeleterious_AC_AN.txt", "/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt")
dSNP_derived <- subset(dSNP_summary, Cat=="Derived_dSNP")

# (take out scaffold seq)
dSNP_Derived_Bins <- apply_SNP_bins(dSNP_derived[-grep("Ha412HOChr00c", dSNP_derived$Chromosome),], 
	10000000, "dSNP")

### merge all classes
All_info1 <- Reduce(function(x, y) merge(x, y, 
	by=c("StartPos", "EndPos", "Chromosome", "Bin"), all=TRUE), Derived_Bins)

### merge (without scaffold, N=325) with dSNP data
All_info <- merge(dSNP_Derived_Bins, 
	All_info1[-grep("Ha412HOChr00c", All_info1$Chromosome),], 
	by=c("StartPos", "EndPos", "Chromosome", "Bin"), all=TRUE)

All_info$StartPos <- as.integer(All_info$StartPos)
All_info$EndPos <- as.integer(All_info$EndPos)

### combine with number of codon info
NumCodons <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/mRNA_bpNumCodonCounts.txt",
	header=FALSE)
colnames(NumCodons) <- c("Chromosome", "StartPos", "EndPos", "feature", "Total_bases", "Total_features", "Num_codons")

All_info2 <- merge(All_info, NumCodons[,c(1:3,7)], by=c("StartPos", "Chromosome"))
All_info2 <- All_info2[order(All_info2$Chromosome, All_info2$StartPos),]

write.table(All_info2, file="/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/Derived_VariantNums.txt",
	row.names=FALSE, quote=FALSE, sep="\t")

```



