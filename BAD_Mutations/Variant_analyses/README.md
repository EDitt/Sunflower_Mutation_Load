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

### Number of derived SNPs for each genotype

First get number of called genotypes for polarized positions ** note: this failed
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

awk 'BEGIN{FS=OFS="\t"}; NR>1 {if ($11!="NA") {print $1,$2}}' ${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt > /scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/PolarizedSNP_Posititions.txt

bcftools view -Oz ${VCF} -R ${OUT_DIR}/IntermediateFiles/PolarizedSNP_Posititions.txt | \
bcftools stats -s - | \
grep "PSC" > ${OUT_DIR}/GenotypeInfo/AllDerived_VariantStats.txt
```

Use scripts to output a table with number of derived variants across all variant classes
```bash
awk 'NR>1 {print $13}' ${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt | sort -u #  annotation classes

# script to return nRefHom, nNonRefHom, nHets, etc. for all genotypes (for all variant classes)
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

See graph of derived dSNPs/sSNPs at: Plots/Germplasm_boxplot.R


### Frequency of derived SNPs for different germplasm groups


```bash
sbatch --export=SAM_INFO='/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/GroupFreqs',\
snps_remove='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/ToRemove_NotDerived.txt' /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Group_Freqs.sh # Submitted batch job 12461287

# this returns tables for all groups of variants with the 
```
test: Submitted batch job 12460630


############# May redo the code below:

#### Private/Shared
I will look at the proportion of shared versus private SNPs for deleterious and synonymous across the different frequency classes as well as plot the frequencies of all SNPs for both groups
First, set up a table with derived frequency information for SNPs among the two different groups as well as variant annotation
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

mkdir -p ${OUT_DIR}/HeteroticGroups

# stats for HA and RHA
group=HA
genotypes=$(awk -v var="$group" -F',' '{if ($11==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO | paste -sd,) # N=131
bcftools view -Ou --samples ${genotypes} ${VCF} |\
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${OUT_DIR}/HeteroticGroups/HA_SNP_info.txt

group=RHA
genotypes=$(awk -v var="$group" -F',' '{if ($11==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO | paste -sd,) # N=104
bcftools view -Ou --samples ${genotypes} ${VCF} |\
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${OUT_DIR}/HeteroticGroups/RHA_SNP_info.txt

# reduce SNP table to just deleterious and synonymous to save memory in R
awk '{print $13}' $SNP_INFO | sort -u
awk '{if ($13=="AllDel" || $13=="Synonymous" || $13=="Variant_type") {print $0}}' $SNP_INFO > ${OUT_DIR}/HeteroticGroups/SNP_Info_Reduced.txt #### <- reduce the one where I took out the derived dSNP issue? (dSNP but not derived relative to debilis)

# combine information into spreadsheet
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SNP_Freq_Groups.R" \
"${OUT_DIR}/HeteroticGroups/HA_SNP_info.txt" \
"HA" \
"${OUT_DIR}/HeteroticGroups/RHA_SNP_info.txt" \
"RHA" \
"${OUT_DIR}/HeteroticGroups/SNP_Info_Reduced.txt" \
"${OUT_DIR}/HeteroticGroups/HA_RHA_Freqs.txt"

```

Use R to bin and plot results
```R
setwd("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups")
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
Heterotic_FreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_Freqs.txt", 
	"HA", "RHA", "_Num_Alt_alleles", "_Num_Ref_alleles", 0.5, 0.1, "MAF")

library(ggplot2)
library(ggpubr)
plot <- function(df) {
  ggplot(data=df, aes(x=Bin, y=PropPrivate, fill=Annotation)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("#FC4E07", "#00AFBB")) +
  theme_minimal() +
  ylab("") + xlab("")
}
p1 <- plot(Heterotic_FreqBins[which(Heterotic_FreqBins$Bin=="(0,0.1]"),])
p2 <- plot(Heterotic_FreqBins[which(Heterotic_FreqBins$Bin!="(0,0.1]"),])
p3 <- annotate_figure(ggarrange(p1,p2, common.legend = TRUE, widths = c(1,3)), 
                bottom = text_grob("Frequency Bin"),
                left=text_grob("Proportion Private", rot=90))
pdf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_Private.pdf")
print(p3)
dev.off()
```

#### Frequency differentiation
Scatterplot of HA/RHA frequencies

Prune VCF with ancestral state calls to obtain derived allele frequency
```bash
#INPUT_VCF
#OUT_PREFIX
```

Joint SFS
```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
snp_freqs <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_Freqs.txt", header=T, sep="\t") # N=914,172
snp_freqs_derived <- subset(snp_freqs, !is.na(Derived_Freq)) # remove snps that are not polarized # N=649,187

# to remove from list (deleterious alleles not derived with respect to H. debilis)
to_remove <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/ToRemove_NotDerived.txt", sep = "\t",
	col.names=c("Chromosome", "Position")) # 6144

snp_freqs_new <- snp_freqs_derived[which(!paste0(snp_freqs_derived$Chromosome, "_", snp_freqs_derived$Position) %in%
	paste0(to_remove$Chromosome, "_", to_remove$Position)),] # 643046

snp_freqs_type <- split(snp_freqs_new, snp_freqs_new$Variant_type)

SFS_plots1 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot1(x, 0.05, "HA_Derived_Freq", "RHA_Derived_Freq", "HA Derived Freq", "RHA Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_JointSFS.pdf")
print(ggarrange(SFS_plots1$Synonymous, SFS_plots1$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

SFS_plots2 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot2(x, 100, "HA_Derived_Freq", "RHA_Derived_Freq", "HA Derived Freq", "RHA Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/HA_RHA_JointSFS2.pdf")
print(ggarrange(SFS_plots2$Synonymous, SFS_plots2$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()
```

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

#### IBS patterns
Scatterplot showing relationship between pairwise IBS for synonymous versus deleterious variants.
Difference between inter and intra heterotic group crosses?

```bash
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

### using the "PCA filter" Plink file which was converted from a VCF with singletons removed (see `${REPO_DIR}/Variant_analyses/PCA.md`)

# synonymous positions
type="Synonymous"
awk -v var="$type" 'BEGIN{OFS=":"}; {if ($13==var) {print $1,$2}}' $SNP_INFO > /scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Synonymous_Pos4Plink.txt

sbatch --export=Input_prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
SNP_List='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Synonymous_Pos4Plink.txt',\
Window_Size='1',\
Step_Size='1',\
Rsquared='0.9',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups',\
output_prefix='Synonymous' ${REPO_DIR}/BAD_Mutations/Variant_analyses/IBS_SNP_Subset.sh # Submitted batch job 8178681 (formerly 8177892)

# deleterious positions
type="AllDel"
awk -v var="$type" 'BEGIN{OFS=":"}; {if ($13==var) {print $1,$2}}' $SNP_INFO > /scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Deleterious_Pos4Plink.txt

sbatch --export=Input_prefix='/scratch/eld72413/SAM_seq/PCA/Sunflower_SAM_SNP_Calling_PCAfilter',\
SNP_List='/scratch/eld72413/SAM_seq/dSNP_results/IntermediateFiles/Deleterious_Pos4Plink.txt',\
Window_Size='1',\
Step_Size='1',\
Rsquared='0.9',\
outputdir='/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups',\
output_prefix='Deleterious' ${REPO_DIR}/BAD_Mutations/Variant_analyses/IBS_SNP_Subset.sh # Submitted batch job 8178683 (formerly 8178631)
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
aggregate(IBS_GroupInfo$IBS, by=list(IBS_GroupInfo$HeteroticGroup1, IBS_GroupInfo$HeteroticGroup2), length) # check


Mod1 <- lm(IBS ~ Variant_type + Cross + Variant_type:Cross,
	data=IBS_GroupInfo[which(IBS_GroupInfo$Cross!="Other"),]) # significant Variant type : cross interaction
drop1(Mod1, test = "F") # F=33.887, p<0.0001

# make wide
IBS_GroupWide <- reshape(IBS_GroupInfo[which(IBS_GroupInfo$Cross!="Other"),],
	idvar=c("Genotype1", "HeteroticGroup1", "Genotype2", "HeteroticGroup2", "Cross"),
	timevar="Variant_type",
	direction="wide")

Mod2 <- lm(IBS.Deleterious ~ IBS.Synonymous + Cross,
	data = IBS_GroupWide)
# slope = 0.7786, within-group cross estimate is positive

## write.table(IBS_GroupWide, file = "/scratch/eld72413/SAM_seq/dSNP_results/HeteroticGroups/IBS_syndel.txt", sep = "\t", quote=FALSE, row.names=FALSE)
```




----- 

Continued with `dSNP_DerAncBINS.R` <- for binning across genome

## Genomic Patterns

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
GFF3=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.gff3
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

awk '{if ($3=="mRNA") {print $1,$4,$5}}' $GFF3 | wc -l # 72995

# use bedtools to find # bp (and # counts) for each mRNA feature in 10Mbp window- divide by 3 to get number of codons
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b $GFF3 -wo \
| awk '{if ($6=="mRNA") {print $0}}' \
| bedtools groupby -g 1,2,3,6 -c 13 -o sum,count \
| awk '{print $1,$2,$3,$4, $5,$6, $5/3}' > ${OUTPUTDIR}/mRNA_bpNumCodonCounts.txt

```

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


#### Older way I binned across genome:
Make windows to bin count of each class of SNPs-
dSNPs, sSNPs, Tolerated SNPs
```bash
# make genome file
cd /scratch/eld72413/SunflowerGenome
awk -v OFS='\t' {'print $1,$2'} "Ha412HOv2.0-20181130.fasta.fai" | head -17 > "GenomeFile.txt"

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
mkdir GenomicBins

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BEDTools/2.30.0-GCC-8.3.0

vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

### dSNPs per window
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Refdeleterious.vcf -c > ${out_dir}/RefDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Altdeleterious.vcf -c > ${out_dir}/AltDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_tolerated.vcf.gz -c > ${out_dir}/Tolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_synonymous.vcf.gz -c > ${out_dir}/Synonymous_10MbCounts.txt

### do the same for the derived alleles
vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/GenomicBins

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedDeleterious.vcf -c > ${out_dir}/RefDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedTolerated.vcf -c > ${out_dir}/RefDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedSynonymous.vcf -c > ${out_dir}/RefDerivedSynonymous_10MbCounts.txt


##
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedDeleterious.vcf -c > ${out_dir}/AltDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedTolerated.vcf -c > ${out_dir}/AltDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedSynonymous.vcf -c > ${out_dir}/AltDerivedSynonymous_10MbCounts.txt

```

