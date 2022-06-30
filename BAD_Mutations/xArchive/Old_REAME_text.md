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




scratch in chunk below (for now)
```bash

###

### updating this
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_RefDerivedDeleterious.vcf |grep "PSC" > Stats/RefDerivedDeleterious_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedDeleterious.vcf | grep "PSC" > Stats/AltDerivedDeleterious_perSampleCounts.txt

bcftools stats -s - SAM_RefDerivedTolerated.vcf | grep "PSC" > Stats/RefDerivedTolerated_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedTolerated.vcf | grep "PSC" > Stats/AltDerivedTolerated_perSampleCounts.txt

bcftools stats -s - SAM_RefDerivedSynonymous.vcf | grep "PSC" > Stats/RefDerivedSynonymous_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedSynonymous.vcf | grep "PSC" > Stats/AltDerivedSynonymous_perSampleCounts.txt

### need total number of genotypes:
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
bcftools stats -s - $vcf | grep "PSC" > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllVariantStats.txt

### need total number of derived genotypes:

vcf=
bcftools stats -s - $vcf | grep "PSC" > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllDerived_VariantStats.txt

module load R/4.0.0-foss-2019b
Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm/Derived_Variant_Numbers.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/Stats" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllVariantStats.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/AlleleNums_Geno.txt"


```

## Code I am no longer using


### Wrangling SNP data to remove variants for which the deleterious allele was not derived relative to outgroup, H. debilis


#### Deleterious SNPs
Compiled Report was used to identify variant classes (see `2.BAD_Mutations/Post_processing.md` and positions files created in `Variant_class_numbers.md`)

First, subset vcf file to make separate vcfs of deleterious (for both reference and alternate alleles)
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

# deleterious in reference
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Reference_DelPositions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Refdeleterious' Subset_vcf.sh # Submitted batch job 4232837

grep -v "#" SAM_Refdeleterious.vcf | wc -l # 11796

# deleterious in alternate
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Alternate_DelPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Altdeleterious' Subset_vcf.sh # Submitted batch job 4232853

grep -v "#" SAM_Altdeleterious.vcf | wc -l # 76016

# use bcftools to get stats for each allele
module load BCFtools/1.13-GCC-8.3.0

# need to bgzip vcf files:
module load HTSlib/1.10.2-GCC-8.3.0
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz

# derived deleterious in alternate
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz

# reference deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Refdeleterious.vcf.gz > AlleleFreqs/SAM_Refdeleterious_AC_AN.txt
# alternate deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Altdeleterious.vcf.gz > AlleleFreqs/SAM_Altdeleterious_AC_AN.txt

# N_ALT : number of alternate alleles
# AC : count of alternate alleles
# AN : number of alleles in called genotypes

### need to calculate frequency of deleterious allele (for deleterious reference it's AN - AC/AN; for deleterious alternate it's AC/AN)
```

Use R to wrangle
```R
#module load R/4.0.0-foss-2019b
#R

dSNPNums <- function (Alt_file, Ref_file, AncestralState_file) {
  Alt_dSNP <- read.table(Alt_file, sep = "\t", header=FALSE)
  colnames(Alt_dSNP) <- c("Chromosome", "Position", "Ref_allele", "Del_allele",
    "Num_Alt_alleles", "Num_alleles")
  Alt_dSNP$dSNP_freq <- Alt_dSNP$Num_Alt_alleles / Alt_dSNP$Num_alleles
  Ref_dSNP <- read.table(Ref_file, sep = "\t", header=FALSE)
  colnames(Ref_dSNP) <- c("Chromosome", "Position", "Del_allele", "Alt_allele",
    "Num_Alt_alleles", "Num_alleles")
  Ref_dSNP$dSNP_freq <- (Ref_dSNP$Num_alleles - Ref_dSNP$Num_Alt_alleles) / Ref_dSNP$Num_alleles
  All_dSNP_freq <- rbind(Ref_dSNP[,c(1:3,7)], Alt_dSNP[,c(1,2,4,7)])
  All_dSNP_freq$Chromosome <- as.factor(All_dSNP_freq$Chromosome)
  All_dSNP_freq <- with(All_dSNP_freq, All_dSNP_freq[order(Chromosome, Position),])
  Polarized <- read.table(AncestralState_file)
  colnames(Polarized) <- c("Chromosome", "Position", "AncestralAllele")
  dSNP_summary_AA <- merge(All_dSNP_freq, Polarized, by=c("Chromosome", "Position"), all.x=TRUE)
  dSNP_summary_AA$Cat <- ifelse(dSNP_summary_AA$AncestralAllele==dSNP_summary_AA$Del_allele, "Ancestral_dSNP",
  ifelse(dSNP_summary_AA$AncestralAllele!=dSNP_summary_AA$Del_allele,"Derived_dSNP", NA))
    return(dSNP_summary_AA)
}

dSNP_summary <- dSNPNums("SAM_Altdeleterious_AC_AN.txt", "SAM_Refdeleterious_AC_AN.txt", "/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt")

#duplicates?
length(unique(paste0(dSNP_summary$Chromosome,"_", dSNP_summary$Position))) # 87794 (out of 87812) 18 duplicated
# how many of these are polarized?
aggregate(dSNP_summary$Position, by=list(dSNP_summary$Cat), length)
# ancestral dSNP= 6141; derived dSNP= 56360 (total 62,501) <- still have 14 duplicate positions

derived_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
Derived_dSNP_histogram <- hist(dSNP_summary[which(dSNP_summary$Cat=="Derived_dSNP"), "dSNP_freq"],
  plot = FALSE, breaks=derived_breaks)
sum(Derived_dSNP_histogram$counts) # 56360

# use "Hist_bins" function from SFS_Info.R
Derived_dSNP <- Hist_bins(dSNP_summary[which(dSNP_summary$Cat=="Derived_dSNP"),],
  derived_breaks, "dSNP_freq", "dSNP")
colnames(Derived_dSNP)[1] <- "breaks"
write.table(Derived_dSNP, file="/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Derived_dSNP_freqbins.txt", sep = "\t", quote=FALSE, row.names=FALSE)

```


Proportions in each frequency class:
```bash
MAF=/scratch/eld72413/SAM_seq/dSNP_results/MAF_Bins.txt
Derived_freq=/scratch/eld72413/SAM_seq/dSNP_results/DerivedFreq_Bins.txt

awk '{if ($4=="AllDel") {print $0}}' $MAF
awk '{if ($4=="NonCoding") {print $0}}' $MAF
awk '{if ($4=="StopLostGained") {print $0}}' $MAF
awk '{if ($4=="Synonymous") {print $0}}' $MAF
awk '{if ($4=="Tolerated") {print $0}}' $MAF

awk '{if ($3=="NonCoding") {print $0}}' $Derived_freq
awk '{if ($3=="StopLostGained") {print $0}}' $Derived_freq
awk '{if ($3=="Synonymous") {print $0}}' $Derived_freq
awk '{if ($3=="Tolerated") {print $0}}' $Derived_freq

awk '{if ($4=="AllDel") {print $0}}' $Derived_freq
```


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