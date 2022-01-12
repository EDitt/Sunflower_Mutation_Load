# Why is there an excess of low frequency variants?

## 1.) Issues with SNP filtering?

### SFS of the raw SNP set compared with the filtered set

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

module load BCFtools/1.13-GCC-8.3.0
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz

raw_vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Sunflower_SAM_SNP_Calling_raw_variants.vcf
# this includes indels...?

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' ${raw_vcf} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Raw_all_alleleFreqInfo.txt

wc -l Raw_all_alleleFreqInfo.txt # 101,509,931
awk '{if ($5==1) {print $0}}' Raw_all_alleleFreqInfo.txt | wc -l # 16,661,402
awk '{if ($7<0.05) {print $0}}' Raw_all_alleleFreqInfo.txt | wc -l # 73,819,616


```

Bin data & plot histogram
```R
# module load R/4.0.0-foss-2019b
# R
library(tidyr)

### Function to make histogram breaks for one dataset
query_to_hist <- function(bcftools_file, hist_breaks) {
	data <- read.table(bcftools_file, sep = "\t", header=FALSE)
	colnames(data) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles", "Alt_Freq")
	data <- separate_rows(data, Alt_Freq, sep = ",")
	data$Alt_Freq <- as.numeric(data$Alt_Freq)
	data$MAF <- ifelse(data$Alt_Freq < 0.5, data$Alt_Freq, 1-data$Alt_Freq)
	y <- hist(data$MAF, plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(hist_breaks[-1], prop)
	colnames(hist_df) <- c("Break", "Proportion")
	return(hist_df)
}

MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

SNP_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt", MAF_breaks)
write.table(SNP_freq_data, "FilteredSNP_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

raw_hist_df <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Raw_all_alleleFreqInfo.txt", MAF_breaks)
write.table(raw_hist_df, "Raw_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

# saved to local computer
library(ggplot2)

raw_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/Raw_Variant_Proportions.txt",
                                                  sep="\t", header = T)
raw_sfs$Type <- "raw"

filtered_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/FilteredSNP_Variant_Proportions.txt",
                                                  sep="\t", header = T)

filtered_sfs$Type <- "filtered"

Both_sfs <- rbind.data.frame(raw_sfs, filtered_sfs)

ggplot(Both_sfs, aes(x=Break, y=Proportion, fill=Type))+ 
  geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)"))
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/RawFiltered.png")
```

Raw dataset also exhibits an excess of rare alleles

## 2.) Related to the pooling scheme used (each line was sequenced from pooled DNA)?

- tried removing singletons and variants that appeared only twice which did not totally remove this issue
- also took out sites with allele count of 1-2:
```bash
# take out sites with allele count of 1-2
awk '{if ($5>2) {print $0}}' $Pos_info | wc -l # 21,699,399
awk '{if ($5>2 && $7<0.05) {print $0}}' $Pos_info | wc -l  # 13,237,478
awk '{if ($5>2) {print $0}}' $Pos_info > All_alleleFreqInfo_noOneorTwo.txt

Rscript "SFS_Info.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo_noOneorTwo.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/MAF_Bins_noOneOrTwo.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/DerivedFreq_Bins_noOneOrTwo.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/FrequencyInfo_noOneOrTwo.RData"
```
Still exhibits an excess of rare alleles

### Separate SFS for homozygous and heterozygous genotypes

Use bcftools to output stats from sites that have no homozygous calls
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.13-GCC-8.3.0
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz

# SFS for only homozygous genotypes
bcftools filter -i 'COUNT(GT="het")=0' --threads 4 $vcf -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Homozygous_FreqInfo.txt
wc -l Homozygous_FreqInfo.txt # 7,308,060

# SFS for sites with at least 1 heterozygote genotype
bcftools filter -i 'COUNT(GT="het")>0' --threads 4 $vcf -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Heterozygous_FreqInfo.txt
wc -l Heterozygous_FreqInfo.txt # 29,821,855

```
```R
# module load R/4.0.0-foss-2019b
# R
library(tidyr)

MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

No_Het_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Homozygous_FreqInfo.txt", MAF_breaks)
write.table(No_Het_freq_data, "NoHet_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

Het_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Heterozygous_FreqInfo.txt", MAF_breaks)
write.table(Het_freq_data, "Het_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

## saved to local computer
setwd("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS")

library(ggplot2)

combine2_sfs <- function(file1, file2, name1, name2) {
  data1 <- read.table(file1, sep="\t", header = T)
  data1$Type <- name1
  data2 <- read.table(file2, sep="\t", header = T)
  data2$Type <- name2
  Combined <- rbind.data.frame(data1, data2)
  return(Combined)
}

graph_folded_sfs <- function(dataframe, x, y, category) {
  p<- ggplot(dataframe, aes(x=x, y=y, fill=category))+ 
    geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
    ylab("Proportion") +
    xlab("Frequency") +
    scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                       labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                  "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)"))
  return(p)
}

HomHet <- combine2_sfs("NoHet_Variant_Proportions.txt", "Het_Variant_Proportions.txt", "homozygous_only", "contains_heterozygotes")

graph_folded_sfs(HomHet, HomHet$Break, HomHet$Proportion, HomHet$Type)
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/HomHet.png")
```

## 3.) Selective sweep or recently expanding population

### Separate SFS for different classes of germplasm
```bash
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv

# non elite varieties- exclude any HA/RHA
genotypes=$(awk -v var="HA" -F',' '{if ($8!~var && $9!="HA412" && $9!="NA" && NR!=1) {print $9}}' $SAM_info | paste -sd,)

bcftools view -Ou --samples ${genotypes} ${vcf} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Non_eliteFreqInfo.txt
# need to remove zeros

awk '{if ($5>0) {print $0}}' Non_eliteFreqInfo.txt > Non_eliteFreqInfoNoZero.txt
wc -l Non_eliteFreqInfoNoZero.txt # 25,031,051
#### will need to recalculate the frequency

# HA only
awk -F',' '{if ($8~/^HA/ && $9!="HA412" && $9!="NA" && NR!=1) {print $0}}' $SAM_info

HA_genotypes=$(awk -F',' '{if ($8~/^HA/ && $9!="HA412" && $9!="NA" && NR!=1) {print $9}}' $SAM_info | paste -sd,)

bcftools view -Ou --samples ${HA_genotypes} ${vcf} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' | \
awk '{if ($5>0) {print $0}}' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/HA_FreqInfo.txt
wc -l HA_FreqInfo # 23,861,044

# RHA only
awk -F',' '{if ($8~/^RHA/ && $9!="HA412" && $9!="NA" && NR!=1) {print $0}}' $SAM_info

RHA_genotypes=$(awk -F',' '{if ($8~/^RHA/ && $9!="HA412" && $9!="NA" && NR!=1) {print $9}}' $SAM_info | paste -sd,)

bcftools view -Ou --samples ${RHA_genotypes} ${vcf} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' | \
awk '{if ($5>0) {print $0}}' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/RHA_FreqInfo.txt
wc -l RHA_FreqInfo.txt # 21,065,435
```


```R
# module load R/4.0.0-foss-2019b
# R
library(tidyr)

MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

# modify function to recalculate alternate allele frequency based on subset
### Function to make histogram breaks for one dataset
query_to_hist <- function(bcftools_file, hist_breaks) {
	data <- read.table(bcftools_file, sep = "\t", header=FALSE)
	colnames(data) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles", "Alt_Freq")
	data$Alt_Freq_Recalc <- data$Num_Alt_alleles / data$Num_alleles
	data$MAF <- ifelse(data$Alt_Freq_Recalc < 0.5, data$Alt_Freq_Recalc, 1-data$Alt_Freq_Recalc)
	y <- hist(data$MAF, plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(hist_breaks[-1], prop)
	colnames(hist_df) <- c("Break", "Proportion")
	return(hist_df)
}

HA_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/HA_FreqInfo.txt", MAF_breaks)
write.table(HA_freq_data, "HA_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

Non_elite_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Non_eliteFreqInfoNoZero.txt", MAF_breaks)
write.table(Non_elite_freq_data, "NonElite_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

RHA_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/RHA_FreqInfo.txt", MAF_breaks)
write.table(RHA_freq_data, "RHA_Variant_Proportions.txt", sep = "\t", quote=FALSE, row.names=FALSE)

### also remove singletons
awk '{if ($5>1) {print $0}}' HA_FreqInfo.txt > HA_FreqInfo_noSingleton.txt # 19,329,156
awk '{if ($5>1) {print $0}}' RHA_FreqInfo.txt > RHA_FreqInfo_noSingleton.txt # 17,259,309
awk '{if ($5>1) {print $0}}' Non_eliteFreqInfoNoZero.txt > Non_eliteFreqInfo_noSingleton.txt # 17,366,354

HA_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/HA_FreqInfo_noSingleton.txt", MAF_breaks)
write.table(HA_freq_data, "HA_Variant_Proportions_noSingle.txt", sep = "\t", quote=FALSE, row.names=FALSE)

Non_elite_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/Non_eliteFreqInfo_noSingleton.txt", MAF_breaks)
write.table(Non_elite_freq_data, "NonElite_Variant_Proportions_noSingle.txt", sep = "\t", quote=FALSE, row.names=FALSE)

RHA_freq_data <- query_to_hist("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SFS_bins_Troubleshoot/RHA_FreqInfo_noSingleton.txt", MAF_breaks)
write.table(RHA_freq_data, "RHA_Variant_Proportions_noSingle.txt", sep = "\t", quote=FALSE, row.names=FALSE)

### saved to local computer
combine3_sfs <- function(file1, file2, file3, name1, name2, name3) {
  data1 <- read.table(file1, sep="\t", header = T)
  data1$Type <- name1
  data2 <- read.table(file2, sep="\t", header = T)
  data2$Type <- name2
  data3 <- read.table(file3, sep="\t", header = T)
  data3$Type <- name3
  Combined <- rbind.data.frame(data1, rbind.data.frame(data2, data3))
  return(Combined)
}

GermplasmGroup <- combine3_sfs("NonElite_Variant_Proportions.txt", "HA_Variant_Proportions.txt", "RHA_Variant_Proportions.txt",
                               "Non-elite", "HA", "RHA")
graph_folded_sfs(GermplasmGroup, GermplasmGroup$Break, GermplasmGroup$Proportion, GermplasmGroup$Type)
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/GermplasmGroup.png")

GermplasmGroup2 <- combine3_sfs("NonElite_Variant_Proportions_noSingle.txt", "HA_Variant_Proportions_noSingle.txt", 
                                "RHA_Variant_Proportions_noSingle.txt",
                               "Non-elite", "HA", "RHA")
graph_folded_sfs(GermplasmGroup2, GermplasmGroup2$Break, GermplasmGroup2$Proportion, GermplasmGroup2$Type)
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/GermplasmGroup_noSingletons.png")
```

In the "non-elite/no singletons", the SFS looks much closer to what would be expected in an ideal population