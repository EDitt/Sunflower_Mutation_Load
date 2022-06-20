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