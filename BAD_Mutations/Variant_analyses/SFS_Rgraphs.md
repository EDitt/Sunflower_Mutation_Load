### Folded site frequency spectrum

```bash
cd /scratch/eld72413/NSFproj/VEP
module load R/3.6.2-foss-2019b
R
```

Function to convert format and identify lowest frequency allele

```R
library(tidyr)

vcftoolsFreqFormat <- function(filename) {
	df <- read.table(filename, sep = "\t", header=TRUE, row.names=NULL)
	colnames(df) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele1Freq", "Allele2Freq")
  	df1 <- df %>% separate(Allele1Freq, c("Base1", "Freq1"), sep=":")
  	df2 <- df1 %>% separate(Allele2Freq, c("Base2", "Freq2"), sep=":")
  	df2$freq <- as.numeric(ifelse(df2$Freq1 < df2$Freq2, df2$Freq1, df2$Freq2))
  	return(df2)
}
```

Get frequency bins to graph on local computer

Missense variants
```R
SAM_MISSENSE_freq <- vcftoolsFreqFormat("Missense/SAM_MISSENSE.frq")
SAM_missense_freqbins <- hist(SAM_MISSENSE_freq$freq, plot=FALSE)
#SAM_missense_probbins <- hist(SAM_MISSENSE_freq$freq, plot=FALSE, freq='FALSE')

WILD_MISSENSE_freq <- vcftoolsFreqFormat("Missense/WILD_MISSENSE.frq")
WILD_missense_freqbins <- hist(WILD_MISSENSE_freq$freq, plot=FALSE)

save(SAM_missense_freqbins, WILD_missense_freqbins, file="MissenseBins.RData")
```

Synonymous variants
```R
SAM_SYNON_freq <- vcftoolsFreqFormat("Synon/SAM_SYNON.frq")
SAM_synon_freqbins <- hist(SAM_SYNON_freq$freq, plot=FALSE)

WILD_SYNON_freq <- vcftoolsFreqFormat("WILD_SYNON.frq")
WILD_synon_freqbins <- hist(WILD_SYNON_freq$freq, plot=FALSE)

save(SAM_synon_freqbins, WILD_synon_freqbins, file="SynonymousBins.RData")
```

Graphs
```R
png("Sam_missense.png")
hist(SAM_MISSENSE_freq$freq)
dev.off()

png("Wild_missense.png")
hist(WILD_MISSENSE_freq$freq)
dev.off()

png("Sam_synonymous.png")
hist(SAM_SYNON_freq$freq)
dev.off()

png("Wild_synonymous.png")
hist(WILD_SYNON_freq$freq)
dev.off()
```

#### Nicer Graphs
Function to convert from hist format

```R
HistBinsFormat <- function(list1, list1_name, list2, list2_name) {
  Freqbins_df <- as.data.frame(cbind(list1$breaks[-1], list1$counts, list2$counts))
  colnames(Freqbins_df) <- c("breaks", list1_name, list2_name)
  sum_list1 <- sum(Freqbins_df[,2])
  sum_list2 <- sum(Freqbins_df[,3])
  Freqbins_df_long <- melt(Freqbins_df, id.vars = c("breaks"))
  Freqbins_df_long$prop <- ifelse(Freqbins_df_long$variable==list1_name, 
                                   Freqbins_df_long$value /  sum_list1,
                                   Freqbins_df_long$value /  sum_list2)
  return(Freqbins_df_long)
}
```

Graphs
```R
load("SynonymousBins.RData")
load("MissenseBins.RData")
library(ggplot2)
library(reshape2)
library(ggthemes)
library(RColorBrewer)

SAM_bins <- HistBinsFormat(SAM_synon_freqbins, "SAM_synon", SAM_missense_freqbins, "SAM_missense")
head(SAM_bins)

q <- ggplot(SAM_bins, aes(x=breaks, y=value, fill=variable))
q2 <- ggplot(SAM_bins, aes(x=breaks, y=prop, fill=variable))
q2 + geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values= c("tomato1", "tomato4")) + theme_minimal() + ylim(0, 0.2)
ggsave("SAM_sfs_newY.png")

Wild_bins <- HistBinsFormat(WILD_synon_freqbins, "Wild_synon", WILD_missense_freqbins, "Wild_missense")
p <- ggplot(Wild_bins, aes(x=breaks, y=value, fill=variable))
p2 <- ggplot(Wild_bins, aes(x=breaks, y=prop, fill=variable))
p2 + geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values= c("skyblue1", "skyblue4")) + theme_minimal()
ggsave("WILD_sfs.png")

Synon_bins <- HistBinsFormat(WILD_synon_freqbins, "Wild_synon", SAM_synon_freqbins, "SAM_synon")
b <- ggplot(Synon_bins, aes(x=breaks, y=prop, fill=variable))
b + geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values= c("#1B9E77", "#D95F02")) + theme_minimal()

Miss_bins <- HistBinsFormat(WILD_missense_freqbins, "Wild_missense", SAM_missense_freqbins, "SAM_missense")
g <- ggplot(Miss_bins, aes(x=breaks, y=prop, fill=variable))
g + geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values= c("#1B9E77", "#D95F02")) + theme_minimal()
```

Testing on subset - 01/15/21
```R
SAM_freq <- vcftoolsFreqFormat("SAM_SNPs_SUBSETFINAL_Biallelic.frq")
SAM_freqbins <- hist(SAM_freq$freq, plot=FALSE)
save(SAM_freqbins, file="SAM_subsetFreqBins.RData")

Freqbins_df <- as.data.frame(cbind(SAM_freqbins$breaks[-1], SAM_freqbins$counts))
colnames(Freqbins_df) <- c("breaks", "counts")
sumcount <- sum(Freqbins_df$counts)
Freqbins_df$prop <- Freqbins_df$counts/sumcount

p <- ggplot(Freqbins_df, aes(x=breaks, y=prop))
p + geom_bar(stat="identity") + theme_minimal()

```

## Deleterious predictions

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
vcftools --vcf SAM_deleterious.vcf --freq --out SAM_deleterious
vcftools --vcf SAM_tolerated.vcf --freq --out SAM_tolerated # redid after removing duplicate positions
vcftools --vcf SAM_synonymous.vcf --freq --out SAM_synonymous

module load R/4.0.0-foss-2019b
R
```

Using vcftoolsFreqFormat function specified above
```R
hist_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

SAM_deleterious_freq <- vcftoolsFreqFormat("SAM_deleterious.frq")
SAM_deleterious_freqbins <- hist(SAM_deleterious_freq$freq, plot=FALSE, breaks=hist_breaks)
save(SAM_deleterious_freqbins, file="DeleteriousBins.RData")

SAM_tolerated_freq <- vcftoolsFreqFormat("SAM_tolerated.frq")
SAM_tolerated_freqbins <- hist(SAM_tolerated_freq$freq, plot=FALSE, breaks=hist_breaks)
save(SAM_tolerated_freqbins, file="ToleratedBins.RData")

SAM_synonymous_freq <- vcftoolsFreqFormat("SAM_synonymous.frq")
SAM_synonymous_freqbins <- hist(SAM_synonymous_freq$freq, plot=FALSE, breaks=hist_breaks)
save(SAM_synonymous_freqbins, file="SynonymousBins.RData")
```

On local computer
```R
setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Results")
load("ToleratedBins.RData")
load("DeleteriousBins.RData")
load("SynonymousBins.RData")

library(ggplot2)
library(reshape2)
library(ggthemes)
library(RColorBrewer)

SAM_bins <- HistBinsFormat(SAM_tolerated_freqbins, "tolerated", 
SAM_deleterious_freqbins, "deleterious")

head(SAM_bins)

Freqbins_df <- as.data.frame(cbind(SAM_synonymous_freqbins$breaks[1:length(SAM_synonymous_freqbins$counts)], SAM_synonymous_freqbins$counts, SAM_tolerated_freqbins$counts, SAM_deleterious_freqbins$counts))

colSums(Freqbins_df) # check numbers

colnames(Freqbins_df) <- c("breaks", "Synonymous", "Tolerated",
                           "Deleterious")
Freqbins_df_long <- melt(Freqbins_df, id.vars = c("breaks"))

Freqbins_df_long$prop <- ifelse(Freqbins_df_long$variable=="Synonymous", 
                            Freqbins_df_long$value / sum(Freqbins_df$Synonymous),
                            ifelse(Freqbins_df_long$variable=="Tolerated",
                                   Freqbins_df_long$value / sum(Freqbins_df$Tolerated),
                            Freqbins_df_long$value / sum(Freqbins_df$Deleterious)))

p <- ggplot(Freqbins_df_long, aes(x=breaks, y=prop, fill=variable))
p + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2")
scale_fill_manual(values= c("#7570B3", "#1B9E77", "#D95F02")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.5)"))
ggsave("SFS_test.eps") # I manually changed the scaling with Rstudio GUI
```