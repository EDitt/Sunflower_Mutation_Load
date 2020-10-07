### Folded site frequency spectrum

```bash
cd /scratch/eld72413/NSFproj/VEP
module load R/3.6.1-foss-2018a-X11-20180131-GACRC
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