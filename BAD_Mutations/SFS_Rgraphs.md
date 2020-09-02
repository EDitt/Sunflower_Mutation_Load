### Folded site frequency spectrum

```bash
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
SAM_MISSENSE_freq <- vcftoolsFreqFormat("SAM_MISSENSE.frq")
SAM_missense_freqbins <- hist(SAM_MISSENSE_freq$freq, plot=FALSE)

WILD_MISSENSE_freq <- vcftoolsFreqFormat("WILD_MISSENSE.frq")
WILD_missense_freqbins <- hist(WILD_MISSENSE_freq$freq, plot=FALSE)

save(SAM_missense_freqbins, WILD_missense_freqbins, file="MissenseBins.RData")
```

Synonymous variants
```R
SAM_SYNON_freq <- vcftoolsFreqFormat("SAM_SYNON.frq")
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
