Analyses to look at the Distribution of dSNPs among germplasm

# Number of variants per SAM line

Which dSNPs are alternate allele state? (from reference)
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load R/4.0.0-foss-2019b
R
```
```R
dsnp_data <- read.table("dsnp_data.table", sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
dsnp_del <- subset(dsnp_data, Result == "Deleterious")

aggregate(dsnp_del$VariantID, by=list(dsnp_del$Refderived, dsnp_del$Altderived), length)
# All 54,445 are in the derived state for the alternate allele only

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Refderived, dsnp_data$Altderived), length)
# Out of all variants, 5 are derived for both reference and alternate allele,
# 284,934 are derived for alternate allele only, 
# 360274 are not derived in either

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Refderived), length) # only 5 of reference are derived?
```
All SNPs annotated as deleterious are in the alternate allele state

Use bcftools to count the number of alternate alleles per genotype

```bash
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_deleterious.vcf > DeleteriousStatsperSample.txt
grep "PSC" DeleteriousStatsperSample.txt > DeleteriousperSampleCounts.txt
```
Saved this file to local computer

```R
setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Results")
allele_counts <- read.table("DeleteriousperSampleCounts.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(allele_counts) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", 
	"nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")
length(allele_counts$sample) #288

key <- read.csv("LineKeywINFO.csv", header=T)

counts_info <- merge(key, allele_counts, by.x="SequenceName", by.y="sample")
length(counts_info$sort) #288

# go by only non-ref homozygous or # of alleles?
counts_info$NumdSNPs <- 2*counts_info$nNonRefHom + counts_info$nHets

```

ANOVA to look at whether categories effect nNonRefHom
Note: I will need to do this with derived dSNPs because right now this pattern likely reflects similarity to the reference genome
```r
library(car)
library(emmeans)
library(ggplot2)
mod1 <- lm(nNonRefHom ~ heterotic_group + Oil_NonOil+
             heterotic_group:Oil_NonOil,
           data=counts_info)
# do I need to use Poisson distribution since I'm using count data?
summary(mod1)
Anova(mod1, test= "Chi") # all categories are significant
hist(resid(mod1)) # normal
hist(counts_info$nNonRefHom) #normal

Mod2 <- glm(nNonRefHom ~ heterotic_group + Oil_NonOil+
             heterotic_group:Oil_NonOil,
                  family=poisson(link = "log"), data = counts_info)
drop1(Mod2, test = "Chi")
Anova(Mod2, test = "F") # all significant

Mod_means <- emmeans(Mod2, ~ heterotic_group | Oil_NonOil, type = "response")
Mod_means
pairs(Mod_means)

# heterotic group?
emmeans(Mod2, ~ heterotic_group, type = "response")

emmip(Mod2, ~ heterotic_group | Oil_NonOil, type = "response")

Mean.df <- as.data.frame(Mod_means)
Mean.df$group <- paste0(Mean.df$heterotic_group,"-", Mean.df$Oil_NonOil)
levels(as.factor(Mean.df$group))
Mean.df_sub <- subset(Mean.df, !is.na(rate))
Mean.df_sub$group <- factor(Mean.df_sub$group, levels=c(
  "HA-INRA", "RHA-INRA", "HA-NonOil", "HA-Oil", "RHA-NonOil",
  "RHA-Oil", "introgressed-NonOil", "introgressed-Oil", "landrace-NonOil",
  "OPV-NonOil", "other-NonOil", "other-Oil"
))
levels(as.factor(Mean.df_sub$group))

p <- ggplot(data=Mean.df_sub, mapping=aes(x=group, y=rate))
p + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=rate-SE, ymax=rate+SE,
                    width=.2))

```