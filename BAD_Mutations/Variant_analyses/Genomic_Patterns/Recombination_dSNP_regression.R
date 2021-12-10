#### Is there a relationship between recombination and # dSNPs/codon?

library(lme4)

######################################
######### SET UP DATAFRAME ###########
######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationBins.RData")
# object: Recomb_binData (Median cM/Mbp summarized per window)


VariantNums <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/Derived_VariantNums.txt",
                          header=T, sep="\t")

Recomb_binData <- lapply(names(Recomb_binData), function(x) {
  Recomb_binData[[x]]$Chromosome <- x; return(Recomb_binData[[x]])
})
Recombination <- do.call("rbind", Recomb_binData)
Recombination$StartPos <- as.integer(gsub('[(]', '', 
                                          sapply(strsplit(as.character(Recombination$Bin), ","), "[", 1)))*1000000

All_info <- merge(VariantNums, Recombination, by=c("Chromosome", "StartPos"))

######################################
########## MODEL dSNP/CODON ##########
######################################

All_info$dSNP_codon <- All_info$Number_dSNP / All_info$Num_codons
hist(All_info$dSNP_codon) # looks pretty normally distributed

Mod1 <- lm(dSNP_codon ~ `Median_cM/Mbp`,
           data = All_info)
drop1(Mod1, test = "Chi")
summary(Mod1) # slightly positive slope

plot(All_info$`Median_cM/Mbp`, All_info$dSNP_codon)
abline(Mod1)

######################################
########## MODEL dSNP/sSNP ###########
######################################

All_info$dSNP_sSNP <- All_info$Number_dSNP / All_info$Number_Synonymous
hist(log2(All_info$dSNP_sSNP+1))
All_info$log2_dSNP_sSNP <- log2(All_info$dSNP_sSNP+1)

hist(All_info$`Median_cM/Mbp`)
hist(log2(All_info$`Median_cM/Mbp`+1))
All_info$log2_cM_Mbp <- log2(All_info$`Median_cM/Mbp`+1)

Mod2 <- lmer(log2_dSNP_sSNP ~ log2_cM_Mbp + (1|Chromosome),
           data = All_info, na.action = na.omit)
Mod2 <- lm(log2_dSNP_sSNP ~ log2_cM_Mbp,
             data = All_info, na.action = na.omit)
drop1(Mod2)
summary(Mod2) # negative slope
hist(resid(Mod2))

plot(All_info$log2_cM_Mbp, All_info$log2_dSNP_sSNP)
abline(Mod2)

cor.test(All_info$log2_cM_Mbp, All_info$log2_dSNP_sSNP, use="pairwise.complete.obs", method="spearman")
# rho= -0.4055, p<0.0001

