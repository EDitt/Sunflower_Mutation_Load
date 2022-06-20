## plots for looking at private/shared alleles among heterotic groups across different frequency bins

library(ggplot2)
library(ggpubr)

source("BAD_Mutations/Variant_analyses/Functions.R")

# plot code formerly in README.md
# continued from Heterotic_AlleleInfo.R

######################################
######### SET UP DATAFRAME ###########
######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/FreqBins.RData")

######################################
#### ACROSS OIL/NON-OIL, ALL SNPS ####
######################################

plot <- function(df) {
  ggplot(data=df, aes(x=Bin, y=PropPrivate, fill=Annotation)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=NumPrivate), position=position_dodge(width=0.9), vjust=1.25, colour="white") +
    scale_fill_manual(values=c(Col_Deleterious, Col_Synonymous),
                      labels = c("Deleterious", "Synonymous")) +
    theme_minimal() +
    ylab("") + xlab("")
}
p1 <- plot(Both_FreqBins[which(Both_FreqBins$Bin=="(0,0.1]"),])
p2 <- plot(Both_FreqBins[which(Both_FreqBins$Bin!="(0,0.1]"),])
p3 <- annotate_figure(ggarrange(p1,p2, common.legend = TRUE, widths = c(1,3)), 
                      bottom = text_grob("Frequency Bin"),
                      left=text_grob("Proportion Private", rot=90))
pdf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/HA_RHA_Private_All.pdf")
print(p3)
dev.off()

######################################
######### OIL ONLY, ALL SNPS #########
######################################

p1 <- plot(Oil_FreqBins[which(Oil_FreqBins$Bin=="(0,0.1]"),])
p2 <- plot(Oil_FreqBins[which(Oil_FreqBins$Bin!="(0,0.1]"),])
p3 <- annotate_figure(ggarrange(p1,p2, common.legend = TRUE, widths = c(1,3)), 
                      bottom = text_grob("Frequency Bin"),
                      left=text_grob("Proportion Private", rot=90))
pdf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/HA_RHA_Oil_Private_All.pdf")
print(p3)
dev.off()

######################################
######## NONOIL ONLY, ALL SNPS #######
######################################

p1 <- plot(NonOil_FreqBins[which(NonOil_FreqBins$Bin=="(0,0.1]"),])
p2 <- plot(NonOil_FreqBins[which(NonOil_FreqBins$Bin!="(0,0.1]"),])
p3 <- annotate_figure(ggarrange(p1,p2, common.legend = TRUE, widths = c(1,3)), 
                      bottom = text_grob("Frequency Bin"),
                      left=text_grob("Proportion Private", rot=90))
pdf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/HA_RHA_NonOil_Private_All.pdf")
print(p3)
dev.off()


######################################
############ DERIVED SNPS ############
######################################

ggplot(data=Both_DerivedFreqBins, aes(x=Bin, y=PropPrivate, fill=Annotation)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=NumPrivate), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_manual(values=c(Col_Deleterious, Col_Synonymous),
                    labels = c("Deleterious", "Synonymous")) +
  theme_minimal() +
  ylab("") + xlab("")

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/HA_RHA_Derived_Private_All.pdf",
       width = 10, height = 5)

