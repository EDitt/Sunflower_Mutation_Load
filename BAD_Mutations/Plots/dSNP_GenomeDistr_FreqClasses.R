### A large graph showing a bargraph of the total number of dSNPs per codon for 10 Mbp regions, and
### the proportion of rare (<10% derived), segregating (10-50% MAF, >90% derived),
### a scatterplot above showing dSNPs/sSNPs for each region

source("BAD_Mutations/Variant_analyses/Functions.R")

library(ggplot2)
library(ggsci)
library(cowplot)
library(magrittr)

######################################
######## READ IN DATAFRAME ###########
######################################

#dSNPNums <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/dSNPperCodonNums.txt",
#                       sep = "\t", header=T)

# put in directory info:
MyDir<-c("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/")

load(paste0(MyDir, "dSNP_genomeGraph.RData"))

levels(as.factor(dSNPNums$FrequencyClass))

# split by chromosome for graphing
DelCod_Chrom <- split(dSNPNums, dSNPNums$Chromosome)

# split dSNP/sSNP dataframe
dSNP_sSNPChrom <- split(dSNP_sSNP, dSNP_sSNP$Chromosome)

######################################
########### COLOR SCHEME #############
######################################

library("scales")
show_col(pal_jco("default")(10))
show_col(pal_npg)

SingletonCol <- "#003C67FF"
LowMAFCol <- "#7AA6DCFF"
HighMAFCol <- "#EFC000FF" 
HighFreqCol <- "#A73030FF"


######################################
############ BAR GRAPHS ##############
######################################


Genome_bargraph <- function(dataset, YMax, HorizIntercept) {
  plot <- ggplot(dataset, aes(x=Mbp, y=NumPerCodon, fill=FrequencyClass, group=FrequencyClass)) +
    geom_bar(stat="identity") +
    #geom_col(width=ColWidth) +
    theme_minimal() +
    scale_fill_manual(
      values= c(HighFreqCol, HighMAFCol,
                LowMAFCol, SingletonCol),
      name = "", labels = c("High Derived Frequency (>0.9)", 
                            "High MAF (0.1-0.5)",
                            "Low MAF (0-0.1, excluding singletons)",
                            "Only in Single Individuals")) +
    ylim(0,YMax) +
    geom_hline(yintercept = HorizIntercept, linetype="dashed")+
  scale_x_continuous(limits = c(0, 230))
  return(plot)
}


dSNPFreqPlots <- lapply(DelCod_Chrom, function(x) {
  Genome_bargraph(x, 0.0022, MeandSNPCod+SDdSNPCod)})

labels <- paste0("Chr ", seq(1,17, by=1))

ggarrange(plotlist = dSNPFreqPlots, labels=labels, common.legend = TRUE)

#ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure5_NumFreqGenome1.pdf")

######################################
####### FACET WRAP ALTERNATIVE #######
######################################

p1 <- Genome_bargraph(dSNPNums, 0.0022, MeandSNPCod+SDdSNPCod) 
p1 + facet_wrap(~Chromosome) ## here the width of the bars are not consistent
p1 + facet_grid(~Chromosome, scales = "free_x", space = "free_x") # on one line
p1 + facet_wrap(~Chromosome, scales = "free_x") # this could be used instead of the scale_x_continuous(limits = c(0, 230)) to allow X axes to vary


######################################
#### ADD 2ND Y-AXIS FOR DSNP/SSNP ####
######################################

# fit dSNP/sSNP between 0.0015 and 0.0020 (diff of 0.0005)
min(dSNP_sSNP$dSNP_sSNP) # 0.04
max(dSNP_sSNP$dSNP_sSNP) # 0.34
# spans 0.3
# divide by 600 and add 0.0015?

Genome_bargraph2 <- function(dataset, YMax, HorizIntercept, dataset2) {
  plot <- ggplot(dataset, aes(x=Mbp, y=NumPerCodon)) +
    geom_bar(stat="identity", aes(fill=FrequencyClass, group=FrequencyClass)) +
    #geom_col(width=ColWidth) +
    theme_minimal() +
    scale_fill_manual(
      values= c(HighFreqCol, HighMAFCol,
                LowMAFCol, SingletonCol),
      name = "", labels = c("High Derived Frequency (>0.9)", 
                            "High MAF (0.1-0.5)",
                            "Low MAF (0-0.1, excluding singletons)",
                            "Only in Single Individuals")) +
    ylim(0,YMax) +
    geom_hline(yintercept = HorizIntercept, linetype="dashed")+
    scale_x_continuous(limits = c(0, 230)) +
    geom_line(data=dataset2, 
              aes(x=Mbp, 
                  y=(dSNP_sSNP/500)+0.0015, linetype="dashed"), 
              show.legend = FALSE)
  return(plot)
}



dSNPFreqPlots2 <- lapply(names(DelCod_Chrom), function(x) {
  Genome_bargraph2(DelCod_Chrom[[x]],
                   0.0022, MeandSNPCod+SDdSNPCod,
                  dSNP_sSNPChrom[[x]]
                  )})

dSNPFreqPlots2[[14]]

ggarrange(plotlist = dSNPFreqPlots2, labels=labels, common.legend = TRUE)




######################################
######################################
######################################  BELOW IS UNDER DEVELOPMENT/SCRATCH
######################################
######################################



######################################
##### ADD 2ND GRAPH W/ dSNP/sSNP? #####
######################################

Marginal_points <- function(plot1, dataset, yMax) {
  plot2 <- axis_canvas(plot1, axis = "x") +
    geom_line(data=dataset, 
              aes(x=Mbp, 
                  y=dSNP_sSNP), linetype="dashed") +
    geom_point(data=dataset, 
               aes(x=Mbp, 
                   y=dSNP_sSNP)) +
    ylab("dSNP/sSNP") 
  #+ylim(0, yMax)
  combined_plot <- insert_xaxis_grob(plot1, plot2, position = "top")
  return(combined_plot)
}


dSNPFreqwMarginPlots <- lapply(names(DelCod_Chrom), function(x) {
  Marginal_points(dSNPFreqPlots[[x]], 
                  dSNP_sSNPChrom[[x]],
                  0.34)})

ggdraw(dSNPFreqwMarginPlots[[3]]) # this looks fine but the grid plots aren't showing up?
plot_grid(plotlist = dSNPFreqwMarginPlots)
plot_test <- lapply(dSNPFreqwMarginPlots, function(x) {ggdraw(x)})
plot_grid(plotlist = plot_test)

#ggarrange(plotlist = dSNPFreqwMarginPlots, labels=labels, common.legend = TRUE)



####
ggplot(data=test1, aes(x=Mbp, y=NumPerCodon)) +
  geom_bar(stat="identity", aes(fill = FrequencyClass)) + 
  scale_fill_manual(values = c(HighFreqCol, HighMAFCol,
                               LowMAFCol, SingletonCol)) +
  theme_minimal() +
  ylim(0,0.0025) +
  geom_point(data=dSNP_sSNPChrom$Ha412HOChr01, 
             aes(x=dSNP_sSNPChrom$Ha412HOChr01$Mbp, 
                 y=dSNP_sSNPChrom$Ha412HOChr01$dSNP_sSNP/50), shape=21) +
  geom_line(data=dSNP_sSNPChrom$Ha412HOChr01, 
            aes(x=dSNP_sSNPChrom$Ha412HOChr01$Mbp, 
                y=dSNP_sSNPChrom$Ha412HOChr01$dSNP_sSNP/50), linetype="dashed")




p1 <- ggplot(data=test1, aes(x=Mbp, y=NumPerCodon)) +
  geom_bar(stat="identity", aes(fill = FrequencyClass)) + 
  scale_fill_manual(values = c(HighFreqCol, HighMAFCol,
                               LowMAFCol, SingletonCol)) +
  theme_minimal() +
  ylim(0,0.002) 

p2 <- axis_canvas(p1, axis = "x", coord_flip = TRUE) +
  geom_line(data=dSNP_sSNPChrom$Ha412HOChr01, 
            aes(x=Mbp, 
                y=dSNP_sSNP), linetype="dashed") +
  geom_point(data=dSNP_sSNPChrom$Ha412HOChr01, 
             aes(x=Mbp, 
                 y=dSNP_sSNP)) +
  coord_flip() 
#+ylab("dSNP/sSNP")

combined_plot <- insert_xaxis_grob(p1, p2, position = "top")
combined_plot %<>% insert_xaxis_grob(., p2, position = "top")

ggdraw(combined_plot)

### scratch
dSNP_sSNPChrom <- split(GenDf_Rel[,c(1,2,13)], GenDf_Rel$Chromosome)
ggplot(data=dSNP_sSNPChrom$Ha412HOChr01, aes(x=Mbp, y=dSNP_sSNP)) +
  geom_point() + geom_line(linetype="dashed") + theme_minimal()


dSNPFreqPlots$Ha412HOChr01 + geom_point(data=dSNP_sSNPChrom$Ha412HOChr01, 
                                        x=dSNP_sSNPChrom$Ha412HOChr01$Mbp, 
                                        y=dSNP_sSNPChrom$Ha412HOChr01$dSNP_sSNP)

#### scratch

p2 <- ggMarginal(p, type = "histogram", fill = "grey")
return(p2)

######################################
############ FACET WRAP ##############        <- is not working with this dataset?
######################################

Newdf <- as.data.frame(GenDfDel_Rel_long[c(1:5)])
p1 <- Genome_bargraph(Newdf, 0.0022, MeandSNPCod+SDdSNPCod) 
p1 + facet_wrap(~Chromosome)

### use geom_col instead of geom_bar?
ggplot(data=GenDfDel_Rel_long, aes(x=Mbp, y=NumPerCodon, fill=FrequencyClass, group = FrequencyClass)) +
  geom_col(width=8) + 
  scale_fill_manual(values = c(HighFreqCol, HighMAFCol,
                               LowMAFCol, SingletonCol)) +
  theme_minimal() +
  facet_wrap(~Chromosome, nrow = 5, ncol = 4, scales = "free") +
  scale_x_continuous() + ylim(0,0.0022)


GenDfDel_Rel_long$ChromNum <- as.numeric(gsub("Ha412HOChr", "", GenDfDel_Rel_long$Chromosome))

MaxLengths <- aggregate(GenDfDel_Rel_long$Mbp, by=list(GenDfDel_Rel_long$Chromosome),
                        max)

#### scratch

### test barplot with chr 1
test1 <- DelCod_Chrom$Ha412HOChr01

test1_long <- reshape(test1[,c(1:5,7)],
                      #varying=list(3:5,7),
                      varying=names(test1[c(3:5,7)]),
                      idvar = c("Chromosome", "Mbp"),
                      times = names(test1[c(3:5,7)]),
                      v.names = "NumPerCodon",
                      direction = "long"
)

ggplot(data=test1, aes(x=Mbp, y=NumPerCodon, fill=FrequencyClass)) +
  geom_bar(stat="identity") + 
  #scale_fill_jco()
  #scale_fill_npg()
  scale_fill_manual(values = c(HighFreqCol, HighMAFCol,
                               LowMAFCol, SingletonCol))


Genome_bargraph(test1, 0.0022, MeandSNPCod+SDdSNPCod)

dSNPFreqPlots$Ha412HOChr01

######################################
####### ALTERNATIVE 2 ######
######################################

dSNPNums$ChromosomeNum <- as.numeric(gsub("Ha412HOChr", "", dSNPNums$Chromosome))
dSNPNums$Chromosome_Position <- as.numeric(paste0(dSNPNums$ChromosomeNum, ".", dSNPNums$Mbp))

ggplot(dSNPNums, aes(x=Chromosome_Position, y=NumPerCodon, fill=FrequencyClass, group=FrequencyClass)) +
  geom_bar(stat="identity") +
  scale_fill_manual(
    values= c(HighFreqCol, HighMAFCol,
              LowMAFCol, SingletonCol),
    name = "", labels = c("High Derived Frequency (>0.9)", 
                          "High MAF (0.1-0.5)",
                          "Low MAF (0-0.1, excluding singletons)",
                          "Only in Single Individuals")) +
  theme_minimal()
