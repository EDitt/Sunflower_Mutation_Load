### Genome patterns

source("BAD_Mutations/Variant_analyses/Functions.R")

library(ggplot2)
library(car)
library(lme4)

######################################
######### SET UP DATAFRAME ###########
######################################

my_files <- list.files(path = "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp", 
                       pattern = ".bed", full.names = TRUE)

GenData <- lapply(my_files, function(x) {read.table(x, header=F, na.strings=c("NaN"))})
names(GenData) <- gsub(".bed", "", my_files)
names(GenData) <- gsub("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/", "", 
                       names(GenData))

# apply common column names
GenData <- lapply(GenData, function(x) {
  colnames(x)[1:3] <- c("Chromosome", "StartPos", "EndPos"); return(x)
})

colnames(GenData$Num_SynonymousNodups_10MbpCounts)[4] <- "NumSynonymous"
colnames(GenData$Num_AllDel_10MbpCounts)[4] <- "NumDeleterious"
colnames(GenData$mRNA_10MbpNumCodonCounts)[5:7] <- c("TotalCodingBases", "TotalmRNA", "Num_Codons")

AllGenData <- merge(GenData$mRNA_10MbpNumCodonCounts,
                    merge(GenData$Num_SynonymousNodups_10MbpCounts, GenData$Num_AllDel_10MbpCounts,
                          by=c("Chromosome", "StartPos", "EndPos")),
                    by=c("Chromosome", "StartPos", "EndPos"))
length(AllGenData$Chromosome) # check -325


######################################
######### CALCULATE METRICS ##########
######################################


AllGenData$dSNP_codon <- AllGenData$NumDeleterious / AllGenData$Num_Codons
AllGenData$sSNP_codon <- AllGenData$NumSynonymous / AllGenData$Num_Codons

AllGenData$dSNP_sSNP <- AllGenData$NumDeleterious / AllGenData$NumSynonymous
AllGenData$dSNPcod_sSNPcod <- AllGenData$dSNP_codon / AllGenData$sSNP_codon

AllGenData$MiddlePos <- ((AllGenData$EndPos - AllGenData$StartPos) / 2) + AllGenData$StartPos
AllGenData$Mbp <- AllGenData$MiddlePos / 1000000


###############################
######## PLOT FUNCTION ########
###############################

Genome_plot <- function(dataset, X_column, y1, y2, y1_col, y2_col, y2scale, yMax, HorizIntercept) {
  plot <- ggplot(dataset, aes(x=X_column, y=y1)) +
    geom_point(fill=y1_col, size=2, shape=21) +
    geom_line(col=y1_col, linetype="dotted") + 
    geom_point(aes(y=y2/y2scale), fill=y2_col, size=2, shape=24) +
    geom_line(aes(y=y2/y2scale), col=y2_col, linetype="dotted") + 
    geom_hline(yintercept = HorizIntercept, linetype="dashed", col=y1_col) +
    theme_minimal() +
    scale_y_continuous(name="dSNP/codon", limits=c(0,yMax),
                       sec.axis = sec_axis(~ .*y2scale, name= "sSNP/codon")) +
    xlab("Mbp")
  return(plot)
}


#ggplot(ChromosomeData$Ha412HOChr01, aes(x=Mbp, y=dSNP_codon)) +
#  geom_point(col=Col_Deleterious, size=2) +
#  geom_line(col=Col_Deleterious, linetype="dotted") + 
#  geom_point(aes(y=sSNP_codon/10), col=Col_Synonymous, size=2) +
#  geom_line(aes(y=sSNP_codon/10), col=Col_Synonymous, linetype="dotted") + 
#  theme_minimal() +
#  scale_y_continuous(sec.axis = sec_axis(~ .*10, name= "sSNP/codon"))

######################################
########### PLOT SNP/CODON ###########
######################################

AllGenData$Chromosome <- factor(AllGenData$Chromosome) # drop unused levels
meanDel <- mean(AllGenData$dSNP_codon)
sdDel <- sd(AllGenData$dSNP_codon)
ChromosomeData <- split(AllGenData, AllGenData$Chromosome)


# test
Genome_plot(ChromosomeData$Ha412HOChr01, ChromosomeData$Ha412HOChr01$Mbp, 
            ChromosomeData$Ha412HOChr01$dSNP_codon, ChromosomeData$Ha412HOChr01$sSNP_codon,
            Col_Deleterious, Col_Synonymous, 10, 0.0025, meanDel+2*sdDel)

GenomePlots <- lapply(ChromosomeData, function(x) {
  Genome_plot(x, x$Mbp, 
              x$dSNP_codon, x$sSNP_codon,
              Col_Deleterious, Col_Synonymous, 10, 0.0025, meanDel+2*sdDel)})

labels <- paste0("Chromosome ", seq(1,17, by=1))

ggarrange(plotlist = GenomePlots, labels=labels) 


## or
Genome_plot(AllGenData, AllGenData$Mbp, AllGenData$dSNP_codon,
            AllGenData$sSNP_codon,
            Col_Deleterious, Col_Synonymous, 10, 0.0025, meanDel+sdDel) +
  facet_wrap(~Chromosome, scales = "free_x") 

######################################
########### PLOT dSNP/sSNP ###########
######################################

Genome_plot2 <- function(dataset, X_column, y, y_col, yMax) {
  plot <- ggplot(dataset, aes(x=X_column, y=y)) +
    geom_point(col=y_col, size=2) +
    geom_line(col=y_col, linetype="dotted") + 
    theme_minimal() +
    ylab("dSNP/sSNP") +
    xlab("Mbp") +
    ylim(0,yMax)
  return(plot)
}

# test
Genome_plot2(ChromosomeData$Ha412HOChr01, ChromosomeData$Ha412HOChr01$Mbp, 
            ChromosomeData$Ha412HOChr01$dSNP_sSNP, "blue", 0.35)

GenomePlots2 <- lapply(ChromosomeData, function(x) {
  Genome_plot2(x, x$Mbp, 
              x$dSNP_sSNP, "blue", 0.36)})

ggarrange(plotlist = GenomePlots2, labels=labels) 

######################################
###### HIGH FREQ SNPS - SETUP ########
######################################

HighFreqNums <- ImportFilesAsList("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp",
                                  "Num_High", "_10MbpCounts.bed", F)


# apply common column names & # merge with codon info
HighFreqNums <- lapply(HighFreqNums, function(x) {
  colnames(x) <- c("Chromosome", "StartPos", "EndPos", "Count"); 
  return(x)
})

HighFreqData <- lapply(HighFreqNums, function(x) {
  merge(x, GenData$mRNA_10MbpNumCodonCounts, by=c("Chromosome", "StartPos", "EndPos"))
})

# add column before 
for (i in seq_along(HighFreqData)) {
  name <- names(HighFreqData[i])
  HighFreqData[[i]]["Count_type"] <- name
  HighFreqData[[i]]["CountPerCodon"] <- HighFreqData[[i]]["Count"] / HighFreqData[[i]]["Num_Codons"]
  HighFreqData[[i]]["Mbp"] <- (HighFreqData[[i]]["StartPos"] + ((HighFreqData[[i]]["EndPos"] - HighFreqData[[i]]["StartPos"]) / 2)) / 1000000
}

######################################
######## HIGH DERIVED FREQ ###########
######################################

HighFreq_der <- do.call("rbind", HighFreqData[1:2])
HighFreq_derWide <- reshape(HighFreq_der[,c(1:4,8:11)], 
                            timevar = "Count_type",
                            idvar = c("Chromosome", "StartPos", "EndPos", "Mbp"),
                            direction = "wide")
mean_HighDerFreqCodon <- mean(HighFreq_derWide$CountPerCodon.DerFreqAllDel) # .000032
sd_HighDerFreqCodon <- sd(HighFreq_derWide$CountPerCodon.DerFreqAllDel)

# split by Chromosome
HighFreq_derChrom <- split(HighFreq_derWide, HighFreq_derWide$Chromosome)


# test
Genome_plot(HighFreq_derChrom$Ha412HOChr01, HighFreq_derChrom$Ha412HOChr01$Mbp, 
            HighFreq_derChrom$Ha412HOChr01$CountPerCodon.DerFreqAllDel, 
            HighFreq_derChrom$Ha412HOChr01$CountPerCodon.DerFreqSynonymousNodups,
            Col_Deleterious, Col_Synonymous, 50, 0.00020, mean_HighDerFreqCodon + sd_HighDerFreqCodon)

GenomeDerFreqPlots <- lapply(HighFreq_derChrom, function(x) {
  Genome_plot(x, x$Mbp, 
              x$CountPerCodon.DerFreqAllDel, 
              x$CountPerCodon.DerFreqSynonymousNodups,
              Col_Deleterious, Col_Synonymous, 50, 0.0002, mean_HighDerFreqCodon + 2*sd_HighDerFreqCodon)})


ggarrange(plotlist = GenomeDerFreqPlots, labels=labels) 

######################################
###### HIGH MINOR ALLELE FREQ ########
######################################

HighMAF <- do.call("rbind", HighFreqData[3:4])


HighMAF_Wide <- reshape(HighMAF[,c(1:4,8:11)], 
                            timevar = "Count_type",
                            idvar = c("Chromosome", "StartPos", "EndPos", "Mbp"),
                            direction = "wide")
mean_HighMAFCodon <- mean(HighMAF_Wide$CountPerCodon.MAFAllDel) # .000156
sd_HighMAFCodon <- sd(HighMAF_Wide$CountPerCodon.MAFAllDel)
max(HighMAF_Wide$CountPerCodon.MAFAllDel) # 0.00038

HighMAF_Wide[which(HighMAF_Wide$CountPerCodon.MAFAllDel> 0.00035),] # Chromosome 10: 4-5,12-14,15-17; & segments on Chr 3,5,6,8

# split by Chromosome
HighMAF_Chrom <- split(HighMAF_Wide, HighMAF_Wide$Chromosome)

# test
Genome_plot(HighMAF_Chrom$Ha412HOChr01, HighMAF_Chrom$Ha412HOChr01$Mbp, 
            HighMAF_Chrom$Ha412HOChr01$CountPerCodon.MAFAllDel, 
            HighMAF_Chrom$Ha412HOChr01$CountPerCodon.MAFSynonymousNodups,
            Col_Deleterious, Col_Synonymous, 20, 0.00050, mean_HighMAFCodon + sd_HighMAFCodon)

GenomeMAFPlots <- lapply(HighMAF_Chrom, function(x) {
  Genome_plot(x, x$Mbp, 
              x$CountPerCodon.MAFAllDel, 
              x$CountPerCodon.MAFSynonymousNodups,
              Col_Deleterious, Col_Synonymous, 20, 0.0005, mean_HighMAFCodon + 2*sd_HighMAFCodon)})


ggarrange(plotlist = GenomeMAFPlots, labels=labels) 







###### SCRATCH
# combine
HighFreq_df <- do.call("rbind", HighFreqData)
HighFreq_df$CountperCodon <- HighFreq_df$Count / HighFreq_df$Num_Codons
hist(HighFreq_df[which(HighFreq_df$Count_type=="DerFreqAllDel"),"CountperCodon"])
hist(HighFreq_df[which(HighFreq_df$Count_type=="DerFreqAllDel"),"Count"])
HighFreq_df[which(HighFreq_df$CountperCodon > 0.00015 & HighFreq_df$Count_type=="DerFreqAllDel"),] # 1 on Chr7
HighFreq_df[which(HighFreq_df$Count > 25 & HighFreq_df$Count_type=="DerFreqAllDel"),] # 1 on Chr 7, 1 on Chr 14
HighFreq_df[which(HighFreq_df$Count > 15 & HighFreq_df$Count_type=="DerFreqAllDel"),]

# split by Chromosome
HighFreq_Chrom <- split(HighFreq_df, HighFreq_df$Chromosome)

# test
test_sub <- subset(HighFreq_Chrom$Ha412HOChr01, Count_type=="DerFreqAllDel")

Genome_plot2(test_sub,
             test_sub$EndPos, test_sub$Count / test_sub$Num_Codons,
            "darkgreen", 0.00007)

