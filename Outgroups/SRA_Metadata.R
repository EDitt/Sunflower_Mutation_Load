setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/HaploblockData")

################################
######### WILD ANNUUS ##########
################################

annuus <- read.csv("WildAnnuus.csv", header=T)
length(annuus$Species) #719
aggregate(annuus$Sample, by=list(annuus$Population), length)
# most populations are represented by 10 individuals, 4 have 11 (ANN_62, ANN_63, ANN_64, ANN71 - Iowa and Kansas, Saskatchewan), 1 has 12 (ANN_70 - South Dakota), 1 has 14 (ANN_65 - Kansas), 1 has 9 (ANN_69 - South Dakota)

# whole populations to download

# ANN_71: Saskatchewan (N=11)
# ANN_70: South Dakota (N=12)
# ANN_65: Kansas (N=14)
# ANN_52: Texas (N=10)

### Combine with SRA Info
SRA <- read.csv("SraRunInfo.csv", header=T)
length(SRA$Run)
aggregate(SRA$Run, by=list(SRA$Submission), length)
# 562, 289, 4

IDs <- annuus$Sample
length(which(SRA$SampleName %in% IDs)) #855
# some individuals represented more than once?
aggregate(SRA$Run, by=list(SRA$SampleName), length)

all <- merge(annuus, SRA, by.x = "Sample", by.y="SampleName", all = TRUE)
str(all$SRA.ID)
str(SRA$BioSample)

aggregate(all$Sample, by=list(all$Population), length)

### Select the 4 populations
Pops <- c("ANN_52", "ANN_65", "ANN_70", "ANN_71")
length(which(all$Population %in% Pops)) # N = 50 (3 more than expected)

PopSubset <- all[which(all$Population %in% Pops),]
length(PopSubset$Sample) #N=50
write.table(PopSubset[,c(1,3,15,16,57)], file="WildAnnuusPops", sep = "\t", quote = FALSE)

### Also select one individual from each population (except pops that are really near each other)
mean(annuus$Reads) # 83.7 M
min(annuus$Reads)  # 35.7 M
max(annuus$Reads)  # 262 M
aggregate(annuus$Reads, by=list(annuus$Population), max)