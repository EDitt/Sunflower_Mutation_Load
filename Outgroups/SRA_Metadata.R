setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/HaploblockData")
setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/HaploblockData")

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
SRA_red <- SRA[,c(1,11,12,21,22,25,26,30)]
length(which(SRA_red$BioSample %in% annuus$SRA.ID)) #855 (all)
length(which(annuus$Sample %in% SRA_red$SampleName)) #613 (out of 720)


#### Another set of SRA info
sra2 <- read.csv("SraRunInfo2.csv", header=T)
SRA_red2 <- sra2[,c(1,11,12,21,22,25,26,30)]
length(which(SRA_red2$BioSample %in% annuus$SRA.ID)) #157 (out of 492)
length(which(annuus$Sample %in% SRA_red2$SampleName)) #100 (out of 720)

## combine SRA dataframes
Both_SRA <- rbind(SRA_red, SRA_red2)
length(which(annuus$Sample %in% Both_SRA$SampleName)) #713 (out of 720)

all <- merge(annuus, Both_SRA, by.x = "Sample", by.y="SampleName", all = TRUE)

aggregate(all$Run, by=list(all$Population), length)
aggregate(all$Run, by=list(all$Population, all$Sample), length)


### Select the 4 populations
Pops <- c("ANN_52", "ANN_65", "ANN_70", "ANN_71")
length(which(all$Population %in% Pops)) # N = 72 

PopSubset <- all[which(all$Population %in% Pops),]
aggregate(PopSubset$Run, by=list(PopSubset$Population, PopSubset$Sample), length) #some individuals represented more than once

write.table(PopSubset[,c(1,3,15,16,22)], file="WildAnnuusPops", sep = "\t", quote = FALSE)
# missing 1 ANN_65 individual

################################
######### ONE PER POP ##########
################################

### Also select one individual from each population (except pops that are really near each other)
mean(annuus$Reads) # 83.7 M
min(annuus$Reads)  # 35.7 M
max(annuus$Reads)  # 262 M
MaxReads <- aggregate(annuus$Reads, by=list(annuus$Population), max)

hist(annuus$Reads)

MaxReadsPop <- annuus[which(annuus$Reads %in% MaxReads$x ),] #N = 71

### Subset to one individual per 4 population
OnePerPop <- all[which(all$Reads %in% MaxReads$x),]
aggregate(OnePerPop$Species, by=list(OnePerPop$Sample), length)

write.table(OnePerPop[,c(1,3,15,16,22)], file="WildAnnuusOnePerPop", sep = "\t", quote = FALSE)

### Sequences that don't have to be merged

filesPerPop <- aggregate(OnePerPop$Species, by=list(OnePerPop$Population), length)

OneFile <- subset(filesPerPop, x == 1) #N=20
OneFile$Group.1

#### Files per subset
### I attempted to download all files, then chose a subset (of the ones successfully downloaded) that covered a geographic range
Subset <- read.table("FirstStab_subset")

filesPerPopSubset <- filesPerPop[which(filesPerPop$Group.1 %in% Subset$V1),]
write.table(filesPerPopSubset, file="WildAnnuusSubset1", sep = "\t", quote = FALSE)




