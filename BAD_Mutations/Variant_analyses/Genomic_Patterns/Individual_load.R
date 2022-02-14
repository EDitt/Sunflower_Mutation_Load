### graphing load per individual and performance (redo 'Genotype_patterns.md')

library(reshape2)
library(car)
library(ggplot2)

######################################
########## SET UP DATAFRAME ##########
######################################

setwd("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns")

dSNP_counts <- read.table("All_dSNP_stats.txt", sep = "\t", header=TRUE,
                            stringsAsFactors = FALSE)
length(dSNP_counts$sample) #288

key <- read.csv("LineKeywINFO.csv", header=T)
length(key$sort) # 290

counts_info <- merge(key, dSNP_counts, by.x="SequenceName", by.y="sample")
length(counts_info$sort) #288

# relativize by called genotypes
counts_info$Rel_nDel <- counts_info$nDeleterious_Hom / counts_info$CalledGenotypes
hist(counts_info$Rel_nDel)

######################################
########## PLOT DISTRIBUTION #########
######################################

hist(counts_info$nDeleterious_Hom)

boxplot(counts_info$nDeleterious_Hom ~ counts_info$Mandel,
        xlab="Mandel labels",
        ylab="Relative number of homozygous deleterious genotypes")

######################################
# RELATIONSHIP W/ PERFORMANCE:HELNUT #
######################################

performance <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/Nut_Pilot/data files/HELNUT16SAMlsmeans_wide.csv", header=T)
length(performance$SAM) # 263

AllData <- merge(performance[,c(2,43,44,53,54,59,60)], counts_info, by.y="SamID", by.x="SAM")
length(AllData$SAM) # 262

# stress / control
AllData$StressControl <- AllData$RGR_Low / AllData$RGR_High



plot(AllData$RGR_High ~ AllData$Rel_nDel)
plot(AllData$RGR_Low ~ AllData$Rel_nDel)
plot(AllData$StressControl ~ AllData$Rel_nDel)

Data_long <- melt(AllData[,c(1:7,17,18,31)], id.vars = c("SAM", "heterotic_group",
                                                         "Oil_NonOil", "Rel_nDel") )
Data_long$variable <- strsplit(as.character(Data_long$variable), "_")
Data_long$PerfMeasure <- as.factor(sapply(Data_long$variable, "[", 1))
Data_long$Condition <- as.factor(sapply(Data_long$variable, "[", 2))

Mod1 <- lm(value ~ Rel_nDel + Condition +
             heterotic_group +
             Rel_nDel:Condition,
           data=Data_long[which(Data_long$PerfMeasure=="RGR"),])
summary(Mod1)
drop1(Mod1, test="Chi") # relative # del x condition p=0.02097
Anova(Mod1) # no main effect of nDel

p <- ggplot(data=Data_long[which(Data_long$PerfMeasure=="RGR"),],
            aes(x=Rel_nDel, y=value, group=Condition))
p + geom_point(aes(color=Condition)) +
  geom_smooth(method='lm', aes(color=Condition)) +
  theme_minimal() +
  ylab("Relative growth rate")

p2 <- ggplot(data=Data_long[which(Data_long$PerfMeasure=="ln.Plant.weight"),],
             aes(x=Rel_nDel, y=value, group=Condition))
p2 + geom_point(aes(color=Condition)) +
  geom_smooth(method='lm', aes(color=Condition)) +
  theme_minimal() +
  ylab("Plant_weight")


######################################
# RELATIONSHIP W/ PERFORMANCE:SALTYHEL #
######################################

performance2 <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/Nut_Pilot/data files/SALTYHEL17_reduced_020518.csv", header=T)
length(performance2$SAM) #263

AllData2 <- merge(performance2[,c(2:6)], counts_info, by.y="SamID", by.x="SAM")
length(AllData2$SAM) # 235

plot(AllData2$ln.Plant.weight_water ~ AllData2$Rel_nDel)
plot(AllData2$ln.Plant.weight_salt ~ AllData2$Rel_nDel)

Data_long2 <- melt(AllData2[,c(1:5,15,16,28)], id.vars = c("SAM", "heterotic_group", "Oil_NonOil", "Rel_nDel") )
Data_long2$variable <- strsplit(as.character(Data_long2$variable), "_")
Data_long2$PerfMeasure <- as.factor(sapply(Data_long2$variable, "[", 1))
Data_long2$Condition <- as.factor(sapply(Data_long2$variable, "[", 2))

Mod_salt <- lm(value ~ Rel_nDel + Condition +
                 heterotic_group +
                 Rel_nDel:Condition,
               data=Data_long2[which(Data_long2$PerfMeasure=="ln.Plant.weight"),])
summary(Mod_salt)
drop1(Mod_salt, test="Chi") # interaction n.s.
Anova(Mod_salt) 

q <- ggplot(data=Data_long2[which(Data_long2$PerfMeasure=="ln.Plant.weight"),],
            aes(x=Rel_nDel, y=value, group=Condition))
q + geom_point(aes(color=Condition)) +
  geom_smooth(method='lm', aes(color=Condition)) +
  theme_minimal() +
  ylab("ln Plant weight")
