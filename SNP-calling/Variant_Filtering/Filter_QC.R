#setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/Variant_graphs")
setwd("/Users/eld72413/Documents/GitHub/Sunflower_Mutation_Load/SNP-calling")
QUAL <- read.table("SAM_SNPs_GTfilteredStatsQUAL.txt", header = FALSE, sep = "\t")

colnames(QUAL)[3:6] <- c("Quality", "Num_SNPs", "Num_ts", "Num_tv")

QUAL$tstv <- QUAL$Num_ts / QUAL$Num_tv

hist(QUAL$tstv)

plot(QUAL$tstv ~ QUAL$Quality, cex=0.5)
abline(v=40, lty=2)
abline(h=1.71)


abline(h=mean(QUAL$tstv))
abline(v=69)
abline(h=mean(QUAL[which(QUAL$Quality>=40),"tstv"]), lty=2)

abline(v=50, lty=3)
abline(h=1.79, lty=2)

plot(QUAL[which(QUAL$Quality >=40),"tstv"] ~ QUAL[which(QUAL$Quality >=40),"Quality"])
plot(QUAL[which(QUAL$Quality >=50),"tstv"] ~ QUAL[which(QUAL$Quality >=50),"Quality"])
plot(QUAL[which(QUAL$Quality >=60),"tstv"] ~ QUAL[which(QUAL$Quality >=60),"Quality"])
plot(QUAL[which(QUAL$Quality >=70),"tstv"] ~ QUAL[which(QUAL$Quality >=70),"Quality"])

sd(QUAL$tstv)
sqrt(length(QUAL$tstv))

sd(QUAL$tstv) / sqrt(length(QUAL$tstv))

mean(QUAL$tstv) - sd(QUAL$tstv)

QUAL_sub <- subset(QUAL, tstv>=1.74)
min(QUAL_sub$Quality) #69

QUAL[which(QUAL$Quality==40),] #ts/tv=1.65
QUAL[which(QUAL$Quality==30),] #ts/tv=1.58

mean(QUAL[which(QUAL$Quality>=40),"tstv"]) #1.79
mean(QUAL[which(QUAL$Quality<40),"tstv"]) 
mean(QUAL[which(QUAL$Quality>=70),"tstv"]) #1.80


##### heterozygous individuals

het <- read.table("Variant_Filtering/SAM_hetIND_GTfiltered.het", header = TRUE, sep = "\t")

hist(het$F) # F is the inbreeding coefficient
# after 2 rounds of selfing, F expected to be 0.75
het[which(het$F < 0),]
het[which(het$F < 0.25),]
het[which(het$F < 0.1),]

het$Num_het <- het$N_SITES - het$O.HOM.
hist(het$Num_het)
