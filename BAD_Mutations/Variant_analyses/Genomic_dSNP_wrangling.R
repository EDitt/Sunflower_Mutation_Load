### Wrangling the number of dSNPs per codon for different frequency classes

source("BAD_Mutations/Variant_analyses/Functions.R")


######################################
######### SET UP DATAFRAME ###########
######################################

GenData <- ImportFilesAsList("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp",
                             ".bed", "", F)
for (i in seq_along(GenData)) {
  name <- names(GenData[i])
  name <- gsub("_10MbpCounts", "", name)
  colnames(GenData[[i]]) <- c("Chromosome", "StartPos", "EndPos", "Count")
  GenData[[i]]["Count_type"] <- name
  #GenData[[i]]["Mbp"] <- (GenData[[i]]["StartPos"] + ((GenData[[i]]["EndPos"] - GenData[[i]]["StartPos"]) / 2))/1000000
  GenData[[i]]["Mbp"] <- GenData[[i]]["StartPos"] /1000000
}

lapply(GenData, function(x) {length(colnames(x))})

colnames(GenData$mRNA_10MbpNumCodonCounts)[5:7] <- c("TotalCodingBases", "TotalmRNA", "Num_Codons")

GenDf <- do.call("rbind", GenData[2:11])

######################################
############ WRANGLING ###############
######################################

# add deleterious vs. synonymous

GenDf$Consequence <- "Fill"
GenDf$Consequence[grep("AllDel",
                       GenDf$Count_type)] <- "Deleterious"
GenDf$Consequence[grep("SynonymousNodups",
                       GenDf$Count_type)] <- "Synonymous"

# check
aggregate(GenDf$EndPos, by=list(GenDf$Count_type, GenDf$Consequence), length)

GenDf$CountCat <- gsub("AllDel", "", GenDf$Count_type)
GenDf$CountCat <- gsub("SynonymousNodups", "", GenDf$CountCat)

GenDf_wide1 <- reshape(GenDf[,c(1,6,4,7,8),],
                       timevar = "CountCat",
                       idvar = c("Chromosome", "Mbp", "Consequence"),
                       direction = "wide")

GenDf_wide1$TotalNum <- GenDf_wide1$Count.Num_HighDerFreq + GenDf_wide1$Count.Num_HighMAF + GenDf_wide1$Count.Num_LowMAF
length(GenDf_wide1[which(GenDf_wide1$TotalNum!=GenDf_wide1$Count.Num_),"Mbp"])
head(GenDf_wide1[which(GenDf_wide1$TotalNum!=GenDf_wide1$Count.Num_),c("TotalNum", "Count.Num_")])
# the difference is the total number includes the ones removed because the ancestral allele != deleterious

# subtract the singletons from the low MAF
GenDf_wide1$Count.Num_LowMAF_minusSingle <- GenDf_wide1$Count.Num_LowMAF - GenDf_wide1$Count.SingleIndiv_Alleles


GenDf_wide2 <- reshape(GenDf_wide1[,c(1:3,5,6,8:10)],
                       timevar = "Consequence",
                       idvar = c("Chromosome", "Mbp"),
                       direction = "wide")

GenDf_wide2$dSNP_sSNP <- GenDf_wide2$TotalNum.Deleterious / GenDf_wide2$TotalNum.Synonymous

# merge with codon info
GenDf_all <- merge(GenDf_wide2, GenData[["mRNA_10MbpNumCodonCounts"]][,c(1,7,9)],
                   by=c("Chromosome", "Mbp"))

# relativize by number of codons
GenDf_Rel <- GenDf_all
GenDf_Rel[3:12] <- lapply(GenDf_Rel[3:12], function(x) {x / GenDf_Rel$Num_Codons})

# reshape back to long (just for deleterious proportional data)
GenDfDel_Rel_long <- reshape(GenDf_Rel[,c(1:5,7)],
                             varying=names(GenDf_Rel[c(3:5,7)]),
                             idvar = c("Chromosome", "Mbp"),
                             times = names(GenDf_Rel[c(3:5,7)]),
                             timevar = c("FrequencyClass"),
                             v.names = "NumPerCodon",
                             direction = "long"
)

#levels(as.factor(GenDfDel_Rel_long$FrequencyClass))
#GenDfDel_Rel_long$FrequencyClass <- factor(GenDfDel_Rel_long$FrequencyClass,
#                                           levels=c("Count.Num_HighDerFreq.Deleterious", "Count.Num_HighMAF.Deleterious",
#                                                    "Count.Num_LowMAF_minusSingle.Deleterious", "Count.SingleIndiv_Alleles.Deleterious"))


######################################
######### SAVE AS A DATAFRAME ########
######################################

write.table(GenDfDel_Rel_long, 
            file = "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/dSNPperCodonNums.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

