# Site Frequency Spectrum for Groups of Germplasm

source("BAD_Mutations/Variant_analyses/Functions.R")

################################
######### SET UP DATA ##########
################################

DirPath=c("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/Group_Freq_bins")
my_files <- list.files(path = DirPath, 
                       pattern = "*.txt", full.names = TRUE)
my_data <- lapply(my_files, read.table, header=T)
names(my_data) <- gsub("\\.txt$", "", my_files)
names(my_data) <- gsub(DirPath, "", names(my_data))
names(my_data) <- gsub("/", "", names(my_data))

levels(my_data$`HA-NonOil_DerivedFreq_Bins`$Annotation)

my_data <- lapply(my_data, function(x) {x$Annotation <- factor(x$Annotation, 
                                  levels=c("NonCodingNodups", "SynonymousNodups",
                                         "Tolerated", "MissenseOther", "AllDel", 
                                          "StartStopLostGainednoDups", "SpliceAcceptorDonorNodups"));
return(x)})


#############################
## PLOT SFS across groups ###
#############################

PlotSFS <- function(df, Annotations, AnnotationCols, AnnotationLabs) {
  df_subset <- subset(df, Annotation %in% Annotations)
  plot <- ggplot(df_subset, aes(x=breaks, y=prop, fill=Annotation)) +
  geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
    scale_fill_manual(values= AnnotationCols,
                      name = "Variant Class", labels = AnnotationLabs) +
    ylab("Proportion") +
    xlab("Frequency") +
    theme(panel.grid.minor = element_blank()) +
    ylim(0,1)
  return(plot)
}


PlotSFS(my_data$`HA-NonOil_DerivedFreq_Bins`, c("SynonymousNodups", "AllDel"),
        c(Col_Synonymous, Col_Deleterious), c("Synonymous", "Deleterious")) # test

Allplots <- lapply(my_data, function(x) {PlotSFS(x, 
                                                 c("SynonymousNodups", "AllDel"),
                                                 c(Col_Synonymous, Col_Deleterious), 
                                                 c("Synonymous", "Deleterious"))})

labels <- gsub("_DerivedFreq_Bins", "", names(Allplots))

ggarrange(plotlist = Allplots, labels=labels) 

##############################
#### COMBINE HA/RHA dSNPS ####
##############################

### merge HA/RHA
Combine_Freqs <- function(dataframe1, dataframe2, annotation, group_name){
  GroupDf <- merge(dataframe1[which(dataframe1$Annotation==annotation),],
                   dataframe2[which(dataframe2$Annotation==annotation),],
        by=c("breaks", "Annotation"))
  GroupDf$count <- GroupDf$count.x + GroupDf$count.y
  GroupDf$prop <- GroupDf$count / sum(GroupDf$count)
  GroupDf$group <- group_name
  return(GroupDf)
}

NonOil_bins <- Combine_Freqs(my_data$`HA-NonOil_DerivedFreq_Bins`,
                                  my_data$`RHA-NonOil_DerivedFreq_Bins`,
                                  "AllDel", "NonOil")
Oil_bins <- Combine_Freqs(my_data$`HA-Oil_DerivedFreq_Bins`,
                          my_data$`RHA-Oil_DerivedFreq_Bins`,
                          "AllDel", "Oil")
my_data$landrace_OPV_DerivedFreq_Bins$group <- "Landrace_OPV"

GroupData <- rbind(my_data$landrace_OPV_DerivedFreq_Bins[which(my_data$landrace_OPV_DerivedFreq_Bins$Annotation=="AllDel"),],
                   rbind(Oil_bins[,c(1,8,7,2,9)],
                         NonOil_bins[,c(1,8,7,2,9)]))

#############################
######## PLOT dSNPS #########
#############################

# 2nd y-axis for lowest freq class to make plot more readable
GroupData$prop2 <- ifelse(GroupData$breaks!=0.05, GroupData$prop * 10,
                          GroupData$prop)

GroupFreqPlot <- function(dataframe, bar_width) {
  ggplot(data=dataframe, aes(x=breaks, y=prop)) +
    geom_bar(aes(fill=group), 
             stat="identity", position = position_dodge(width=bar_width)) +
    scale_fill_manual(values= c(landrace_OPV, HA_NonOil, HA_Oil),
                        name = "Germplasm_Group", 
                      labels = c("Landrace/OPV", "NonOil", "Oil")) +
      ylab("Proportion") +
      xlab("Frequency") +
      theme(panel.grid.minor = element_blank()) +
    #geom_vline(xintercept = 0.075, linetype="dashed") +
    ylim(0,0.1) +
    #scale_y_continuous(
    #  name="Proportion",
    #  sec.axis = sec_axis(trans = ~./10, name="Proportion")
    #) +
    theme_minimal()
}

GroupFreqPlot(GroupData, 0.04)

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/dSNP_SFS_1.pdf")
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/dSNP_SFS_2.pdf")


pdf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/dSNP_SFS.pdf")
print(GroupFreqPlot(GroupData, 0.04))
dev.off()

##############################
#### COMBINE HA/RHA sSNPS ####
##############################

NonOil_bins_sSNP <- Combine_Freqs(my_data$`HA-NonOil_DerivedFreq_Bins`,
                             my_data$`RHA-NonOil_DerivedFreq_Bins`,
                             "SynonymousNodups", "NonOil")
Oil_bins_sSNP <- Combine_Freqs(my_data$`HA-Oil_DerivedFreq_Bins`,
                          my_data$`RHA-Oil_DerivedFreq_Bins`,
                          "SynonymousNodups", "Oil")

GroupData_sSNP <- rbind(my_data$landrace_OPV_DerivedFreq_Bins[which(my_data$landrace_OPV_DerivedFreq_Bins$Annotation=="SynonymousNodups"),],
                   rbind(Oil_bins_sSNP[,c(1,8,7,2,9)],
                         NonOil_bins_sSNP[,c(1,8,7,2,9)]))

ggplot(data=GroupData_sSNP, aes(x=breaks, y=prop)) +
  geom_bar(aes(fill=group), stat="identity", position = position_dodge()) +
  scale_fill_manual(values= c(landrace_OPV, HA_NonOil, HA_Oil),
                    name = "Germplasm_Group", 
                    labels = c("Landrace/OPV", "NonOil", "Oil")) +
  ylab("Proportion") +
  xlab("Frequency") +
  theme(panel.grid.minor = element_blank()) +
  #ylim(0,0.15) +
  theme_minimal()

##################################
# dSNP/sSNPS ACROSS FREQ CLASSES #
##################################

dSNP_sSNP_bins <- merge(GroupData_sSNP,
                        GroupData,
                        by=c("breaks", "group"))

dSNP_sSNP_bins$ratio <- dSNP_sSNP_bins$count.y / dSNP_sSNP_bins$count.x

dSNP_sSNP_ratioPlot <-
  ggplot(data=dSNP_sSNP_bins, aes(x=breaks, y=ratio)) +
    geom_bar(aes(fill=group), stat="identity", position = position_dodge()) +
    scale_fill_manual(values= c(landrace_OPV, HA_NonOil, HA_Oil),
                     name = "Germplasm_Group", 
                     labels = c("Landrace/OPV", "NonOil", "Oil")) +
    xlab("Frequency") +
    ylab("Number dSNP/sSNP") +
    theme(panel.grid.minor = element_blank()) +
    theme_minimal()

pdf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/dSNP_sSNP_GroupFreqs.pdf")
print(dSNP_sSNP_ratioPlot)
dev.off()
