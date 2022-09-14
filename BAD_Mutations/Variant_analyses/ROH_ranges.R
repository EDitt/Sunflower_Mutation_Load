### code to plot ROH hotspots across germplasm groups

source("BAD_Mutations/Variant_analyses/Functions.R")

######################################
######### SET UP DATAFRAME ###########
######################################

GroupROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/ROH_GroupHotpots.txt",
                       header=F, 
                       col.names=c("Chromosome", "StartPos", "EndPos",
                                   "Group", "Num_Individuals"))

GroupROH$Mbp <- round((GroupROH$StartPos / 1000000) + 1, digits = 0)



######################################
### UPPER/LOWER QUANTILES PER GROUP ###
######################################

Group_ROHList <- split(GroupROH, GroupROH$Group)

# maximum
lapply(Group_ROHList, function(x) {max(x$Num_Individuals)})
# HA-Nonoil-50, HA-Oil-77, introgressed-15, landrace/OPV-11, RHA-NonOil-19, RHA-Oil-77

# total number of individuals per group
aggregate(key$sort, by=list(key$Groups), length)
# HA-NonOil- 55, HA-Oil-90, introgressed-18, landrace/OPV-13, RHA-NonOil-21, RHA-Oil-91

NumInd <- list("HA-NonOil" = 55,
               "HA-Oil" = 90,
               "introgressed" = 18,
               "landrace_OPV" = 13,
               "RHA-NonOil" = 21,
               "RHA-Oil" = 91)
# 2 HA-Oils aren't used in dataset but are in spreadsheet

# proportion
Group_ROHList <- lapply(names(Group_ROHList), function(x) {
  Group_ROHList[[x]]["Proportion"] <- Group_ROHList[[x]]["Num_Individuals"] / NumInd[[x]];
  return(Group_ROHList[[x]])
})

# add column for upper and lower quantiles

Group_ROHList <- lapply(Group_ROHList, function(x) {
  AddUpperLowerQuantile(x, "Num_Individuals", 0.975, 0.025, "97.5%", "2.5%")
}
)

# recombine lists
Group_ROH_new <- do.call("rbind", Group_ROHList)

# to visualize
ggplot(Group_ROH_new[which(Group_ROH_new$Group!="introgressed"),], 
       aes(x=Mbp, y=Proportion)) +
  geom_smooth(aes(color=Group), method="loess") +
  geom_point(data=Group_ROH_new[which(Group_ROH_new$PointCol=="97.5%" &
                                        Group_ROH_new$Group!="introgressed"),],
             aes(x=Mbp, y=Proportion, color=Group), alpha=0.6) +
  facet_wrap(~Chromosome, 
             scales = "free_x")

ggplot(Group_ROH_new[which(Group_ROH_new$Group!="introgressed"),], 
       aes(x=Mbp, y=Proportion)) +
  geom_smooth(aes(color=Group), method="loess") +
  geom_point(data=Group_ROH_new[which(Group_ROH_new$PointCol!="none" &
                                        Group_ROH_new$Group!="introgressed"),],
             aes(x=Mbp, y=Proportion, color=Group), alpha=0.6) +
  facet_wrap(~Chromosome, 
             scales = "free_x")


######################################
########## SHARED HOTSPOTS ###########
######################################

aggregate(Group_ROH_new$Chromosome, by=list(Group_ROH_new$Group,
                                            Group_ROH_new$PointCol), length)

# long to wide to get hotspots shared across groups
Group_ROH_wide <- reshape(Group_ROH_new[,c(1:4,8)],
                          idvar=c("Chromosome", "StartPos", "EndPos"),
                          timevar = "Group",
                          direction = "wide")

Group_ROH_wide_hotspot <- subset(Group_ROH_wide, `PointCol.HA-NonOil`=="97.5%" &
                                   `PointCol.HA-Oil`=="97.5%" &
                                   `PointCol.RHA-NonOil`=="97.5%" &
                                   `PointCol.RHA-Oil`=="97.5%")

Group_ROH_wide_hotspot$Category <- ifelse(Group_ROH_wide_hotspot$PointCol.landrace_OPV=="97.5%",
                                          "all", "elite")

aggregate(Group_ROH_wide_hotspot$Chromosome, by=list(Group_ROH_wide_hotspot$Category), length)
Group_ROH_wide_hotspot$Chromosome <- factor(Group_ROH_wide_hotspot$Chromosome) # reduce levels

# order
Group_ROH_wide_hotspot <- Group_ROH_wide_hotspot[order(Group_ROH_wide_hotspot$Chromosome,
                                                         Group_ROH_wide_hotspot$StartPos),]
# make a list with ranges:
hotspot_split <- split(Group_ROH_wide_hotspot, Group_ROH_wide_hotspot$Chromosome) # seven chromosomes
# 3,4, ie chromosome 9,12 has mixed elite/all

hotspot_range_list <- lapply(hotspot_split, function(x) {MakeROHSegments(x)})
hotspot_range <- do.call("rbind", hotspot_range_list)

write.table(hotspot_range, "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Shared_hotspot_range.txt",
            quote=FALSE, sep = "\t", row.names=FALSE)

hotspot_range$length <- hotspot_range$EndPos - hotspot_range$StartPos
hist(hotspot_range$length)
hotspot_range[which(hotspot_range$length==max(hotspot_range$length)),] # chromosome 12: 88,000,000 92,750,000

############################################
# MAKE CONTINUOUS SEGMENTS for OUTSIDE ROH #
###########################################

OutsideROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/NonHotspot_range.bed")
colnames(OutsideROH) <- c("Chromosome", "StartPos", "EndPos")

# order
OutsideROH <- OutsideROH[order(OutsideROH$Chromosome,
                               OutsideROH$StartPos),]
# make a list with ranges:
OutsideROH_split <- split(OutsideROH, OutsideROH$Chromosome)

OutsideROH_split_list <- lapply(OutsideROH_split, function(x) {MakeROHSegments(x)})
OutsideROH_range <- do.call("rbind", OutsideROH_split_list)

write.table(OutsideROH_range, "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Outside_hotspot_range.txt",
            quote=FALSE, sep = "\t", row.names=FALSE)

######################################
############# HA vs. RHA #############
######################################

HetGrp_ROH <- CombineGroups_ROH(Group_ROHList, 1, 2, 5, 6,
                                "HA", "RHA", 145, 112)

HetGrp_ROH_wide <- reshape(HetGrp_ROH[,c(1:3,6,8)],
                          idvar=c("Chromosome", "StartPos", "EndPos"),
                          timevar = "Group",
                          direction = "wide")
HetGrp_ROH_wide$Category <- ifelse(HetGrp_ROH_wide$PointCol.HA=="97.5%" &
                                     HetGrp_ROH_wide$PointCol.RHA=="97.5%",
                                   "shared", ifelse(HetGrp_ROH_wide$PointCol.HA=="97.5%" &
                                                      HetGrp_ROH_wide$PointCol.RHA!="97.5%",
                                                    "HA_hotspot",
                                                    ifelse(HetGrp_ROH_wide$PointCol.HA!="97.5%" &
                                                             HetGrp_ROH_wide$PointCol.RHA=="97.5%",
                                                           "RHA_hotspot", "none")))

aggregate(HetGrp_ROH_wide$Chromosome, by=list(HetGrp_ROH_wide$Category), length)
HetGrp_ROH_wide[which(HetGrp_ROH_wide$Category=="shared"),] # more than in the previous 'shared'

Het_Grp_Diff <- subset(HetGrp_ROH_wide, Category=="HA_hotspot" | Category=="RHA_hotspot")

# order
Het_Grp_Diff <- Het_Grp_Diff[order(Het_Grp_Diff$Chromosome,
                                   Het_Grp_Diff$StartPos),]
Het_Grp_Diff$Chromosome <- factor(Het_Grp_Diff$Chromosome)

# make a list with ranges:
Het_split <- split(Het_Grp_Diff, list(Het_Grp_Diff$Category,
                                   Het_Grp_Diff$Chromosome), drop=TRUE)

Het_range_list <- lapply(Het_split, function(x) {MakeROHSegments(x)})

RHA_lists <- names(Het_range_list[grep("RHA", names(Het_range_list))])
RHA_ranges <- do.call("rbind", Het_range_list[RHA_lists])
RHA_ranges$Het_Group <- "RHA"

HA_lists <- names(Het_range_list[grep("^HA", names(Het_range_list))])
HA_ranges <- do.call("rbind", Het_range_list[HA_lists])
HA_ranges$Het_Group <- "HA"

Het_ranges <- rbind(HA_ranges, RHA_ranges)

write.table(Het_ranges, "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Het_group_ranges.txt",
            quote=FALSE, sep = "\t", row.names=FALSE)

Het_ranges$MbpStart <- round((Het_ranges$StartPos / 1000000) + 1, digits = 0)
Het_ranges$MbpEnd <- round((Het_ranges$EndPos / 1000000) + 1, digits = 0)
Het_ranges$length <- Het_ranges$EndPos - Het_ranges$StartPos
Het_ranges[which(Het_ranges$length==max(Het_ranges$length)),] # chromosome 10: 50,250,000 - 54,000,000
hist(Het_ranges$length)

Het_ranges[which(Het_ranges$length>2000000),] # chromosome 10 and chromosome 4
Het_ranges[which(Het_ranges$length>1500000),] # chromosome 9

# to visualize
ggplot(HetGrp_ROH, 
       aes(x=Mbp, y=Proportion)) +
  geom_smooth(aes(color=Group), method="loess") +
  geom_point(data=HetGrp_ROH[which(HetGrp_ROH$PointCol=="97.5%"),],
             aes(x=Mbp, y=Proportion, color=Group), alpha=0.5) +
  facet_wrap(~Chromosome, 
             scales = "free_x") +
  geom_segment(data=Het_ranges, aes(x=MbpStart, y=0.5, xend=MbpEnd, yend=0.5, color=Het_Group)) +
  theme_minimal()

######################################
########### OIL vs. NONOIL ###########
######################################

OilGrp_ROH <- CombineGroups_ROH(Group_ROHList, 2, 6, 1, 5,
                                "Oil", "Non-Oil", 181, 76)

OilGrp_ROH_wide <- reshape(OilGrp_ROH[,c(1:3,6,8)],
                           idvar=c("Chromosome", "StartPos", "EndPos"),
                           timevar = "Group",
                           direction = "wide")

OilGrp_ROH_wide$Category <- ifelse(OilGrp_ROH_wide$PointCol.Oil=="97.5%" &
                                     OilGrp_ROH_wide$`PointCol.Non-Oil`=="97.5%",
                                   "shared", ifelse(OilGrp_ROH_wide$PointCol.Oil=="97.5%" &
                                                      OilGrp_ROH_wide$`PointCol.Non-Oil`!="97.5%",
                                                    "Oil_hotspot",
                                                    ifelse(OilGrp_ROH_wide$PointCol.Oil!="97.5%" &
                                                             OilGrp_ROH_wide$`PointCol.Non-Oil`=="97.5%",
                                                           "NonOil_hotspot", "none")))

aggregate(OilGrp_ROH_wide$Chromosome, by=list(OilGrp_ROH_wide$Category), length)

Oil_Grp_Diff <- subset(OilGrp_ROH_wide, Category=="Oil_hotspot" | Category=="NonOil_hotspot")

# order
Oil_Grp_Diff <- Oil_Grp_Diff[order(Oil_Grp_Diff$Chromosome,
                                   Oil_Grp_Diff$StartPos),]
Oil_Grp_Diff$Chromosome <- factor(Oil_Grp_Diff$Chromosome)

# make a list with ranges:
Oil_split <- split(Oil_Grp_Diff, list(Oil_Grp_Diff$Category,
                                      Oil_Grp_Diff$Chromosome), drop=TRUE)

Oil_range_list <- lapply(Oil_split, function(x) {MakeROHSegments(x)})

NonOil_lists <- names(Oil_range_list[grep("NonOil", names(Oil_range_list))])
NonOil_ranges <- do.call("rbind", Oil_range_list[NonOil_lists])
NonOil_ranges$Oil_Group <- "NonOil"

Oil_lists <- names(Oil_range_list[grep("^Oil", names(Oil_range_list))])
Oil_ranges <- do.call("rbind", Oil_range_list[Oil_lists])
Oil_ranges$Oil_Group <- "Oil"

Oil_group_ranges <- rbind(Oil_ranges, NonOil_ranges)

write.table(Oil_group_ranges, "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Oil_group_ranges.txt",
            quote=FALSE, sep = "\t", row.names=FALSE)

Oil_group_ranges$MbpStart <- round((Oil_group_ranges $StartPos / 1000000) + 1, digits = 0)
Oil_group_ranges$MbpEnd <- round((Oil_group_ranges$EndPos / 1000000) + 1, digits = 0)
Oil_group_ranges$length <- Oil_group_ranges$EndPos - Oil_group_ranges$StartPos
Oil_group_ranges[which(Oil_group_ranges$length==max(Oil_group_ranges$length)),] # chromosome 8: 76,000,000-77,500,000
hist(Oil_group_ranges$length)

Oil_group_ranges[which(Oil_group_ranges$length>1200000),] # mostly all oil hotspots: Chr 2, Chr8, Chr10, Chr11, Chr16, Chr17
# one NonOil hotspot: Chr2

ggplot(OilGrp_ROH, 
       aes(x=Mbp, y=Proportion)) +
  geom_smooth(aes(color=Group), method="loess") +
  geom_point(data=OilGrp_ROH[which(OilGrp_ROH$PointCol=="97.5%"),],
             aes(x=Mbp, y=Proportion, color=Group), alpha=0.5) +
  facet_wrap(~Chromosome, 
             scales = "free_x") +
  geom_segment(data=Oil_group_ranges, aes(x=MbpStart, y=0.5, xend=MbpEnd, yend=0.5, color=Oil_Group)) +
  theme_minimal()

######################################
############# FREQ DIST #############
######################################



#################################
########### FUNCTIONS ###########
#################################

AddUpperLowerQuantile <- function(dataframe, ColName, Upper_quantile, Lower_quantile,
                                  UpperName, LowerName){
  Upperquantile <- quantile(dataframe[,ColName], Upper_quantile, na.rm = T)
  Lowerquantile <- quantile(dataframe[,ColName], Lower_quantile, na.rm = T)
  dataframe$PointCol <- ifelse(dataframe[,ColName]>Upperquantile,
                               UpperName,
                               ifelse(dataframe[,ColName]<Lowerquantile,
                                      LowerName, 
                                      "none"))
  return(dataframe)
}

# function to combine HA/RHA or Oil/NonOil
CombineGroups_ROH <- function(List, Group1a, Group1b, Group2a, Group2b,
                              Group1_name, Group2_name,
                              Group1_NumberInd, Group2_NumberInd) {
  Group1 <- merge(List[[Group1a]][,c(1:6)],
                  List[[Group1b]][,c(1:6)],
                  by=c("Chromosome", "StartPos", "EndPos", "Mbp"))
  Group1$TotalNum <- Group1$Num_Individuals.x + Group1$Num_Individuals.y
  Group1 <- AddUpperLowerQuantile(Group1, "TotalNum", 0.975, 0.025, "97.5%", "2.5%")
  Group1$Proportion <- Group1$TotalNum / Group1_NumberInd
  Group1$Group <- Group1_name
  Group2 <- merge(List[[Group2a]][,c(1:6)],
                  List[[Group2b]][,c(1:6)],
                  by=c("Chromosome", "StartPos", "EndPos", "Mbp"))
  Group2$TotalNum <- Group2$Num_Individuals.x + Group2$Num_Individuals.y
  Group2 <- AddUpperLowerQuantile(Group2, "TotalNum", 0.975, 0.025, "97.5%", "2.5%")
  Group2$Proportion <- Group2$TotalNum / Group2_NumberInd
  Group2$Group <- Group2_name
  Combined <- rbind(Group1[,c("Chromosome", "StartPos", "EndPos", "Mbp",
                              "TotalNum", "PointCol", "Proportion", "Group")],
                    Group2[,c("Chromosome", "StartPos", "EndPos", "Mbp",
                              "TotalNum", "PointCol", "Proportion", "Group")])
  return(Combined)
}

## function to make ranges
MakeROHSegments <- function(dataframe) {
  dataframe <- dataframe[order(dataframe$Chromosome, dataframe$StartPos),]
  dataframe$PrevPos <- "fill"
  dataframe$NextPos <- "fill"
  dataframe$PrevPos[1] <- "start"
  if(length(dataframe$Chromosome)==1){
    return(dataframe[,c(1:3)])} else{
      dataframe$PrevPos[2:length(dataframe$PrevPos)] <- dataframe$EndPos[1:length(dataframe$EndPos)-1]
      dataframe$NextPos[length(dataframe$NextPos)] <- "end"
    dataframe$NextPos[1:length(dataframe$Next)-1] <- dataframe$StartPos[2:length(dataframe$StartPos)] 
    dataframe$section <- ifelse(dataframe$PrevPos==dataframe$StartPos &
                                 dataframe$NextPos!=dataframe$EndPos,
                                "end",
                                ifelse(dataframe$PrevPos!=dataframe$StartPos &
                                        dataframe$NextPos==dataframe$EndPos,
                                     "start",
                                     ifelse(dataframe$PrevPos!=dataframe$StartPos &
                                              dataframe$NextPos!=dataframe$EndPos,
                                            "single", "middle")))
    dataframe_start <- subset(dataframe, section=="start")
    dataframe_end <- subset(dataframe, section=="end")
    dataframe_segments1 <- cbind(dataframe_start[,c(1:2)],
                               dataframe_end[,c(3)])
    colnames(dataframe_segments1)[3] <- "EndPos"
    dataframe_segments2 <- subset(dataframe, section=="single")
    dataframe_segments <- rbind(dataframe_segments1, 
                              dataframe_segments2[,c(1:3)])
  return(dataframe_segments)
    }
}

  
### SCRATCH
test1 <- hotspot_split$Ha412HOChr09

test1$PrevPos <- "fill"
test1$NextPos <- "fill"

test1$PrevPos[1] <- "start"
test1$PrevPos[2:length(test1$PrevPos)] <- test1$EndPos[1:length(test1$EndPos)-1]
test1$NextPos[length(test1$NextPos)] <- "end"

test1$NextPos[1:length(test1$Next)-1] <- test1$StartPos[2:length(test1$StartPos)] 

test1$section <- ifelse(test1$PrevPos==test1$StartPos &
                           test1$NextPos!=test1$EndPos,
                        "end",
                        ifelse(test1$PrevPos!=test1$StartPos &
                                  test1$NextPos==test1$EndPos,
                               "start",
                               ifelse(test1$PrevPos!=test1$StartPos &
                                        test1$NextPos!=test1$EndPos,
                                      "single", "middle")))

test1_start <- subset(test1, section=="start")
test1_end <- subset(test1, section=="end")
test1_segments1 <- cbind(test1_start[,c(1:2)],
                         test1_end[,c(3)])
colnames(test1_segments1)[3] <- "EndPos"
test1_segments2 <- subset(test1, section=="single")
test1_segments <- rbind(test1_segments1, 
                            test1_segments2[,c(1:3)])

test1 <- MakeROHSegments(test1)


test2 <- hotspot_split[[5]]

test2 <- MakeROHSegments(test2)
test3 <- MakeROHSegments(hotspot_split[[3]])
test4 <- MakeROHSegments(hotspot_split[[4]])

test4[,c(1:3,13)]
test4_subset_start <- subset(test4, section=="start")
test4_subset_end <- subset(test4, section=="end")

test4_segments1 <- cbind(test4_subset_start[,c(1:2)],
                         test4_subset_end[,c(3)])
colnames(test4_segments1)[3] <- "EndPos"
test4_segments <- rbind(test4_segments1, 
                        subset(test4, section=="single"))

test5 <- MakeROHSegments(hotspot_split[[4]])

test4_subset1_wide <- reshape(test4_subset1[,c(1:3,13)],
                              idvar = "Chromosome",
                              timevar = "section",
                              direction = "wide"
)
MakeROHSegments(hotspot_split[[3]])