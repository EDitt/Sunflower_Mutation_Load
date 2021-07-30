### combining sequencing reads + coverage info into a table

ImportData <- function (DirPath, suffix_pattern) {
  my_files <- list.files(path = DirPath, pattern = suffix_pattern, full.names = TRUE)
  my_data <- lapply(my_files, function(x) {read.table(x, header=TRUE, sep = "\t")})
  names(my_data) <- gsub("\\.txt$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  return(my_data)
}

#my_files <- list.files(path = "/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Manuscript/coverage_data",
 #                      pattern = "*.txt", full.names = TRUE)
#my_data <- lapply(my_files, function(x) {read.table(x, header=TRUE, sep = "\t")})

Cov <- ImportData("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Manuscript/coverage_data/",
                  "*.txt")
str(Cov)

All_Cov <- do.call("rbind", Cov)
length(All_Cov$Sample.name) #290
length(unique(All_Cov$Sample.name)) #288
All_Cov[duplicated(All_Cov$Sample.name),] #PPN008, PPN042
All_Cov[which(All_Cov$Sample.name=="PPN008" | All_Cov$Sample.name == "PPN042"),] # re-ran these 2 with diff encoding specified
which(All_Cov$Sample.name=="PPN008" | All_Cov$Sample.name == "PPN042") # will remove rows 162 + 164
All_Cov <- All_Cov[-c(162,164),]

Mapping <- ImportData("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Manuscript/coverage_data/MappingStats",
                      "*.txt")
str(Mapping)

All_Mapping <- do.call("rbind", Mapping)
length(All_Mapping$Sample.name) #291
length(unique(All_Mapping$Sample.name)) #288
All_Mapping[duplicated(All_Mapping$Sample.name),] #PPN082, PPN042, PPN008
All_Mapping[which(All_Mapping$Sample.name=="PPN082" | All_Mapping$Sample.name == "PPN042" |
                    All_Mapping$Sample.name=="PPN008"),]
which(All_Mapping$Sample.name=="PPN082" | All_Mapping$Sample.name == "PPN042" |
        All_Mapping$Sample.name=="PPN008") # remove rows 61,164,166
All_Mapping <- All_Mapping[-c(61,164,166),]

# merge datasets
AllData <- merge(All_Cov, All_Mapping, by="Sample.name")
length(AllData$Sample.name) #288
write.csv(AllData, file="SNP-calling/MappingData.csv", row.names = FALSE)
### combine with line info
SAM <- read.csv("SNP-calling/All_SAM_Info.csv", header=T)

SAM_data <- merge(SAM, AllData, by.x = "Line", by.y = "Sample.name")
length(SAM_data$Line) #288
colnames(SAM_data)[20] <- "Mean_Coverage"

### save to Dropbox
write.csv(SAM_data[,c(1,3:4,6:9,23,24,20)], 
          file = "/Users/emilydittmar/Dropbox/Sunflower_dSNP/Figures_Tables/TableS1.csv",
          row.names = FALSE)
  


