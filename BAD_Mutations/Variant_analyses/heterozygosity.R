### plotting proportion of heterozygotes for different MAF classes

# srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
source("BAD_Mutations/Variant_analyses/Functions.R")

######################################
######## WRANGLE GATK OUTPUT #########
######################################

het_data <- ImportFilesAsList("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/Heterozygosity", "_Het_table.txt", "All", 
                              TRUE)

het_data <- lapply(het_data, function(x) {PropHet_MAF(x)})

het_lists <- lapply(het_data, function(x) {
  het_MAF <- split(x, x$MAFbin)
})

HetBreaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25)

Het_bins <- lapply(names(het_lists), function(x) {
  lapply(het_lists[[x]], function(xi) {Hist_bins(xi, HetBreaks, "PropHets", x)})
})

Hetbins_combined <- lapply(names(Het_bins[[1]]), function(x) {
  rbind(Het_bins[[1]][[x]], Het_bins[[2]][[x]])
})
names(Hetbins_combined) <- names(Het_bins[[1]])

# fix column names
Hetbins_combined <- lapply(Hetbins_combined, function(x) {
  colnames(x)[c(1,3)] <- c("breaks", "count"); return(x)
})

Hetbins_combined <- lapply(names(Hetbins_combined), function(x) {
  Hetbins_combined[[x]][,"MAF_bin"] <- x; return(Hetbins_combined[[x]])
})

Hetbin_data <- do.call("rbind", Hetbins_combined)

write.table(Hetbin_data, file="/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/Heterozygosity/Hetbin_data.txt",
            row.names=FALSE, quote=FALSE, sep = "\t")

######################################
######### SET UP DATAFRAME 2 ##########
######################################

########### save a dataframe for looking at heterozygosity across the genome

## add columns and rbind and to use for looking at heterozygosity across the genome (for bedtools can do chr, pos-1, pos, and prophets, tab-separated)

for (i in seq_along(het_data)) {
  name <- names(het_data[i])
  het_data[[i]]["Variant_type"] <- name
}

het_genomeData <- do.call("rbind", het_data)


write.table(het_genomeData, file="/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/Heterozygosity/HetGenome_data.txt",
            row.names=FALSE, quote=FALSE, sep = "\t")

########################################
# PLOT PROP. HETS FOR DIFF MAF CLASSES #
########################################

hetbins <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterozygosity/Hetbin_data.txt",
                      header = T, sep = "\t")
levels(hetbins$count)

# make a second axis to see the smaller proportions
hetbins$ForSecondAxis <- ifelse(hetbins$breaks=="0.05",
                                hetbins$prop,
                                hetbins$prop*3)

ggplot(data=hetbins, aes(x=breaks, y=ForSecondAxis)) +
         geom_bar(stat="identity", position = position_dodge(), 
                  aes(fill=count)) +
  scale_fill_manual(name = "",
                    values=c(Col_Deleterious, Col_Synonymous),
                    labels = c("Deleterious", "Synonymous")) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
                                        labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                                   "[0.20,0.25)")) +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~ .*0.3334)) +
         facet_wrap(~MAF_bin, scales = "free",
                    nrow=5) +
  geom_vline(xintercept=0.075, linetype = "dashed") +
  ylab("") + xlab("")

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/heterozygosity_MAFclasses.pdf",
       width = 10, height = 15)





