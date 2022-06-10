# Calculate Fst across genome

### Make lists of HA and RHA / Oil and NonOil lines
```bash
mkdir /scratch/eld72413/SAM_seq/Fst

# uploaded file: LineKeysINFO.txt to directory

SAM_info=/scratch/eld72413/SAM_seq/Fst/LineKeywINFO.csv

# HA lines
awk -F',' '{if ($11=="HA" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_info | wc -l # 131
awk -F',' '{if ($11=="HA" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_info > /scratch/eld72413/SAM_seq/Fst/HA_lines.txt

# RHA lines
awk -F',' '{if ($11=="RHA" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_info | wc -l # 104
awk -F',' '{if ($11=="RHA" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_info > /scratch/eld72413/SAM_seq/Fst/RHA_lines.txt

# Oil lines
awk -F',' '{if ($12=="Oil" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_info | wc -l #151
awk -F',' '{if ($12=="Oil" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_info > /scratch/eld72413/SAM_seq/Fst/Oil_lines.txt

# NonOil lines
awk -F',' '{if ($12=="NonOil" && $9!="NA" && $9!="HA412") {print $0}}' $SAM_info | wc -l #108
awk -F',' '{if ($12=="NonOil" && $9!="NA" && $9!="HA412") {print $9}}' $SAM_info > /scratch/eld72413/SAM_seq/Fst/NonOil_lines.txt
```

### Calculate Fst for each variant (parallelize among chromosomes)

HA/RHA
```bash
sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
LIST1='/scratch/eld72413/SAM_seq/Fst/HA_lines.txt',\
LIST2='/scratch/eld72413/SAM_seq/Fst/RHA_lines.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/Fst/HA_RHA/HA_RHA' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Fst_calc.sh # Submitted batch job 8084025

awk '{if ($3 > 0.1) {print $0}}' HA_RHA_Ha412HOChr01.weir.fst | wc -l # 385,459 (out of 1,847,871)
awk '{if ($3 > 0.1) {print $0}}' HA_RHA_Ha412HOChr10.weir.fst | wc -l # 786,185 (out of 1,802,564)
```

Oil/NonOil
```bash
sbatch --export=GenomeFile='/scratch/eld72413/SunflowerGenome/GenomeFile.txt',\
VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',\
LIST1='/scratch/eld72413/SAM_seq/Fst/Oil_lines.txt',\
LIST2='/scratch/eld72413/SAM_seq/Fst/NonOil_lines.txt',\
OUT_PREFIX='/scratch/eld72413/SAM_seq/Fst/Oil_NonOil/Oil_NonOil' \
/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Fst_calc.sh # Submitted batch job 8093662
```

### Plot with R
```R
# srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# module load R/4.0.0-foss-2019b
# R

library(ggplot2)
library(ggpubr)

##### all chromosomes on one plot?

ImportFilesAsList <- function (Dir, Suffix, Prefix) {
	my_files <- list.files(path = Dir, pattern = Suffix, full.names = TRUE)
  	my_data <- lapply(my_files, function(x) {read.table(x, header=T, na.strings=c("NaN"))})
  	names(my_data) <- gsub(Suffix, "", my_files)
  	names(my_data) <- gsub(Dir, "", names(my_data))
  	names(my_data) <- gsub("/", "", names(my_data))
  	names(my_data) <- gsub(Prefix, "", names(my_data))
  	return(my_data)
  }

Fst_all <- ImportFilesAsList("/scratch/eld72413/SAM_seq/Fst/HA_RHA", "*.weir.fst", "HA_RHA_")

PlotFst <- function (dataframe, quantile_threshold, Fst_plot_threshold) {
	threshold <- quantile(dataframe$WEIR_AND_COCKERHAM_FST, c(quantile_threshold), na.rm = T)
	dataframe$Position <- dataframe$POS / 1000000
	#dataframe$outlier <- ifelse(dataframe$WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background")
	outliers <- subset(dataframe, WEIR_AND_COCKERHAM_FST > 0.4)
	#dataframe_sub <- subset(dataframe, WEIR_AND_COCKERHAM_FST > Fst_plot_threshold)
	plot <- ggplot(data=dataframe[!is.na(dataframe$WEIR_AND_COCKERHAM_FST),], aes(x=Position, y=WEIR_AND_COCKERHAM_FST)) +
    #geom_point(aes(color=outlier), shape = 16) +
    #scale_color_manual(values=c("grey", "red")) +
    geom_smooth(color="black") +
    geom_point(data = outliers, aes(x=Position, y=WEIR_AND_COCKERHAM_FST), color = "red", shape = 16) +
    ylab("Fst") + xlab("Position (Mbp)") +
    ylim(0,1) +
    theme_minimal()
    return (plot)
}


### Fst_plot_threshold to save filesize and memory for plotting
Fst_plots <- lapply(Fst_all, function(x) {PlotFst(x, 0.9995, 0.1)})

labels <- paste0("Chromosome ", seq(1,17, by=1))

pdf("Fst_AllVariants4.pdf")
print(ggarrange(plotlist = Fst_plots, labels=labels, 
	#ncol = 2, nrow = 9, 
	#label.y = 0.95, hjust = -1, 
	#font.label = list(size=12, face="bold"), 
	common.legend = TRUE))
dev.off()





#### scratch

# test
plot_test <- PlotFst(Fst_all$Ha412HOChr01, 0.9995, 0.1)

pdf("CHR1_test.pdf")
print(plot_test)
dev.off()



	Filename <- paste0(Dir, "/", Prefix, Chrom_Name, ".weir.fst")
	fst <- read.table(Filename, header=T, na.strings = c("NaN"))
	threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, c(threshold), na.rm = T)
	fst$outlier <- ifelse(fst$WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background")
	plot <- ggplot(data=fst[!is.na(fst$WEIR_AND_COCKERHAM_FST),], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) +
    geom_point(aes(color=outlier)) +
    scale_color_manual(values=c("grey", "red")) +
    geom_smooth(color="black") +
    theme_minimal()
    return (plot)
}


  my_files <- list.files(path = DirPath, pattern = "*.weir.fst", full.names = TRUE)
  my_data <- lapply(my_files, read.table)
  names(my_data) <- gsub("*.weir.fst", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  names(my_data) <- gsub("HA_RHA_", "", names(my_data))
  return(my_data)
}

###
ImportTxts <- function (DirPath) {
  my_files <- list.files(path = DirPath, pattern = "*.weir.fst", full.names = TRUE)
  my_data <- lapply(my_files, read.table)
  names(my_data) <- gsub("*.weir.fst", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  names(my_data) <- gsub("HA_RHA_", "", names(my_data))
  return(my_data)
}


```


Code I used to plot 1 chromosome
```R
GraphFst <- function (Dir, Chrom_Name, Prefix, threshold) {
	Filename <- paste0(Dir, "/", Prefix, Chrom_Name, ".weir.fst")
	fst <- read.table(Filename, header=T, na.strings = c("NaN"))
	threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, c(threshold), na.rm = T)
	fst$outlier <- ifelse(fst$WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background")
	plot <- ggplot(data=fst[!is.na(fst$WEIR_AND_COCKERHAM_FST),], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) +
    geom_point(aes(color=outlier)) +
    scale_color_manual(values=c("grey", "red")) +
    geom_smooth(color="black") +
    theme_minimal()
    return (plot)
}

plot_test <- GraphFst("/scratch/eld72413/SAM_seq/Fst/HA_RHA", "Ha412HOChr01", "HA_RHA_", 0.9995)
# ggsave("/scratch/eld72413/SAM_seq/Fst/HA_RHA/Chrom1_test.png", plot_test) error when I did this in the terminal window on linux

pdf("myplot.pdf")
print(plot_test)
dev.off()

plot_chr10 <- GraphFst("/scratch/eld72413/SAM_seq/Fst/HA_RHA", "Ha412HOChr10", "HA_RHA_", 0.9995)

pdf("chr10")
print(plot_chr10)
dev.off()


### to plot all points
PlotFst <- function (dataframe, quantile_threshold, Fst_plot_threshold) {
	threshold <- quantile(dataframe$WEIR_AND_COCKERHAM_FST, c(quantile_threshold), na.rm = T)
	dataframe$outlier <- ifelse(dataframe$WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background")
	dataframe$Position <- dataframe$POS / 1000000
	dataframe_sub <- subset(dataframe, WEIR_AND_COCKERHAM_FST > Fst_plot_threshold)
	plot <- ggplot(data=dataframe_sub[!is.na(dataframe_sub$WEIR_AND_COCKERHAM_FST),], aes(x=Position, y=WEIR_AND_COCKERHAM_FST)) +
    geom_point(aes(color=outlier), shape = 16) +
    scale_color_manual(values=c("grey", "red")) +
    geom_smooth(color="black") +
    ylab("Fst") + xlab("Position (Mbp)") +
    theme_minimal()
    return (plot)
}
```