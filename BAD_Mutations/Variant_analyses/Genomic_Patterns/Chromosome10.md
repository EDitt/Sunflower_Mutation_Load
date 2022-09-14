### Generate a haplotype graphic for Chromosome 10

Put all Chromosome 10 ROH SNPs into 1 folder
```bash
awk '{print $3}' /scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/sample_name_key.txt > ${OUT_DIR}/IntermediateFiles/SAM_list.txt

ROH_SNP_DIR=/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/intermediates
Chrom10_SNPs=/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/ROH_SNPs


while read line; do
	echo "printing out chromosome 10 SNPs for $line"
	awk '{if ($1=="Ha412HOChr10") {print $0}}' ${ROH_SNP_DIR}/ROH_${line}_SNPstats.txt > ${Chrom10_SNPs}/ROH_${line}_SNPstats_Chr10.txt
done < ${OUT_DIR}/IntermediateFiles/SAM_list.txt
# no ROH on chromosome 10: Hopi, PPN022, PPN038, PPN042, PPN043, PPN071, PPN073, PPN100, PPN103, PPN114, PPN135, PPN149, PPN190, PPN213, PPN239, PPN286 (N=16)


srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
R
```

Combine using R
```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")

Chrom10_rohSNPs <- ImportFilesAsList("/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/ROH_SNPs", "_SNPstats_Chr10.txt", "ROH_", FALSE)

ColLengths <- lapply(Chrom10_rohSNPs, function(x) {length(colnames(x))})
which(ColLengths < 7) # the 
Chrom10_rohSNPs[which(ColLengths == 7)]

Chrom_rohSNPs_subset <- Chrom10_rohSNPs[which(ColLengths == 7)] # N=272

for (i in seq_along(Chrom_rohSNPs_subset)) {
    name <- names(Chrom_rohSNPs_subset[i])
    colnames(Chrom_rohSNPs_subset[[i]]) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele", "Num_Alt_alleles", "Num_alleles", "Alt_Freq")
    Chrom_rohSNPs_subset[[i]]["Genotype"] <- name
  	}

Chrom_rohSNPs_subset <- lapply(Chrom_rohSNPs_subset, function(x) {
	x$geno <- ifelse(x$Num_Alt_alleles==2,
	"Alternate", 
	ifelse(x$Num_Alt_alleles==0, 
	"Reference", "Heterozygote")); return(x)
	})

Chrom10_rohSNPs_df <- do.call("rbind", Chrom_rohSNPs_subset) # 91,115,579 without reducing

write.table(Chrom10_rohSNPs_df, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/Chr10ROH_allSNPs.txt',
	quote=FALSE, row.names=FALSE, sep="\t")

# subset to test on local computer
Chrom10_rohSNPs_df_TEST <- Chrom10_rohSNPs_df[sample(nrow(Chrom10_rohSNPs_df), 100000),]
write.table(Chrom10_rohSNPs_df_TEST, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/Chr10ROH_allSNPsTEST.txt',
	quote=FALSE, row.names=FALSE, sep="\t")


#######
# make segments to reduce data size
Chrom_rohSNPs_toPlot <- lapply(Chrom_rohSNPs_subset, function(x) {
	MakeROHSegments(x)
	})

# dataframes for segments and points separate:
Chrom_rohSNPs_segments <- lapply(Chrom_rohSNPs_toPlot, function(x) {
	x[[1]]
	})
Chrom10_rohSegments_df <- do.call("rbind", Chrom_rohSNPs_segments) # 4,320,185
write.table(Chrom10_rohSegments_df, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/ROH_segments.txt',
	quote=FALSE, row.names=FALSE, sep="\t")


Chrom_rohSNPs_points <- lapply(Chrom_rohSNPs_toPlot, function(x) {
	x[[2]]
	})
Chrom10_rohPoints_df <- do.call("rbind", Chrom_rohSNPs_points) # 4,090,959

write.table(Chrom10_rohPoints_df, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/ROH_points.txt',
	quote=FALSE, row.names=FALSE, sep="\t")


# small test df
Chrom10_rohSegments_df_TEST <- Chrom10_rohSegments_df[sample(nrow(Chrom10_rohSegments_df), 100000),]

write.table(Chrom10_rohSegments_df_TEST, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/Chr10ROH_segmentsTEST.txt',
	quote=FALSE, row.names=FALSE, sep="\t")


Chrom10_rohPoints_df_TEST <- Chrom10_rohPoints_df[sample(nrow(Chrom10_rohPoints_df), 100000),]

write.table(Chrom10_rohSegments_df_TEST, 
	file='/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10/Chr10ROH_segmentsTEST.txt',
	quote=FALSE, row.names=FALSE, sep="\t")

```

Function to make continuous segments
```R
MakeROHSegments <- function(dataframe) {
dataframe <- dataframe[order(dataframe$Position),]
dataframe$PrevGeno <- "fill"
dataframe$NextGeno <- "fill"
dataframe$PrevGeno[1] <- "start"
dataframe$NextGeno[length(dataframe$NextGeno)] <- "end"
dataframe$PrevGeno[2:length(dataframe$PrevGeno)] <- dataframe$geno[1:length(dataframe$geno)-1]
dataframe$NextGeno[1:length(dataframe$NextGeno)-1] <- dataframe$geno[2:length(dataframe$geno)] 
dataframe$section <- ifelse(dataframe$PrevGeno==dataframe$geno &
                            dataframe$NextGeno!=dataframe$geno,
                                "end",
                                ifelse(dataframe$PrevGeno!=dataframe$geno &
                                        dataframe$NextGeno==dataframe$geno,
                                     "start",
                                     ifelse(dataframe$PrevGeno!=dataframe$geno &
                                              dataframe$NextGeno!=dataframe$geno,
                                            "single", "middle")))
dataframe_start <- subset(dataframe, section=="start")
dataframe_end <- subset(dataframe, section=="end")
dataframe_segments <- cbind(dataframe_start[,c("Chromosome", "Position", "Genotype", "geno", "section")],
                               dataframe_end[,c("Position", "geno", "section")])
dataframe_points <- subset(dataframe, section=="single")
return(list(dataframe_segments, dataframe_points[,c("Chromosome", "Position", "Genotype", "geno", "section")]))
}

test2 <- Chrom_rohSNPs_subset[[3]]
test2_b <- MakeROHSegments(Chrom_rohSNPs_subset)

Chrom_rohSNPs_toPlot <- lapply(Chrom_rohSNPs_subset, function(x) {
	MakeROHSegments(Chrom_rohSNPs_subset)
	})
```







#### scratch below
For each genotype, get ref/alt SNP for each position

```bash
### all the *_SNPstat.txt have alleles at ref/alt positions
# /scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/intermediates


awk '{if ($1=="Ha412HOChr10") {print $0}}' $genotype_file | head -25
```

# test on one genotype

```bash
chr10_dir=/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Chr10
genotype_file=/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/intermediates/ROH_SF_33_SNPstats.txt
roh_file=/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/intermediates/ROH_SF_33.bed

# chromosome, pos 1 (minus 1), pos 2, kb, # of codons?
awk '{if ($1=="Ha412HOChr10") {print $0}}' $roh_file | head
awk '{if ($1=="Ha412HOChr10") {print $0}}' $roh_file | wc -l # 42

awk '{if ($1=="Ha412HOChr10") {print $1"\t"$2+1}}' $roh_file | head # starting positions

# genotype file has chromosome, snp, ref allele, alt allele, # alt alleles, # alleles
awk '{if ($1=="Ha412HOChr10") {print $0}}' $genotype_file | head -25

awk '{if ($1=="Ha412HOChr10") {print $0}}' ROH_PPN001_SNPstats.txt | head -25
```










Use: /scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/GroupFreqs/${Group}_SNP_info.txt