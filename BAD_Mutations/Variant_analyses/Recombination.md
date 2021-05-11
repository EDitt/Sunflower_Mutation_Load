# Recombination across genome

Using John Bowers genetic map information, I used the file I manipulated for the truth SNPs: `SNP_Genetic_Map_Unique.txt` for the cM distance of the markers. This file contains 6984 markers that mapped uniquely. 
- The first column is the linkage group (chromosome number), second column is the locus name, third column is the distance in cM.

I also used the vcf file that I created with SNPutils- `MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf` where I mapped the SNPs to the new genome build. (N=6523)
- Here the first column is the chromosome, 2nd is the position in bp, third is Locus name


```bash
GeneticMap=/scratch/eld72413/SNParray/SNP_Genetic_Map_Unique.txt
RemappedVCF=/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf

cd /scratch/eld72413/SAM_seq/Recombination

grep -v "#" $RemappedVCF | awk '{print $1,"\t",$2,"\t",$3}' > SNParray_BPpositions.txt

awk '{print $1,"\t",$2,"\t",$3}' $GeneticMap > SNParray_cMpositions.txt


# I need to change the chromosome names to match:

awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' # need to remove the four contigs
awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' | cat -n > ChromNames.txt 


while read line; do
	OldChromName=$(echo $line | awk '{print $1}')
	NewChromName=$(echo $line | awk '{print $2}')
	echo "$OldChromName to $NewChromName"
	sed -i 's/'^"${OldChromName}"'\b/'"${NewChromName}"'/g' SNParray_cMpositions.txt
done < ChromNames.txt

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load R/4.0.0-foss-
R
```

Use R to merge datasets + graph
```R
cM_pos <- read.table("SNParray_cMpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(cM_pos) <- c("Chromosome", "Locus", "cM")
length(cM_pos$Locus) # 6984

bp_pos <- read.table("SNParray_BPpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(bp_pos) <- c("Chromosome", "bp", "Locus")
length(bp_pos$Locus) # 6524

length(which(cM_pos$Chromosome %in% bp_pos$Chromosome))
length(which(cM_pos$Locus %in% bp_pos$Locus)) #6524
## there are spaces in the Locus names
cM_pos$Locus <- gsub(" ","",cM_pos$Locus)
bp_pos$Locus <- gsub(" ","",bp_pos$Locus)

Recombination <- merge(cM_pos, bp_pos, by=c("Chromosome", "Locus"))
length(Recombination$Locus) #5958
save(Recombination, file="Recombination.RData")
```

On local computer: