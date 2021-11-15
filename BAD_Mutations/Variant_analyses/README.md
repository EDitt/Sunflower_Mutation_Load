# dSNP Analysis

## Navigation: Jump to Section

- [Subset VCF](#subset-vcf-by-consequence)
- [Site Frequency Spectra](#site-frequency-spectra)
- [Polarize SNPs](#polarize-snps)
- [Genomic Patterns](#genomic-patterns)
	- [Recombination across genome](#recombination-across-genome)
	- [Binning SNP classes](#binning-snp-classes)
- 

---

## Subset VCF by Consequence
Compiled Report was used to identify variant classes (see `2.BAD_Mutations/Post_processing.md`)

Created VCFs of the variants in the different allele classes (see `Subset_VCFbyConsequence.md`)

##### Numbers in various allele classes (excluding duplicates)
Need to parse VeP table


## Site Frequency Spectra
All dSNPs
```bash

module load BCFtools/1.13-GCC-8.3.0

# reference deleterious
bcftools query -f '%CHROM\t%POS\t%AC\t%AN\n' SAM_Refdeleterious.vcf.gz > AlleleFreqs/SAM_Refdeleterious_AC_AN.txt
# alternate deleterious
bcftools query -f '%CHROM\t%POS\t%AC\t%AN\n' SAM_Altdeleterious.vcf.gz > AlleleFreqs/SAM_Altdeleterious_AC_AN.txt

# N_ALT : number of alternate alleles
# AC : count of alternate alleles
# AN : number of alleles in called genotypes

### need to calculate frequency of deleterious allele (for deleterious reference it's AN - AC/AN; for deleterious alternate it's AC/AN)
# because dSNPs are polarized to some degree, I cannot compare full set to tolerated or synonymous 
```

All SNPs across the three classes of variants: deleterious, tolerated, synonymous
```bash

```

## Polarize SNPs
Created ancestral fasta file with ANGSDS. See: ANGSD directory

See `Polarize_SNPs.md` for process of creating ancestral state table



## Genomic Patterns

### Recombination across genome

Using John Bowers genetic map information, I used the file I manipulated for the truth SNPs: `SNP_Genetic_Map_Unique.txt` for the cM distance of the markers. This file contains 6984 markers that mapped uniquely. 
- The first column is the linkage group (chromosome number), second column is the locus name, third column is the distance in cM.

I also used the vcf file that I created with SNPutils- `MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf` where I mapped the SNPs to the new genome build. (N=6523)
- Here the first column is the chromosome, 2nd is the position in bp, third is Locus name

##### Create Recombination Data Files
```bash
GeneticMap=/scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
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
```

To graph Recombination patterns across genome, I used R (see Genomic_Patterns/Recombination.R script).

### Binning SNP classes across Genome
Make windows to bin count of each class of SNPs-
dSNPs, sSNPs, Tolerated SNPs
```bash
# make genome file
cd /scratch/eld72413/SunflowerGenome
awk -v OFS='\t' {'print $1,$2'} "Ha412HOv2.0-20181130.fasta.fai" | head -17 > "GenomeFile.txt"

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
mkdir GenomicBins

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BEDTools/2.30.0-GCC-8.3.0

vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

### dSNPs per window
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Refdeleterious.vcf -c > ${out_dir}/RefDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Altdeleterious.vcf -c > ${out_dir}/AltDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_tolerated.vcf.gz -c > ${out_dir}/Tolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_synonymous.vcf.gz -c > ${out_dir}/Synonymous_10MbCounts.txt

### do the same for the derived alleles
vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/GenomicBins

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedDeleterious.vcf -c > ${out_dir}/RefDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedTolerated.vcf -c > ${out_dir}/RefDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedSynonymous.vcf -c > ${out_dir}/RefDerivedSynonymous_10MbCounts.txt


##
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedDeleterious.vcf -c > ${out_dir}/AltDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedTolerated.vcf -c > ${out_dir}/AltDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedSynonymous.vcf -c > ${out_dir}/AltDerivedSynonymous_10MbCounts.txt

```


