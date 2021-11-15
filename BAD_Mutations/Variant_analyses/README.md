# dSNP Analysis

## Subset VCF by Consequence
Compiled Report was used to identify variant classes (see 2.BAD_Mutations/Post_processing.md)

Created VCFs of the variants in the different allele classes (see Subset_VCFbyConsequence.md)

## Genomic Patterns

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


