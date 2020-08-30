## Using ANGSD-wrapper to run Pop Gen analyses: https://github.com/ANGSD-wrapper/angsd-wrapper  

#### Pre-processing

Sequence data was aligned to the HA412v.2 genome - see SNP-callling for H. annuus data and Outgroups for other Helianthus spp.  
All data were processed using the "SAM_Processing" handler from sequence handling.  
The resulting .bam files were then used in 'Realigner_Target_Creator' and 'Indel_Realigner' using sequence handling and GATK v. 3.8-1

#### Setup

To install dependencies, first needed to load gsl module
```bash
module load GSL/2.6-GCC-8.3.0

# other dependencies needed:
module load SAMtools/1.10-GCC-8.3.0
module load gnuplot/5.2.2-foss-2018a
```

Interval List
Format interval list for chromosomal sequence for ANGSD

(below is when testing)
```bash
INTERVALS20k=/scratch/eld72413/SAM_seq/results2/VCF_results_new/N_Intervals/INTERVALS_20k_atNs.bed

awk '{print $1":"$2"-"$3}' $INTERVALS20k > Chromsome_regions.txt

# for testing (first 5 chromosomes):
head -49 Chromsome_regions.txt > Chromosome_regionsTest.txt

```

Interval list for chromosomal sequence for ANGSD
It is highly recommended to use an intervals file, as ANGSD is computationally expensive. I will use bedtools to choose random intervals for this. (later I might pipe into `bedtools shuffle` to exclude strings of N's)  
The genome file here is chromosomal sequence only
```bash
FASTA_INDEX=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta.fai

# genome file for bedtools that contains only chromosomal sequence
head -17 $FASTA_INDEX | awk -v OFS='\t' '{print $1,$2}' > ChromosomeLengths.txt

# use bedtools random
module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

bedtools random -l 100000 -n 250 -seed 56 -g ChromosomeLengths.txt | sort -V | awk '{print $1":"$2"-"$3}' > Random250x100k_regions.txt
```

Then I was able to clone the repository (latest commit `6d10630`) and install dependencies

#### 1. Ancestral Sequence

```bash
./angsd-wrapper Ancestral /home/eld72413/DelMut/Sunflower_Mutation_Load/ANGSD/ConfigFiles/Ancestral_Sequence_Config 
```

#### 2. Inbreeding
Tested with subset of samples (Group 7, N=9) using the first 5 chromosomes as a regions file
Using interactive job
```bash
qsub -I -q s_interq -l walltime=12:00:00 -l nodes=1:ppn=8 -l mem=22gb
```

Subset 1 - 138 samples
```bash
find $(pwd -P) -name "*.bam" | sort -V > Subset1_BamRealigned.txt

# some samples needed to be re-indexed
parallel/20200422-GCCcore-8.3.0
module load SAMtools/1.10-GCC-8.3.0

find -maxdepth 1 -name "*.bam" | parallel samtools index {}


for f in `find -maxdepth 1 -name "*.bam"`; do
	echo $f
	samtools index $f
done

```

```bash
./angsd-wrapper Inbreeding /scratch/eld72413/NSFproj/ANGSD_FILES/Inbreeding_Coefficients_Config
```
later used the ANGSD_Job.sh script to run ^ the above

#### 3. SFS
```bash
./angsd-wrapper SFS /scratch/eld72413/NSFproj/ANGSD_FILES/Site_Frequency_Spectrum_Config
```

### A smaller subset of genome for NSF presentation:
```bash
module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

bedtools random -l 10000 -n 150 -seed 65 -g ChromosomeLengths.txt | sort -V | awk '{print $1":"$2"-"$3}' > GenomeSubset/Random150x10k_regions.txt
```