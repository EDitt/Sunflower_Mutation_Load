## Using ANGSD-wrapper to run Pop Gen analyses: https://github.com/ANGSD-wrapper/angsd-wrapper  


#### Setup

To install dependencies, first needed to load gsl module
```bash
module load GSL/2.6-GCC-8.3.0

# other dependencies needed:
module load SAMtools/1.10-GCC-8.3.0
module load gnuplot/5.2.2-foss-2018a
```

Interval List
Format interval list for chromosomal sequence for ANGSD]
```bash
INTERVALS20k=/scratch/eld72413/SAM_seq/results2/VCF_results_new/N_Intervals/INTERVALS_20k_atNs.bed

awk '{print $1":"$2"-"$3}' $INTERVALS20k > Chromsome_regions.txt

# for testing (first 5 chromosomes):
head -49 Chromsome_regions.txt > Chromosome_regionsTest.txt

```

Then I was able to clone the repository (latest commit `6d10630`) and install dependencies

#### 1. Ancestral Sequence

```bash
./angsd-wrapper Ancestral /home/eld72413/DelMut/Sunflower_Mutation_Load/ANGSD/ConfigFiles/Ancestral_Sequence_Config 
```
```bash
./angsd-wrapper Inbreeding /scratch/eld72413/NSFproj/ANGSD_FILES/Inbreeding_Coefficients_Config
```