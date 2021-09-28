# Shared vs. private dSNPs for the different heterotic groups

## Subset VCF files for different heterotic groups among Oil lines
I'm using the reference and alternate dSNP vcf files separately. This only matters for variants that are *fixed* within a heterotic group since I'm (for now) only requiring 1 allele for each heterotic group

```bash
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv
awk -v var="RHA-Oil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #66
awk -v var="HA-Oil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #63

# RHA Oil, N= 66
sbatch --export=group='RHA-Oil',alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf',ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf',OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups' VCF_SubsetByGenotype.sh # Submitted batch job 4335543

# HA Oil, N= 63
sbatch --export=group='HA-Oil',alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf',ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf',OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups' VCF_SubsetByGenotype.sh # Submitted batch job 4335545
```

### Subset shared/private dSNPs for HA-Oil and RHA-Oil
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups
HA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AllDel.vcf.gz
RHA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/RHA-Oil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Oil_intersections ${HA_oil} ${RHA_oil}
cd ${OUTPUTDIR}/Oil_intersections
grep -v "#" 0000.vcf | wc -l # 10,505 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 10,599 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 29,582 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 29,582 RHA-oil shared by both
```

### Any fixed dSNPs?
```bash
for file in ${OUTPUTDIR}/Oil_intersections/*.vcf;
 do bcftools stats ${file} > ${file}_stats
done

# private to HA-oil - 1 over 98%; 4536 singletons
# private to RHA-oil - 1 over 35% (highest allele frequency bin shown); 4373 singletons
# shared - 4 over 99%; 3763 singletons
```

## Non-oil lines
```bash
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv
awk -v var="RHA-NonOil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #25
awk -v var="HA-NonOil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #51

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm

# RHA NonOil, N= 25
sbatch --export=group='RHA-NonOil',alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf',ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf',OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups' VCF_SubsetByGenotype.sh # Submitted batch job 4335807
# 16982 alt deleterious, 11578 ref deleterious

# HA NonOil, N= 51
sbatch --export=group='HA-NonOil',alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf',ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf',OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups' VCF_SubsetByGenotype.sh # Submitted batch job 4335808
# 28740 alt deleterious, 11767 ref deleterious

### Subset shared/private dSNPs for HA-Oil and RHA-Oil

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups
HA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-NonOil_AllDel.vcf.gz
RHA_Nonoil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/RHA-NonOil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/NonOil_intersections ${HA_NonOil} ${RHA_Nonoil}
cd ${OUTPUTDIR}/NonOil_intersections
grep -v "#" 0000.vcf | wc -l # 15,502 (private to HA-Nonoil)
grep -v "#" 0001.vcf | wc -l # 3,555 (private to RHA-Nonoil)
grep -v "#" 0002.vcf | wc -l # 25,005 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 25,005 RHA-oil shared by both

for file in ${OUTPUTDIR}/NonOil_intersections/*.vcf;
 do bcftools stats ${file} > ${file}_stats
done

# private to HA-Nonoil - 1 over 98%; 5331 singletons (a few more high frequency variants than I've seen in other groups)
# private to RHA-Nonoil - 2 over 89% (highest allele frequency bin shown); 1821 singletons
# shared - 9 over 99%; 1768 singletons
```

## Compare to Synonymous Variants

##### Oil lines
```bash
module load BCFtools/1.10.2-GCC-8.3.0
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv
OUTPUT_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups
synon_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_synonymous.vcf

group="HA-Oil"
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)

# minor allele count is greater than 0
bcftools view -Ou --samples ${genotypes} ${synon_vcf} | bcftools filter -i 'MAC > 0' -Oz -o ${OUTPUT_DIR}/HA-Oil_Synon.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/HA-Oil_Synon.vcf.gz # needs to be indexed for bcftools
bcftools view -Ov ${OUTPUT_DIR}/HA-Oil_Synon.vcf.gz | grep -v "#" | wc -l # 427,141

group="RHA-Oil"
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)
bcftools view -Ou --samples ${genotypes} ${synon_vcf} | bcftools filter -i 'MAC > 0' -Oz -o ${OUTPUT_DIR}/RHA-Oil_Synon.vcf.gz
bcftools view -Ov ${OUTPUT_DIR}/RHA-Oil_Synon.vcf.gz | grep -v "#" | wc -l # 415,510
tabix -p vcf ${OUTPUT_DIR}/RHA-Oil_Synon.vcf.gz

bcftools isec -p /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/Oil_intersections/synonymous \
${OUTPUT_DIR}/HA-Oil_Synon.vcf.gz \
${OUTPUT_DIR}/RHA-Oil_Synon.vcf.gz

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/Oil_intersections/synonymous

grep -v "#" 0000.vcf | wc -l # 95,910 (private to HA-Oil)
grep -v "#" 0001.vcf | wc -l # 84,279 (private to RHA-Oil)
grep -v "#" 0002.vcf | wc -l # 331,231 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 331,231
```

# check to make sure filtering not biased by the different criteria
```bash

HA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AllDel.vcf.gz
RHA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/RHA-Oil_AllDel.vcf.gz

bcftools filter -i 'MAC > 0' -Ou ${HA_oil} | bcftools view -Ov | grep -v "#" | wc -l # 36,048

bcftools view -Ov ${HA_oil} | grep -v "#" | wc -l # 40,087
# fixed differences?
```


# Derived Shared vs. private dSNPs for the different heterotic groups







######### SCRATCH USED FOR DEVELOPING SCRIPT
### Setup
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
```

### Groupings
```bash
# test code in script with RHA oil
group="RHA-Oil"
alt_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf
ref_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf
OUTPUT_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups

awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #66

# Number of alt alleles: # 28,487
# Number of ref alleles: # 11,694

# concatenated file:
grep -v "#" RHA-Oil_AllDel.vcf | wc -l # 40,181
```

# compare number same/different
```bash
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/intersect_tests
HA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AllDel.vcf
RHA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/RHA-Oil_AllDel.vcf

bgzip -c ${HA_oil} > ${HA_oil}.gz
tabix -p vcf ${HA_oil}.gz
bgzip -c ${RHA_oil} > ${RHA_oil}.gz
tabix -p vcf ${RHA_oil}.gz

### needs to be compressed!
# -c all means that it considers positions identical regarless of whether ALT allele matches or not
# should not matter here because all sites biallelic
bcftools isec -p ${OUTPUTDIR}/Oil ${HA_oil}.gz ${RHA_oil}.gz

grep -v "#" 0000.vcf | wc -l # 10,505 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 10,599 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 29,582 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 29,582 RHA-oil shared by both

# totals 50,686
# HA Oil: 40,087 (agrees with these numbers)
# RHA Oil: 40,181 (agrees with these numbers)
# different with the -c all flag?
bcftools isec -p ${OUTPUTDIR}/Oil2 -c all ${HA_oil}.gz ${RHA_oil}.gz

### numbers exactly the same
```

#### HA-oil
```bash
# start with HA-oil:

SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv

# last two lines are not included in vcf
group="HA-Oil"

# need comma separated list of samples
# need to remove the 2 samples that are not in the VCF
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)
awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l # 63

# for dSNPs that are the alternate allele
alt_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf
grep -v "#" $alt_vcf | wc -l # 76,016

OUTPUT_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups

module load BCFtools/1.10.2-GCC-8.3.0

bcftools query -l ${alt_vcf}

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="alt") > 0' -o ${OUTPUT_DIR}/${group}_AltDel.vcf

grep -v "#" HA-Oil_AltDel.vcf | wc -l #28,311

# does this number differ if I use GT="AA" (require homozygous)?

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="AA") > 0' -o ${OUTPUT_DIR}/${group}_AltDel_hom.vcf

grep -v "#" HA-Oil_AltDel_hom.vcf | wc -l #15,979

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="AA" | GT="het") > 0' -o ${OUTPUT_DIR}/${group}_AltDel2.vcf

grep -v "#" HA-Oil_AltDel2.vcf | wc -l # 28,311


# for dSNPs that are the reference allele
ref_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf
grep -v "#" $ref_vcf | wc -l # 11,796

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="ref") > 0' -o ${OUTPUT_DIR}/${group}_RefDel.vcf
grep -v "#" HA-Oil_RefDel.vcf | wc -l # 11,774

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="RR") > 0' -o ${OUTPUT_DIR}/${group}_RefDel_hom.vcf
grep -v "#" HA-Oil_RefDel_hom.vcf | wc -l # 11,774

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="ref"| GT="het") > 0' -o ${OUTPUT_DIR}/${group}_RefDel2.vcf
grep -v "#" HA-Oil_RefDel2.vcf | wc -l # 11,776 - need this for the reference alleles (if I'm allowing hets)


### combine the files:
bcftools concat ${OUTPUT_DIR}/${group}_AltDel.vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf -o ${OUTPUT_DIR}/${group}_AllDel.vcf
# Checking the headers and starting positions of 2 files
# Concatenating /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AltDel.vcf	0.315904 seconds
# Concatenating /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_RefDel2.vcf
#The chromosome block Ha412HOChr00c00004 is not contiguous, consider running with -a.

grep -v "#" HA-Oil_AllDel.vcf | wc -l # 28,080 (expecting 40,087, difference of 12,007)

bcftools concat -a ${OUTPUT_DIR}/${group}_AltDel.vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf -o ${OUTPUT_DIR}/${group}_AllDel.vcf
# to use -a flag, files need to be gzipped?

bgzip -c ${OUTPUT_DIR}/${group}_AltDel.vcf > ${OUTPUT_DIR}/${group}_AltDel.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_AltDel.vcf.gz
bgzip -c ${OUTPUT_DIR}/${group}_RefDel2.vcf > ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz

# also needs index
bcftools concat -a ${OUTPUT_DIR}/${group}_AltDel.vcf.gz ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz -Ov -o ${OUTPUT_DIR}/${group}_AllDel.vcf
grep -v "#" HA-Oil_AllDel.vcf | wc -l # 40,087 <- this gives the right number!

# also test with MergeVcfs (Picard)
module load picard/2.21.6-Java-11
java -jar /apps/eb/picard/2.21.6-Java-11/picard.jar MergeVcfs \
          I=${OUTPUT_DIR}/${group}_AltDel.vcf \
          I=${OUTPUT_DIR}/${group}_RefDel2.vcf \
          O=output_variants.vcf
grep -v "#" output_variants.vcf | wc -l #0

# count variants in gzipped vcf?
bcftools view -Ov ${OUTPUT_DIR}/${group}_AltDel.vcf.gz | grep -v "#" | wc -l # 28,311

#### tests
bcftools concat $alt_vcf $ref_vcf -o ${OUTPUT_DIR}/AllDel_TEST.vcf
bcftools concat $alt_vcf $alt_vcf -o ${OUTPUT_DIR}/AllDel_TEST2.vcf
grep -v "#" AllDel_TEST.vcf | wc -l # 75,982 (should be 87,812) (difference of 11,830)
grep -v "#" AllDel_TEST2.vcf | wc -l # 75,982
cat $alt_vcf $ref_vcf > ${OUTPUT_DIR}/AllDel_TEST3.vcf
grep -v "#" AllDel_TEST3.vcf | wc -l # 87,812
cat $alt_vcf $alt_vcf > ${OUTPUT_DIR}/AllDel_TEST4.vcf
grep -v "#" AllDel_TEST4.vcf | wc -l  # 152,032

### why is the concatenated vcf for HA oil not the same as the sum of the 2?

bcftools query -l ${OUTPUT_DIR}/${group}_AltDel.vcf | wc -l # 63
bcftools query -l ${OUTPUT_DIR}/${group}_RefDel2.vcf | wc -l # 63


##############
# for all HAs
awk -F',' '{if ($8~"HA" && $8!~"RHA") {print $9}}' $SAM_info

# for all RHAs
awk -F',' '{if ($8~"RHA") {print $9}}' $SAM_info





bcftools filter -i 'AC > 0'  # what is the difference?


#GT = "AA" if I want homozygous only




```