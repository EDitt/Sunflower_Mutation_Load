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
why are the numbers different?

# Derived Shared vs. private (derived) dSNPs for the different heterotic groups

### HA & RHA Oil

##### dSNPs
```bash
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv
awk -v var="RHA-Oil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #66
awk -v var="HA-Oil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #63

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm

# RHA Oil, N= 66
sbatch --export=group='RHA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedDeleterious.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedDeleterious.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/dSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4973906

# HA Oil, N= 63
sbatch --export=group='HA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedDeleterious.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedDeleterious.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/dSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4973985

### Subset shared/private dSNPs for HA-Oil and RHA-Oil
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.10.2-GCC-8.3.0
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/dSNPs
HA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/dSNPs/HA-Oil_AllDel.vcf.gz
RHA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/dSNPs/RHA-Oil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Oil_intersections ${HA_Oil} ${RHA_Oil}
cd ${OUTPUTDIR}/Oil_intersections
grep -v "#" 0000.vcf | wc -l # 7,307 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 7,366 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 16,164 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 16,164 RHA-oil shared by both

for file in ${OUTPUTDIR}/Oil_intersections/*.vcf;
 do bcftools stats ${file} > ${file}_stats
done
```
The group RHA-Oil has at least one deleterious alternate allele at 18584 sites
The group RHA-Oil has at least one deleterious reference allele at 4946 sites
Checking the headers and starting positions of 2 files

The group HA-Oil has at least one deleterious alternate allele at 18469 sites
The group HA-Oil has at least one deleterious reference allele at 5003 sites
Checking the headers and starting positions of 2 files

##### Tolerated
```bash
# RHA Oil, N= 66
sbatch --export=group='RHA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedTolerated.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedTolerated.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/Tolerated' \
VCF_SubsetByGenotype.sh # Submitted batch job 4973986

# HA Oil, N= 63
sbatch --export=group='HA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedTolerated.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedTolerated.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/Tolerated' \
VCF_SubsetByGenotype.sh # Submitted batch job 4973987

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/Tolerated
HA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/Tolerated/HA-Oil_AllDel.vcf.gz
RHA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/Tolerated/RHA-Oil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Oil_intersections ${HA_Oil} ${RHA_Oil}
cd ${OUTPUTDIR}/Oil_intersections
grep -v "#" 0000.vcf | wc -l # 38,674 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 37,233 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 173,100 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 173,100
```
The group RHA-Oil has at least one deleterious alternate allele at 122843 sites
The group RHA-Oil has at least one deleterious reference allele at 87490 sites
Checking the headers and starting positions of 2 files

The group HA-Oil has at least one deleterious alternate allele at 123837 sites
The group HA-Oil has at least one deleterious reference allele at 87937 sites
Checking the headers and starting positions of 2 files

##### Synonymous
```bash
# RHA Oil, N= 66
sbatch --export=group='RHA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedSynonymous.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedSynonymous.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/sSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4979178

# HA Oil, N= 63
sbatch --export=group='HA-Oil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedSynonymous.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedSynonymous.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/sSNPs' \
VCF_SubsetByGenotype.sh  # Submitted batch job 4979179

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/sSNPs
HA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/sSNPs/HA-Oil_AllDel.vcf.gz
RHA_Oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/sSNPs/RHA-Oil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Oil_intersections ${HA_Oil} ${RHA_Oil}
cd ${OUTPUTDIR}/Oil_intersections
grep -v "#" 0000.vcf | wc -l # 56,395 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 48,430 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 304,177 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 304177
```
The group RHA-Oil has at least one deleterious alternate allele at 194308 sites
The group RHA-Oil has at least one deleterious reference allele at 158299 sites
Checking the headers and starting positions of 2 files

The group HA-Oil has at least one deleterious alternate allele at 201333 sites
The group HA-Oil has at least one deleterious reference allele at 159239 sites
Checking the headers and starting positions of 2 files

### HA & RHA Non-Oil

##### dSNPs
```bash
SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv
awk -v var="RHA-NonOil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #25
awk -v var="HA-NonOil" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #51

cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm

# RHA NonOil, N= 25
sbatch --export=group='RHA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedDeleterious.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedDeleterious.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/dSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4976054

# HA NonOil, N= 51
sbatch --export=group='HA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedDeleterious.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedDeleterious.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/dSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4976055

### Subset shared/private dSNPs for HA-NonOil and RHA-NonOil
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.10.2-GCC-8.3.0
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/dSNPs
HA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/dSNPs/HA-NonOil_AllDel.vcf.gz
RHA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/dSNPs/RHA-NonOil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Non-Oil_intersections ${HA_NonOil} ${RHA_NonOil}
cd ${OUTPUTDIR}/Non-Oil_intersections
grep -v "#" 0000.vcf | wc -l # 10,637 (private to HA-Nonoil)
grep -v "#" 0001.vcf | wc -l # 2,457 (private to RHA-Nonoil)
grep -v "#" 0002.vcf | wc -l # 13,128 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 13,128 RHA-oil shared by both

for file in ${OUTPUTDIR}/Non-Oil_intersections/*.vcf;
 do bcftools stats ${file} > ${file}_stats
done
```
The group RHA-NonOil has at least one deleterious alternate allele at 10703 sites
The group RHA-NonOil has at least one deleterious reference allele at 4882 sites
Checking the headers and starting positions of 2 files

The group HA-NonOil has at least one deleterious alternate allele at 18769 sites
The group HA-NonOil has at least one deleterious reference allele at 4996 sites
Checking the headers and starting positions of 2 files

##### Tolerated
```bash
# RHA NonOil, N= 25
sbatch --export=group='RHA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedTolerated.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedTolerated.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/Tolerated' \
VCF_SubsetByGenotype.sh # Submitted batch job 4976056

# HA NonOil, N= 51
sbatch --export=group='HA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedTolerated.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedTolerated.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/Tolerated' \
VCF_SubsetByGenotype.sh # Submitted batch job 4976057

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/Tolerated
HA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/Tolerated/HA-NonOil_AllDel.vcf.gz
RHA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/Tolerated/RHA-NonOil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Non-Oil_intersections ${HA_NonOil} ${RHA_NonOil}
cd ${OUTPUTDIR}/Non-Oil_intersections
grep -v "#" 0000.vcf | wc -l # 58,674 (private to HA-Nonoil)
grep -v "#" 0001.vcf | wc -l # 11,967 (private to RHA-Nonoil)
grep -v "#" 0002.vcf | wc -l # 152,524 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 152,524
```

The group RHA-NonOil has at least one deleterious alternate allele at 77876 sites
The group RHA-NonOil has at least one deleterious reference allele at 86615 sites
Checking the headers and starting positions of 2 files

The group HA-NonOil has at least one deleterious alternate allele at 123328 sites
The group HA-NonOil has at least one deleterious reference allele at 87870 sites
Checking the headers and starting positions of 2 files

##### Synonymous
```bash
# RHA NonOil, N= 25
sbatch --export=group='RHA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedSynonymous.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedSynonymous.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/sSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4979204

# HA NonOil, N= 51
sbatch --export=group='HA-NonOil',\
alt_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_AltDerivedSynonymous.vcf',\
ref_vcf='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/SAM_RefDerivedSynonymous.vcf',\
OUTPUT_DIR='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/sSNPs' \
VCF_SubsetByGenotype.sh # Submitted batch job 4979207

OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/sSNPs
HA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/sSNPs/HA-NonOil_AllDel.vcf.gz
RHA_NonOil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/HeteroticGroupOverlap/NonOil/sSNPs/RHA-NonOil_AllDel.vcf.gz

bcftools isec -p ${OUTPUTDIR}/Non-Oil_intersections ${HA_NonOil} ${RHA_NonOil}
cd ${OUTPUTDIR}/Non-Oil_intersections
grep -v "#" 0000.vcf | wc -l # 83,441 (private to HA-Nonoil)
grep -v "#" 0001.vcf | wc -l # 16,088 (private to RHA-Nonoil)
grep -v "#" 0002.vcf | wc -l # 273,197 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 273,197
```

The group RHA-NonOil has at least one deleterious alternate allele at 132545 sites
The group RHA-NonOil has at least one deleterious reference allele at 156740 sites
Checking the headers and starting positions of 2 files

The group HA-NonOil has at least one deleterious alternate allele at 197523 sites
The group HA-NonOil has at least one deleterious reference allele at 159115 sites
Checking the headers and starting positions of 2 files