# Numbers in each variant class

### Numbers in various allele classes (excluding duplicates)
- The number of deleterious and tolerated sites were identified and duplicates excluded from variant lists (see `2.BAD_Mutations/Post_processing.md`)

After removing duplicate sites:
- Total Deleterious: 87,812
	- Reference allele deleterious: 11,796
	- Alternate allele deleterious: 76,016
- Total Tolerated: 553,720

- Also want to identify- (see below for code used to parse numbers)
	- the number of premature stop codons: 19,081
	- stop codon lost: 2,241
	- total number of missense variants: 708,022 (708,311- 289 represented as early stop codons or stop lost)
	- total number of synonymous variants: (excluding duplicates within and across more severe consequences): 826,378
	- total number of intergenic variants: 24,148,083


### VeP annotation classes represented:

```bash
# VeP categories:
grep -v "#" $VEP | awk '{print $7}' | sort -u

```

In order of severity (according to VeP):
1.) Splice donor/acceptor variant
		splice_acceptor_variant
		splice_donor_variant
2.) Stop/Start Lost/Gained
		start_lost
		start_lost,splice_region_variant
		stop_gained
		stop_gained,splice_region_variant
		stop_lost
		stop_lost,splice_region_variant
3.) Missense Variant
		missense_variant
4.) Synonymous (includes splice region variant and stop retained variants)
		missense_variant,splice_region_variant (because this wasn't in the group tested by BAD_Mutations)
		stop_retained_variant
		synonymous_variant
		splice_region_variant,stop_retained_variant
		splice_region_variant,synonymous_variant
		splice_region_variant,3_prime_UTR_variant
		splice_region_variant,5_prime_UTR_variant
		splice_region_variant,intron_variant
5.) Non-coding (includes 5' UTR, 3' UTR, intron, upstream, downstream, intergenic)
		3_prime_UTR_variant
		5_prime_UTR_variant
		downstream_gene_variant
		intergenic_variant
		intron_variant
		upstream_gene_variant


## Parse VeP

```bash
DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs

cd $DIR

wc -l Reference_DelPositions.txt # 11796
wc -l Alternate_DelPositionsNoDups.txt # 76016
wc -l ToleratedPositionsNoDups.txt # 553720
#wc -l Synonymous_positionsNoDups.txt # 827102

VEP=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
grep -v "#" $VEP | wc -l # 43,019,114
grep -v "#" $VEP | awk '{print $2}' | sort -u | wc -l # 37,120,112 (why diff. from VCF: #37,129,915?)
# total of 5,899,002 duplicates

# how many categories represented?
awk '{print $7}' $VEP | sort -u

### premature stop codon:
grep -v "#" $VEP | awk '{if ($7~"stop_gained") {print $2}}' | wc -l # 19,101
grep -v "#" $VEP | awk '{if ($7~"stop_gained") {print $2}}' | sort -u | wc -l # 19,081
grep -v "#" $VEP | awk '{if ($7~"stop_gained") {print $2}}' | sort -u > All_StopGained.txt

### loss of stop codon:
grep -v "#" $VEP | awk '{if ($7~"stop_lost") {print $2}}' | wc -l # 2,241
grep -v "#" $VEP | awk '{if ($7~"stop_lost") {print $2}}' | sort -u | wc -l # 2,241
grep -v "#" $VEP | awk '{if ($7~"stop_lost") {print $2}}' | sort -u > All_StopLost.txt

### overlap?:
sort All_StopGained.txt All_StopLost.txt | uniq -d | wc -l # 1

### how many overlap predicted deleterious missense mutations?
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | awk 'BEGIN {OFS=":"; print $0}{$1=$1}1' | wc -l # 87,813
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | awk 'BEGIN {OFS=":"; print $0}{$1=$1}1' | sort -u | wc -l # 87,795
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | awk 'BEGIN {OFS=":"; print $0}{$1=$1}1' | sort -u > All_dSNP_Positions.txt
sort All_dSNP_Positions.txt All_StopGained.txt | uniq -d | wc -l # 34
sort All_dSNP_Positions.txt All_StopLost.txt | uniq -d | wc -l # 6

### all missense
grep -v "#" $VEP | awk '{if ($7~"missense_variant") {print $2}}' | wc -l # 712,631
grep -v "#" $VEP | awk '{if ($7~"missense_variant") {print $2}}' | sort -u | wc -l # 708,311 (previously 699,805 until ~)
grep -v "#" $VEP | awk '{if ($7~"missense_variant") {print $2}}' | sort -u > All_missense_variant.txt
# 4,270 duplicates within missense
# 553,720 tolerated & 87,812 deleterious = 641,532 (remainder due to not being able to align?)

sort All_missense_variant.txt All_StopGained.txt | uniq -d | wc -l # 232 (34 represented as deleterious or tolerated)
sort All_missense_variant.txt All_StopLost.txt | uniq -d | wc -l # 57 (6 represented as deleterious or tolerated)
# subtract 232+57 = 289 from total missense
uniq -d <(sort <(sort -u All_missense_variant.txt) <(sort -u All_StopLost.txt All_StopGained.txt)) | wc -l # 289

### synonymous (some synonymous also annotated as 'splice region variant' but going by most severe consequence, thus using ==)
grep -v "#" $VEP | awk '{if ($7=="synonymous_variant") {print $2}}' | wc -l # 835,622 (846,651 using ~)
grep -v "#" $VEP | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u | wc -l # 832,090 (843,074 using ~)
# 3,577 duplicates within synonymous

### Number of synonymous mutations annotated as missense OR as premature stop codon/stop lost
grep -v "#" $VEP | awk '{if ($7=="synonymous_variant") {print $2}}' | sort -u > All_SynonymousList.txt

sort All_missense_variant.txt All_SynonymousList.txt | uniq -d | wc -l # 5384
sort All_StopGained.txt All_SynonymousList.txt | uniq -d | wc -l # 274
# both?
uniq -d <(sort <(sort -u All_SynonymousList.txt) <(sort -u All_missense_variant.txt All_StopGained.txt)) | wc -l #5658

sort All_StopLost.txt All_SynonymousList.txt | uniq -d | wc -l # 54
uniq -d <(sort <(sort -u All_SynonymousList.txt) <(sort -u All_missense_variant.txt All_StopLost.txt)) | wc -l # 5438
uniq -d <(sort <(sort -u All_SynonymousList.txt) <(sort -u All_StopGained.txt All_StopLost.txt)) | wc -l # 328

### out of 832,090 unique synonymous variants:
# 5384 annotated as missense
# 274 annotated as stop gained
# 54 annotated as stop lost
### after removing 5,712 variants annotated in a more severe variant class, there is: 826,378 synonymous variants

```

### Final position lists to use for subsetting VCF:

##### Synonymous
```bash
## make a list of missense and stop codon gained/lost
uniq <(sort <(sort -u All_missense_variant.txt) <(sort -u All_StopLost.txt All_StopGained.txt)) | wc -l # 729,343
# total unique missense: 708,022, + 19081 (premature stop) + 2241 (stop lost) - 1 (annotated as stop gained & lost)

uniq <(sort <(sort -u All_missense_variant.txt) <(sort -u All_StopLost.txt All_StopGained.txt)) > Missense_StopLostGained.txt

# check- synonymous sites to remove from list:
sort All_SynonymousList.txt Missense_StopLostGained.txt | uniq -d | wc -l # 5712 synonymous mutations also annotated as missense or stop lost/gained
# 5384 annotated as missense, 274 annotated as stop gained, 54 annotated as stop lost

# Synonymous that don't include missense or stop lost/gained (output with tab delimiter)
comm -23 --check-order All_SynonymousList.txt Missense_StopLostGained.txt | wc -l # 826,378

comm -23 --check-order All_SynonymousList.txt Missense_StopLostGained.txt | awk '{$1=$1}1' FS=':' OFS='\t' > FinalPositionFiles/Synonymous_positionsNoDups.txt

```

##### Stop Gained/Lost
```bash
sort -u All_StopGained.txt All_StopLost.txt | wc -l # 21,321 (1 overlap that is removed)
sort -u All_StopGained.txt All_StopLost.txt > All_StopGainedLost.txt
# remove dSNP positions (already accounted for, N=40):
comm -23 --check-order All_StopGainedLost.txt All_dSNP_Positions.txt | wc -l # 21,281
comm -23 --check-order All_StopGainedLost.txt All_dSNP_Positions.txt | awk '{$1=$1}1' FS=':' OFS='\t' > FinalPositionFiles/StopLostGained_positionsNoDups.txt

### overlap with tolerated?
cat ToleratedPositionsNoDups.txt| awk '{$1=$1}1' FS='\t' OFS=':' | sort -u | comm -23 --check-order All_StopGainedLost.txt - | wc -l # 21,093
cat ToleratedPositionsNoDups.txt| awk '{$1=$1}1' FS='\t' OFS=':' | sort - All_StopGainedLost.txt | uniq -d | wc -l # 228

```

##### Deleterious & Tolerated
```bash
# deleterious
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | wc -l
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | sort -u | wc -l # 87794
cat Reference_DelPositions.txt Alternate_DelPositionsNoDups.txt | sort -u > FinalPositionFiles/All_DelPositions.txt

# tolerated
### remove Stop Lost/Gained
cat ToleratedPositionsNoDups.txt| awk '{$1=$1}1' FS='\t' OFS=':'  | sort -u | comm -23 --check-order - All_StopGainedLost.txt | wc -l # 553,492

cat ToleratedPositionsNoDups.txt| awk '{$1=$1}1' FS='\t' OFS=':'  | sort -u | comm -23 --check-order - All_StopGainedLost.txt | awk '{$1=$1}1' FS=':' OFS='\t' > FinalPositionFiles/ToleratedPositionsNoDups.txt

```

##### Non-coding
Will consider categories: intergenic, downstream_gene_variant, intron_variant, upstream_gene_variant (not including splice variants)
Intergenic: 24,148,083
Intron: 2,819,301
Up or Downstream gene variant: 10,164,676
Total=  (37,132,060 - duplicates)
```bash
### Number of intergenic mutations annotated as synonymous, missense, or premature stop
grep -v "#" $VEP | awk '{if ($7=="intergenic_variant") {print $2}}' | wc -l # 24,148,083
grep -v "#" $VEP | awk '{if ($7=="intergenic_variant") {print $2}}' | sort -u | wc -l # 24,148,083
grep -v "#" $VEP | awk '{if ($7=="intergenic_variant") {print $2}}' | sort -u > All_Intergenic.txt

sort All_Intergenic.txt All_SynonymousList.txt | uniq -d | wc -l # 0
uniq -d <(sort <(sort -u All_Intergenic.txt) <(sort -u All_missense_variant.txt All_StopGained.txt)) | wc -l #0


grep -v "#" $VEP | awk '{if ($7=="intron_variant") {print $2}}' | wc -l # 2,873,055
grep -v "#" $VEP | awk '{if ($7=="intron_variant") {print $2}}' | sort -u | wc -l # 2,819,301
grep -v "#" $VEP | awk '{if ($7=="intron_variant") {print $2}}' | sort -u > All_Intron.txt

grep -v "#" $VEP | awk '{if ($7=="upstream_gene_variant" || $7=="downstream_gene_variant") {print $2}}' | wc -l # 13,646,713
grep -v "#" $VEP | awk '{if ($7=="upstream_gene_variant" || $7=="downstream_gene_variant") {print $2}}' | sort -u | wc -l # 10,164,676
grep -v "#" $VEP | awk '{if ($7=="upstream_gene_variant" || $7=="downstream_gene_variant") {print $2}}' | sort -u > All_UpDownStream

uniq -d <(sort <(sort -u All_Intergenic.txt) <(sort -u All_Intron.txt All_UpDownStream)) | wc -l # 0
sort All_Intron.txt All_UpDownStream | uniq -d | wc -l # 1,155,817

uniq <(sort <(sort -u All_Intergenic.txt) <(sort -u All_Intron.txt All_UpDownStream)) > All_Noncoding.txt
wc -l All_Noncoding.txt # 35,976,243 (all duplicates in intron-up/down stream)

comm -23 --check-order All_Noncoding.txt All_SynonymousList.txt | wc -l # 35,575,304
comm -23 --check-order All_Noncoding.txt Missense_StopLostGained.txt | wc -l # 35,618,092

# exclude total of everything other than non-coding
sort -u All_SynonymousList.txt Missense_StopLostGained.txt | wc -l # 1,555,721
sort -u All_SynonymousList.txt Missense_StopLostGained.txt | comm -23 --check-order All_Noncoding.txt - | wc -l # 35,219,747
sort -u All_SynonymousList.txt Missense_StopLostGained.txt | comm -23 --check-order All_Noncoding.txt - | awk '{$1=$1}1' FS=':' OFS='\t' > FinalPositionFiles/NonCodingPositionsNoDups.txt
```

