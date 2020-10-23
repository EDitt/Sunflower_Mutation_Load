#PBS -S /bin/bash
#PBS -q batch
#PBS -N VariantQC
#PBS -l nodes=1:ppn=6
#PBS -l walltime=22:00:00
#PBS -l mem=22gb

#PBS -M dittmare@gmail.com
#PBS -m abe

module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
GATK_JAR=/usr/local/apps/eb/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk

module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1

INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120/Sunflower_SAM_SNP_Calling_snps.filtered.vcf
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120/QC
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed


num_sites=$(grep -v "#" "${INPUT_VCF}" | wc -l)
echo "Number of sites is ${num_sites}"

cd $OUTPUT_DIR
vcftools --TsTv-by-qual ${INPUT_VCF}

echo "Using GATK's variants to table function"
gatk VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F TYPE -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F TRANSITION -F HET \
     -GF GQ -GF DP \
     -O "${OUTPUT_DIR}/FilteredVariants.table"
echo "done"
