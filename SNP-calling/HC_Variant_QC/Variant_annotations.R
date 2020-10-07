# Plot distribution of annotation across raw SNPs

# compare to plots shown here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
# hard filtering recommendations listed here are meant to be lenient so may not be appropriate 
# for filtering for a HC subset

#########################
######## SETUP ##########
#########################

library(ggplot2)
library(gridExtra)

VCF <- read.csv('SNP-calling/HC_Variant_QC/RawVariants.table', 
                header = T, na.strings=c("","NA"), sep = "\t")
dim(VCF) # 13,952,515 x 12

length(which(!is.na(VCF$GQ))) # 0 NA for GQ annotation

head(VCF)

#########################
####### FIGURES #########
#########################

### QUAL by depth
# normalize quality by depth coverage
# looks similar to the example in the GATK link above
# generic recommendation is to filter out QD values below 2

hist(VCF$QD)

QD <- ggplot(VCF, aes(x=QD)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) + 
  geom_vline(xintercept=2, size=0.7, lty=2) +
  xlab("Quality by Depth")
  
  

### Fisher Strand
## measures strand bias
## Fisher's Exact Test of raw counts
## recommendation to fail variants with a FS value greater than 60
hist(VCF$FS)
VCF$log_FS <- log(VCF$FS)
  
FS <- ggplot(VCF, aes(x=FS)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
    geom_vline(xintercept=60, size=0.7, lty=2) + 
  scale_x_log10() +
  xlab("Fisher Strand (log-scaled)")
  #geom_vline(xintercept=0, size=0.7, lty=1) + 
# Removed 2347067 rows containing non-finite values (stat_density)
## with no strand bias, FS=0 (this results in infinite values for log)

### Strand Odds Ratio (SOR)
## Like FS, but takes into account ratios of reads that cover both alleles
## recommendation to fail variants with a FS value greater than 3

SOR <- ggplot(VCF, aes(x=SOR)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  geom_vline(xintercept=3, size=0.7, lty=2) +
  xlab("Strand Odds Ratio")

### Mapping Quality (MQ)
## the root mean square mapping quality over all the reads at a site.
## instead of just looking at the average mapping quality, this
## takes into account standard deviation
## When mapping qualities are good, MQ should be around 60
## filtering recommendations filter out a variant less than 40 (could also do 50?)

hist(VCF$MQ)
max(na.omit(VCF$MQ)) # 303.92

MQ <- ggplot(VCF, aes(x=MQ)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  geom_vline(xintercept=40, size=0.7, lty=2) + xlim(0,100) + 
  geom_vline(xintercept = 60, lty=3) +
  xlab("Mapping Quality")

### Mapping Quality Rank Sum Test (MQRankSum)
## A positive value means the mapping qualities of the read supporting the alternate
## allele are higher than those supporting the reference allele; a negative value
## indicates the mapping qualities of the reference allele are higher than those 
## supporting the alternate allele. A value close to zero is best
## Hard filter threshold = -12.5

hist(VCF$MQRankSum)
min(na.omit(VCF$MQRankSum)) # -17.91
max(na.omit(VCF$MQRankSum)) # 6.92

MQRS <- ggplot(VCF, aes(x=MQRankSum)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  geom_vline(xintercept=-12.5, size=0.7)

### more stringent threshold?
MQRSb <- ggplot(VCF, aes(x=MQRankSum)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  xlim(-5,5) + geom_vline(xintercept=-2, size=0.7, lty=2) +
  xlab("Mapping Quality Rank Sum Test")
#+ geom_vline(xintercept=-1.5, size=0.7, lty=2)
### change threshold to 2?

### ReadPosRankSumTest (ReadPosRankSum)
## compares whether the positions of the reference and alternate alleles are different within the reads
## an allele only near the ends of reads is indicative of error
## negative value = alternate allele is found at the ends of reads more often than the reference allele
## positive = reference allele is found at the ends of reads more often than the alternate
## zero means little difference between positions of the reference and alternate alleles in the reads
## hard filter threshold to remove variants with a value less than -8.0
hist(VCF$ReadPosRankSum)

RPRS <- ggplot(VCF, aes(x=ReadPosRankSum)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  geom_vline(xintercept=-8, size=0.7, lty=2)

# more stringent threshold
RPRSb <- ggplot(VCF, aes(x=ReadPosRankSum)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  xlim(-2.5,2.5) + 
  #geom_vline(xintercept=-1.5, size=0.7, lty=2) +
  geom_vline(xintercept=-2, size=0.7, lty=2) +
  xlab("Read Position Rank Sum Test")
### change threshold to 2?

### Depth (DP) (cutoff of 5)
DP <- ggplot(VCF, aes(x=DP)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) + 
  xlim(0,10000) +
  geom_vline(xintercept=5, size=0.7, lty=2) +
  xlab("Depth")

ggplot(VCF, aes(x=DP)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) + 
  xlim(0,1000) +
  geom_vline(xintercept=5, size=0.7, lty=2)


hist(VCF$DP)
max(VCF$DP) # 351952
min(VCF$DP) # 1
mean(VCF$DP) # 1955.92

### QUAL score (cutoff of 40?)

hist(VCF$QUAL)
QUAL <- ggplot(VCF, aes(x=QUAL)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  xlim(0,1000) +
  geom_vline(xintercept=40, size=0.7, lty=2) +
  xlab("Quality Score")

#########################
######## COMBINE ########
#########################

grid.arrange(QD, FS, SOR, MQ, MQRSb, RPRSb, DP, QUAL, nrow=4)
#ggsave("SNP-calling/HC_Variant_QC/RawVariants_annotations", device = "pdf")

#########################
### NUM AFTER FILTERS ###
#########################

# stringent filtering
VCF$filter <- ifelse(VCF$QUAL < 40, "fail",
                     ifelse(VCF$QD < 2, "fail",
                            ifelse(VCF$DP < 5, "fail",
                                   ifelse(VCF$MQ < 40, "fail",
                                          ifelse(VCF$MQRankSum < -2, "fail",
                                                 ifelse(VCF$FS > 60, "fail",
                                                        ifelse(VCF$ReadPosRankSum < -2, "fail",
                                                               ifelse(VCF$SOR > 3, "fail",
                                                                      "pass"))))))))
aggregate(VCF$POS, by=list(VCF$filter), length)
#   fail 2,925,325
#    pass 1,498,681

#when MQ threshold was 50
aggregate(VCF$POS, by=list(VCF$filter), length)
#   fail 3,777,465
#    pass  954,985
length(VCF$POS) #5,247,197

# more lenient filtering
VCF$filter <- ifelse(VCF$QUAL < 40, "fail",
                     ifelse(VCF$QD < 2, "fail",
                            ifelse(VCF$DP < 5, "fail",
                                   ifelse(VCF$MQ < 40, "fail",
                                          ifelse(VCF$MQRankSum < -12.5, "fail",
                                                 ifelse(VCF$FS > 60, "fail",
                                                        ifelse(VCF$ReadPosRankSum < -8, "fail",
                                                               ifelse(VCF$SOR > 3, "fail",
                                                                      "pass"))))))))

aggregate(VCF$POS, by=list(VCF$filter), length)
#   fail  2,130,417
#   pass  2,285,817

#########################
### DIST AFTER FILTERS ###
#########################


Pass <- '#A9E2E4'
Fail <- '#F4CCCA'

QD <- ggplot(VCF, aes(x=QD)) + geom_density(fill=c(alpha('#A9E2E4', 0.3))) +
  geom_vline(xintercept=2, size=0.7)


ggplot(data=subset(VCF, !is.na(filter)), aes(x=QD, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  geom_vline(xintercept=2, size=0.7)

ggplot(data=subset(VCF, !is.na(filter)), aes(x=MQRankSum, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  xlim(-10,5)

ggplot(data=subset(VCF, !is.na(filter)), aes(x=FS, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  xlim(0,50)

ggplot(data=subset(VCF, !is.na(filter)), aes(x=SOR, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  xlim(0,10)

ggplot(data=subset(VCF, !is.na(filter)), aes(x=ReadPosRankSum, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  xlim(-5,5)

ggplot(data=subset(VCF, !is.na(filter)), aes(x=MQ, group=filter)) + 
  geom_density(aes(fill=filter)) +
  scale_fill_manual(values=c(alpha(Fail, 0.6), alpha(Pass, 0.6))) +
  xlim(0,100) +
  geom_vline(xintercept=60, size=0.7, lty=2)

