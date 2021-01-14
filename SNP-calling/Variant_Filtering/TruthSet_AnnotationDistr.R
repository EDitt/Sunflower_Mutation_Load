
library(tidyr)

### Looking at annotation distributions of concordant distribution

# 5583 - 5028 = 555 (10%) were removed from my filters

setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Variant_graphs")

ConcAnnAll <- read.table("ConcordantVariants.table", header = T, na.strings=c("","NA"), sep = "\t")

ConcAnnSite <- ConcAnnAll[,c(1:11)]
ConcAnnGQ <- ConcAnnAll[,c(12:299)]
ConcAnnDP <- read.table("ConcordantVariantsDP.table", header = T, na.strings=c("","NA"), sep = "\t")

######### SITE ANNOTATIONS #############

# 1.) How many sites failed filtering?

aggregate(ConcAnnSite$POS, by=list(ConcAnnSite$FILTER), length) # All passed filtering

# 2.) More than 20% no-call, low GQ, low depth

LowGQ <- apply(ConcAnnGQ, 1, function(x) {length(which(x < 6))})
HighDP <- apply(ConcAnnDP[,c(3:290)], 1, function(x) {length(which(x > 50))})

hist(LowGQ)
hist(HighDP)

ConcAnnSite$LowGQ <- LowGQ
ConcAnnSite$HighDP <- HighDP
length(which(HighDP > 0)) #73

### total samples that would be excluded
ConcAnnSite$NCalledAfterFilt <- ConcAnnSite$NCALLED - (ConcAnnSite$LowGQ + ConcAnnSite$HighDP)
hist(ConcAnnSite$NCalledAfterFilt)

### how many are excluded? (< 230 or 80% of samples )
length(which(ConcAnnSite$NCalledAfterFilt < 230)) #225 (4%)
hist(ConcAnnSite$NCalledAfterFilt)

### histogram of GQ values
GQ_all <- gather(ConcAnnGQ)
hist(GQ_all$value)
GQ_summary <- boxplot(GQ_all$value, plot=FALSE)
hist(GQ_summary$out) # lowest value is 65 (no low outliers)
quantile(GQ_all$value, c(.02, .05, .10), na.rm=TRUE) #10% is 9, 5% is 6, 1% is 0, 2% is 3

### histogram of DP values
DP_all <- gather(ConcAnnDP[,c(3:290)])
hist(DP_all$value)
hist(DP_all[which(DP_all$value < 50), "value"])
quantile(DP_all$value, c(.5, .90, .99, .999), na.rm=TRUE) # 50% is 7, 90% is 15, 99% is 23, 99.9% is 36

length(which(ConcAnnSite$NCALLED < 230)) #122 (2%)
length(which(ConcAnnSite$NCALLED < 250)) #191 (3.4%)
length(which(ConcAnnSite$NCALLED < 259)) #242

hist(ConcAnnSite$NCALLED)
boxplot(ConcAnnSite$NCALLED)

### most of the truth sites that were removed were due to not having enough called samples

# 3.) QUAL > 40 Filter

hist(ConcAnnSite$QUAL)
length(which(ConcAnnSite$QUAL < 40)) #3 (.05%)
length(which(!ConcAnnSite$QUAL > 40))

# 4.) Heterozygous filtering (sites with more than 20% heterozygotes)

ConcAnnSite$PropHet <- ConcAnnSite$HET / ConcAnnSite$NCALLED
hist(ConcAnnSite$PropHet)

length(which(!ConcAnnSite$PropHet < 0.2)) #154

# 5.) High ExcessHet values (>5)

hist(ConcAnnSite$ExcessHet)
length(which(ConcAnnSite$ExcessHet > 5)) #95 (1.7%)

# 6.) Multi-allelic sites

aggregate(ConcAnnSite$CHROM, by=list(ConcAnnSite$MULTI.ALLELIC), length) # 251 (4.5%)


#### Summary
### 5583 - 5028 = 555 lost in my annotations-

# 225 lost due to high (more than 20%) no-calls, low GQ, high DP 
#       most of this was due to no-calls- 242 had less than 90% no calls
# 3 removed due to QUAL < 40 filter
# 154 had higher than 20% heterozygosity
# 95 had high ExcessHet
# (why does this only total 477?)

# 251 are multi-allelic