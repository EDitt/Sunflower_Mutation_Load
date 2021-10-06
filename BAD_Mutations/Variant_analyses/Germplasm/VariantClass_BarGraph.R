
# graphing private allele classes across heterotic groups
setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups")


private <- read.csv("Derived.csv", header=T)
names(private[c(4,6,8)])

private_long <- reshape(private, varying=list(names(private[c(4,6,8)])),
                        direction = "long",
                        times = names(private[c(4,6,8)]))

names(private_long)[7] <- "Proportion"

library(ggplot2)
p <- ggplot(data=private_long, aes(x=time, y=Proportion, fill=Category))
p + geom_bar(stat="identity") + 
  facet_wrap(~ Type) +
  theme_bw() +
  scale_x_discrete(labels = c("dSNPs", "Tolerated", "sSNPs"))

chisq.test(private[which(private$Type=="Oil","dSNPs")], 
           private[which(private$Type=="Oil","Tolerated")])

