### pie chart for Sunflower dSNP project

slices <- c(35219747, 826378, 620210, 87812, 21281)
annotation <- c("Non-coding", "Synonymous", "Nonsynonymous- tolerated or unknown", "Nonsynonymous- deleterious",
                "Stop Lost/Gained")
pct <- round(slices/sum(slices)*100, digit = 2)

labels <- (paste0(annotation, " ", pct, "%"))

pie(slices, labels=labels, col=rainbow(length(labels)))

# as function
piechart_data <- function(slices, labels) {
  pct <- round(slices/sum(slices)*100, digit = 2)
  New_labels <- (paste0(labels, " ", pct, "%"))
  return(pie(slices, labels=New_labels, col=rainbow(length(New_labels))))
}

piechart_data(slices, annotation) ## coding regions are too 

### non-coding, synonymous, non-synonymous

piechart_data(c(35219747, 826378, 708022), 
              c("Non-coding", "Synonymous", "Nonsynonymous"))

### coding regions only

piechart_data(c(826378, 620210, 87812, 21281),
              c("Synonymous", "Nonsynonymous- tolerated or unknown", "Nonsynonymous- deleterious",
                "Stop Lost/Gained"))