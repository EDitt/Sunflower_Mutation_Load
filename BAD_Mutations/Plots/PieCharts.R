## Pie Charts showing the proportion of variants in different annotation classes


source("BAD_Mutations/Variant_analyses/Functions.R")

#######################################
########## CODING/NON-CODING ##########
#######################################

### non-coding v. coding

piechart_data(c(35453937, 1666175), 
              c("Non-coding", "Coding"),
              c("grey", "white"))


pdf(NULL)
dev.control(displaylist="enable")
piechart_data(c(35453937, 1666175), 
              c("Non-coding", "Coding"),
              c("grey", "white"))
p1 <- recordPlot()
invisible(dev.off())

grid::grid.newpage()
p1

setEPS
postscript("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Pie1.eps", 
           height=10, width=10)
p1
dev.off()

#######################################
######### CODING ANNOTATIONS ##########
#######################################

CodingCatNums <- c(927677, 553720, 66750, 87794, 23383, 6851)
CodingCats <- c("Synonymous", "Nonsynonymous - tolerated",
                "Nonsynonymous- unalignable",
                "Nonsynonymous- deleterious",
                "Stop Lost/Gained",
                "Splice Variant")
ColorCats <- c(Col_Synonymous, Col_Tolerated,
               "#66A61E",
               Col_Deleterious,
               Col_StopStart,
               Col_Splice)

piechart_data(CodingCatNums, 
              CodingCats,
              ColorCats)

pdf(NULL)
dev.control(displaylist="enable")
piechart_data(CodingCatNums, 
              CodingCats,
              ColorCats)
p2 <- recordPlot()
invisible(dev.off())

grid::grid.newpage()
p2

setEPS
postscript("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Pie2.eps", 
           height=10, width=10)
p2
dev.off()

############################
######### ARRANGE ##########
############################


plot_grid(p1,p2, nrow=1)