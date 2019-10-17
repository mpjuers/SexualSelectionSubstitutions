# Generates candidate loci using given FDR.

fdr <- 0.20
nSNPs <- 2.4e6
data <- read.csv("../../Data/matingBehaviorSNPs.csv", header=TRUE, stringsAsFactors=FALSE)
data$P.value <- as.numeric(gsub('\\*', '', as.character(data$P.value)))
data$bhCrit <- 1:nrow(data) / nSNPs * fdr
data$SNP <- with(data, P.value < bhCrit)
locationsOfInterest <- data$Genomic.Location[data$SNP == TRUE]
write.table(
    locationsOfInterest, "../../Data/locationsOfInterestCourtship.txt",
    sep="\n", row.names=FALSE, col.names=FALSE
    )
