# Generates candidate loci using given FDR.
args <- commandArgs()

fdr <- 0.20
nSNPs <- as.numeric(args[1])
data <- read.csv(args[2], header=TRUE, stringsAsFactors=FALSE)
data$P.value <- as.numeric(gsub('\\*', '', as.character(data$P.value)))
data$bhCrit <- 1:nrow(data) / nSNPs * fdr
data$SNP <- with(data, P.value < bhCrit)
locationsOfInterest <- data$Genomic.Location[data$SNP == TRUE]
write.table(
    locationsOfInterest,
    paste("Data/Interest/", tools::file_path_sans_ext(basename(args[2])), ".interest.txt", sep=''),
    sep="\n", row.names=FALSE, col.names=FALSE
    )
