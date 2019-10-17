# Generates candidate loci using given FDR.
args <- commandArgs(trailingOnly=TRUE)

fdr <- 0.20
nSNPs <- 2.4e6
for (arg in args) {
    data <- read.csv(arg, header=TRUE, stringsAsFactors=FALSE)
    data$P.value <- as.numeric(gsub('\\*', '', as.character(data$P.value)))
    data$bhCrit <- 1:nrow(data) / nSNPs * fdr
    data$SNP <- with(data, P.value < bhCrit)
    locationsOfInterest <- data$Genomic.Location[data$SNP == TRUE]
    write.table(
        locationsOfInterest,
        paste("Data/Interest/", tools::file_path_sans_ext(basename(arg)), ".interest.txt", sep=''),
        sep="\n", row.names=FALSE, col.names=FALSE
        )
}
rm(data)
