library(ggplot2)
library(reshape2)
setwd(dir = "C:/Users/sebsd/Desktop/GWAS")

# Set file names (adjust path if needed)
grm_file <- "merged_new.grm.gz"
id_file <- "merged_new.grm.id"

# Load IDs
ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE)
colnames(ids) <- c("FID", "IID")

# Load GRM
grm_data <- read.table(gzfile(grm_file), header = FALSE)
colnames(grm_data) <- c("ID1", "ID2", "N_SNPs", "GRM")

# Filter out self comparisons if needed (ID1 == ID2)
grm_offdiag <- subset(grm_data, ID1 != ID2)

# Summary stats of relatedness (GRM values)
summary(grm_offdiag$GRM)
hist(grm_offdiag$GRM, breaks=50, main="Histogram of genetic relatedness (off-diagonal)", xlab="Genetic Relatedness")
