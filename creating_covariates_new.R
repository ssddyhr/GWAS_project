library(readr)
library(dplyr)
library(data.table)
setwd("C:/Users/sebsd/Desktop/GWAS/new_data")
# Load adjusted height data
height_data <- read_table("height.txt", col_names = TRUE)
colnames(height_data) <- c("IID", "height")
height_data <- height_data %>%
  distinct(IID, .keep_all = TRUE) %>%
  mutate(IID = as.character(IID))

# Load .fam file
fam_data <- read_table("merged_new.fam", col_names = FALSE)
colnames(fam_data) <- c("FID", "IID", "father", "mother", "sex", "phenotype")
fam_data <- fam_data %>% select(-phenotype) %>%
  mutate(FID = as.character(FID),
         IID = as.character(IID))

# Merge height with fam data
initial_merged_data <- inner_join(height_data, fam_data, by = "IID")

# Load metadata and extract chip info
metadata <- read_delim("metadata.txt", delim = "\t")
metadata_chip <- metadata %>%
  mutate(user_numeric = as.numeric(user)) %>%
  filter(!is.na(user_numeric)) %>%
  mutate(IID = as.character(user_numeric)) %>%
  select(IID, chip)

# Merge chip info
combined_pheno_covar_data <- inner_join(initial_merged_data, metadata_chip, by = "IID")

# Load PLINK-generated PCA file
pca_data <- read_table2("merged_new_pca.eigenvec", col_names = FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_data) - 2)))
pca_data_top3 <- pca_data %>% select(FID, IID, PC1, PC2, PC3)

# Merge PCA data with phenotype and covariates
plink_data_ready <- inner_join(combined_pheno_covar_data, pca_data_top3, by = c("FID", "IID")) %>%
  select(FID, IID, height, SEX = sex, chip, PC1, PC2, PC3)

# Ensure correct column types
plink_data_clean <- plink_data_ready %>%
  mutate(FID = as.character(FID),
         IID = as.character(IID),
         SEX = as.numeric(SEX),
         chip = as.character(chip),
         PC1 = as.numeric(PC1),
         PC2 = as.numeric(PC2),
         PC3 = as.numeric(PC3))

# Create and save phenotype file
phenotype_data <- plink_data_clean %>%
  select(FID, IID, height)
fwrite(phenotype_data, "plink_phenotype.txt", sep = " ", col.names = TRUE, na = "NA", quote = FALSE)

# Create and save a single covariate file
covariate_data_for_plink <- plink_data_clean %>%
  select(FID, IID, PC1, PC2, PC3, SEX, chip)
fwrite(covariate_data_for_plink, "plink_covariates.txt", sep = " ", col.names = TRUE, na = "NA", quote = FALSE)