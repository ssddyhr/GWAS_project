library(qqman)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)

setwd("C:/Users/sebsd/Desktop/GWAS/new_data")
# Set the correct path to your PLINK GLM output file
plink_glm_file_path <- "C:/Users/sebsd/Desktop/GWAS/new_data/assoc_results.assoc.linear"

# Read the PLINK GLM output file
gwas_results <- fread(plink_glm_file_path,
                      header = TRUE,
                      sep = "auto",
                      na.strings = "NA",
                      data.table = FALSE
)

# Prepare the data for qqman
gwas_results_clean <- gwas_results %>%
  filter(TEST == "ADD") %>%
  filter(!is.na(P), P > 0)

# Rename columns to match qqman's requirements
colnames(gwas_results_clean)[colnames(gwas_results_clean) == "#CHROM"] <- "CHR"
colnames(gwas_results_clean)[colnames(gwas_results_clean) == "POS"] <- "BP"
colnames(gwas_results_clean)[colnames(gwas_results_clean) == "ID"] <- "SNP"

# Convert non-numeric chromosome names to numeric
gwas_results_clean$CHR[gwas_results_clean$CHR == "X"] <- 23
gwas_results_clean$CHR[gwas_results_clean$CHR == "Y"] <- 24
gwas_results_clean$CHR[gwas_results_clean$CHR == "MT"] <- 25

gwas_results_clean$CHR <- as.numeric(gwas_results_clean$CHR)

# Remove rows where CHR became NA or explicitly filter out MT (chromosome 25)
gwas_results_clean <- gwas_results_clean %>%
  filter(CHR != 25) %>%
  drop_na(CHR)

# Create the Manhattan Plot
cat("\nGenerating Manhattan plot...\n")
manhattan(gwas_results_clean,
          main = "Manhattan Plot for Height Association",
          col = c("grey", "skyblue"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          cex.axis = 0.8,
          cex.lab = 1.2,
          chrlabs = c(1:22, "X", "Y")
)

# Create the QQ Plot
cat("\nGenerating QQ plot...\n")
qq(gwas_results_clean$P,
   main = "QQ Plot for Height Association",
   cex.axis = 0.8,
   cex.lab = 1.2
)

cat("\nPlotting complete.\n")

# Identify Genome-Wide Significant Variants
genome_wide_significant_variants <- gwas_results_clean %>%
  filter(P < 5e-8) %>%
  arrange(P)

if (nrow(genome_wide_significant_variants) > 0) {
  print(genome_wide_significant_variants)
} else {
  cat("No genome-wide significant variants found (P < 5e-8) across the entire genome.\n")
  
  # Show highly suggestive variants if no genome-wide significant hits
  suggestive_variants <- gwas_results_clean %>%
    filter(P < 1e-5) %>%
    arrange(P)
  
  if (nrow(suggestive_variants) > 0) {
    print(head(suggestive_variants, 20))
  } else {
    cat("No highly suggestive variants found (P < 1e-5) either.\n")
  }
}

##### conditional for rs6033553 #####

# --- Configuration ---
conditional_glm_file_path <- "C:/Users/sebsd/Desktop/GWAS/new_data/conditional_height_association_rs6033553_fixed.height.glm.linear"
conditioned_snp_id <- "rs6033553"
genomic_lambda <- 1.00272 # From PLINK log: "--adjust: Genomic inflation est. lambda (based on median chisq) = 1.00272."

# --- Load Conditional Results ---
cat("Loading conditional GWAS results...\n")
conditional_results <- fread(conditional_glm_file_path, header = TRUE, sep = "auto", na.strings = "NA", data.table = FALSE)

# Clean and format columns for plotting
conditional_results_clean <- conditional_results %>%
  filter(TEST == "ADD") %>%
  filter(!is.na(P), P > 0) %>%
  rename(SNP = "ID", CHR = "#CHROM", BP = "POS", P_unadjusted = "P")

# Apply Genomic Control adjustment
conditional_results_clean <- conditional_results_clean %>%
  mutate(
    CHISQ_UNADJ = qchisq(1 - P_unadjusted, df = 1),
    CHISQ_GC = CHISQ_UNADJ / genomic_lambda,
    P = pchisq(CHISQ_GC, df = 1, lower.tail = FALSE)
  ) %>%
  mutate(P = ifelse(P == 0, .Machine$double.xmin, P)) %>%
  select(SNP, CHR, BP, P) # Keep only necessary columns for plotting

# Convert chromosome names to numeric and remove MT
conditional_results_clean$CHR[conditional_results_clean$CHR == "X"] <- 23
conditional_results_clean$CHR[conditional_results_clean$CHR == "Y"] <- 24
conditional_results_clean$CHR[conditional_results_clean$CHR == "MT"] <- 25 # Temporarily map MT to 25
conditional_results_clean <- conditional_results_clean %>%
  mutate(CHR = as.numeric(CHR)) %>%
  filter(CHR != 25) %>% # Filter out MT (now numeric 25)
  drop_na(CHR)

# --- Display result for the conditioned SNP ---
rs6033553_conditional_info <- conditional_results_clean %>%
  filter(SNP == conditioned_snp_id)

cat(paste0("\n--- Result for ", conditioned_snp_id, " in Conditional Analysis (GC-adjusted P) ---\n"))
if (nrow(rs6033553_conditional_info) > 0) {
  print(rs6033553_conditional_info)
} else {
  cat(paste0(conditioned_snp_id, " not found in conditional results (expected after conditioning).\n"))
}

# --- Visualize: Manhattan Plot using qqman ---
cat("\nGenerating Manhattan plot using qqman...\n")
manhattan(
  conditional_results_clean,
  main = paste0("Conditional Manhattan Plot for Height (Conditioned on ", conditioned_snp_id, ")"),
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  col = c("grey", "skyblue"),
  genomewide = -log10(5e-8),
  suggestiveline = -log10(1e-5),
  highlight = NULL,
  logp = TRUE,
  cex.axis = 0.8,
  chrlabs = c(1:22, "X", "Y") # Exclude MT from labels
)

# --- Save the plot to a file ---
plot_filename <- "conditional_manhattan_plot_height_qqman.png"
png(filename = plot_filename, width = 1800, height = 700, res = 100)
manhattan(
  conditional_results_clean,
  main = paste0("Conditional Manhattan Plot for Height (Conditioned on ", conditioned_snp_id, ")"),
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  col = c("grey", "skyblue"),
  genomewide = -log10(5e-8),
  suggestiveline = -log10(1e-5),
  highlight = NULL,
  logp = TRUE,
  cex.axis = 0.8,
  chrlabs = c(1:22, "X", "Y") 
)
dev.off()

cat(paste0("\nConditional Manhattan plot generated and saved as '", plot_filename, "'.\n"))
install.packages("writexl")
library(writexl)

write_xlsx(suggestive_variants,"suggetivevariants")

