# GWAS-workflow
## Splitting up individual chips
Splitting up chips
First in the filtering process we subset out data based on the chips, making b-files based on our original data for each. We made a R-script that splits our metadata.txt into 5 different cohorts based on the chips used, with a separate file for no chip. And then using the command
## Pr SNP QC
```
plink --bfile gwas_data --keep HTS_iSelect_HD.txt --make-bed --out HTS_iSelect_HD_subset.
```
### Genotyping efficiency / call rate (missingness)
PLINK removes all SNPs where more than 5% of individuals are missing a genotype for that SNP. Because the missingness was inflated (0.5 ‘ish), because we had a combined b-file to begin with, and when i tried to filter by missing individuals, it was way too high and would remove all the individuals. So this way I don't remove too many individuals. Even though the literature says to do it by individual first, and the remove individual SNP’s. Now i can use the --mind flag to remove by individual. I do this though the flag --geno.
```
 plink --bfile OmniExpress_subset --geno 0.05 --make-bed --out OmniExpress_geno

```
### HWE

### MAF

## Per Individual QC
### Sex check

# Litterature
[Plink, v1.9](https://www.cog-genomics.org/plink/1.9/)
https://pmc.ncbi.nlm.nih.gov/articles/PMC3066182/


