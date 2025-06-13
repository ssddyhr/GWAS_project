# QC
## Splitting up individual chips
Splitting up chips
First in the filtering process we subset out data based on the chips, making b-files based on our original data for each. We made a R-script that splits our metadata.txt into 5 different cohorts based on the chips used, with a separate file for no chip. And then using the command
```
plink --bfile gwas_data --keep HTS_iSelect_HD.txt --make-bed --out HTS_iSelect_HD_subset.
```
I've chosen to only analyse the ones with the height phenotype from my phenotype. It contains 1106 individuals.

Using the command.
```
plink --bfile OmniExpress_plus --keep keep.txt --pheno height.txt --make-bed --out OmniExpress_plus_height
```

## Sex check

Imputing sex did not work, as it did not infer any problematic individuals. It's best to be more conservative, so remove individuals with problematic sex instead.

```
plink --bfile OmniExpress_plus_height --impute-sex --make-bed --out OmniExpress_plus_imputed
```

First i do a sexcheck.
```
plink --bfile OmniExpress_subset  --check-sex 0.2 0.8 --out OmniExpress_subset
```
Then i filter and remove. 
```
grep -v "OK" HTS_iSelect_HD_sex.sexcheck > HTS_wrongsex.txt
```
```
plink --bfile nochip_subset --remove nochib_wrongsex.txt --make-bed --out nochip_sexflt
```

## Pr SNP QC

### Genotyping efficiency / call rate (missingness)
PLINK removes all SNPs where more than 5% of individuals are missing a genotype for that SNP. Because the missingness was inflated (0.5 ‘ish), because we had a combined b-file to begin with, and when i tried to filter by missing individuals, it was way too high and would remove all the individuals. So this way I don't remove too many individuals. Even though the literature says to do it by individual first, and the remove individual SNP’s. Now i can use the --mind flag to remove by individual. I do this though the flag --geno. The subset with no assigned chip, removed many of the SNPs, which makes sense, since the might in reality be from different SNPs. so the best option might just be to remove them.

```
plink --bfile HTS_iSelect_HD_sexflt --geno 0.05 --make-bed --out HTS_iSelect_HD_sexflt_geno
```

### identification of SNPs demonstrating a significant deviation from Hardy-Weinberg equilibrium (HWE) and removal. 

```
plink --bfile HTS_iSelect_HD_geno --hwe 1e-5 --make-bed --out HTS_iSelect_HD_HWE
```

### MAF
The final step when conducting QC is to remove all SNPs with a very low MAF. Typically a MAF threshold of 1-2% is applied.
```
plink --bfile OmniExpress_sexflt --geno 0.05 --maf 0.01 --make-bed --out OmniExpress_MAF
```
Or in one command.
```
plink --bfile OmniExpress_sexflt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out OmniExpress_sexflt
```
## Per Individual QC

### Sample call rate (missingness)

```
plink --bfile HTS_iSelect_HD_filtered --mind 0.03 --make-bed --out HTS_iSelect_HD_filtered
```

### Heterozygosity (HE)
```
plink --bfile GWA-QC --missing --out GWA-QC
```
Creates GWA-QC.imiss and GWA-QC.lmiss.
```
plink --bfile GWA-QC --het --out GWA-QC
```
This command will create the file GWA-data.het, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of non-missing genotypes [N(NM)] per individual.
From that i can calculate the observed heterozygosity rate per individual, and then set a threshhold rate that is more than 3 s.d
Then i remove the remaining individuals.

```
plink --bfile GWA --remove wrong_het_missing_values.txt --make-bed --out GWA

```

### Relatednesss (Identification of duplicated or related individuals)

```
plink --bfile GWA --indep-pairwise 500kb 5 0.2 --out GWA

plink --bfile GWA --extract GWA.prune.in --genome --min 0.185 --out GWA
ls -lht
```
```

plink --bfile GWA --remove wrong_ibd.txt --make-bed --out GWA
```
## Merging files

```
plink --bfile HTS_iSelect_HD_idb \
      --merge-list mergelist.txt \
      --make-bed \
      --out /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/merged_new
```
Had problems with my covariates, so ended up switching to plink2, which worked.

## PCA
Creating PC's
```
plink2 --bfile merged_new \
       --pca 10 \
       --out merged_new_pca
```

## Association study

```
plink2 --bfile /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/merged_new \
    --glm no-x-sex \
    --pheno /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/plink_phenotype_new.txt \
    --pheno-name height \
    --covar /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/plink_covariates_new.txt \
    --covar-name SEX,chip,PC1,PC2,PC3 \
    --out height_association_filtered
```
### Conditioning on the effect of a specific SNP (rs6033553).
```
plink2 --bfile /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/merged_new \
    --glm no-x-sex \
    --pheno /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/plink_phenotype_new.txt \
    --pheno-name height \
    --covar /home/animaldyhr/populationgenomics/students/animaldyhr/project_f/newdata/plink_covariates_new.txt \
    --covar-name SEX,chip,PC1,PC2,PC3 \
    --condition rs6033553 \
    --out conditional_height_association_rs6033553_fixed
```

# Litterature
[Plink, v1.9](https://www.cog-genomics.org/plink/1.9/)
[PMC3066182](https://pmc.ncbi.nlm.nih.gov/articles/PMC3066182/)
[PMC3025522](https://pmc.ncbi.nlm.nih.gov/articles/PMC3025522/)
