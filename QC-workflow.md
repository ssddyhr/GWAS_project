# QC
## Splitting up individual chips
Splitting up chips
First in the filtering process we subset out data based on the chips, making b-files based on our original data for each. We made a R-script that splits our metadata.txt into 5 different cohorts based on the chips used, with a separate file for no chip. And then using the command
```
plink --bfile gwas_data --keep HTS_iSelect_HD.txt --make-bed --out HTS_iSelect_HD_subset.
```
I've chosen to only analyse the ones with the height phenotype from my phenotype. It contains 1106 individuals.

People left: 
275 people remaining HTS

128 people remaining for Illumina

399 people (180 males, 183 females, 36 ambiguous) loaded from .fam.
Ambiguous sex IDs written to nochip_height.nosex .
232 remaining

291 people (105 males, 121 females, 65 ambiguous) loaded from .fam.
Ambiguous sex IDs written to OmniExpress_height.nosex .
170 phenotype values present after --pheno.
--keep: 170 people remaining.

omni express plus
484 people (304 males, 180 females) loaded from .fam.
266 phenotype values present after --pheno.
--keep: 266 people remaining.

Using the command.
```
plink --bfile OmniExpress_plus --keep keep.txt --pheno height_adjusted.txt --make-bed --out OmniExpress_plus_height
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

Warning: --hwe observation counts vary by more than 10%, due to the X
chromosome.  You may want to use a less stringent --hwe p-value threshold for X
chromosome variants.
So i can either split them up either into autosome and sex-chromsome, then applying the HWE on the autosome, and a less stringent sex-chromosome. Or maybe just dont apply a HWE Removal on the sex-chromsomes at all. Afterwards i would then merge it.

I split into autosome.

```
plink --bfile HTS_iSelect_HD_sexflt --autosome --make-bed --out HTS_iSelect_HD_sexflt_auto

```
and X-chromosome.
```
plink --bfile HTS_iSelect_HD_sexflt --chr X --make-bed --out HTS_iSelect_HD_sexflt_chrX
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
Creates GWA-QC.imiss and GWA-QC.lmiss. The fourth column in the file GWA-data.imiss (N_MISS) is the number of missing SNPs and  (F_MISS) is the proportion of missing SNPs per individual.

```
plink --bfile GWA-QC --het --out GWA-QC
```
This command will create the file GWA-data.het, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of non-missing genotypes [N(NM)] per individual.
From that i can calculate the observed heterozygosity rate per individual, and then set a threshhold rate that is more than 3 s.d. from the mean, and saving them to a txt, using R. Then i remove the remaining individuals.

```
plink --bfile GWA-QC --remove wrong_het_missing_values.txt --make-bed --out GWA-QC

```

### Relatednesss (Identification of duplicated or related individuals)

To identify duplicated or related individuals we will calculate the identity by descent (IBD) matrix. This works best if it is done on a set of non-correlated SNPs. So first we will “prune” the data and create a list of SNPs where no pair (within a given genomic interval) has an r2 value greater than a given threshold, typically chosen to be 0.2. This can be done by the indep-pairwise command, using 500kb as window size and 5 variants as step size:

```
plink --bfile GWA-QC --indep-pairwise 500kb 5 0.2 --out GWA-QC

plink --bfile GWA-data_sflt_mflt_hflt --extract GWA-data_sflt_mflt_hflt.prune.in --genome --min 0.185 --out GWA-data_sflt_mflt_hflt
ls -lht
```
Use R to make these new 
```

plink --bfile GWA-data_sflt_mflt_hflt --remove wrong_ibd.txt --make-bed --out GWA-data_sflt_mflt_hflt_iflt
```

# Stratification an PCA's

```
plink --bfile /home/animaldyhr/populationgenomics/students/animaldyhr/project/merged_all_chips_qc --pca --out merged_all_chips_qc
```

# Assosiation study

```
plink --bfile /home/animaldyhr/populationgenomics/students/animaldyhr/project/merged_all_chips_qc \
    --pheno formatted_height.txt \
    --linear \
    --pheno-name height \
    --covar all_covariates.txt \
    --covar-name chip,SEX,PC1,PC2,PC3 \
    --dummy-coding \
    --out height_association_by_chip_pc_gender
```
"
# Litterature
[Plink, v1.9](https://www.cog-genomics.org/plink/1.9/)
[PMC3066182](https://pmc.ncbi.nlm.nih.gov/articles/PMC3066182/)
[PMC3025522](https://pmc.ncbi.nlm.nih.gov/articles/PMC3025522/)


