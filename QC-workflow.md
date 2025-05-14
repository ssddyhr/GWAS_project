# GWAS-workflow
## Splitting up individual chips
Splitting up chips
First in the filtering process we subset out data based on the chips, making b-files based on our original data for each. We made a R-script that splits our metadata.txt into 5 different cohorts based on the chips used, with a separate file for no chip. And then using command

plink --bfile gwas_data --keep HTS_iSelect_HD.txt --make-bed --out HTS_iSelect_HD_subset.

## Genotyping efficiency / call rate (missingness)

