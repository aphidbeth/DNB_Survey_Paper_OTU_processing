# DNB_Survey_Paper_OTU_processing
Processing of PineBiome OTU tables to extract needle cast pathogen abundance for DNB survey paper

## Input datasets
Minimally cleaned OTU tables & sequence and taxonomy tables from https://doi.org/10.5281/zenodo.20179422 

## Methodology
An overview of the processing steps are described in Blast_verification_method.pdf

## Output datasets

The core output dataset is a table of relative abundances of needle pathogens that have been blast verified: 
FinalRelativeAbundances.csv. This dataset allows comparison of abundances of needle pathogens across tree genotypes. 

Intermediate data files are also available including:
- The updated phyloseq objects with the manually checked needle pathogens are listed under Data/Phyloseq_objects

- Relative abundances for needle pathogens across three timepoints recorded in Data/T<>_relative_abundance_nc_taxa.tsv

Note FinalRelativeAbundances.csv is the same as T3_relative_abundance_nc_taxa.tsv but with some changes to the decimal place rounding, removal of non-tree samples and renaming of the trees by the general pine trial genotypes rather than PineBiome sample codes. 
