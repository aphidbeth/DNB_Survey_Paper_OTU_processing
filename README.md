# DNB_Survey_Paper_OTU_processing
Processing of PineBiome OTU tables to extract needle cast pathogen abundance for DNB survey paper

## Input datasets
Minimally cleaned OTU tables & sequence and taxonomy tables from https://doi.org/10.5281/zenodo.20179422 

## Methodology
An overview of the processing steps are described in Blast_verification_method.pdf

## Output datasets
The updated phyloseq objects with the manually checked needle pathogens are listed under Data/Phyloseq_objects

Relative abundances for needle pathogens are recorded in Data/T<>_relative_abundance_nc_taxa.tsv

FinalRelativeAbundances.csv is the same as T3_relative_abundance_nc_taxa.tsv but with the general pine trial genotypes rather than PineBiome sample codes. 
