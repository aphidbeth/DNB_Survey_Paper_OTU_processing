# Prerocessing of OTU tables for needle cast pathogen analysis

# Written by Beth Moore 10/02/2026

# This script has following steps: 
# 1) Convert OTU tables to relative abundance
# 2) Create fasta files containing the raw sequences for all the OTUs of potential
#   needle cast pathogens
# 3) *** RUN MASSBLASTER = outwith R in the plutoF platform ***
# 4) Read in and parse the massblaster results htmls
# 5) Identify agreement/conflict and create a list of pathogen OTUs to rename
# 6) Output the renamed pathogen OTUs


#===============================================================================
# SETUP

rm(list=ls())

library("tidyverse"); packageVersion("tidyverse")
library("qiime2R"); packageVersion("qiime2R")
library("phyloseq"); packageVersion("phyloseq")
library("microViz"); packageVersion("microviz")
library("microbiome"); packageVersion("microbiome")
library("xml2"); packageVersion("xml2")
library("rgbif"); packageVersion("rgbif")
library( "rvest"); packageVersion( "rvest")
library("DECIPHER"); packageVersion("DECIPHER")

# > library("tidyverse"); packageVersion("tidyverse")
# [1] ‘2.0.0’
# > library("qiime2R"); packageVersion("qiime2R")
# [1] ‘0.99.6’
# > library("phyloseq"); packageVersion("phyloseq")
# [1] ‘1.50.0’
# > library("microViz"); packageVersion("microviz")
# [1] ‘0.12.5’
# > library("microbiome"); packageVersion("microbiome")
# [1] ‘1.28.0’
# > library("xml2"); packageVersion("xml2")
# [1] ‘1.3.8’
# > library("rgbif"); packageVersion("rgbif")
# [1] ‘3.8.4’
# > library( "rvest"); packageVersion( "rvest")
# [1] ‘1.0.4’
# > library("DECIPHER"); packageVersion("DECIPHER")
# [1] ‘3.2.0’

# R.version
# _                                
# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          4.2                              
# year           2024                             
# month          10                               
# day            31                               
# svn rev        87279                            
# language       R                                
# version.string R version 4.4.2 (2024-10-31 ucrt)
# nickname       Pile of Leaves  


# Import OTU tables
# These can be downloaded from zenodo repository: https://doi.org/10.5281/zenodo.20179422 

T1_min_clean <- readRDS("~/OTU_tables/T1_ITS_min_clean_OTU.rds") # Edit filepaths as needed
T2_min_clean <- readRDS("~/OTU_tables/T2_ITS_min_clean_OTU.rds")
T3_min_clean <- readRDS("~/OTU_tables/T3_ITS_min_clean_OTU.rds")

# Check these are phyloseq objects
T1_min_clean
T2_min_clean
T3_min_clean

# Import the sequence tables 
# These can be downloaded from zenodo repository: https://doi.org/10.5281/zenodo.20179422 

T1_seqs <- read.table("~/seqs-and-taxonomy_T1_ITS.tsv", header = TRUE) # Edit filepaths as needed
T2_seqs <- read.table("~/seqs-and-taxonomy_T2_ITS.tsv", header = TRUE)
T3_seqs <- read.table("~/seqs-and-taxonomy_T3_ITS.tsv", header = TRUE)

# Merge into a single sequence file for easy use
sequences <- rbind(T1_seqs, T2_seqs, T3_seqs) %>% unique()

str(sequences)

#===============================================================================
# 1) CONVERT TO RELATIVE ABUNDANCE

# Note to Annika: This is different from the first iteration as I do not rarefy 
# instead just converting to relative abundance
T1_min_clean_RA <- transform_sample_counts(T1_min_clean, function(x) x / sum(x))
T2_min_clean_RA <- transform_sample_counts(T2_min_clean, function(x) x / sum(x))
T3_min_clean_RA <- transform_sample_counts(T3_min_clean, function(x) x / sum(x))

# check 
min(otu_table(T1_min_clean_RA))
max(otu_table(T1_min_clean_RA))
min(otu_table(T2_min_clean_RA))
max(otu_table(T2_min_clean_RA))
min(otu_table(T3_min_clean_RA))
max(otu_table(T3_min_clean_RA))

#===============================================================================
# EXTRACT NEEDLE CAST PATHOGENS

# List of pathogens requring verification provided by Annika:
# Dothistroma septosporum
# Lophodermium seditiosum
# Lophodermella sulcigena
# Coleosporium tussilaginis
# Lophodermella conjuncta
# Cyclaneusma spp

# Additional species checked & had some presence so also assessed:
# Gremmeniella abietina
# Gremmenia infestans
# Pachyramichloridium pini
# Rhizosphaera spp
# Sirococcus piceicola

# Extract all OTUs at genus level for these taxa:
Genus_of_interest <- c("Dothistroma",
                       "Lophodermium",
                       "Lophodermella",
                       "Coleosporium",
                       "Cyclaneusma",
                       "Gremmeniella",
                       "Gremmenia",
                       "Pachyramichloridium",
                       "Rhizosphaera",
                       "Sirococcus")

T1_taxa <- tax_table(T1_min_clean_RA) %>% as.data.frame() %>% mutate(id = rownames(.))
T2_taxa <- tax_table(T2_min_clean_RA) %>% as.data.frame() %>% mutate(id= rownames(.))
T3_taxa <- tax_table(T3_min_clean_RA) %>% as.data.frame() %>% mutate(id = rownames(.))

all_taxa <- rbind(T1_taxa, T2_taxa, T3_taxa) %>% unique()
str(all_taxa)

# Split into taxa tables per genus
Dothistroma_sklearn_tt <- all_taxa %>% filter(Genus %in% "Dothistroma") %>% mutate(OTU_string = rownames(.))
Lophodermium_sklearn_tt <- all_taxa %>% filter(Genus %in% "Lophodermium")  %>% mutate(OTU_string = rownames(.))
Lophodermella_sklearn_tt <- all_taxa %>% filter(Genus %in% "Lophodermella")  %>% mutate(OTU_string = rownames(.))
Coleosporium_sklearn_tt <- all_taxa %>% filter(Genus %in% "Coleosporium")  %>% mutate(OTU_string = rownames(.))
Cyclaneusma_sklearn_tt <- all_taxa %>% filter(Genus %in% "Cyclaneusma")  %>% mutate(OTU_string = rownames(.))
Gremmeniella_sklearn_tt <- all_taxa %>% filter(Genus %in% "Gremmeniella")  %>% mutate(OTU_string = rownames(.))
Gremmenia_sklearn_tt <- all_taxa %>% filter(Genus %in%  "Gremmenia") %>% mutate(OTU_string = rownames(.))
Pachyramichloridium_sklearn_tt <- all_taxa %>% filter(Genus %in% "Pachyramichloridium") %>% mutate(OTU_string = rownames(.))
Rhizosphaera_sklearn_tt <- all_taxa %>% filter(Genus %in% "Rhizosphaera") %>% mutate(OTU_string = rownames(.))
Sirococcus_sklearn_tt <- all_taxa %>% filter(Genus %in% "Sirococcus") %>% mutate(OTU_string = rownames(.))


# Write fasta files per genus:

# write_nc_fasta <- function(nc_genus, sequences_file, tax_table) {
#   
#   select_taxa <- filter(tax_table, Genus %in% nc_genus)
#   
#   select_seqs <- filter(sequences, id %in% select_taxa$id)
#   
#   genus_fasta <- as.vector(rbind(paste0(">", select_seqs$id), select_seqs$Sequence))
#   
#   outfile <- paste0("Data/Potential_Needle_Cast_Seqs/", nc_genus, ".fasta")
#   
#   message("Writing fasta file: ", outfile)
#   
#   writeLines(genus_fasta, outfile)
# }
# 
# 
# # Apply to all the potential needle cast Genus
# lapply(Genus_of_interest, write_nc_fasta, sequences = sequences, tax_table = all_taxa)

#===============================================================================
# 3) MASSBLASTER

# This is run on the PlutoF platform: https://doi.org/10.4137/EBO.S627 
# Within the massBLASTer module, BLAST+ 2.13.0 was run on the UNITE fungi and INSD 
# fungal databases. 

# Files are searchable under DNB_survey_paper_{Genus}.fasta
# html results were downloaded and are parsed below.
# MassBlaster html consist of 

#===============================================================================
# 4) IMPORT AND PARSE MASSBLASTER RESULTS

source("Scripts/parse_massblaster_function.R")

Coleosporium_massBLASTer_tt<- massBLASTer_to_tax_table("Data/massBLASTer_results/Coleosporium_massBLASTer.html")  
Cyclaneusma_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Cyclaneusma_massBLASTer.html")  
Dothistroma_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Dothistroma_massBLASTer.html")  
Gremmenia_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Gremmenia_massBLASTer.html")  
Gremmeniella_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Gremmeniella_massBLASTer.html")
Lophodermella_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Lophodermella_massBLASTer.html")
Lophodermium_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Lophodermium_massBLASTer.html")
Pachyramichloridium_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Pachyramichloridium_massBLASTer.html")
Rhizosphaera_massBLASTer_tt <- massBLASTer_to_tax_table("Data/massBLASTer_results/Rhizosphaera_massBLASTer.html")

#===============================================================================
# 5) ASSESS CONFLICT WITH THE TAXA ASSIGNED BY SKLEARN

source("Scripts/compare_taxonomic_classification_function.R")

Colesporium_res <- get_taxonomic_agreement(Coleosporium_massBLASTer_tt, Coleosporium_sklearn_tt)
Cyclaneusma_res <- get_taxonomic_agreement(Cyclaneusma_massBLASTer_tt, Cyclaneusma_sklearn_tt)
Dothistroma_res <- get_taxonomic_agreement(Dothistroma_massBLASTer_tt, Dothistroma_sklearn_tt)
Gremmenia_res <- get_taxonomic_agreement(Gremmenia_massBLASTer_tt, Gremmenia_sklearn_tt)
Gremmeniella_res <- get_taxonomic_agreement(Gremmeniella_massBLASTer_tt, Gremmeniella_sklearn_tt)
Lophodermella_res <- get_taxonomic_agreement(Lophodermella_massBLASTer_tt, Lophodermella_sklearn_tt)
Lophodermium_res <- get_taxonomic_agreement(Lophodermium_massBLASTer_tt, Lophodermium_sklearn_tt)
Pachyramichloridium_res <- get_taxonomic_agreement(Pachyramichloridium_massBLASTer_tt, Pachyramichloridium_sklearn_tt)
Rhizosphaera_res <- get_taxonomic_agreement(Rhizosphaera_massBLASTer_tt, Rhizosphaera_sklearn_tt)

#===============================================================================
# 6) EXPORT THOSE THAT REQUIRE MANUAL VERIFICATION

all_nc_results <- rbind(Colesporium_res,
                        Cyclaneusma_res,
                        Dothistroma_res,
                        Gremmenia_res,
                        Gremmeniella_res,
                        Lophodermella_res,
                        Lophodermium_res,
                        Pachyramichloridium_res,
                        Rhizosphaera_res)

for_verification <- filter(all_nc_results, consensus_assignment == "Requires manual verification")

# Write to a csv file with the following columns:
  df_for_verification <- data.frame(
    OTU_string = for_verification$OTU_string,
    massblaster_assignment = for_verification$df1_assignment,
    sklearn_assignment = for_verification$df2_assignment,
    agreement_level = for_verification$agreement_level,
    raw_seq = sequences$Sequence[match(for_verification$OTU_string, sequences$id)], # from here onwards are blank for filling in
    synonyms =  rep("", length(for_verification$OTU_string)), 
    massBLASTer_accession = rep("", length(for_verification$OTU_string)),
    sklearn_accession = rep("", length(for_verification$OTU_string)),
    massBLASTer_accession_culture = rep("", length(for_verification$OTU_string)),
    sklearn_accession_culture = rep("", length(for_verification$OTU_string)),
    massBLASTer_reference =  rep("", length(for_verification$OTU_string)),
    sklearn_reference =  rep("", length(for_verification$OTU_string)),
    manual_assignment = rep("", length(for_verification$OTU_string)),
    justification = rep("", length(for_verification$OTU_string))
    
  )
  
write.csv(df_for_verification, "Data/Needle_cast_OTUs_for_manual_verification_EMPTY.csv")  

#===============================================================================
# MANUAL VERIFICATION STEPS:

# 1) If sequences were not assigned to the species level there is little that can
# be done to decide accuracy of different approaches. In this case check taxonomy
# for errors and assign highest matching taxonomic levels between the two assignments
# 2) Check species assignments for synonyms. If this is the cause of the conflict rename to the 
# most recent accepted species name.
# 3) Find 

# For step 3 I have written a function that compares the two accession sequences 
# and outputs an alignment of ITS2 region which can then be viewed using the 
#  interactive viewer from DECIPHER to help with decision making:
source("Scripts/Prep_alignment_function.R")

##----INTERACTIVE VERSIONS-----------------------------------------------------------------
# Used for looking at each OTU with Doth 
subset_df <- filter(sequences, id == "c5f4cc5548cbc1cbe8a7fcc3e525d9be")
raw_seq <- subset_df$Sequence
prep_alignment("DQ926955.1", "AY808302.2", raw_seq) %>% BrowseSeqs()
##-----------------------------------------------------------------------------------------


# Comparison of Coleosporium_tussilaginis and Coleosporium_clematidis ITS2 seqs from culture collections:
subset_df <- filter(df_for_verification, OTU_string == "b393a92eb102f3acd011a341080fac70")
raw_seq <- subset_df$raw_seq
prep_alignment("MZ159705.1", "OM200310.1", raw_seq) %>% BrowseSeqs()
# Very similar across the ITS2 region but raw seq more similar to D. septosporum


# Comparison of Dothistroma septosporum & Dothistroma pini ITS2 seqs from culture collections:
subset_df <- filter(df_for_verification, OTU_string == "efe3427678ed5f964b38717745e611e3")
raw_seq <- subset_df$raw_seq
prep_alignment("DQ926955.1", "AY808302.2", raw_seq) %>% BrowseSeqs()
# Very similar across the ITS2 region but raw seq more similar to D. septosporum

# Comparison of Pseudocercospora_dodonaeae & Dothistroma pini from culture collections:
subset_df <- filter(df_for_verification, OTU_string == "27c45c93626d5f7a03b951a9f3459a2d")
raw_seq <- subset_df$raw_seq
prep_alignment("NR_177115.1", "AY808302.2", raw_seq) %>% BrowseSeqs()
# Raw seq more similar to D. pini

# Comparison of Grovesiella abieticola	Gremmeniella abietina:
subset_df <- filter(df_for_verification, OTU_string == "aeb85a7e2c5fa2de1714c35cd2ff5b6a")
raw_seq <- subset_df$raw_seq
prep_alignment("MN082655", "MH857809", raw_seq) %>% BrowseSeqs()
# Hard to tell which is the better alignment -> might be best to calculate a distance

# Comparison of Lophodermella_sulcigena	Lophodermella_arcuata:
subset_df <- filter(df_for_verification, OTU_string == "52477eaf55c92e5fb5ba12b789c74586")
raw_seq <- subset_df$raw_seq
prep_alignment("OQ288965", "MT906333", raw_seq) %>% BrowseSeqs()
# More similar to Lophodermella_sulcigena

# Comparison of L. pinastri and L. congigeum
subset_df <- filter(df_for_verification, OTU_string == "2a94d86c4065746da315c7c5474be208")
raw_seq <- subset_df$raw_seq
prep_alignment("MH856647","AY247751", raw_seq) %>% BrowseSeqs()
# Hard to tell for some of the OTUs but most look better aligned to pinastri and there
# is a fair amount of sequence variance between the two accessions

# Comparison of Lophodermium_actinothyrium and Lophodermium_herbarum
subset_df <- filter(df_for_verification, OTU_string == "2a94d86c4065746da315c7c5474be208")
raw_seq <- subset_df$raw_seq
prep_alignment("AY100663.1","MK584985.1", raw_seq) %>% BrowseSeqs()
# Two sequences are more similar than the raw seq to either -> relabel to lophodermium

#===============================================================================
# IMPORT AND INTEGRATE MANUAL REASSIGNMENT

verified_OTUs <- read.csv("Data/Needle_cast_OTUs_for_manual_verification_FILLED.csv")

# Replace the taxa name with the manually verified name in the needle cast df
all_nc_results$concensus_taxa <- ifelse(all_nc_results$OTU_string %in% verified_OTUs$OTU_string, 
                                        verified_OTUs$manual_assignment[match(all_nc_results$OTU_string, verified_OTUs$OTU_string)], 
                                        all_nc_results$consensus_assignment)

# Update the taxonomy with the new assignments
# !!! Note this could be done better: I just reassign at the species level here

match_OTU_to_concensus_taxa <- function(phyloseq_object, all_nc_results_df) {
  taxa_phyloseq <- as.data.frame(tax_table(phyloseq_object))
  taxa_phyloseq$concensus_taxa <- ifelse(rownames(taxa_phyloseq) %in% all_nc_results_df$OTU_string, all_nc_results_df$concensus_taxa[match(rownames(taxa_phyloseq), all_nc_results_df$OTU_string)], taxa_phyloseq$Species)  
  tax_table(phyloseq_object) <- as.matrix(taxa_phyloseq)
  return(phyloseq_object)
}

# Also manually replace the consensus for other Dothistroma hits to Dothistroma_sp given the lack of variation between the ITS2 regions
replace_doth <- function(phyloseq_object) {
  taxa_phyloseq <- as.data.frame(tax_table(phyloseq_object))
  taxa_phyloseq$concensus_taxa[taxa_phyloseq$Genus == "Dothistroma"] <- "Dothistroma_sp"  
  tax_table(phyloseq_object) <- as.matrix(taxa_phyloseq)
  return(phyloseq_object)
}


T1_min_clean_RA_verified <- match_OTU_to_concensus_taxa(T1_min_clean_RA, all_nc_results) %>%
                            replace_doth(.)
T2_min_clean_RA_verified <- match_OTU_to_concensus_taxa(T2_min_clean_RA, all_nc_results) %>%
                            replace_doth(.)
T3_min_clean_RA_verified <- match_OTU_to_concensus_taxa(T3_min_clean_RA, all_nc_results) %>%
                            replace_doth(.)

# Quick check of phyloseqs taxa:
T1_min_clean_RA_verified %>% subset_taxa(Genus == "Dothistroma") %>% plot_bar(fill = "concensus_taxa")

# Also add the concensus column to the sequences df:

sequences_mod <- sequences %>%
  left_join(all_nc_results %>% select(OTU_string, concensus_taxa), by = c("id" = "OTU_string"))

#==========================================================================================================
# SUMMARISE proportion of needle cast pathogens by species:

concensus_taxa_of_interest <- c("Dothistroma_sp",
                                "Lophodermium_seditiosum",
                                "Lophodermium_conigenum",
                                "Lophodermella_sulcigena",
                                "Lophodermella_conjuncta",
                                "Coleosporium_tussilaginis",
                                "Cyclaneusma_minus",
                                "Gremmeniella_abietina",
                                "Gremmenia_infestans",
                                "Pachyramichloridium_pini",
                                "Rhizosphaera_kalkhoffii",
                                "Sirococcus_piceicola")

rel_abund_concensus_taxa <- function(phyloseq_object, concensus_taxa_of_interest) {
  # Subset the phyloseq object to only include the taxa of interest
  phyloseq_subset <- prune_taxa(tax_table(phyloseq_object)[, "concensus_taxa"] %in% concensus_taxa_of_interest, phyloseq_object)
  
  # Phyloseq is already in relative abundance, so we can directly melt the subsetted phyloseq object to get the relative abundance of the taxa of interest in each sample
  rel_abund_df <- psmelt(phyloseq_subset)
  
  # However the consensus taxa have multiple OTUs assigned to them, so we need to sum the relative abundance of all OTUs assigned to each consensus taxa in each sample
  rel_abund_df <- rel_abund_df %>%
    group_by(Sample, concensus_taxa) %>%
    summarise(rel_abundance = sum(Abundance)) %>%
    ungroup() %>%
    pivot_wider(names_from = concensus_taxa, values_from = rel_abundance, values_fill = 0)
  
  return(rel_abund_df)
} 

T1_nc_rel_abund <- rel_abund_concensus_taxa(T1_min_clean_RA_verified, concensus_taxa_of_interest)
T2_nc_rel_abund <- rel_abund_concensus_taxa(T2_min_clean_RA_verified, concensus_taxa_of_interest)
T3_nc_rel_abund <- rel_abund_concensus_taxa(T3_min_clean_RA_verified, concensus_taxa_of_interest)

#===============================================================================
# OUTPUT DATA OBJECTS WITH MANUALLY VERIFIED NEEDLE CAST PATHOGEN OTUS

write_rds(T1_min_clean_RA_verified, "Data/Phyloseq_objects/T1_min_clean_nc_verified.rds")
write_rds(T2_min_clean_RA_verified, "Data/Phyloseq_objects/T2_min_clean_nc_verified.rds")
write_rds(T3_min_clean_RA_verified, "Data/Phyloseq_objects/T3_min_clean_nc_verified.rds")

write.table(sequences_mod, "Data/OTU_tables_and_sequences/sequences_with_consensus_nc_taxa.tsv", quote = FALSE)

write_OTU_table <- function(timepoint, amplicon, suffix, physeq_obj, output_path) {
  # Output full OTU table
  OTU_df <- as.data.frame(otu_table(physeq_obj))
  OTU_list <- data.frame(OTU_string = rownames(OTU_df))
  final_OTU <- cbind(OTU_list, OTU_df)
  
  # Write to folder
  write.table(final_OTU, 
              file.path(output_path, paste0(timepoint, "_", amplicon, suffix, "_OTU_table.tsv")),
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t")
}

# Run for each dataset
write_OTU_table("T1", "ITS", "min_clean",T1_min_clean_RA_verified, "Data/OTU_tables_and_sequences/")
write_OTU_table("T2", "ITS", "min_clean", T2_min_clean_RA_verified, "Data/OTU_tables_and_sequences/")
write_OTU_table("T3", "ITS", "min_clean", T3_min_clean_RA_verified, "Data/OTU_tables_and_sequences/")

# Write relative abundance tables
write.table(T1_nc_rel_abund, "Data/T1_relative_abundance_nc_taxa.tsv", quote = FALSE)
write.table(T2_nc_rel_abund, "Data/T2_relative_abundance_nc_taxa.tsv", quote = FALSE)
write.table(T3_nc_rel_abund, "Data/T3_relative_abundance_nc_taxa.tsv", quote = FALSE)


#===============================================================================
# Check abundance of Dothistroma reads: 

# filter the original phyloseq object to only include the OTUs of interest
prune_taxa("27c45c93626d5f7a03b951a9f3459a2d", T1_min_clean) %>% plot_bar(fill = "OTU")
prune_taxa("34a36fd08aa3c27686c0cb1568914881", T1_min_clean) %>% plot_bar(fill = "OTU")
prune_taxa("8e052f8bc08076b9ac7f22fa0c4575b9", T1_min_clean) %>% plot_bar(fill = "OTU")
prune_taxa("8e052f8bc08076b9ac7f22fa0c4575b9", T2_min_clean) %>% plot_bar(fill = "OTU")
prune_taxa("8e052f8bc08076b9ac7f22fa0c4575b9", T3_min_clean) %>% plot_bar(fill = "OTU")
prune_taxa("94ebf675ee5e022772ce135751db8bb5", T3_min_clean) %>% plot_bar(fill = "OTU")

#===============================================================================
# Final Dataframe curation for use in molecular validation

# The final dataframe underwent some manual curation in excel before being save as a
# csv. This took the Data/T3_relative_abundance_nc_taxa.tsv and:
# - Removed of non tree samples
# - Rounding the number of decimal places
# - Converted the sample names to just keep the 4 digit identifiers
# - Reordered the columns
# - Changed "Dothistroma_sp" to "Dothistroma_septosporum"

# Output saved as FinalRelativeAbundances.csv

