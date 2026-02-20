# Function for converting MassBlaster html outputs into a taxa table

# Written by Beth Moore with assistance from chatgpt/copilot for sections
# 11/02/2026


# Input: 
  # Standard html output from MassBlaster in regular mode where each block
  # starts with the Query x of Y: OTU id

# Output: 
  # A taxa table with the following columns:
    # id: OTU id
    # Kingdom
    # Phylum
    # Class
    # Order
    # Family
    # Species

# Notes
  # This approach goes through a series of steps to reduce the top alignments to
  # single consensus alignment based on the score evalue and percentage identity.
  

massBLASTer_to_tax_table <- function(massblaster_html){
    
    if (!requireNamespace("xml2", quietly = TRUE)) {
      stop("Package 'xml2' must be installed to use this function.")
    }
    
    if (!requireNamespace("rgbif", quietly = TRUE)) {
      stop("Package 'rgbif' must be installed to use this function.")
    }
    
    if (!requireNamespace("rvest", quietly = TRUE)) {
      stop("Package 'rvest' must be installed to use this function.")
    }
  
  #=====================================================
  # STEP 1: Parse massblaster html
  #=====================================================
  
    page <- read_html(massblaster_html)
    
    # --- 1. Find all query header nodes ---
    query_nodes <- xml_find_all(page, "//b[contains(text(),'Query ')]")
    
    # --- 2. Extract full text from parent node to get query ID and length ---
    query_full_texts <- sapply(query_nodes, function(node) {
      xml_text(xml_parent(node))
    })
    
    # --- 3. Extract query ID and length ---
    query_ids <- str_extract(query_full_texts, "(?<=:\\s)[^(]+") %>% str_trim()
    query_lengths <- str_extract(query_full_texts, "\\((\\d+)\\s*bp") %>% str_extract("\\d+") %>% as.integer()
    
    results_list <- vector("list", length(query_nodes))
    
    # --- 4. Loop over each query ---
    for(i in seq_along(query_nodes)) {
      
      # Find the first table immediately after the query node
      tbl_node <- xml_find_first(query_nodes[[i]], "following-sibling::table[1]")
      
      # Fallback if no sibling table
      if(is.na(tbl_node)) {
        tbl_node <- xml_find_first(query_nodes[[i]], "following::table[1]")
      }
      
      # Skip if no table found
      if(is.na(tbl_node)) {
        warning(paste("No table found for query", query_ids[i]))
        next
      }
      
      # Convert HTML table to data frame
      df <- suppressWarnings(html_table(tbl_node, fill = TRUE))
      
      # Replace NAs with quotations so they don't get dropped
      df[is.na(df)] <- ""
      
      # Standard column names
      std_names <- c("Reference","SH_0.5","SH_1.0","SH_1.5",
                     "Taxon","Score","Evalue","Prcnt","MisM",
                     "Qstart","Qend","Rstart","Rend")
      names(df) <- std_names[1:ncol(df)]
      
      # Convert all columns to character to avoid type issues
      df <- mutate_all(df, as.character)
      
      # Add query info
      df$query_id <- query_ids[i]
      df$query_length <- query_lengths[i]
      
      results_list[[i]] <- df
    }
    
    # Combine all results into a single table
    results_df <- bind_rows(results_list)

    print(colnames(results_df))

    #==============================================================
    # STEP 2: Fill in the higher taxonomic levels for the blast hits
    #================================================================
    
    # The blast assignments do not classify to higher taxonomic levels 
    # but matching at higher taxonomic levels is important for the descision 
    # process. To find the higher level taxonomy we use name_backbone_verbose 
    # function from the gbifR package. 
  
    # Function to pull the taxonomy for each hit
    fill_taxonomy <- function(name) {
      
      res <- name_backbone_verbose(name)
      
      # Get the two parts of the output list
      res_dat <- res$data
      res_alt <- res$alternatives
      
      # Check if res$data or res$alternatives is empty or missing
      if (is.null(res_dat) || nrow(res_dat) == 0) {
        warning(paste("No data found for:", name))
        res_dat <- data.frame(Tax_df_name = "Data_df", kingdom = "NA") 
      }
      
      if (is.null(res_alt) || nrow(res_alt) == 0) {
        warning(paste("No alternatives found for:", name))
        res_alt <- data.frame(Tax_df_name = "Alternative_df", kingdom = "NA")  
      }
      
      colnames_all <- c(colnames(res_dat), colnames(res_alt)) %>% unique()
      
      missing_columns_data <- setdiff(colnames_all, colnames(res_dat))
      missing_columns_alternatives <- setdiff(colnames_all, colnames(res_alt))
      
      # Add missing columns to the data frames before merging
      for (col in missing_columns_data) {
        res_dat[[col]] <- NA
      }
      
      # Add missing columns to the data frames before merging
      for (col in missing_columns_alternatives){
        res_alt[[col]] <- NA
      }
      
      # Reorder the columns into the same order
      res_dat <- select(res_dat, all_of(colnames_all))
      res_alt <- select(res_alt, all_of(colnames_all))
      
      Options <- rbind(res_dat,
                       res_alt)
      
      # Filter to fungal options
      
      fungalOptions <- dplyr::filter(Options, kingdom == "Fungi")
      print(fungalOptions)
      
      # If still empty, fallback to the first alternative
      if (nrow(fungalOptions) == 0) {
        warning(paste("No fungal taxonomy found for:", name))
        fungalOptions <- Options
      }
      
      custom_order <- c("EXACT",
                        "VARIANT",
                        "HIGHERRANK",
                        "NONE")
      
      # Convert the charVec column to a factor with the custom order
      fungalOptions$matchType <- factor(fungalOptions$matchType, levels = custom_order)
      
      # Sort by confidence DESC
      fungalOptions <- fungalOptions[order(-fungalOptions$confidence, fungalOptions$matchType), ]
      
      # TKE THE TOP MATCH
      fungalOptions <- fungalOptions[1,]
      
      # Select just the tax columns we are interested in
      tax_cols <- c("kingdom","phylum","class","order","family","genus","species")
      subset_cols <- dplyr::select(fungalOptions, dplyr::any_of(tax_cols))
      
      # Add original input name
      subset_cols$Taxon <- name
      
      return(subset_cols)
    }

    # Make a vector of the blast hits
    taxa_vector <- results_df$Taxon %>% unique()
    
    # Apply fill_taxonomy to each element
    taxonomy_list <- lapply(taxa_vector, fill_taxonomy)
    
    # Ensure all columns are present for merging
    all_columns <- unique(unlist(lapply(taxonomy_list, colnames)))
    
    taxonomy_list <- lapply(taxonomy_list, function(df) {
      missing_cols <- setdiff(all_columns, colnames(df))
      if(length(missing_cols) > 0) {
        # Add missing columns with NA values
        df[missing_cols] <- NA
      }
      return(df)
    })
    
    # Merge into a single dataframe
    tax_lookup_df <- do.call(rbind, taxonomy_list)
    tax_lookup_df <- tax_lookup_df[!duplicated(tax_lookup_df), ] # make sure no duplicates
    
    # Merge these with the massblaster hits dataframe
    results_df_w_tax <-  left_join(results_df, tax_lookup_df, by = "Taxon")

  #================================================================
  # STEP 3: Select "the best" hit per query for massBLASTer results
  #================================================================
  
    # As we want to prefer lower level assignments we need to add a vector that captures
    # the level of classification
    # Lets add a column that summarises the taxonomic level to which a hit is classified:
    results_df_w_tax$ClassVec <- rep("Unclassified", times = nrow(results_df_w_tax))
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$kingdom))]  <- "Kingdom"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$phylum))]   <- "Phylum"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$class))]   <- "Class"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$order))]  <- "Order"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$family))] <- "Family"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$genus))] <- "Genus"
    results_df_w_tax$ClassVec[!(is.na(results_df_w_tax$species))] <- "Species"
    
    results_df_w_tax$ClassVec <- factor(results_df_w_tax$ClassVec, 
                                            levels = c("Species", "Genus", "Family", "Class", "Phylum", "Order", "Kingdom" ))
    
    # Lets filter out the assignments that are below the threshold we require:
    # e value -20 => from massblaster output we just get evalues of 0
    # query coverage > 80%
    # query length > 100bp
    # percentage identity >97.5 
    
    results_df_w_tax[,c("Evalue", "Prcnt", "Qend", "query_length")] <- lapply(results_df_w_tax[,c("Evalue", "Prcnt", "Qend", "query_length")], as.numeric ) # make columns numeric for filtering
    
    results_df_w_tax$query_cov <- (results_df_w_tax$Qend/results_df_w_tax$query_length)*100
      
    results_df_w_tax_f <- results_df_w_tax %>% 
      filter(query_length > 100, 
             query_cov >80,
             Evalue == 0 ,
             Prcnt > 97.5 )
    
    str(results_df_w_tax)
    str(results_df_w_tax_f)
    
    # Ensure we keep a value in for all queries by checking those dropped in filtering
    # and adding blank rows for them.
    
    df1 <- results_df_w_tax$query_id %>% unique() %>% length()
    df2 <- results_df_w_tax_f$query_id %>% unique() %>% length() 
    if (df1 != df2) {
      print("Warning: filtering removed some entire queries. Adding in blank classifications for these queries")
      diff_queries <- setdiff(unique(results_df_w_tax$query_id), unique(results_df_w_tax_f$query_id))
      # for each query add back into the main df with all columns NA with the exception of the query_id and taxon 
      # which we set to "No_blast_hits_passed_threshold"
      
      # build rows for missing queries
      blank_template <- results_df_w_tax_f[0, ] 
      add_back_rows <- lapply(diff_queries, function(q) {
        row <- blank_template[1, ]         # make one blank row
        row[,] <- NA                       # set everything to NA
        row$query_id <- q                  # fill query_id
        row$Taxon <- "No_blast_hits_passed_threshold"
        row
      }) %>% dplyr::bind_rows()
      
      # append back into filtered dataframe
      results_df_w_tax_f <- dplyr::bind_rows(results_df_w_tax_f, add_back_rows)
    }
    
    # Order the dataframe by highest
    # score then by lowest classification
    results_df_w_tax_f_ordered <- results_df_w_tax_f[
      order(results_df_w_tax_f$query_id,
            -as.numeric(results_df_w_tax_f$Score),
            results_df_w_tax_f$ClassVec),
    ]
    
  # Print all the taxons found for each query
  
  tally_hits_with_print <- function(df) {
    df %>%
      count(query_id, Taxon, name = "n_hits") %>%
      arrange(query_id, desc(n_hits), Taxon) %>%
      group_by(query_id) %>%
      # print one line per query (side-effect), then keep the rows
      mutate(
        .printed = {
          line <- paste0(query_id[1], ": ",
                         paste0(Taxon, "=", n_hits, collapse = "; "))
          message(line)
          NA
        }
      ) %>%
      select(-.printed) %>%
      slice_head(n=1) %>%
      dplyr::rename(MostFreqTaxon = Taxon) %>%
      ungroup()
  }

  # Choose the first line for each query

  top_tally <- tally_hits_with_print(results_df_w_tax_f_ordered)
  
  massblaster_per_query <- results_df_w_tax_f_ordered %>%
    group_by(query_id)%>%
    slice_head(n = 1) %>%
    ungroup()
 
  merged_approaches <- merge(top_tally, massblaster_per_query) %>% 
    select(query_id, MostFreqTaxon, Taxon) %>% dplyr::rename(HighestScoreTaxon = Taxon)

  print("\nMerged Approaches:\n")
  print(merged_approaches)
  print("\n proceeding with Highest Scoring Taxon\n")
  
  # Rename the massblaster cols to prevent column conflicts when merging:
  massblaster_per_query <- massblaster_per_query %>% dplyr::rename(OTU_string = query_id,
                                                            Kingdom = kingdom,
                                                            Phylum = phylum,
                                                            Class = class,
                                                            Order = order,
                                                            Family = family,
                                                            Genus = genus,
                                                            Species = species)
  massblaster_per_query$Species<- str_replace_all(massblaster_per_query$Species, " ", "_")
  
  #Output the dataframe
  return(massblaster_per_query)
}


 


  
  
