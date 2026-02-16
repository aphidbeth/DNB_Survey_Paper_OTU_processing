# Function for comparison of taxonomic conflicts between two approaches

# Written by Beth Moore with assistance from chatgpt/copilot for sections
# 11/02/2026

# Input: Two taxa tables with the following columns
  # id: OTU id
  # Kingdom
  # Phylum
  # Class
  # Order
  # Family
  # Species

# Output: A dataframe identifying agreement level between the two approaches

get_taxonomic_agreement <- function(df1, df2) {
  
  ranks <- c("Kingdom", "Phylum", "Class",
            "Order", "Family", "Genus", "Species")
  
  otu_col = "OTU_string"
  
  # ----Verify otu are same in both dataframes ----------------------------
  if (!(otu_col %in% names(df1))) stop("OTU column not found in df1")
  if (!(otu_col %in% names(df2))) stop("OTU column not found in df2")
  
  otu1 <- df1[[otu_col]]
  otu2 <- df2[[otu_col]] 
  
  if (!setequal(otu1, otu2)) {
    stop("df1 and df2 do not contain identical OTU_string values")
  }
  
  df2 <- df2[match(otu1, df2[[otu_col]]), ] # make sure they are in the same order
  
  # verify that order matches exactly now
  if (!all(df1[[otu_col]] == df2[[otu_col]])) {
    stop("Internal OTU_string ordering mismatch after alignment")
  }
  
  # ----  Check required rank columns exist ---
  if (!all(ranks %in% colnames(df1))) stop("Not all rank columns found in df1")
  if (!all(ranks %in% colnames(df2))) stop("Not all rank columns found in df2")
  
  # ---- Compare df1 and df2 across ranks -----

  normalize_chars <- function(x) {
    x <- trimws(x)
    x
  } # trims whitespace
  
  df1_cmp <- as.data.frame(lapply(df1[ranks], normalize_chars), stringsAsFactors = FALSE)
  df2_cmp <- as.data.frame(lapply(df2[ranks], normalize_chars), stringsAsFactors = FALSE)
  
  comp_df <- df1_cmp == df2_cmp
  
  # Convert NAs (from NA == value or NA == NA) into FALSE
  comp_df[is.na(comp_df)] <- FALSE
  
  # Ensure it's a matrix for apply()
  comparison_matrix <- as.matrix(comp_df)

  
  # ---- Determine deepest matching rank per row ----
  agreement <- apply(comparison_matrix, 1, function(x) {
    if (!any(x)) return(NA)
    ranks[max(which(x))]
  })
  
  # --- Determine the level of classification for each approach------
  
  # Deepest (rightmost) non-NA, non-empty rank
  get_class_level <- function(row) {
    # TRUE for filled ranks
    filled <- row != "" & !is.na(row)
    if (!any(filled)) return(NA)
    ranks[max(which(filled))]
  }
  
  get_class_assignment <- function(row) {
    # TRUE for filled ranks
    filled <- row != "" & !is.na(row)
    if (!any(filled)) return(NA)
    row[max(which(filled))]
  }
  
  
  df1_class_level <- apply(df1[ranks], 1, get_class_level)
  df2_class_level <- apply(df2[ranks], 1, get_class_level)
  
  df1_class_assign <- apply(df1[ranks], 1, get_class_assignment)
  df2_class_assign  <- apply(df2[ranks], 1, get_class_assignment)
  
  
  #--- Create a df of the results -----------
  df_out <- data.frame(OTU_string = df1[[otu_col]],
                       df1_class_level =  df1_class_level,
                       df2_class_level = df2_class_level,
                       df1_assignment = df1_class_assign,
                       df2_assignment = df2_class_assign,
                       agreement_level = agreement)
  
  #---- Assess if manual verification required ----
  
  rank_order <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
  df_out <- df_out %>%
    mutate(
      df1_class_level = factor(df1_class_level, levels = rank_order, ordered = TRUE),
      df2_class_level = factor(df2_class_level, levels = rank_order, ordered = TRUE),
      
      consensus_assignment = case_when(
        
        # 1️⃣ One level is NA → accept the other
        is.na(df1_class_level) & !is.na(df2_class_level) ~ df2_assignment,
        is.na(df2_class_level) & !is.na(df1_class_level) ~ df1_assignment,
        
        # 2️⃣ Same level & identical assignment → accept
        df1_class_level == df2_class_level &
          df1_assignment == df2_assignment ~ df1_assignment,
        
        # 3️⃣ Different levels → check if shallower = agreement
        df1_class_level != df2_class_level ~ {
          
          deeper_is_df1 <- as.numeric(df1_class_level) > as.numeric(df2_class_level)
          
          ifelse(
            deeper_is_df1,
            if_else(df2_class_level == agreement_level, df1_assignment, "Requires manual verification"),
            if_else(df1_class_level == agreement_level, df2_assignment, "Requires manual verification")
          )
        },
        
        # 4️⃣ Default
        TRUE ~ "Requires manual verification"
      )
    )
  
  return(df_out)
}



  