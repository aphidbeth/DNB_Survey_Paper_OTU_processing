# Function for getting a global alignment between two sequences from genbank

# Written by Beth Moore with assistance from chatgpt/copilot for sections
# 13/02/2026

# Input:
  # asc1 = genbank accession number of the first sequence
  # asc2 = genbank accession number of the second sequence
  # raw_seq = string of the raw sequence for the given OTU of interest
  # motif_left = "GCATCGATGAAGAACGCAGC"  # default to ITS2 => conserved end of 5.8S rRNA,
  # motif_right = "GACCTG" # default to ITS2 => conserved start of 28S rRNA 

# Output:
  # An alignment for visual comparison with BrowseSeqs()


prep_alignment <- function(asc1, asc2, raw_seq,
                                       motif_left = "GCATCGATGAAGAACGCAGC" ,
                                       motif_right = "GACCTG" ) {

  
  # List required packages
  required_packages <- c("ape", "DECIPHER")
  
  # Check if each package is installed, stop with error if not
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it first."))
    }
  }
  
  # Load the packages
  library(ape)       # for read.GenBank(), phylogenetic utilities
  library(DECIPHER)  # for DNAStringSet and sequence handling
  
  # Download sequences
  seq1 <- read.GenBank(asc1, as.character = T)
  seq2 <- read.GenBank(asc2, as.character = T)
  
  # Collapse into a single string
  
  seq1_flat <-unlist(seq1) %>%  paste0(., collapse = "")
  seq2_flat <-unlist(seq2) %>%  paste0(., collapse = "")
  
  full_seqs <- DNAStringSet(c(seq1_flat, seq2_flat, raw_seq))
  names(full_seqs) <- c("Seq1", "Seq2", "Raw_seq")
  
  # Extract the sequence are of interest

  # Find motifs
  left_hits  <- vmatchPattern(motif_left, full_seqs)
  right_hits <- vmatchPattern(motif_right, full_seqs)
  
  
  # Function to extract ITS2 between motifs
  # seq   : a Biostrings sequence (e.g., DNAString)
  # left  : match/range object for the left motif (e.g., output from matchPattern)
  # right : match/range object for the right motif
  extract_ITS2 <- function(seq, left, right) {
    # Helper to test presence of at least one match
    has_left  <- !is.null(left)  && length(left)  > 0
    has_right <- !is.null(right) && length(right) > 0
    
    if (has_left && has_right) {
      # Take the first occurrence of each motif
      s_idx <- end(left)[1] + 1
      e_idx <- start(right)[1] - 1
      
      # Guard against overlaps or invalid bounds
      if (s_idx > e_idx) {
        warning("Left and right motifs overlap or are in an unexpected order; returning full sequence.")
        return(seq)
      }
      return(subseq(seq, start = s_idx, end = e_idx))
      
    } else if (has_left && !has_right) {
      # Return from end of left motif to sequence end
      warning("End motif not found; returning sequence from left motif to end.")
      s_idx <- end(left)[1] + 1
      # Ensure start is within bounds
      s_idx <- max(1, min(s_idx, length(seq)))
      return(subseq(seq, start = s_idx))
      
    } else if (!has_left && has_right) {
      # Return from sequence start to start of right motif - 1
      warning("Start motif not found; returning sequence from start to right motif.")
      e_idx <- start(right)[1] - 1
      # Ensure end is within bounds
      e_idx <- max(1, min(e_idx, length(seq)))
      return(subseq(seq, end = e_idx))
      
    } else {
      # Neither motif found: return full sequence
      warning("Neither start nor end motifs found; returning full sequence.")
      return(seq)
    }
  }
  
  # Apply to all sequences
  its2_seqs <- DNAStringSet(mapply(extract_ITS2, full_seqs, left_hits, right_hits))
  names(its2_seqs) <- c(asc1, asc2, "raw_seq")
  
  # -----------------------------
  # Align ITS2 sequences
  # -----------------------------
  alignment <- AlignSeqs(its2_seqs, iterations = 10)  # pairwise alignment for two sequences
  
  return(alignment)
  
} 

