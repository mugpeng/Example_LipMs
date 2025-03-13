# Source the main pipeline functions
#!/usr/bin/env Rscript

# LiP-MS Analysis Pipeline
# This script performs statistical analysis of LiP-MS data following the specified workflow:
# 1. Load and preprocess MaxQuant data
# 2. Perform statistical testing to identify significantly changed peptides and proteins
# 3. Filter protein abundance changes using suitable cutoffs
# 4. Normalize peptide-level LiP data using protein abundance changes
# 5. Filter normalized peptide-level abundance changes
# 6. Visualize results with volcano plots

# Set parameters
log2FC_cutoff <- 1  # Log2 fold change cutoff (4-fold)
# q_value_cutoff <- 0.01  # q-value cutoff
p_value_cutoff <- .01

# Preprocess raw_peptide
preprocess_raw_peptide <- function(raw_peptide){
  message("Preprocess MaxQuant data...")
  raw_peptide$Accession_seq <- paste0(raw_peptide$`Master Protein Accessions`, "-", 
                                      raw_peptide$Sequence)
  # Count occurrences of each ID
  id_counts <- table(raw_peptide$Accession_seq)
  
  # Get IDs that appear more than once
  duplicate_ids <- names(id_counts[id_counts > 1])
  
  # Create a new vector to store modified IDs
  modified_sequence <- raw_peptide$Sequence
  
  # Process only the duplicated IDs
  for (dup_id in duplicate_ids) {
    # Find positions of this duplicate ID
    positions <- which(raw_peptide$Accession_seq == dup_id)
    seq <- raw_peptide[raw_peptide$Accession_seq %in% dup_id,]$Sequence %>% unique()
    # Add suffixes (_2, _3, etc.) to all but the first occurrence
    for (i in 1:length(positions)) {
      modified_sequence[positions[i]] <- paste0(seq, "_", i)
    }
  }
  raw_peptide$precursor <- modified_sequence
  raw_peptide
}

# Function to load MaxQuant data
transfer_maxquant_data <- function(raw_peptide) {
  message("Transfer MaxQuant data...")
  # preprocess pipetide data
  peptides <- raw_peptide %>%
    # Convert to long format for analysis
    pivot_longer(
      cols = contains("Abundances "),
      names_to = "sample",
      values_to = "intensity"
    ) %>%
    # Clean sample names
    mutate(
      sample = gsub("Abundances (Normalized): ", "", sample, fixed = T),
      # Log2 transform intensities
      intensity_log2 = log2(intensity),
      Accession = `Master Protein Accessions`
    ) %>% select(c(Accession, Sequence, precursor, Accession_seq,
                   sample, intensity, intensity_log2))
  # Extract condition information from sample names
  # Assuming sample names follow a pattern like "Condition_Replicate"
  peptides <- peptides %>%
    mutate(condition = gsub("^F\\d+:\\s*(PK-)?(CTR|1\\.5|5).*$", "\\2", sample))
  if(any(grepl("PK", peptides$sample))){
    peptides_pk = peptides[grepl("PK", peptides$sample),]
    peptides_tripson = peptides[!grepl("PK", peptides$sample),]
    return(list(peptides_pk = peptides_pk, peptides_tripson = peptides_tripson))
  }
  return(peptides)
}

# Clean all NA precursor/Sequence
clean_peptides_na <- function(peptides){
  peptides_split <- split(peptides, peptides$precursor)
  a <- lapply(peptides_split, function(x){
    all(is.na(x[["intensity_log2"]]))
  }) 
  a <- do.call(c,a)
  peptides[peptides$precursor %in% names(a)[!a],]
}

# Example usage with custom parameters
run_lip_ms_analysis <- function(peptides_list, 
                                output_dir = "Output/results",
                                log2FC_cutoff = 1,
                                p_value_cutoff = .01,
                                only_correct_sig_protein = T,
                                normalize_method = "satterthwaite",
                                use_gene_id = T,
                                modifications_df = NULL) {
  # Validate inputs
  if(missing(peptides_list)) {
    stop("peptides_list is required")
  }

  # Define output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Load data
  message("Performing data clean...")
  peptides_pk <- peptides_list$peptides_pk
  peptides_ctrl <- peptides_list$peptides_tripson
  
  # Clean peptides
  peptides_pk_clean <- clean_peptides_na(peptides_pk)
  peptides_ctrl_clean <- clean_peptides_na(peptides_ctrl)
  
  # Quality control
  message("Performing quality control...")
  
  # Data completeness QC using our custom function with condition coloring
  qc_completeness_pk <- qc_data_completeness(
    data = peptides_pk_clean,
    sample = sample,
    grouping = Sequence,
    intensity = intensity_log2,
    condition = condition
  )
  ggsave(file.path(output_dir, "pk_qc_data_completeness.pdf"), qc_completeness_pk, width = 10, height = 8)
  qc_completeness_ctrl <- qc_data_completeness(
    data = peptides_ctrl_clean,
    sample = sample,
    grouping = Sequence,
    intensity = intensity_log2,
    condition = condition
  )
  ggsave(file.path(output_dir, "ctrl_qc_data_completeness.pdf"), qc_completeness_ctrl, width = 10, height = 8)
  
  # PCA analysis
  qc_pca_plot_pk <- qc_pca(
    data = peptides_pk_clean,
    sample = sample,
    grouping = Sequence,
    intensity = intensity_log2,
    condition = condition
  )
  ggsave(file.path(output_dir, "pk_qc_pca.pdf"), qc_pca_plot_pk, width = 10, height = 8)
  qc_pca_plot_ctrl <- qc_pca(
    data = peptides_ctrl_clean,
    sample = sample,
    grouping = Sequence,
    intensity = intensity_log2,
    condition = condition
  )
  ggsave(file.path(output_dir, "ctrl_qc_pca.pdf"), qc_pca_plot_ctrl, width = 10, height = 8)
  
  # Calculate protein abundance for ctrl
  message("Calculating protein abundance...")
  protein_abundance <- calculate_protein_abundance(
    data = peptides_ctrl_clean,
    sample = sample,
    protein_id  = Accession,
    intensity_log2 = intensity_log2,
    min_n_peptides = 1,
    # condition = condition,
    peptide = Sequence,
    precursor = precursor,
    method = "sum",  # Can also use "mean", "median", or "top_n"
    retain_columns = c(condition)
  )
  
  # Assign missingness for proteins for statistical analysis
  message("Assigning missingness for protein...")
  proteins_with_missingness <- assign_missingness(
    data = protein_abundance,
    sample = sample,
    condition = condition,
    grouping = Accession,
    intensity = intensity_log2,
    completeness_MAR = 0.7,
    completeness_MNAR = 0.5,
    ref_condition = "CTR",
  )
  
  # Impute missingness proteins data
  proteins_imputed <- impute(
    data = proteins_with_missingness,
    sample = sample,
    grouping = Accession,
    intensity_log2 = intensity_log2,
    condition = condition,
    comparison = comparison,
    missingness = missingness,
    method = "ludovic",
    # retain_columns = c(protein, peptide_intensity)
  )
  
  # Perform statistical testing for proteins
  message("Performing statistical testing for proteins...")
  protein_diff <- calculate_diff_abundance(
    data = proteins_imputed,
    sample = sample,
    condition = condition,
    grouping = Accession,
    intensity_log2 = imputed_intensity,
    missingness = missingness,
    comparison = comparison,
    method = "t-test"  # Can also use "moderated_t-test" or "proDA"
  )
  
  # Filter significant protein changes
  message("Filtering significant protein changes...")
  significant_proteins <- protein_diff %>%
    filter(abs(diff) >= log2FC_cutoff, pval < p_value_cutoff)
  
  # Save significant proteins
  write_csv(significant_proteins, file.path(output_dir, "significant_proteins.csv"))
  
  # Assign missingness for peptides for statistical analysis
  message("Assigning missingness for peptides...")
  peptides_with_missingness <- assign_missingness(
    data = peptides_pk_clean,
    sample = sample,
    condition = condition,
    grouping = precursor,
    intensity = intensity_log2,
    ref_condition = "CTR",
    completeness_MAR = 0.7,
    completeness_MNAR = 0.5,
    retain_columns = c("Accession")
  )
  
  # Impute missingness peptides data
  peptides_imputed <- impute(
    data = peptides_with_missingness,
    sample = sample,
    grouping = precursor,
    intensity_log2 = intensity_log2,
    condition = condition,
    comparison = comparison,
    missingness = missingness,
    method = "ludovic",
    retain_columns = c("Accession")
  )
  
  # Perform statistical testing for peptides
  message("Performing statistical testing for peptides...")
  peptide_diff <- calculate_diff_abundance(
    data = peptides_imputed,
    sample = sample,
    condition = condition,
    grouping = precursor,
    intensity_log2 = intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "t-test",
    retain_columns = c("Accession")
  )
  
  # Normalize peptide-level LiP data using protein abundance changes
  message("Normalizing peptide-level LiP data...")

  # Correct LiP data for protein abundance changes using the original method
  # Only correct significant proteins
  peptide_diff_normalized <- correct_lip_for_abundance(
    lip_data = peptide_diff,
    trp_data = protein_diff,
    protein_id = Accession,
    grouping = precursor,
    comparison = comparison,
    diff = diff,
    std_error = std_error,
    only_correct_sig_protein = only_correct_sig_protein,
    # retain_columns = c("missingness"),
    method = "satterthwaite"
  )
  peptide_diff_normalized$Accession_seq_id <- paste0(peptide_diff_normalized$Accession, "_", peptide_diff_normalized$precursor)
  
  # Add modifications info
  if(!is.null(modifications_df)){
    peptide_diff_normalized <- merge(peptide_diff_normalized,
                                 modifications_df,
                                 by = "precursor", all.x = T)
    peptide_diff_normalized$Accession_seq_id <- paste0(
      peptide_diff_normalized$Accession_seq_id, "_", peptide_diff_normalized$Modifications
    )
    peptide_diff_normalized$Accession_seq_id <- gsub("_NA", "", peptide_diff_normalized$Accession_seq_id)
  }
  
  # Add Gene id
  if(use_gene_id){
    peptide_diff_normalized <- merge(peptide_diff_normalized,
                                 gene_protein_df,
                                 by = "Accession", all.x = T)
    peptide_diff_normalized <- peptide_diff_normalized %>%
      mutate(Accession_seq_id = case_when(
        !is.na(GeneName) ~ str_replace(Accession_seq_id, Accession, GeneName),
        TRUE ~ Accession_seq_id
      ))
  }
  
  # Filter normalized peptide-level abundance changes
  message("Filtering significant peptide changes...")
  significant_peptides <- peptide_diff_normalized %>%
    filter(abs(adj_diff) >= log2FC_cutoff, pval < p_value_cutoff)
  
  # Save significant peptides
  write_csv(significant_peptides, file.path(output_dir, "significant_peptides.csv"))
  
  # Create volcano plot for peptides
  message("Creating volcano plot...")
  volcano <- volcano_plot(
    data = peptide_diff_normalized,
    grouping = Accession_seq_id,
    log2FC = adj_diff,
    significance = pval,
    method = "significant",
    split_by = comparison,
    plot_ncol = 2,
    significance_cutoff = p_value_cutoff,
    log2FC_cutoff = log2FC_cutoff,
    show_counts = T,
    label_top_n = 5
  )
  # Save volcano plot
  ggsave(file.path(output_dir, "lip_ms_volcano_plot.pdf"), volcano, width = 10, height = 8)
  
  # Add Group info
  peptide_diff_normalized <- peptide_diff_normalized %>%
    dplyr::mutate(regulation = dplyr::case_when(
      (adj_diff > log2FC_cutoff) & (pval < p_value_cutoff) ~ "Up",
      (adj_diff < -log2FC_cutoff) & (pval < p_value_cutoff) ~ "Down",
      TRUE ~ "Not Significant"
    ))
  
  # Summarize results
  message("Analysis complete!")
  message(paste("Found", nrow(significant_proteins), "significantly changed proteins"))
  message(paste("Found", nrow(significant_peptides), "significantly changed peptides after normalization"))
  message("Results saved to the output directory:", output_dir)
  return(peptide_diff_normalized)
}

