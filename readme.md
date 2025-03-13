# LipMS Analysis Pipeline Documentation

## Overview

The LipMS Analysis Pipeline is a comprehensive R-based workflow for analyzing Limited Proteolysis Mass Spectrometry (LiP-MS) data. This pipeline processes raw MaxQuant data to identify significant changes in protein abundance and structural alterations through peptide-level analysis.



## Major content

- It processes LipMs data following the instruction from "Measuring protein structural changes on a proteome-wide scale using limited proteolysiscoupled mass spectrometry":

```
(iv) Perform two separate analyses. (i) Analyze control samples treated with trypsin only as in a standard quantitative proteomics pipeline to obtain protein abundance changes. (ii) Analyze LiP samples by treating every peptide as an independent entity to obtain peptide-level abundance changes. (v) Perform statistical testing to identify significantly changed peptides and proteins with the R package SafeQuant. This enables calculation of median abundance changes and associated P values corrected for multiple testing (q values). Abundance changes from Progenesis (or MaxQuant or alternative) can be imported into SafeQuant. (vi) Filter the protein abundance changes using suitable cutoffs. We recommend using a log2 abundance change cutoff of twofold and q values <0.01. (vii) Use significant protein abundance changes as normalization factors for peptide-level LiP data. This is achieved by dividing peptide-level abundance changes across samples by the (significant) abundance change of the respective protein. For proteins that do not significantly change abundance, use a normalization factor of 1 (i.e., no correction). (viii) Filter normalized peptide-level abundance changes from LiP data using suitable cutoffs. We recommend using a log2 abundance change cutoff of twofold and q values <0.01. (ix) Plot the results in the form of a volcano plot, representing peptide abundance changes versus the associated q values. Proteins associated with peptides that significantly change abundance in the LiP samples are deemed structurally variant. These peptides identify the specific protein region undergoing the structural change.
```

- make some revises to some script to better fit the pipeline, and modify the volcano script for better visualization. 
- Add some new functions for data processing.



## Usage Example

```r
# Load gene-protein mapping data
gene_protein_df <- fread("Input/Human_protein_gene_id.csv")
colnames(gene_protein_df) <- c("Accession", "GeneName")
gene_protein_df$`UniProtKB-ID` <- NULL

# Load and preprocess raw data
raw_peptide <- fread("Input/lipMs_test-Peng-0313.csv")
raw_peptide2 <- preprocess_raw_peptide(raw_peptide)

# Extract modifications information
modifications_df <- unique(
  raw_peptide2[,c("precursor", "Modifications")]
) %>% na.omit()

# Process data
peptides_list <- transfer_maxquant_data(raw_peptide2)

# Run the complete analysis
peptides_diff <- run_lip_ms_analysis(peptides_list,
                    modifications_df = modifications_df)
```

## 



directly run `run.Rmd` which majorly use function from `example_analysis_function.R`.



## Input Data Requirements

1. Peptide data with or without Proteinase K (PK) treatment
2. Peptide abundances data (typically labeled with "Abundances" or "Intensity")
3. Sequence and Protein Accessions information
4. Optional: Gene names and modification information



## Workflow

The LipMS analysis follows this workflow:

1. **Data Loading and Preprocessing**
   - Load raw MaxQuant data
   - Preprocess peptide data to handle duplicates
   - Convert to long format for analysis
   - Log2 transform intensities

2. **Quality Control**
   - Data completeness assessment
   - Principal Component Analysis (PCA)
   - Visualization of QC metrics

3. **Protein Abundance Calculation**
   - Calculate protein abundance from peptide data
   - Assign missingness for statistical analysis
   - Impute missing values

4. **Statistical Testing**
   - Perform statistical testing for proteins
   - Filter significant protein changes
   - Assign missingness for peptides
   - Perform statistical testing for peptides

5. **Normalization and Correction**
   - Normalize peptide-level LiP data using protein abundance changes
   - Correct for protein abundance changes

6. **Results Filtering and Visualization**
   - Filter significant peptide changes
   - Create volcano plots
   - Add regulation information (Up/Down/Not Significant)

7. **Results Export**
   - Save significant proteins and peptides to CSV files
   - Generate visualization outputs as PDF files

## Key Functions

- `preprocess_raw_peptide()`: Prepares raw peptide data for analysis
- `transfer_maxquant_data()`: Converts MaxQuant data to the required format
- `clean_peptides_na()`: Removes entries with all NA values
- `run_lip_ms_analysis()`: Main function that executes the complete workflow



## Parameters

The main analysis function `run_lip_ms_analysis()` accepts the following parameters:

- `peptides_list`: Required. List containing processed peptide data
- `output_dir`: Output directory for results (default: "Output/results")
- `log2FC_cutoff`: Log2 fold change cutoff for significance (default: 1)
- `p_value_cutoff`: P-value cutoff for significance (default: 0.01)
- `only_correct_sig_protein`: Whether to only correct significant proteins (default: TRUE)
- `normalize_method`: Method for normalization (default: "satterthwaite")
- `use_gene_id`: Whether to use gene IDs (default: TRUE)
- `modifications_df`: Optional dataframe with modification information

## Output Files

The pipeline generates the following output files in the specified output directory:

- `pk_qc_data_completeness.pdf`: QC plot for PK data completeness
- `ctrl_qc_data_completeness.pdf`: QC plot for control data completeness
- `pk_qc_pca.pdf`: PCA plot for PK data
- `ctrl_qc_pca.pdf`: PCA plot for control data
- `significant_proteins.csv`: List of significantly changed proteins
- `significant_peptides.csv`: List of significantly changed peptides after normalization
- `lip_ms_volcano_plot.pdf`: Volcano plot of peptide changes

