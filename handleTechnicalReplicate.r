library(Seurat)
library(Matrix)

# Helper function to sum two technical replicates
sum_technical_reps <- function(mat1, mat2) {
  # Find all unique barcodes across both runs
  all_barcodes <- union(colnames(mat1), colnames(mat2))
  all_genes <- rownames(mat1) # Assuming genes are the same
  
  # Create a template matrix of zeros
  # We use sparse matrices to save memory
  combined_mat <- Matrix(0, nrow = length(all_genes), ncol = length(all_barcodes), 
                         dimnames = list(all_genes, all_barcodes), sparse = TRUE)
  
  # Add counts from run 1 and run 2
  combined_mat[, colnames(mat1)] <- combined_mat[, colnames(mat1)] + mat1
  combined_mat[, colnames(mat2)] <- combined_mat[, colnames(mat2)] + mat2
  
  return(combined_mat)
}

# Apply to your data
hom_total.data <- sum_technical_reps(hom1.data, hom2.data)
het1_total.data <- sum_technical_reps(het1.data, het2.data)
het2_total.data <- sum_technical_reps(het02_1.data, het02_2.data)
wt_total.data <- sum_technical_reps(wt1.data, wt2.data)


# Create Seurat objects
# Create individual objects
hom_obj <- CreateSeuratObject(counts = hom_total.data, project = "Homozygous")
het1_obj <- CreateSeuratObject(counts = het1_total.data, project = "Heterozygous_1")
het2_obj <- CreateSeuratObject(counts = het2_total.data, project = "Heterozygous_2")
wt_obj <- CreateSeuratObject(counts = wt_total.data, project = "WT")

# Add metadata for genotype
hom_obj$genotype <- "HOM"
het1_obj$genotype <- "HET"
het2_obj$genotype <- "HET"
wt_obj$genotype <- "WT"