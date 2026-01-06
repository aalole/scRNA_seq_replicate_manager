# scRNA_seq_replicate_manager

This code will efficiently handle merging two biological replicate RNA feature barcode matrices into 1 single matrices for n samples, resolving the pain point of replicate-induced batch defect.

# Step 1: Align and Sum the Matrices

Because we are using .h5 files from Cell Ranger, the gene lists should be identical. However, the number of barcodes (cells) might vary slightly between runs if the sequencing depth allowed Cell Ranger to "call" a cell in one run but not the other.

# Step 2: Create Seurat Objects

Once you have the summed matrices, you can create your Seurat objects for each genotype.
