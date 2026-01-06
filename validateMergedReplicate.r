# Pick one genotype to test (e.g., WT)
total_sum_merged <- sum(wt_total.data)
total_sum_individual <- sum(wt1.data) + sum(wt2.data)

# This should return TRUE
print(total_sum_merged == total_sum_individual)

# Check cell counts
message("Cells in Run 1: ", ncol(wt1.data))
message("Cells in Run 2: ", ncol(wt2.data))
message("Cells in Merged: ", ncol(wt_total.data)) 
# Note: Merged should be <= (Run1 + Run2) because many barcodes are shared.

# Create temporary objects for the individual runs just for comparison
wt1_obj <- CreateSeuratObject(counts = wt1.data, project = "Run1")
wt2_obj <- CreateSeuratObject(counts = wt2.data, project = "Run2")
wt_merged_obj <- CreateSeuratObject(counts = wt_total.data, project = "Merged")

# Combine for plotting
comparison <- merge(wt1_obj, y = c(wt2_obj, wt_merged_obj))

# Visualize the improvement in depth
VlnPlot(comparison, features = c("nCount_RNA", "nFeature_RNA"), group.by = "orig.ident")

# Find a common barcode
common <- intersect(colnames(wt1.data), colnames(wt2.data))
cell_id <- common[1] # Look at the first shared cell

# Plot Gene Counts: Run 1 vs Run 2
plot(as.numeric(wt1.data[, cell_id]), 
     as.numeric(wt2.data[, cell_id]),
     main = paste("Correlation for Cell:", cell_id),
     xlab = "Run 1 Counts", ylab = "Run 2 Counts",
     pch = 16, col = rgb(0,0,0,0.2))
abline(0, 1, col = "red") # Perfect 1:1 line