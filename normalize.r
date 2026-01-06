# 1. Run SCTransform
# We regress out 'percent.mt' to ensure clusters are driven by biology, not cell stress.
# This step may take a few minutes.
obj <- SCTransform(filtered_merged, 
                   vars.to.regress = "percent.mt", 
                   verbose = FALSE)

# 2. Run PCA (Principal Component Analysis)
# PCA compresses the 3,000 variable genes into "dimensions" of variance.
obj <- RunPCA(obj, verbose = FALSE)


# Visualizing the PCA 'Elbow'
ElbowPlot(obj, ndims = 50)

# 3. Run UMAP and find neighbors
# We use 30 dims based on the Elbow Plot
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

# 4. Plot by Genotype
DimPlot(obj, reduction = "umap", group.by = "orig.ident") + 
  ggtitle("UMAP by Genotype")

# 5. Plot by Cluster
DimPlot(obj, reduction = "umap", label = TRUE) + 
  ggtitle("Initial Clustering")