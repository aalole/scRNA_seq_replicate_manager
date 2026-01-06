# 1. Create Seurat objects with the 'min.cells' filter
hom_obj  <- CreateSeuratObject(counts = hom_total.data, project = "HOM", min.cells = 3)
het1_obj <- CreateSeuratObject(counts = het1_total.data, project = "HET1", min.cells = 3)
het2_obj <- CreateSeuratObject(counts = het2_total.data, project = "HET2", min.cells = 3)
wt_obj   <- CreateSeuratObject(counts = wt_total.data, project = "WT", min.cells = 3)

# 2. Add Mitochondrial % 
# Use "^MT-" for Human or "^mt-" for Mouse
hom_obj[["percent.mt"]]  <- PercentageFeatureSet(hom_obj, pattern = "^mt-")
het1_obj[["percent.mt"]] <- PercentageFeatureSet(het1_obj, pattern = "^mt-")
het2_obj[["percent.mt"]] <- PercentageFeatureSet(het2_obj, pattern = "^mt-")
wt_obj[["percent.mt"]]   <- PercentageFeatureSet(wt_obj, pattern = "^mt-")

# 3. Merge for visualization
# This allows us to see all samples side-by-side in one plot
raw_merged <- merge(wt_obj, y = c(het1_obj, het2_obj, hom_obj), 
                    add.cell.ids = c("WT", "HET1", "HET2", "HOM"))


##### QC Visualization #####


# Violin Plot: Distribution of Genes, UMIs, and MT%
VlnPlot(raw_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter Plot: Relationship between UMIs and MT%
# High UMIs with high MT% often indicates broken cells/low quality
plot1 <- FeatureScatter(raw_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


########################### Next Steps ###########################
# Apply filtering
filtered_merged <- subset(raw_merged, 
                          subset = nFeature_RNA > 100 & 
                                   percent.mt < 10)

# Print comparison of cell counts
message("Cells before filtering: ", ncol(raw_merged))
message("Cells after filtering: ", ncol(filtered_merged))



####################### hpost filrtering QC plot#########################
# Post-filter Violin Plot
VlnPlot(filtered_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Post-filter Scatter Plot
plot3 <- FeatureScatter(filtered_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(filtered_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4