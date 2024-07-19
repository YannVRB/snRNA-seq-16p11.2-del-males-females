library(Seurat)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(writexl)
library(DESeq2)

#Load the CellRanger outputs of individual sample
WTmale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample795/outs/filtered_feature_bc_matrix")
WTmale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample809/outs/filtered_feature_bc_matrix")
hetmale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample796/outs/filtered_feature_bc_matrix")
hetmale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample807/outs/filtered_feature_bc_matrix")
WTfemale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample798/outs/filtered_feature_bc_matrix")
WTfemale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample803/outs/filtered_feature_bc_matrix")
hetfemale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample806/outs/filtered_feature_bc_matrix")
hetfemale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample811/outs/filtered_feature_bc_matrix")

#Initialize the Seurat object with the raw (non-normalized data).
WTmale1 <- CreateSeuratObject(counts = WTmale1.data, project = "WT_male", min.cells = 3, min.features = 200)
WTmale2 <- CreateSeuratObject(counts = WTmale2.data, project = "WT_male", min.cells = 3, min.features = 200)
hetmale1 <- CreateSeuratObject(counts = hetmale1.data, project = "WT_male", min.cells = 3, min.features = 200)
hetmale2 <- CreateSeuratObject(counts = hetmale2.data, project = "WT_male", min.cells = 3, min.features = 200)
WTfemale1 <- CreateSeuratObject(counts = WTfemale1.data, project = "WT_female", min.cells = 3, min.features = 200)
WTfemale2 <- CreateSeuratObject(counts = WTfemale2.data, project = "WT_female", min.cells = 3, min.features = 200)
hetfemale1 <- CreateSeuratObject(counts = hetfemale1.data, project = "het_female", min.cells = 3, min.features = 200)
hetfemale2 <- CreateSeuratObject(counts = hetfemale2.data, project = "het_female", min.cells = 3, min.features = 200)

#Standard pre-processing workflow
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats.
WTmale1[["percent.mt"]] <- PercentageFeatureSet(WTmale1, pattern = "^mt-")
WTmale2[["percent.mt"]] <- PercentageFeatureSet(WTmale2, pattern = "^mt-")
hetmale1[["percent.mt"]] <- PercentageFeatureSet(hetmale1, pattern = "^mt-")
hetmale2[["percent.mt"]] <- PercentageFeatureSet(hetmale2, pattern = "^mt-")
WTfemale1[["percent.mt"]] <- PercentageFeatureSet(WTfemale1, pattern = "^mt-")
WTfemale2[["percent.mt"]] <- PercentageFeatureSet(WTfemale2, pattern = "^mt-")
hetfemale1[["percent.mt"]] <- PercentageFeatureSet(hetfemale1, pattern = "^mt-")
hetfemale2[["percent.mt"]] <- PercentageFeatureSet(hetfemale2, pattern = "^mt-")

#Visualize QC metrics as a violin plot
VlnPlot(WTmale1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(WTmale2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hetmale1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hetmale2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(WTfemale1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(WTfemale2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hetfemale1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hetfemale2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Use QC metrics to filter cells.
WTmale1 <- subset(WTmale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTmale2 <- subset(WTmale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetmale1 <- subset(hetmale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetmale2 <- subset(hetmale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTfemale1 <- subset(WTfemale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTfemale2 <- subset(WTfemale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetfemale1 <- subset(hetfemale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetfemale2 <- subset(hetfemale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Setup conditions in metadata "dataset". This will help us later when all the samples are integreated into one Seurat object.
WTmale1$dataset <- "WTmale"
WTmale2$dataset <- "WTmale"
WTfemale1$dataset <- "WTfemale"
WTfemale2$dataset <- "WTfemale"
hetmale1$dataset <- "hetmale"
hetmale2$dataset <- "hetmale"
hetfemale1$dataset <- "hetfemale"
hetfemale2$dataset <- "hetfemale"

#Normalizing the data.
WTmale1 <- NormalizeData(WTmale1, verbose = FALSE)
WTmale2 <- NormalizeData(WTmale2, verbose = FALSE)
hetmale1 <- NormalizeData(hetmale1, verbose = FALSE)
hetmale2 <- NormalizeData(hetmale2, verbose = FALSE)
WTfemale1 <- NormalizeData(WTfemale1, verbose = FALSE)
WTfemale2 <- NormalizeData(WTfemale2, verbose = FALSE)
hetfemale1 <- NormalizeData(hetfemale1, verbose = FALSE)
hetfemale2 <- NormalizeData(hetfemale2, verbose = FALSE)

#Identification of highly variable features (feature selection).
WTmale1 <- FindVariableFeatures(WTmale1, selection.method = "vst", nfeatures = 2000)
WTmale2 <- FindVariableFeatures(WTmale2, selection.method = "vst", nfeatures = 2000)
hetmale1 <- FindVariableFeatures(hetmale1, selection.method = "vst", nfeatures = 2000)
hetmale2 <- FindVariableFeatures(hetmale2, selection.method = "vst", nfeatures = 2000)
WTfemale1 <- FindVariableFeatures(WTfemale1, selection.method = "vst", nfeatures = 2000)
WTfemale2 <- FindVariableFeatures(WTfemale2, selection.method = "vst", nfeatures = 2000)
hetfemale1 <- FindVariableFeatures(hetfemale1, selection.method = "vst", nfeatures = 2000)
hetfemale2 <- FindVariableFeatures(hetfemale2, selection.method = "vst", nfeatures = 2000)

#Integration.
#Find integration anchors.
anchors <- FindIntegrationAnchors(object.list = list(WTmale1, WTmale2, WTfemale1, WTfemale2, hetmale1, hetmale2, hetfemale1, hetfemale2), dims = 1:20)
#Perform the integration
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

#Scaling the data.
combined <- ScaleData(combined, verbose = FALSE)

#Perform linear dimensional reduction.
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

#Cluster the cells.
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP).
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

#Create a new UMAP using the integrated assay.
DimPlot(combined, reduction = "umap", assay = "integrated")

#Create a new UMAP using the integrated assay and split between the 4 datasets.
DimPlot(combined, reduction = "umap", split.by = "dataset", assay = "integrated")

#Compute celltype-specific biomarkers in each individual cluster against the rest of the cells then compare top significant biomarkers against litterature.
cluster1.markers <- FindMarkers(combined, ident.1 = "1")

#Dot plot (bubble plot) for celltype-specific biomarkers.
Idents(combined) <- factor(Idents(combined), levels = c("iSPNs", "dSPNs", "eSPNs", "Pax6", "Lamp5", "Vip", "Sst_Chodl", "Sst_Pvalb", "L2/3 IT CTX", "L4 IT CTX", "L5 IT CTX", "L6 IT CTX", "L2/3 IT ENTI", "L5 PT CTX", "Car3", "L5/6 NP CTX", "L6 CT CTX", "L6b CT CTX", "Endothelial", "Oligodendrocytes", "Cholinergic interneurons", "Astrocytes"))
markers.to.plot <- c("Penk", "Tac1", "Casz1", "Gabrg1", "Kit", "Vip", "Nos1", "Pvalb", "Stard8", "Adam33", "Dkkl1", "C1ql3", "Ndst4", "Hcn1", "Car3", "Fgd6", "Col5a1", "Ctgf", "Golim4", "Ptgds", "Ecel1", "Slc1a3")
DotPlot(combined, features = markers.to.plot, cols = c("lightblue", "orange", "darkblue", "darkorange"), dot.scale = 8, split.by = "dataset") +  RotatedAxis()

#Setup Seurat object fo differential analyses.
combined$celltype.dataset <- paste(Idents(combined), combined$dataset, sep = "_")
combined$celltype <- Idents(combined)
Idents(combined) <- "celltype.dataset"

#Identify differential expressed genes between 16p11.2 del/+ male and wildtype male.
dSPNs.DEGs.male <- FindMarkers(combined, ident.1 = "dSPNs_hetmale", ident.2 = "dSPNs_WTmale", verbose = FALSE)


### Heatmap of genes inside the 16p11.2 deletion
library(ComplexHeatmap)
library(circlize)

# Create a vector of genes located inside the 16p11.2 deletion
genes_of_interest <- c("Coro1a", "Mapk3", "Gm9967", "Ypel3", "Tbx6", "Ppp4c", "Aldoa", "Fam57b",
                       "Gm15676", "Doc2a", "Ino80e", "Hirip3", "Taok2", "Tmem219", "Kctd13",
                       "Gm21984", "Asphd1", "Sez6l2", "D830044I16Rik", "Cdipt", "Mvp", "Prrt2",
                       "Gm42742", "Maz", "Kif22", "RP23-142A14.9", "AI467606", "Qprt", "Cd2bp2",
                       "Tbc1d10b", "Mylpf", "Spn")

# Filter dataframes to include only the genes of interest
dSPNsmale_filtered <- dSPNsmale[rownames(dSPNs.DEGs.male) %in% genes_of_interest, , drop = FALSE]
iSPNsmale_filtered <- iSPNsmale[rownames(iSPNs.DEGs.male) %in% genes_of_interest, , drop = FALSE]
dSPNsfemale_filtered <- dSPNsfemale[rownames(dSPNs.DEGs.female) %in% genes_of_interest, , drop = FALSE]
iSPNsfemale_filtered <- iSPNsfemale[rownames(iSPNs.DEGs.female) %in% genes_of_interest, , drop = FALSE]

# Ensure the rows are in the same order for each dataframe
dSPNsmale_ordered <- dSPNsmale_filtered[match(genes_of_interest, rownames(dSPNsmale_filtered)), , drop = FALSE]
iSPNsmale_ordered <- iSPNsmale_filtered[match(genes_of_interest, rownames(iSPNsmale_filtered)), , drop = FALSE]
dSPNsfemale_ordered <- dSPNsfemale_filtered[match(genes_of_interest, rownames(dSPNsfemale_filtered)), , drop = FALSE]
iSPNsfemale_ordered <- iSPNsfemale_filtered[match(genes_of_interest, rownames(iSPNsfemale_filtered)), , drop = FALSE]

# Combine the ordered dataframes into a single dataframe
combined_df <- cbind(dSPNsmale_ordered[, "avg_log2FC", drop = FALSE],
                     iSPNsmale_ordered[, "avg_log2FC", drop = FALSE],
                     dSPNsfemale_ordered[, "avg_log2FC", drop = FALSE],
                     iSPNsfemale_ordered[, "avg_log2FC", drop = FALSE])

# Rename columns for clarity
colnames(combined_df) <- c("dSPNsmale", "iSPNsmale", "dSPNsfemale", "iSPNsfemale")

# Define the rownames to be removed
rows_to_remove <- c("NA", "NA.1", "NA.2", "NA.3", "NA.4", "NA.5", "NA.6", "NA.7")

# Remove the specified rows
combined_df <- combined_df[!rownames(combined_df) %in% rows_to_remove, ]


# Convert the combined dataframe to a matrix
combined_matrix <- as.matrix(combined_df)

# Define the color scale
col_fun <- colorRamp2(c(min(combined_matrix), max(combined_matrix)), c("blue", "white"))

# Sort the matrix rows based on the minimum value in each row
sorted_indices <- order(apply(combined_matrix, 1, min))
sorted_matrix <- combined_matrix[sorted_indices, ]

# Create the heatmap with the specified color gradient and clustering
heatmap <- Heatmap(sorted_matrix, name = "avg_logFC",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_row_names = TRUE, show_column_names = TRUE,
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 8),
                   col = col_fun)

draw(heatmap, merge_legends = TRUE)

# Save as PDF
pdf("heatmap.pdf", width = 3, height = 4)  # Adjust width and height as needed
draw(heatmap, merge_legends = TRUE)
dev.off()


## Pseudobulk analysis of D1 and D2 SPNs for quadrant plot
combined$ClusterNames <- combined@active.ident
bulk <- AggregateExpression(combined, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("ClusterNames", "sample")
# Get the cell names
cell_names <- Cells(bulk)

# Identify cells that contain "D2-SPNs" in their name
d1_spns_cells <- grep("eSPNs", cell_names, value = TRUE)


# Subset the Seurat object to keep only these cells
bulk_d1_spns <- subset(bulk, cells = d1_spns_cells)

# Extract the count data
count_data <- GetAssayData(bulk_d1_spns, slot = "counts")

# Create metadata dataframe
metadata <- data.frame(cell = colnames(count_data))
metadata$condition <- ifelse(grepl("WTfemale", metadata$cell), "WTfemale", 
                             ifelse(grepl("hetfemale", metadata$cell), "hetfemale", 
                                    ifelse(grepl("WTmale", metadata$cell), "WTmale", "hetmale")))

# Ensure that cell names are the row names of metadata
row.names(metadata) <- metadata$cell
metadata$cell <- NULL

# Ensure that the count data columns and metadata rows match
if (all(colnames(count_data) == rownames(metadata))) {
  # Create the DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
  
  # Filter low count genes (optional but recommended)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run the DESeq2 analysis
  dds <- DESeq(dds)
  
  # Define the comparisons
  comparisons <- list(
    c("hetmale", "WTmale"),
    c("hetfemale", "WTfemale"),
    c("WTfemale", "WTmale")
  )
  
  # Define the file path
  output_file <- "DESeq2_results_comparisons.xlsx"
  
  # Initialize the workbook
  wb <- createWorkbook()
  
  for (comparison in comparisons) {
    # Extract the results for the current comparison
    res <- results(dds, contrast = c("condition", comparison[1], comparison[2]))
    
    # Order results by adjusted p-value
    res <- res[order(res$padj), ]
    
    # Convert DESeq2 results to a data frame
    res_df <- as.data.frame(res)
    
    # Add a new sheet with the comparison name
    sheet_name <- paste(comparison[1], "vs", comparison[2], sep = "_")
    addWorksheet(wb, sheet_name)
    
    # Write the data frame to the worksheet
    writeData(wb, sheet = sheet_name, x = res_df, rowNames = TRUE)
  }
  
  # Save the workbook to the file
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("Results have been successfully exported to", output_file)
} else {
  cat("The count data columns and metadata rows do not match.")
}

# Define the path to the Excel file
input_file <- "DESeq2_results_comparisons.xlsx"

# Read the results for the specified comparisons
results_WTfemale_vs_WTmale <- read_excel(input_file, sheet = "WTfemale_vs_WTmale")
results_hetmale_vs_WTmale <- read_excel(input_file, sheet = "hetmale_vs_WTmale")

# Ensure the rownames are the gene names
rownames(results_WTfemale_vs_WTmale) <- results_WTfemale_vs_WTmale$...1
rownames(results_hetmale_vs_WTmale) <- results_hetmale_vs_WTmale$...1

# Extract log2 fold changes
log2FC_WTfemale_vs_WTmale <- results_WTfemale_vs_WTmale$log2FoldChange
log2FC_hetmale_vs_WTmale <- results_hetmale_vs_WTmale$log2FoldChange

# Create a data frame with the log2 fold changes
fold_changes <- data.frame(
  gene = rownames(results_WTfemale_vs_WTmale),
  log2FC_WTfemale_vs_WTmale = log2FC_WTfemale_vs_WTmale,
  log2FC_hetmale_vs_WTmale = log2FC_hetmale_vs_WTmale
)

# List of specified genes
specified_genes_d1 <- c(
  "Rnf121", "Gm26825", "Malat1", "Kcnq1ot1", "Gm10600", "Grin2b", "Gm42418", "Nisch", "Ank3", "Ttc3", "Trank1", 
  "Ddx5", "Crocc", "Ssbp4", "Atp6v0b", "Aldoa", "Gm48715", "Cdk9", "Gm26724", "Rab11fip2", "Meg3", "Ppp3ca", "Eif1", 
  "Gsdme", "Elavl3", "Rbm39", "Dock3", "Kcnq2", "Plekha6", "Syt7", "Gm16105", "Sf3b3", "Cacnb4", "Tspoap1", 
  "Sgcz", "Gabarapl2", "Grik2", "AI314180", "Ypel3", "Kif1b", "Cacna1a", "Fam193b", "Camk2b", "Ryr3", "Ppargc1b", 
  "Tmem181a", "Grin1", "Acp1", "Srcin1", "Golga4", "Negr1", "Tmem131", "Gm26699", "Kmt2e", "Atp9b", "Rpl21", 
  "Klhl33", "Pcdhgc4", "Gm27032", "Ttc19", "Hirip3", "Sptbn4", "Ahi1", "Lrrc7", "Gm16152", "Tshz1", 
  "Il1rapl1", "Zfp445", "Pebp1", "Gm10848", "Tra2a", "Actn1", "Aste1", "Timp4", "Ppfia4", "Zfp846", "Emsy", "Pitpnm1", 
  "Zfp644", "Trim35", "Gm42722", "Adcy5", "4833420G17Rik", "Cpne5", "Pcsk2", "Mpped2", "Prrc2b", "Prrc2a", "Adgrl1", 
  "Zc3h11a", "Snhg11", "Clip4", "Strn", "Safb2", "Ankrd17", "AY036118", "Unc13a", "Mirg"
)

# Subset the fold changes data frame
subset_fold_changes <- fold_changes[fold_changes$gene %in% specified_genes_d1, ]

# Scatter plot
p1 <- ggplot(subset_fold_changes, aes(x = log2FC_WTfemale_vs_WTmale, y = log2FC_hetmale_vs_WTmale)) +
  geom_point() +
  geom_text_repel(aes(label = gene)) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  stat_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson", label.x = -4, label.y = 4) +
  labs(title = "D1-SPN DEGs between 16p males vs wt females",
       x = "Pseudobulk log2 fold change (WTfemale vs WTmale)",
       y = "Pseudobulk log2 fold change (hetmale vs WTmale)") +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Save the plot as a PDF
ggsave("D1-SPN DEGs between 16p males vs wt females.pdf", plot = p1, width = 8, height = 6)
# Save the subset_fold_changes dataframe to an Excel file
write_xlsx(subset_fold_changes, "D1-SPN DEGs between 16p males vs wt females subset_fold_changes.xlsx")

# Identify cells that contain "D2-SPNs" in their name
d2_spns_cells <- grep("D2-SPNs", cell_names, value = TRUE)


# Subset the Seurat object to keep only these cells
bulk_d2_spns <- subset(bulk, cells = d2_spns_cells)

# Extract the count data
count_data <- GetAssayData(bulk_d2_spns, slot = "counts")

# Create metadata dataframe
metadata <- data.frame(cell = colnames(count_data))
metadata$condition <- ifelse(grepl("WTfemale", metadata$cell), "WTfemale", 
                             ifelse(grepl("hetfemale", metadata$cell), "hetfemale", 
                                    ifelse(grepl("WTmale", metadata$cell), "WTmale", "hetmale")))

# Ensure that cell names are the row names of metadata
row.names(metadata) <- metadata$cell
metadata$cell <- NULL

# Ensure that the count data columns and metadata rows match
if (all(colnames(count_data) == rownames(metadata))) {
  # Create the DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
  
  # Filter low count genes (optional but recommended)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run the DESeq2 analysis
  dds <- DESeq(dds)
  
  # Define the comparisons
  comparisons <- list(
    c("hetmale", "WTmale"),
    c("hetfemale", "WTfemale"),
    c("WTfemale", "WTmale")
  )
  
  # Define the file path
  output_file <- "DESeq2_results_comparisons.xlsx"
  
  # Initialize the workbook
  wb <- createWorkbook()
  
  for (comparison in comparisons) {
    # Extract the results for the current comparison
    res <- results(dds, contrast = c("condition", comparison[1], comparison[2]))
    
    # Order results by adjusted p-value
    res <- res[order(res$padj), ]
    
    # Convert DESeq2 results to a data frame
    res_df <- as.data.frame(res)
    
    # Add a new sheet with the comparison name
    sheet_name <- paste(comparison[1], "vs", comparison[2], sep = "_")
    addWorksheet(wb, sheet_name)
    
    # Write the data frame to the worksheet
    writeData(wb, sheet = sheet_name, x = res_df, rowNames = TRUE)
  }
  
  # Save the workbook to the file
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("Results have been successfully exported to", output_file)
} else {
  cat("The count data columns and metadata rows do not match.")
}

# Define the path to the Excel file
input_file <- "DESeq2_results_comparisons.xlsx"

# Read the results for the specified comparisons
results_WTfemale_vs_WTmale <- read_excel(input_file, sheet = "WTfemale_vs_WTmale")
results_hetmale_vs_WTmale <- read_excel(input_file, sheet = "hetmale_vs_WTmale")

# Ensure the rownames are the gene names
rownames(results_WTfemale_vs_WTmale) <- results_WTfemale_vs_WTmale$...1
rownames(results_hetmale_vs_WTmale) <- results_hetmale_vs_WTmale$...1

# Extract log2 fold changes
log2FC_WTfemale_vs_WTmale <- results_WTfemale_vs_WTmale$log2FoldChange
log2FC_hetmale_vs_WTmale <- results_hetmale_vs_WTmale$log2FoldChange

# Create a data frame with the log2 fold changes
fold_changes <- data.frame(
  gene = rownames(results_WTfemale_vs_WTmale),
  log2FC_WTfemale_vs_WTmale = log2FC_WTfemale_vs_WTmale,
  log2FC_hetmale_vs_WTmale = log2FC_hetmale_vs_WTmale
)

# List of specified genes
specified_genes_d2 <- c("Gsdme", "Gm42418", "Grin2b", "Malat1", "Gm26825", "Gm16105", "Rnf121", "Pcsk2", "Srcin1", 
                     "Ttc3", "Ube3a", "Atp6v0b", "Trim35", "Kcnq1ot1", "Sgcz", "Acp1", "Ypel3", "Sez6", "Ddx5", 
                     "Cbl", "Rpl21", "Tmem158", "Celf2", "Kalrn", "Ank3", "Trim2", "Opcml", "Penk", "Rbm39", 
                     "Tmem181a", "Unc79", "Il1rapl1", "Traip", "Gm26699", "Gm10600", "4833420G17Rik", "Chd3", 
                     "Itsn2", "Timp4", "Ryr3", "Snhg11", "Aste1", "Ftl1-ps1", "Ttc19", "Snhg12", "Dgkh", 
                     "2900055J20Rik", "Fam193b", "Strip2", "Golga4", "Adgrb1", "Prpf4b", "Ssbp4", "3100002H09Rik", 
                     "Rab11fip2", "Gm26724", "Unc13a", "Gm27032", "Gm16152", "Gm48715", "Epb41", "Baz2a", "Actn1", 
                     "Snx29", "Ddx50", "Ablim2", "Pde4a", "Rb1cc1", "Aldoa", "Gm5136", "Maz", "Hirip3", "Trank1", 
                     "Gabarapl2", "Pitpnm1", "Atxn2", "Herc1", "Camk2b", "Ppp3ca", "Negr1", "Ylpm1", "Tspoap1", 
                     "A230057D06Rik", "Crocc", "Adgrl1", "Elavl3", "Mprip", "Lrp1", "Eml5", "Adcy5", "Mtcl1", 
                     "Hnrnpu", "Celf3", "Pebp1", "Rasgef1a", "Eif1", "Ubb", "Ddx17", "Sez6l2", "Trrap", "Zfp445", 
                     "Fam228a", "Gm26917", "Gm15594", "Rasgrp2", "Usp7", "Slc16a6", "Zfp950", "Lage3", "Prrc2b", "Pde1b", "Kmt2e", "1110059E24Rik", "Clip4")

# Subset the fold changes data frame
subset_fold_changes <- fold_changes[fold_changes$gene %in% specified_genes_d2, ]

# Scatter plot
p2 <- ggplot(subset_fold_changes, aes(x = log2FC_WTfemale_vs_WTmale, y = log2FC_hetmale_vs_WTmale)) +
  geom_point() +
  geom_text_repel(aes(label = gene)) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  stat_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson", label.x = -4, label.y = 4) +
  labs(title = "D2-SPN DEGs between 16p males vs wt females",
       x = "Pseudobulk log2 fold change (WTfemale vs WTmale)",
       y = "Pseudobulklog2 fold change (hetmale vs WTmale)") +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Save the plot as a PDF
ggsave("D2-SPN DEGs between 16p males vs wt females.pdf", plot = p2, width = 8, height = 6)
# Save the subset_fold_changes dataframe to an Excel file
write_xlsx(subset_fold_changes, "D2-SPN DEGs between 16p males vs wt females subset_fold_changes.xlsx")
