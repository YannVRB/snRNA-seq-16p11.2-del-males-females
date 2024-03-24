library(Seurat)

WTmale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample795/outs/filtered_feature_bc_matrix")
WTmale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample809/outs/filtered_feature_bc_matrix")
hetmale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample796/outs/filtered_feature_bc_matrix")
hetmale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_males/Sample807/outs/filtered_feature_bc_matrix")
WTfemale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample798/outs/filtered_feature_bc_matrix")
WTfemale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample803/outs/filtered_feature_bc_matrix")
hetfemale1.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample806/outs/filtered_feature_bc_matrix")
hetfemale2.data <- Read10X(data.dir = "/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11.2_females/Sample811/outs/filtered_feature_bc_matrix")

WTmale1 <- CreateSeuratObject(counts = WTmale1.data, project = "WT_male", min.cells = 3, min.features = 200)
WTmale2 <- CreateSeuratObject(counts = WTmale2.data, project = "WT_male", min.cells = 3, min.features = 200)
hetmale1 <- CreateSeuratObject(counts = hetmale1.data, project = "WT_male", min.cells = 3, min.features = 200)
hetmale2 <- CreateSeuratObject(counts = hetmale2.data, project = "WT_male", min.cells = 3, min.features = 200)
WTfemale1 <- CreateSeuratObject(counts = WTfemale1.data, project = "WT_female", min.cells = 3, min.features = 200)
WTfemale2 <- CreateSeuratObject(counts = WTfemale2.data, project = "WT_female", min.cells = 3, min.features = 200)
hetfemale1 <- CreateSeuratObject(counts = hetfemale1.data, project = "het_female", min.cells = 3, min.features = 200)
hetfemale2 <- CreateSeuratObject(counts = hetfemale2.data, project = "het_female", min.cells = 3, min.features = 200)

WTmale1[["percent.mt"]] <- PercentageFeatureSet(WTmale1, pattern = "^mt-")
WTmale2[["percent.mt"]] <- PercentageFeatureSet(WTmale2, pattern = "^mt-")
hetmale1[["percent.mt"]] <- PercentageFeatureSet(hetmale1, pattern = "^mt-")
hetmale2[["percent.mt"]] <- PercentageFeatureSet(hetmale2, pattern = "^mt-")
WTfemale1[["percent.mt"]] <- PercentageFeatureSet(WTfemale1, pattern = "^mt-")
WTfemale2[["percent.mt"]] <- PercentageFeatureSet(WTfemale2, pattern = "^mt-")
hetfemale1[["percent.mt"]] <- PercentageFeatureSet(hetfemale1, pattern = "^mt-")
hetfemale2[["percent.mt"]] <- PercentageFeatureSet(hetfemale2, pattern = "^mt-")

WTmale1 <- subset(WTmale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTmale2 <- subset(WTmale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetmale1 <- subset(hetmale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetmale2 <- subset(hetmale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTfemale1 <- subset(WTfemale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WTfemale2 <- subset(WTfemale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetfemale1 <- subset(hetfemale1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hetfemale2 <- subset(hetfemale2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

WTmale1$dataset <- "WTmale"
WTmale2$dataset <- "WTmale"
WTfemale1$dataset <- "WTfemale"
WTfemale2$dataset <- "WTfemale"
hetmale1$dataset <- "hetmale"
hetmale2$dataset <- "hetmale"
hetfemale1$dataset <- "hetfemale"
hetfemale2$dataset <- "hetfemale"

WTmale1 <- NormalizeData(WTmale1, verbose = FALSE)
WTmale2 <- NormalizeData(WTmale2, verbose = FALSE)
hetmale1 <- NormalizeData(hetmale1, verbose = FALSE)
hetmale2 <- NormalizeData(hetmale2, verbose = FALSE)
WTfemale1 <- NormalizeData(WTfemale1, verbose = FALSE)
WTfemale2 <- NormalizeData(WTfemale2, verbose = FALSE)
hetfemale1 <- NormalizeData(hetfemale1, verbose = FALSE)
hetfemale2 <- NormalizeData(hetfemale2, verbose = FALSE)

WTmale1 <- FindVariableFeatures(WTmale1, selection.method = "vst", nfeatures = 2000)
WTmale2 <- FindVariableFeatures(WTmale2, selection.method = "vst", nfeatures = 2000)
hetmale1 <- FindVariableFeatures(hetmale1, selection.method = "vst", nfeatures = 2000)
hetmale2 <- FindVariableFeatures(hetmale2, selection.method = "vst", nfeatures = 2000)
WTfemale1 <- FindVariableFeatures(WTfemale1, selection.method = "vst", nfeatures = 2000)
WTfemale2 <- FindVariableFeatures(WTfemale2, selection.method = "vst", nfeatures = 2000)
hetfemale1 <- FindVariableFeatures(hetfemale1, selection.method = "vst", nfeatures = 2000)
hetfemale2 <- FindVariableFeatures(hetfemale2, selection.method = "vst", nfeatures = 2000)

anchors <- FindIntegrationAnchors(object.list = list(WTmale1, WTmale2, WTfemale1, WTfemale2, hetmale1, hetmale2, hetfemale1, hetfemale2), dims = 1:20)

combined <- IntegrateData(anchorset = anchors, dims = 1:20)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined, reduction = "umap", split.by = "dataset")

#Compute celltype-specific biomarkers in each individual cluster agaisnt the rest of the cells then compare top significant biomarkers against litterature
cluster1.markers <- FindMarkers(combined, ident.1 = "1")

#Dot plot (bubble plot) for celltype-specific biomarkers
Idents(combined) <- factor(Idents(combined), levels = c("iSPNs", "dSPNs", "eSPNs", "Pax6", "Lamp5", "Vip", "Sst_Chodl", "Sst_Pvalb", "L2/3 IT CTX", "L4 IT CTX", "L5 IT CTX", "L6 IT CTX", "L2/3 IT ENTI", "L5 PT CTX", "Car3", "L5/6 NP CTX", "L6 CT CTX", "L6b CT CTX", "Endothelial", "Oligodendrocytes", "Cholinergic interneurons", "Astrocytes"))
markers.to.plot <- c("Penk", "Tac1", "Casz1", "Gabrg1", "Kit", "Vip", "Nos1", "Pvalb", "Stard8", "Adam33", "Dkkl1", "C1ql3", "Ndst4", "Hcn1", "Car3", "Fgd6", "Col5a1", "Ctgf", "Golim4", "Ptgds", "Ecel1", "Slc1a3")
DotPlot(combined, features = markers.to.plot, cols = c("lightblue", "orange", "darkblue", "darkorange"), dot.scale = 8, split.by = "dataset") +  RotatedAxis()
#Setup Seurat object fo differential analyses
combined$celltype.dataset <- paste(Idents(combined), combined$dataset, sep = "_")
combined$celltype <- Idents(combined)
Idents(combined) <- "celltype.dataset"