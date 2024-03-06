# scRNA-seq data obtained from: Werba et al. Nat Comm. 2023 (https://www.nature.com/articles/s41467-023-36296-4)
# Werba et al. deposited data downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205013

########################################################################################################################
# PREPARE R ENVIRONMENT
########################################################################################################################
# Libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
set.seed(1)



########################################################################################################################
# LOAD SINGLE CELL DATA
########################################################################################################################

############################################################
# Load metadata and scRNA-seq data

# Annotation file for samples (primary or metastasis sample type)
annot_sample <- read.table("Metadata/sample_info.txt", sep = '\t', header = T, quote = "", comment.char="")

# Load in all samples
sample_list <- list()
for(sample_name in list.files("Data/")){
    sample_list[[sample_name]] <- Read10X(data.dir = paste("Data/",sample_name, sep = ""))
    rm(sample_name)
}

# Create Seurat objects
for(patient in names(sample_list)){
    sample_list[[patient]] <- CreateSeuratObject(counts = sample_list[[patient]], project = patient)
    rm(patient)
}






########################################################################################################################
# SEURAT ANALYSIS: FILTER ALL SAMPLES
########################################################################################################################

############################################################
# Set quality cutoffs

# Quality cutoffs
nReads = 1500     # same as Werba et al.
nGenes = 500      # same as Werba et al.
pct.mt = 15       # same as Werba et al.
pct.erythroid = 1 # same as Werba et al.



############################################################
# Add % mitochondrial and erythroid reads, then filter by quality

for(sample in names(sample_list)){
    sample_name = sample
    sample <- sample_list[[sample]]
    
    cells_before = length(Cells(sample))
    
    sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
    sample <- PercentageFeatureSet(sample, pattern = "^HBA1$|^HBA2$|^HBB$|^HBM$|^ALAS2$", assay = "RNA", col.name = "percent.erythroid")
    
    sample <- subset(x = sample, 
                     subset = 
                         nCount_RNA > nReads & 
                         nFeature_RNA > nGenes &
                         percent.mt < pct.mt &
                         percent.erythroid <= pct.erythroid
    )
    
    cells_after = length(Cells(sample))
    
    print(paste(sample_name,": ", cells_before, " to ", cells_after, sep = ""))
    
    sample_list[[sample_name]] <- sample
    rm(sample,sample_name, cells_before, cells_after)
}



############################################################
# Split liver mets and primary tumors into separate lists

sample_list_all <- sample_list
sample_list_primary <- sample_list_all[annot_sample$Patient[annot_sample$Source=="Pancreas"]]
sample_list_liver <- sample_list_all[annot_sample$Patient[annot_sample$Source=="Liver"]]






########################################################################################################################
# SEURAT ANALYSIS: PRIMARY TUMORS
########################################################################################################################

# Set variables to primary tumors
PlotPrefix <- "Werba2023_PrimaryTumors"  # Plots will be printed with this prefix
sample_list <- sample_list_primary



############################################################
# SCTransform, Integrate, UMAP, and Clusters

# Run SCTransform on each sample
sample_list <- lapply(X = sample_list, FUN = function(x){
    x <- SCTransform(x, vst.flavor = "v2", verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)

# Run PCA (for RPCA integration) on each sample
sample_list <- lapply(X = sample_list, FUN = function(x){
    x <- RunPCA(x, verbose = TRUE)
})

# RPCA integration (fast)
sample <- FindIntegrationAnchors(object.list = sample_list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
sample <- IntegrateData(anchorset = sample, normalization.method = "SCT")

# PCA, UMAP, clustering
sample <- RunPCA(sample, verbose = FALSE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30, verbose = FALSE)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
sample <- FindClusters(sample, resolution = 0.7)

rm(sample_list,features)



############################################################
# Set assay to RNA and normalize

DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample)
sample <- ScaleData(sample, features = rownames(sample))



############################################################
# Find cluster markers

sample.allMarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



############################################################
# Collect cell type info based on significant marker genes

# Cell type markers from Fig 1C in Werba et al. (2023) (see Werba et al. methods)
celltype_list <- list(
    "T/NK" = c("CD3E"),
    "Epithelial" = c("KRT19"),
    "Endothelial" = c("VWF","PECAM1"),
    "Myeloid" = c("CD68"),
    "B/Plasma" = c("CD79A"),
    "Mesenchyme" = c("DCN"),
    "Mast" = c("KIT"),
    "Prolif. Epithelial" = c("KRT19","MKI67"),
    "Prolif. Lymphoid" = c("CD3E","MKI67"),
    "Prolif. Myeloid" = c("CD68","MKI67")
)

# Make list of Seurat clusters that will be annotated with cell type
celltype_tally <- list()
for(clust in unique(as.character(sample$seurat_clusters))){
    celltype_tally[[clust]] = character()
    rm(clust)
}

# Parse marker gene results (p-adj <= 0.05) for cell type markers in 'celltype_list'.
# Annotate clusters in 'celltype_tally' with cell types based on matching marker gene enrichment
for(cell_type in names(celltype_list)){
    clusters = c()
    for(gene in celltype_list[[cell_type]]){
        clusters = c(clusters, as.character(sample.allMarkers$cluster[sample.allMarkers$gene==gene & sample.allMarkers$p_val_adj<=0.05]))
    }
    
    if(length(celltype_list[[cell_type]])>1){
        clusters = as.data.frame(table(clusters))
        clusters <- clusters$clusters[clusters$Freq==length(celltype_list[[cell_type]])]
    }
    
    clusters <- as.character(clusters)
    
    if(!(is.null(clusters))){
        for(clust in clusters){
            celltype_tally[[clust]] = c(celltype_tally[[clust]],cell_type)
        }
    }
    rm(cell_type,clusters,gene,clust)
}

# Show results of celltype tally
print(celltype_tally)



############################################################
# Refine cell type identification manually based on marker gene results

# Results of 'celltype_tally' list:
#   - Clusters 16,18,19,29,30 had no 'celltype_list' marker gene match
#   - Cluster 33 are doublets (express both epithelial and endothelial markers)

# Manually annotate clusters 16,18,19,29,30, and 33 (markers obtained from Werba et al. (2023) Methods section)
# These markers describe granular cell types in Werba et al. (2023), but here we label cells as broad cell types
celltype_tally[["16"]] = "T/NK"       # cluster 16 has NK-cell markers GNLY and NCAM1
celltype_tally[["18"]] = "Myeloid"    # cluster 18 has cDC2 markers CD1C and FCER1A
celltype_tally[["19"]] = "Mesenchyme" # cluster 19 has Pericyte markers ACTA2, RGS5, CSPG4, NOTCH3
celltype_tally[["29"]] = "Mesenchyme" # cluster 29 has Schwann markers SOX10, S100B, PLP1
celltype_tally[["30"]] = "Mesenchyme" # cluster 30 has pDC markers LILRA4, PLD4
celltype_tally[["33"]] = "Doublets"   # cluster 33 has Epithelial-Endothelial dual expression



############################################################
# Annotate cell types and remove doublets

# Prepare metadata data frame to add to Seurat object
celltype <- data.frame(cluster = sample$seurat_clusters)
celltype$cluster = as.character(celltype$cluster)

# Use celltype_tally info to add cell type results to data frame
for(clust in names(celltype_tally)){
    if(length(celltype_tally[[clust]])>0){
        celltype$cluster[celltype$cluster==clust] = celltype_tally[[clust]][length(celltype_tally[[clust]])]
    }
    
    rm(clust)
}


# Add cell type annotation to Seurat object
sample <- AddMetaData(sample, celltype, "Celltype")

# Remove doublets
sample <- sample[,sample$Celltype!="Doublets"]


# Plot UMAP annotated by cell type
# Plot with various annotation settings
PlotSuffix = "UMAP_CellType"
pdf(width = 6, height = 6, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
DimPlot(sample, label = F, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + NoLegend() + xlab("") + ylab("")
DimPlot(sample, label = F, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + xlab("") + ylab("")
DimPlot(sample, label = T, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + NoLegend() + xlab("") + ylab("")
dev.off()



############################################################
# Dot plots of markers by cell type

Idents(sample) <- "Celltype"
levels(sample) <- rev(c("Epithelial","Prolif. Epithelial","Mesenchyme","Endothelial","Myeloid","Prolif. Myeloid","T/NK","Prolif. Lymphoid","B/Plasma","Mast"))

PlotSuffix = "DotPlot_CellType"
pdf(width = 6.5, height = 5, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
features = c("CD3E","KRT19","VWF","CD68","MKI67","CD79A","DCN","KIT")
DotPlot(sample, features = features, scale = T) + RotatedAxis()
dev.off()



############################################################
# Module Scores: pORG and pSUB gene lists

features <- list(
    pORG = c("TNFSF18","CST2","RPRD1B","BARD1","SORD2P","ZNF217","CHMP4B","SLFN13","MMEL1","SPATA24","OSER1","PHACTR3","PKIB","ARHGAP12","RBM41","MPZL3","GDE1","FGD6","PLAAT2","STAU1","ZG16B","RIOK3","ALOX5","CDS1","PIK3CB","HPGD","HIST1H1T","GID8","NUDT21","PAK1","CYP4F3","ARL17A","PLAAT3","SRP9","ANKRD18B","CCL20","CFB","SPOPL","CTNNBL1","NFE2L3","ACSL5","CIAO2A","SCOC","SAA4","PPP1R3C","WFDC13","FUNDC1","ADNP","DLEU7","NME7","EFNB2","LRRIQ1","MAP3K2","LNPK","ZDHHC13"),
    pSUB = c("RHCG","KRT17","GJB6","P2RY2","S100A3","C16orf74","NXPH4","HES2","PPP1R14B","DNER","FSCN1","TMPRSS11E","KRT6A","KRT7","GJB2","GSDMC","SEMA3F","KRT6C","BANF1","RHOV","PC","TNFRSF6B","PLXNA1","DENND2C","FHOD3","HAS3","PLAU","NRP2","KYNU","PORCN","SYT12","LVRN","TBC1D2","LPAR3","RAC3","SLC16A3","KYAT1","GPR153","SLCO1C1","PHKA1","PLEKHG5","NIBAN2","CNTNAP3","STRIP2","CALB2","USP31","EGFR","LRRC8A","FGFBP1","CSNK1E","COL7A1"))

sample <- AddModuleScore(
    object = sample,
    features = features,
    ctrl = 100,
    name = names(features),
    assay = "RNA"
)


# Plot module scores on UMAP
# Plot with various annotation settings
PlotSuffix = "ModuleScore_UMAP_pORG_pSUB"
pdf(width = 7, height = 7, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
FeaturePlot(sample, "pORG1", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("") + ggtitle("Module Score: pORG")
FeaturePlot(sample, "pORG1", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend() + xlab("") + ylab("") + ggtitle("Module Score: pORG")
FeaturePlot(sample, "pSUB2", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("") + ggtitle("Module Score: pSUB")
FeaturePlot(sample, "pSUB2", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend() + xlab("") + ylab("") + ggtitle("Module Score: pSUB")
dev.off()



############################################################
# Save Seurat object

saveRDS(sample, file = paste(getwd(), "/", PlotPrefix, "_SeuratObject.rds", sep = ""))






########################################################################################################################
# SEURAT ANALYSIS: LIVER METASTASES
########################################################################################################################

# Set variables to liver mets
PlotPrefix <- "Werba2023_LiverMets"
sample_list <- sample_list_liver



############################################################
# SCTransform, Integrate, UMAP, and Clusters

# Run SCTransform on each sample
sample_list <- lapply(X = sample_list, FUN = function(x){
    x <- SCTransform(x, vst.flavor = "v2", verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)

# Run PCA (for RPCA integration) on each sample
sample_list <- lapply(X = sample_list, FUN = function(x){
    x <- RunPCA(x, verbose = TRUE)
})

# RPCA integration (fast)
sample <- FindIntegrationAnchors(object.list = sample_list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
sample <- IntegrateData(anchorset = sample, normalization.method = "SCT")

# PCA, UMAP, clustering
sample <- RunPCA(sample, verbose = FALSE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30, verbose = FALSE)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
sample <- FindClusters(sample, resolution = 0.7)

rm(sample_list,features)



############################################################
# Set assay to RNA and normalize

DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample)
sample <- ScaleData(sample, features = rownames(sample))



############################################################
# Find cluster markers

sample.allMarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



############################################################
# Collect cell type info based on significant marker genes

# Cell type markers from Fig 1C in Werba et al. (2023) (see Werba et al. methods)
celltype_list <- list(
    "T/NK" = c("CD3E"),
    "Epithelial" = c("KRT19"),
    "Endothelial" = c("VWF","PECAM1"),
    "Myeloid" = c("CD68"),
    "B/Plasma" = c("CD79A"),
    "Mesenchyme" = c("DCN"),
    "Mast" = c("KIT"),
    "Prolif. Epithelial" = c("KRT19","MKI67"),
    "Prolif. Lymphoid" = c("CD3E","MKI67"),
    "Prolif. Myeloid" = c("CD68","MKI67")
)

# Make list of Seurat clusters that will be annotated with cell type
celltype_tally <- list()
for(clust in unique(as.character(sample$seurat_clusters))){
    celltype_tally[[clust]] = character()
    rm(clust)
}

# Parse marker gene results (p-adj <= 0.05) for cell type markers in 'celltype_list'.
# Annotate clusters in 'celltype_tally' with cell types based on matching marker gene enrichment
for(cell_type in names(celltype_list)){
    clusters = c()
    for(gene in celltype_list[[cell_type]]){
        clusters = c(clusters, as.character(sample.allMarkers$cluster[sample.allMarkers$gene==gene & sample.allMarkers$p_val_adj<=0.05]))
    }
    
    if(length(celltype_list[[cell_type]])>1){
        clusters = as.data.frame(table(clusters))
        clusters <- clusters$clusters[clusters$Freq==length(celltype_list[[cell_type]])]
    }
    
    clusters <- as.character(clusters)
    
    if(!(is.null(clusters))){
        for(clust in clusters){
            celltype_tally[[clust]] = c(celltype_tally[[clust]],cell_type)
        }
    }
    rm(cell_type,clusters,gene,clust)
}

# Show results of celltype tally
print(celltype_tally)



############################################################
# Refine cell type identification manually based on marker gene results

# Results of 'celltype_tally' list:
#   - Clusters 11,13,15,16,20,23,29 had no 'celltype_list' marker gene match
#   - Clusters 22 and 32 are doublets (express both epithelial and myeloid markers)

# Manually annotate clusters 16,18,19,29,30, and 33 (markers obtained from Werba et al. (2023) Methods section)
# These markers describe granular cell types in Werba et al. (2023), but here we label cells as broad cell types
celltype_tally[["11"]] = "Epithelial" # cluster 11 has basal-like marker KRT17 (Werba et al. 2023 figure 2E)
celltype_tally[["13"]] = "T/NK"       # cluster 13 has NK-cell markers GNLY and NCAM1
celltype_tally[["15"]] = "Epithelial" # cluster 15 has basal-like marker KRT17 (Werba et al. 2023 figure 2E)
celltype_tally[["16"]] = "T/NK"       # cluster 16 has NK cell markers GZMH and GNLY
celltype_tally[["20"]] = "Hepatocyte" # cluster 20 has hepatocyte markers ALB and HP (Aizarani N et al. 2019 Nature)
celltype_tally[["23"]] = "Myeloid"    # cluster 23 has cDC2 markers CD1C and FCER1A
celltype_tally[["29"]] = "T/NK"       # cluster 29 has T cell marker IL7R

celltype_tally[["22"]] = "Doublets"   # cluster 33 has Epithelial-Myeloid dual expression
celltype_tally[["32"]] = "Doublets"   # cluster 33 has Epithelial-Myeloid dual expression



############################################################
# Annotate cell types and remove doublets

# Prepare metadata data frame to add to Seurat object
celltype <- data.frame(cluster = sample$seurat_clusters)
celltype$cluster = as.character(celltype$cluster)

# Use celltype_tally info to add cell type results to data frame
for(clust in names(celltype_tally)){
    if(length(celltype_tally[[clust]])>0){
        celltype$cluster[celltype$cluster==clust] = celltype_tally[[clust]][length(celltype_tally[[clust]])]
    }
    
    rm(clust)
}


# Add cell type annotation to Seurat object
sample <- AddMetaData(sample, celltype, "Celltype")

# Remove doublets
sample <- sample[,sample$Celltype!="Doublets"]


# Plot UMAP annotated by cell type
# Plot with various annotation settings
PlotSuffix = "UMAP_CellType"
pdf(width = 6, height = 6, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
DimPlot(sample, label = F, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + NoLegend() + xlab("") + ylab("")
DimPlot(sample, label = F, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + xlab("") + ylab("")
DimPlot(sample, label = T, group.by = "Celltype", raster = T, raster.dpi = c(2048,2048), pt.size = 2.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Cell Type") + NoLegend() + xlab("") + ylab("")
dev.off()



############################################################
# Dot plots of markers by cell type

Idents(sample) <- "Celltype"
levels(sample) <- rev(c("Epithelial","Prolif. Epithelial","Hepatocyte","Mesenchyme","Endothelial","Myeloid","Prolif. Myeloid","T/NK","Prolif. Lymphoid","B/Plasma"))

PlotSuffix = "DotPlot_CellType"
pdf(width = 6.5, height = 5, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
features = c("CD3E","KRT19","ALB","VWF","CD68","MKI67","CD79A","DCN")
DotPlot(sample, features = features, scale = T) + RotatedAxis()
dev.off()



############################################################
# Module Scores: pORG and pSUB gene lists

features <- list(
    pORG = c("TNFSF18","CST2","RPRD1B","BARD1","SORD2P","ZNF217","CHMP4B","SLFN13","MMEL1","SPATA24","OSER1","PHACTR3","PKIB","ARHGAP12","RBM41","MPZL3","GDE1","FGD6","PLAAT2","STAU1","ZG16B","RIOK3","ALOX5","CDS1","PIK3CB","HPGD","HIST1H1T","GID8","NUDT21","PAK1","CYP4F3","ARL17A","PLAAT3","SRP9","ANKRD18B","CCL20","CFB","SPOPL","CTNNBL1","NFE2L3","ACSL5","CIAO2A","SCOC","SAA4","PPP1R3C","WFDC13","FUNDC1","ADNP","DLEU7","NME7","EFNB2","LRRIQ1","MAP3K2","LNPK","ZDHHC13"),
    pSUB = c("RHCG","KRT17","GJB6","P2RY2","S100A3","C16orf74","NXPH4","HES2","PPP1R14B","DNER","FSCN1","TMPRSS11E","KRT6A","KRT7","GJB2","GSDMC","SEMA3F","KRT6C","BANF1","RHOV","PC","TNFRSF6B","PLXNA1","DENND2C","FHOD3","HAS3","PLAU","NRP2","KYNU","PORCN","SYT12","LVRN","TBC1D2","LPAR3","RAC3","SLC16A3","KYAT1","GPR153","SLCO1C1","PHKA1","PLEKHG5","NIBAN2","CNTNAP3","STRIP2","CALB2","USP31","EGFR","LRRC8A","FGFBP1","CSNK1E","COL7A1"))

sample <- AddModuleScore(
    object = sample,
    features = features,
    ctrl = 100,
    name = names(features),
    assay = "RNA"
)


# Plot module scores on UMAP
# Plot with various annotation settings
PlotSuffix = "ModuleScore_UMAP_pORG_pSUB"
pdf(width = 7, height = 7, file = paste(PlotPrefix, "_", PlotSuffix, ".pdf", sep = ""))
FeaturePlot(sample, "pORG1", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("") + ggtitle("Module Score: pORG")
FeaturePlot(sample, "pORG1", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend() + xlab("") + ylab("") + ggtitle("Module Score: pORG")
FeaturePlot(sample, "pSUB2", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("") + ggtitle("Module Score: pSUB")
FeaturePlot(sample, "pSUB2", order = T, pt.size = 2.5, raster = T, raster.dpi = c(2048,2048)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend() + xlab("") + ylab("") + ggtitle("Module Score: pSUB")
dev.off()



############################################################
# Save Seurat object

saveRDS(sample, file = paste(getwd(), "/", PlotPrefix, "_SeuratObject.rds", sep = ""))


