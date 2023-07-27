# SC_RNA_Repo_Data_2

## About the Dataset: 
This is a new data set with 8 samples. 7 (SRR9897621, SRR14615558, SRR9897623, SRR9897624, SRR9897625, SRR12539462, SRR12539463) are BLCA samples, 1 () is a normal bladder mucosa sample. Among them, SRR14615558 is a paracancerous tissue sample from BC4, i.e, SRR9897624.   

## Reading the Files and Preparing Seurat Object: 
```R
#________________________Reading the files______________________#
# create list of samples
samples_data_2 <- list.files("/Users/andrew/Documents/Saindhabi/Data_2/")

# read files into Seurat objects
for (file in samples_data_2){
  print(paste0(file))
  seurat_data_2 <- Read10X(data.dir <- paste0("/Users/andrew/Documents/Saindhabi/Data_2/", file))
  seurat_obj_2 <- CreateSeuratObject(counts = seurat_data_2, 
                                     min.features = 100, 
                                     project = file)
  assign(file, seurat_obj_2)
}

# now merging all objects into one Seurat obj 
# There are 7 tumor samples and 1 normal (SRR14615558) sample in this cohort
merged_seurat_data_2 <- merge(x = SRR9897621,
                              y = c(SRR9897623,
                                    SRR9897624,
                                    SRR9897625,
                                    SRR12539462,
                                    SRR12539463,
                                    SRR14615558,
                                    SRR12603780),
                              add.cell.id = samples_data_2)
```
## Computing Percentage Mitochondrial Ratio: 
```
# Compute percent mito ratio
merged_seurat_data_2$mitoRatio <- PercentageFeatureSet(object = merged_seurat_data_2, pattern = "^MT-")
merged_seurat_data_2$mitoRatio <- merged_seurat_data_2@meta.data$mitoRatio / 100

# adding cell column
merged_seurat_data_2$cells <- rownames(merged_seurat_data_2@meta.data)
# re-setting the rownames
rownames(merged_seurat_data_2@meta.data) <- merged_seurat_data_2@meta.data$cells
```
## Filteration:
```
# Filteration
filtered_seurat_data_2 <- subset(merged_seurat_data_2, 
                          subset= nCount_RNA >= 1000 &
                                  nFeature_RNA <= 6000 & 
                                  mitoRatio < 0.10)
```
## Normalization and Integration: 
```
#________________________Integration using Harmony____________________________#
#integration using harmony need sevral steps to be undertaken:

# Perform log-normalization and feature selection, as well as SCT normalization on global object
merged_seurat_data_2 <- filtered_seurat_data_2 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("mitoRatio", "orig.ident"))

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merged_seurat_data_2 <- RunPCA(merged_seurat_data_2, assay = "SCT", npcs = 50)
```
```
#PC_ 1 
#Positive:  DCN, CFD, LUM, MT2A, MGP, GSN, COL1A2, FBLN1, CCL2, SFRP2 
#           IGFBP6, MFAP4, COL1A1, SOD3, COL6A2, COL3A1, CCN1, SERPINF1, C11orf96, IGFBP7 
#           VIM, CCDC80, PTGDS, CCN2, SPARC, PI16, RARRES2, MT1E, NNMT, MEG3 
#Negative:  SRGN, HLA-DRA, CXCR4, CD74, LYZ, SAMSN1, RGS1, TYROBP, CCL4, HLA-DPB1 
#           CD52, PTPRC, CD69, GPR183, CCL5, HLA-DRB1, HLA-DPA1, FCER1G, AREG, LAPTM5 
#           SIPA1L1, HLA-DQA1, CCL3, KLRB1, BCL2A1, STAT4, EREG, SYTL3, RHOH, IL1B 

#PC_ 2 
#Positive:  ADIRF, SPINK1, TAGLN, MYL9, RGS5, S100A6, IGFBP7, ACTA2, FABP4, CSTB 
#           PSCA, MT1M, KRT13, KRT19, CLDN4, MUSTN1, UCA1, DHRS2, S100P, MT1A 
#           KRT17, RPL41, CD24, SNCG, ID1, GDF15, NDUFA4L2, UPK2, FABP5, MYH11 
#Negative:  HLA-DRA, TYROBP, CD74, C1QB, C1QA, HLA-DPB1, HLA-DPA1, C1QC, CCL3, FCER1G 
#           HLA-DRB1, LYZ, APOE, AIF1, CCL4, CCL3L1, SRGN, CCL4L2, HLA-DQA1, RGS1 
#           HLA-DQB1, IL1B, MS4A6A, GPR183, IFI30, SPP1, CD14, DOCK4, ALOX5AP, LST1 

#PC_ 3 
#Positive:  ACKR1, ADAMTS9, CCL14, CLDN5, ZNF385D, MCTP1, TNFRSF6B, AQP1, IGFBP7, IL6 
#           SPARCL1, PLVAP, MT2A, VWF, SELE, CSF3, SLCO2A1, FLT1, TCF4, TLL1 
#           PECAM1, ECSCR, RAMP3, SPRY1, ADGRL4, CALCRL, ICAM1, SERPINA3, PCAT19, CNKSR3 
#Negative:  DCN, CFD, CXCR4, LUM, CCL5, CD52, COL1A2, CD69, CCL4, RGS1 
#           FBLN1, SFRP2, TRAC, CD3D, NKG7, KLRB1, PTPRC, FYN, GZMK, TRBC1 
#           CD2, GSN, GZMA, RHOH, COL1A1, IGFBP6, CD7, MGP, CD247, IL32 

#PC_ 4 
#Positive:  CFD, DCN, ACKR1, CCL14, GSN, CLDN5, LUM, ZNF385D, TNFRSF6B, MCTP1 
#           MGP, FBLN1, AQP1, CSF3, SFRP2, PLVAP, SERPINA3, SELE, VWF, ITM2A 
#           PI16, CCL5, CCN1, IGFBP6, TLL1, SLCO2A1, CXCR4, PECAM1, ADGRL4, ECSCR 
#Negative:  TAGLN, MYL9, RGS5, ACTA2, MT1M, MT1A, TPM2, MT2A, C11orf96, IGFBP7 
#           CRISPLD2, MUSTN1, MYH11, PPP1R14A, FRZB, ADAMTS4, MT1E, NDUFA4L2, MT1X, BGN 
#           GJA4, PTP4A3, CALD1, CACNA1C, ADRA2A, HEYL, WFDC1, PTPRG, MYLK, COX4I2 

#PC_ 5 
#Positive:  CCL5, KLRB1, TRAC, GZMA, CD3D, CCL4, IL32, CD52, NKG7, CD2 
#           LINC01871, IGFBP7, TRBC1, TAGLN, CD7, FYN, CXCR4, RGS1, PTPRC, MT2A 
#           MYL9, ADAMTS9, CD69, TRBC2, GZMB, RGS5, CD3E, GZMK, ACTA2, SYTL3 
#Negative:  LYZ, SPINK1, HLA-DRA, FTH1, PTGDS, FTL, TYROBP, CD74, CXCL1, EREG 
#           CXCL8, ADIRF, CSTB, PSCA, FCER1G, IL1B, IGHM, CD79A, CCL11, C1QA 
#           UPK2, S100A6, SPRR3, S100A8, APOE, S100A9, IGFBP5, IGKC, C1QB, SFRP2 
```
```
merged_seurat_data_2 <- RunTSNE(merged_seurat_data_2, assay = "SCT", npcs = 50)
```
```
# Integration
#install.packages("harmony")

library(harmony)

harmonized_seurat_data_2 <- RunHarmony(merged_seurat_data_2,
                                       group.by.vars = "orig.ident", 
                                       reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

harmonized_seurat_data_2 <- RunUMAP(harmonized_seurat_data_2, reduction = "harmony", assay = "SCT", dims = 1:40)
#harmonized_seurat_data_2 <- RunUTSNE(harmonized_seurat_data_2, reduction = "harmony", assay = "SCT", dims = 1:40)
```
## Clustering: 
```
#________________________Cluster identification and Inspect the effects of Harmony batch removal ____________#

# to set reduction to harmony and find the clusters
harmonized_seurat_data_2 <- FindNeighbors(object = harmonized_seurat_data_2, reduction = "harmony")
harmonized_seurat_data_2 <- FindClusters(harmonized_seurat_data_2, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8))

# visualization
Idents(harmonized_seurat_data_2) <- harmonized_seurat_data_2@meta.data$SCT_snn_res.0.1

# color cells based on the sample name
# Plot UMAP 
png(filename = "harmony_UMAP_y_sample_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(harmonized_seurat_data_2,
        group.by = "orig.ident",
        reduction = "umap")
dev.off()
```
![harmony_UMAP_y_sample_data_2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/d340fe8a-bc71-4cd6-8b1c-6b4c431d3704)

## Clusters with Labels: 
```
#________________________SuperCluster Identification____________#

png(filename = "harmony_umap_cluster_with_label_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(harmonized_seurat_data_2,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()
```
![harmony_umap_cluster_with_label_data_2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/4897f9fc-329e-4053-8d99-10220c2e41c2)

## Marker Identification: 
```
# let's visualize cells expressing supercluster markers:
# CD31: PECAM1
markers <- c("EPCAM", "PECAM1", "COL1A1", "PDGFRA", "RGS5", "CD79A", "LYZ", "CD3D", "TPSAB1")

png(filename = "umap_superCluster_cells_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = harmonized_seurat_data_2,
            features = markers,
            order = TRUE,
            min.cutoff = "q10",
            reduction = "umap",
            label = TRUE,
            repel = TRUE)

dev.off()
```
![umap_superCluster_cells_data_2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/961b7603-e86e-4883-bf6f-bcb5f12a831c)

### Visualising All Markers:
```
#______________________________ All markers________________________________
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers_data_2 <- FindAllMarkers(object = harmonized_seurat_data_2, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 

saveRDS(markers_data_2, "harmony_markers_data_2.RDS")

# mutate the markers dataframe
# Extract top 10 markers per cluster

top10_data_2 <- markers_data_2 %>%
  mutate(delta_pct = (pct.1 - pct.2)) %>%
  #filter(avg_log2FC > 1.5) %>%  # only keep rows where avg_log2FC > 1.5
  group_by(cluster) %>%
  top_n(n = 10, wt = delta_pct)

data.table::fwrite(top10_data_2, "harmony_blca_top10_all_markers_data_2.csv")

# Visualization of top markers in each cluster:
# Cluster markers:
cluster_markers_10_data_2 <- top10_data_2 %>% 
  group_by(cluster) %>% 
  summarize(genes = paste(gene, collapse = ","))
data.table::fwrite(cluster_markers_10_data_2, "cluster_markers_10_data_2.csv")

# feature plot for top markers
plotList = list()

for(cluster in 1:nrow(cluster_markers_10_data_2)){
  mkr = unlist(strsplit(cluster_markers_10_data_2$genes[cluster], ","))
  plotList[[cluster]] = FeaturePlot(object = harmonized_seurat_data_2,
                                    features = mkr,
                                    order = TRUE,
                                    min.cutoff = "q10",
                                    reduction = "umap",
                                    label = TRUE,
                                    repel = TRUE)
}
```
## Iterate over All Clusters:
```
# Cluster 0
png(filename = "harmony_blca_clsuter_markers_cluster0_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
plotList[[1]]
dev.off()
```
[Cluster0](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122897224)
[Cluster1](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122898703)
[Cluster2]
[Cluster3]
[Cluster4]
[Cluster5]
[Cluster6]
[Cluster7]
[Cluster8]
[Cluster9]
[Cluster10]
[Cluster11]
[Cluster12]

```
# LYZ cells
png(filename = "LYZ_harmony_blca_clsuter_marker_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = harmonized_seurat_data_2,
            features = c("LYZ"),
            order = TRUE,
            min.cutoff = "q10",
            reduction = "umap",
            label = TRUE,
            repel = TRUE)
dev.off()
```
![LYZ_harmony_blca_clsuter_marker_data_2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/3c177009-e5b0-4b93-bcce-4deb32e0385a)

## Renaming the Clusters:
```
# renaming clusters

# Rename all identities
harmonized_seurat_data_2 <- RenameIdents(object = harmonized_seurat_data_2, 
                                  "0" = "Luminal cells",
                                  "1" = "Basal cells", 
                                  "2" = "Luminal markers, hormonal metabolism related genes & cell adhesion",
                                  "3" = "Immune related cells",
                                  "4" = "Mesenchymal cells",
                                  "5" = "Squamous-like cells",
                                  "6" = "T-lymphocytes",
                                  "7" = "Epithelial cells - Luminal-like",
                                  "8" = "B-lymphocytes",
                                  "9" = "Heterogenous stromal cells",
                                  "10" = "Urothelial or transitional cell-like",
                                  "11" = "Heterogenous cancer cells",
                                  "12" = "Luminal-like cells")

# Plot the UMAP withy new labells
png(filename = "harmont_blca_umap_with_label_data_2.png", width = 16, height = 8.135, units = "in", res = 600)
DimPlot(object = harmonized_seurat_data_2, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
dev.off()
```
![harmont_blca_umap_with_label_data_2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/0fea68cb-a5a9-4f74-a93f-9339884b4261)
