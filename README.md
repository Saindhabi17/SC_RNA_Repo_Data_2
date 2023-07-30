# SC_RNA_Repo_Data_2

## About the Dataset: 
This is a new data set with 8 samples. 7 (SRR9897621, SRR14615558, SRR9897623, SRR9897624, SRR9897625, SRR12539462, SRR12539463) are BLCA samples, 1 (SRR12603780) is a normal bladder mucosa sample. Among them, SRR14615558 is a paracancerous tissue sample from BC4, i.e, SRR9897624.   

## Analysing only the tumours:

### Reading Files & Preparing Seurat Object: 
Only the tumour samples are taken in account here. There are 6 samples now. 

```R
#________________________Reading the files______________________#
# create list of samples
samples_data_2_n <- list.files("/Users/andrew/Documents/Saindhabi/Data_2/")

# read files into Seurat objects
for (file in samples_data_2_n){
  print(paste0(file))
  seurat_data_2_n <- Read10X(data.dir <- paste0("/Users/andrew/Documents/Saindhabi/Data_2/", file))
  seurat_obj_2_n <- CreateSeuratObject(counts = seurat_data_2_n, 
                                     min.features = 100, 
                                     project = file)
  assign(file, seurat_obj_2_n)
}

# updated sample name 
samples_2_n <- samples_data_2_n[-c(3,4)]

# now merging all objects into one Seurat obj 
# There are 7 tumor samples and 1 normal (SRR14615558) sample in this cohort
# Here we are working with 6 tumor samples for the time being.
merged_seurat_data_2_n <- merge(x = SRR9897621,
                              y = c(SRR9897623,
                                    SRR9897624,
                                    SRR9897625,
                                    SRR12539462,
                                    SRR12539463),
                              add.cell.id = samples_2_n)
```
## Quality Control
First, I have explored the metadata of the merged Seurat file. 
```R
# Explore merged metadata
View(merged_seurat_data_2_n@meta.data)
```
There are 3 columns in the merged meta data. They are-

1. orig.ident: The first column contains the sample identity as known. By default it shows the value provided for the project argument when loading in the data.
2. nCount_RNA: This column represents the number of UMIs per cell. UMI (unique molecular identifiers) is used to determine whether a read is a biological or technical duplicate (PCR duplicate). There can be 2 types of duplicates - Biological Duplicates - Reads with different UMIs mapping to the same transcript derived from different molecules, and Technical Duplicates - Reads with the same UMI originated from the same molecule. For the biological duplicates each read should be counted where for the technical duplicates reads should be counted as a single one.
3. nFeature_RNA: This column represents the number of genes detected per cell.

### Recommended Features to Add to the Metadata
1. Novelty Score: It is the number of genes detected per UMI. More genes detected per UMI, more complex the data will be.
2. Mitochondrial Ratio: This metric will give us a percentage of cell reads originating from the mitochondrial genes (coming from dying cells).
```R
#Add number of genes per UMI for each cell to metadata
merged_seurat_data_2_n$log10GenesPerUMI <- log10(merged_seurat_data_2_n$nFeature_RNA) / log10(merged_seurat_data_2_n$nCount_RNA)

# Compute percent mito ratio
merged_seurat_data_2_n$mitoRatio <- PercentageFeatureSet(object = merged_seurat_data_2_n, pattern = "^MT-")
merged_seurat_data_2_n$mitoRatio <- merged_seurat_data_2_n@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata_2_n <- merged_seurat_data_2_n@meta.data
# Add cell IDs to metadata
metadata_2_n$cells <- rownames(metadata_2_n)

# adding sample type to metadata. The orginal file could be download from SRA explorer
SampleType_n <- c("BC1", "BC3", "BC4", "BC5", "BC6", "BC7")
names(SampleType_n) <- c("SRR9897621", "SRR9897623", "SRR9897624", "SRR9897625", "SRR12539462", "SRR12539463")

metadata_2_n$sampleType <- stringr::str_replace_all(metadata_2_n$orig.ident, SampleType_n)

# Rename columns
metadata_2_n <- metadata_2_n %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA,
                sample = sampleType)

# Add metadata back to Seurat object
merged_seurat_data_2_n@meta.data <- metadata_2_n

# Create .RData object to load at any time
save(merged_seurat_data_2_n, file="merged_filtered_seurat_data_2_n.RData")
```

## Visualization

### Cell counts per sample 
```R
# Visualize the number of cell counts per sample
png(filename = "cell_counts_before_QC.png", width = 16, height = 8.135, units = "in", res = 300)
bqcc <- metadata_2_n %>% 
  ggplot(aes(x=seq_folder, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells before QC")
dev.off()
```
### UMI per cell
Typically, we expect the UMI counts per cell to be higher than 500, which is the lower limit of the expected range. If the UMI counts range between 500-1000, the data is still usable, but deeper sequencing may have been beneficial for these cells.

```R
# Visualize the number UMIs/transcripts per cell
png(filename = "UMI_per_Transcript.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_2_n %>% 
  ggplot(aes(x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 1000) +
  labs(fill = "Sample")
dev.off()
```
![UMI_per_Transcript](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/054f68fa-d16d-4253-9194-6580cc236b89)

From the plots, it is clear that the cells have way more than 1000 UMI.

### Genes Detected Per Cell
In scRNA-seq, the number of genes detected per cell is a crucial quality metric that we expect to be similar to the UMI detection, albeit slightly lower.

For high-quality data, the proportional histogram of genes detected per cell should show a single peak that represents encapsulated cells. However, if there is a small shoulder or a bimodal distribution to the left of the main peak, this could indicate a few things. It could be due to some failed cells or biologically different cell types, such as quiescent cell populations or less complex cells of interest. For instance, larger cells or different cell types may have higher gene counts.
```R
# Visualize the distribution of genes detected per cell via histogram
png(filename = "Genes_detected_per_cell.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_2_n %>% 
  ggplot(aes(x=nGene, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 500) +
  labs(fill = "Sample")
dev.off()
```
![Genes_detected_per_cell](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/5c4230bc-ccbb-4ba3-8e66-9d87aa715325)

### Novelty Score
The novelty score, computed as the ratio of nGenes over nUMI, measures the complexity of RNA species in each cell. A low number of genes detected in a cell with many captured transcripts (high nUMI) indicates low complexity or novelty. This could be due to an artifact, contamination, or represent a specific cell type (e.g. red blood cells). A good quality cell typically has a novelty score above 0.80.
```R
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
png(filename = "Novelty_score.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_2_n %>%
  ggplot(aes(x=log10GenesPerUMI, fill=sample)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~seq_folder) +
  xlab("Novelty Score") +
  geom_vline(xintercept = 0.8)
dev.off()
```
![Novelty_score](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/69083263-55c0-4d58-8f29-26586ae894ba)

### Mitochondrial Gene Expression Detected per Cell
High level of expression from mitochondria indicate dying or dead cells. Basically poor quality samples are those that exceed 0.2 mitochondria ratio mark.
```R
# Visualize the distribution of mitochondrial gene expression detected per cell
png(filename = "Mito_expression.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_2_n %>%
  ggplot(aes(x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) + 
  theme_classic() +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 0.2)
dev.off()
```
![Mito_expression](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/44017c29-62d1-4854-b473-b658b9a870eb)

```R
# Joint Filtering- nUMI, nGene and mitoRatio 
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(filename = "Joint_filtering.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_2_n %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~seq_folder)
dev.off()
```
![Joint_filtering](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/00703859-e9cf-48d2-98b0-20ec869f32c5)

There are samples that shows high-quality cells ; high nUMI, high nGene, low number of cells with high mitoRatio and also there are some samples that would clearely benfit from filtering, as they have low quality cells. We expect to see that dying cells to show high level of mitoRatio and low nUMI and nGene.

Basically, it is not uncommon to observe cells with high numbers of UMIs and nGene with, but also high mitoRatio. These cells may be stressed or damaged, but they could also represent a heterogeneous population of cells with distinct metabolic states.

To investigate the potential cause of high mitochondrial expression ratios, it is important to examine the expression of specific mitochondrial genes and compare them to other genes in the cell. If the expression of mitochondrial genes is elevated relative to other genes, this could suggest mitochondrial dysfunction. Additionally, examining the expression of other stress or damage markers, such as heat shock proteins or cell cycle genes, can also provide insight into the health and state of the cell.
# Filtering

## The Cell-level Filtering
-nUMI > 500, -nGene > 250 & < 6000, -log10GenesPerUMI or Novelty Score > 0.8, -mitoRatio < 0.10
```R
# Cell level filtering 
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat_data_2_n <- subset(merged_seurat_data_2_n, 
                          subset= nUMI >= 500 &
                            nGene >= 250 &
                            nGene <= 6000 & 
                            log10GenesPerUMI > 0.80 & 
                            mitoRatio < 0.10)
```
## The Gene-level Filtering
Keeping only genes which are expressed in 100 or more cells (usually this is 10)
```R
# Gene level filtering:
# Extract counts
counts_2_n <- GetAssayData(object = filtered_seurat_data_2_n, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero_2_n <- counts_2_n > 0
# Sums all TRUE values and returns TRUE if more than 100 TRUE values per gene
keep_genes_2_n <- Matrix::rowSums(nonzero_2_n) >= 100

# Only keeping those genes expressed in more than 100 cells
filtered_counts_2_n <- counts_2_n[keep_genes_2_n, ]
# Reassign to filtered Seurat object
filtered_seurat_data_2_n <- CreateSeuratObject(filtered_counts_2_n, meta.data = filtered_seurat_data_2_n@meta.data)

# Create .RData object to load at any time
save(filtered_seurat_data_2_n, file="seurat_filtered_data_2_n.RData")
```

## Re assess QC Metric:
```R
# Save filtered subset to new metadata
metadata_clean_2_n <- filtered_seurat_data_2_n@meta.data

# to see drop in filtering cells:

met_before_2_n <- data.frame(unclass(table(metadata_2_n$seq_folder)))
met_before_2_n$QCgroup <- "before"
met_before_2_n$cell<- rownames(met_before_2_n)
names(met_before_2_n)[1] <- 'count'

met_after_2_n <- data.frame(unclass(table(metadata_clean_2_n$seq_folder)))
met_after_2_n$QCgroup <- "after"
met_after_2_n$cell<- rownames(met_after_2_n)
names(met_after_2_n)[1] <- 'count'
# count
cell_count_2_n <- data.frame(rbind(met_before_2_n, met_after_2_n))
```
## Visualization of Re-assessment:
```R
# visualization :
png(filename = "nCells_before_after.png", width = 16, height = 8.135, units = "in", res = 300)
cell_count_2_n %>% ggplot(aes(x=cell, y=count, fill=QCgroup)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_manual(values = c("#808080", "#FFBF00")) +
  xlab("samples") +
  ggtitle("nCells count before and after QC")
dev.off()
```
![nCells_before_after](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/f53a9aaa-9f0e-485f-a2f6-e77d5d730026)

```R
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(filename = "correlation.png", width = 16, height = 8.135, units = "in", res = 300)
metadata_clean_2_n %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~seq_folder)
dev.off()
```
![correlation](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/9a62138d-6d73-4fda-88f6-a01841f5b329)

# Normalization and Regressing Out Unwanted Variation
The ultimate goal is to define clusters of cells and identify cell types in the samples. To achieve this, there are several steps. The 1st of them is normalization and regressing out the unwanted variation.

The first step is identifying unwanted variability by exploring data and covariates such as cell cycle and mitochondrial gene expression. Both biological source of variation (e.g. effect of cell cycle on transcriptome) and technical source should be explored and account for.

Next we have to normalize and remove unwanted variability using Seurat's SCTransform function. The normalization step is necessary to make expression counts comparable across genes and/or samples. The counts of mapped reads for each gene is proportional to the expression of RNA (“interesting”) in addition to many other factors (“uninteresting” such as sequencing depth and gene length).

Normalization is the process of adjusting raw count values to account for the “uninteresting” factors. For simplicity , normalization is assumed as two step process: scaling and transforming. In scaling the goal is to multiply each UMI count by a cell specific factor to get all cells to have the same UMI counts.For transformation simple approaches like log-transformation showed to be not that useful, especially in the case of genes with high expression but showing decent performance for low/intreemediate expressed genes. So we cannot treat all genes the same. The proposed solution for data transformation is Pearson residuals (inmplemented in Seurat's SCTransform function), which applies a gene-specific weight to each measurement based on the evidence of non-uniform expression across cells. This weight is higher for genes expressed in a smaller fraction of cells, making it useful for detecting rare cell populations. The weight takes into account not just the expression level but also the distribution of expression.

## Exploring sources of unwanted variation
Here, the goal is to evaluate the effects of cell cycle and mitochondrial expression -

## For Cell Cycle: 
1. First scoring the cells for cell cycle genes
2. Then determining whether the cell cycle is a major source of variation in our data set using PCA.

```R
# Sources of Unwanted Variation:

# Normalizing the counts
# This normalization method is solely for the purpose of exploring the sources of variation in our data.
seurat_phase_data_2_n <- NormalizeData(filtered_seurat_data_2_n, normalization.method = "LogNormalize", scale.factor = 10000)

# Loading cell cycle markers
load("/Users/andrew/Downloads/cycle.rda")

# Scoring cells for cell cycle
seurat_phase_data_2_n <- CellCycleScoring(seurat_phase_data_2_n, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# Viewing cell cycle scores and phases assigned to cells                                 
View(seurat_phase_data_2_n@meta.data) 
table(seurat_phase_data_2_n$Phase)
```
## Cells in different cell cycle phases -
1. G1 : 52468
2. G2M : 15399
3. S : 31209
From this, it is clear most of the cells are in G1 and S, which makes sense.
```R
# Identifying the most variable genes and scaling them
seurat_phase_data_2_n <- FindVariableFeatures(seurat_phase_data_2_n, 
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = TRUE)

# Identifying the 10 most highly variable genes
top10_data_2_n <- head(VariableFeatures(seurat_phase_data_2_n), 10)

# plotting variable features with and without labels

plot1_2_n <- VariableFeaturePlot(seurat_phase_data_2_n)
plot2_2_n <- LabelPoints(plot = plot1_2_n, points = top10_data_2_n, repel = TRUE)
dev.new()
library(ggplot2)
ggsave('plot1_2_n.png', plot1_2_n)
ggsave('plot2_2_n.png', plot2_2_n)
print(plot2_2_n)
print(plot1_2_n)
```
```R
# Checking quartile values for mitoRatio, we will use this variable later to mitigate unwanted source of variation in dataset
summary(seurat_phase_data_2_n@meta.data$mitoRatio)

# Turning mitoRatio into categorical factor vector based on quartile values
seurat_phase_data_2_n@meta.data$mitoFr <- cut(seurat_phase_data_2_n@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.015, 0.025, 0.045, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

# Scaling the counts
# This step is essential for PCA , clustering and heatmap generation
seurat_phase_data_2_n <- ScaleData(seurat_phase_data_2_n)
saveRDS(seurat_phase_data_2_n, "seurat_phase_data_2_n.rds")

# Performing PCA
seurat_phase_data_2_n <- RunPCA(seurat_phase_data_2_n)
```
```R
# PC_ 1 
# Positive:  ADIRF, CSTB, IGFBP3, UCA1, CD24, FTH1, KRT20, SPINK1, PCDH7, FHL2 
#            ID3, GDF15, LINC01285, BRI3, CRTAC1, GCLC, KRT8, IGFBP2, IGFBP5, TMEM97 
#            MID1, TRIM31, CAMK2N1, SCHLAP1, HSD17B2, SPOCD1, KRT18, GADD45A, C4orf3, C5orf17 
# Negative:  S100A8, S100A9, MDM2, MT1X, CRIP1, RAB21, IGKC, DMKN, S100A2, DEFB1 
#            LINC02484, S100A7, CPM, H19, KRT6A, IGLC2, KRT13, IGHG1, LY6D, DSG3 
#            CTAG2, IGLC1, AC018816.1, MT2A, LGALS7, SERPINB3, LINC02154, S100A14, SERPINB4, IGHG3 
# PC_ 2 
# Positive:  MDM2, DMKN, S100A2, S100A14, S100A8, S100A9, RAB21, S100A16, CRABP2, MT1X 
#            IFI27, DEFB1, PTN, LINC02484, KRT13, IGKC, KRT6A, HMGCS1, AC018816.1, LY6D 
#            S100A7, YWHAZ, H19, LINC00355, PPDPF, DSG3, MSMO1, CPM, MUC1, CDK4 
# Negative:  SRGN, RGS1, HLA-DPB1, VIM, HLA-DPA1, ALOX5AP, C1QC, HLA-DQB1, LAPTM5, C1QA 
#            CORO1A, C1QB, TYROBP, SAMSN1, CD74, ARHGAP15, HLA-DQA1, HLA-DRB1, FCER1G, CCL4 
#            MS4A6A, PTPRC, CXCR4, HLA-DRA, DOCK4, CD53, GPR183, FCGR2A, CELF2, ITGB2 
# PC_ 3 
# Positive:  LGALS1, NUPR1, PTN, IFI27, FAM25A, GPX3, UPK1B, CRCT1, PRAP1, DMKN 
#            MDM2, BOK-AS1, DIRC3-AS1, KRT20, TYROBP, S100A8, C1QB, C1QA, KRT33B, APOE 
#            S100A14, PLA2G2A, AC073114.1, HMGCS1, KRT80, FCER1G, C1QC, ALOX5AP, RGS1, GPNMB 
# Negative:  HS6ST3, SPINK1, EEF1A1, RBMS3, CYP24A1, FABP4, RPL7, TSHZ2, BCAM, SLC7A11 
#           NRXN3, CES1, DPP10, SCN11A, RHEX, RPS3A, SLIT3, SCGB3A1, EVA1C, LINC02506 
#           RACK1, KCTD16, PWRN1, MPPED2, SNHG29, AL109930.1, RPL7A, CNTN3, PELI2, NSG1 
# PC_ 4 
# Positive:  IFI27, SCN11A, SERPINE1, MT1X, CCND1, S100A8, RAB21, TIMP3, ELOVL5, SLC35F1 
#            PHLDA2, DMKN, KYNU, DEFB1, ZFAND2A, FABP5, MDM2, IFITM3, LTO1, HMGCS1 
#            BMP2, CDK4, LDLR, DNAJB1, WFDC2, KRT13, LINC00511, FDPS, IGKC, LEAP2 
# Negative:  TFF2, S100A4, TFF1, C10orf99, UPK2, PSORS1C2, FXYD4, NDUFA4L2, CKB, UPK1A 
#            RNASET2, UPK3A, PNCK, LINC02672, IGFBP2, CSTB, GDF15, ADIRF, FTH1, TFF3 
#            SNCG, PDE10A, CXCL2, PVALB, TBC1D30, SPINK1, AL158206.1, CXCL8, ZNF350-AS1, CADM1 
# PC_ 5 
# Positive:  LRP1B, KRT13, C5orf17, RBMS3, TSHZ2, ERBB4, RIMS2, EVA1C, SLC7A11, NRXN3 
#            FABP4, SLC35F1, AL589693.1, LINC02506, HS6ST3, DPP10, CES1, AL109930.1, SNHG29, KCTD16 
#            CYP24A1, NKAIN2, RPL12, DLC1, SLIT3, PDE7B, MPPED2, TRAC, DIRC3-AS1, PELI2 
# Negative:  CXCL8, SULT1E1, GSTM3, GDF15, LTO1, LEAP2, PLAU, INSIG1, BPGM, CD74 
#            TMEM97, TTTY14, FTH1, CXCL1, ISG15, UCA1, HMGCS2, CA2, CITED4, HLA-DRB5 
#            ABTB2, RNASET2, FBLN1, ALDOC, BMP2, GPC5, LCN2, EDN1, HLA-DRA, IFI27
```
```R
# Plotting the PCA colored by cell cycle phase
no_split_2_n <- DimPlot(seurat_phase_data_2_n,
                    reduction = "pca",
                    group.by= "Phase")

no_split_2_n 

with_split_2_n <- DimPlot(seurat_phase_data_2_n,
                      reduction = "pca",
                      group.by= "Phase",
                      split.by= "Phase")

with_split_2_n

no_split_2_n + with_split_2_n
```
<img width="958" alt="PCA_dta_2_n" src="https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/92453ace-d508-4dda-9366-892b240f35fe">

```R
# Mitochondrial Expression:
# Plotting the PCA colored by mitochondrial expression
no_split_2_n_mito <- DimPlot(seurat_phase_data_2_n,
                    reduction = "pca",
                    group.by= "mitoFr")
with_split_2_n_mito <- DimPlot(seurat_phase_data_2_n,
                      reduction = "pca",
                      group.by= "mitoFr",
                      split.by= "mitoFr")
no_split_2_n_mito + with_split_2_n_mito
```
![Mito_Data_2_n](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/assets/133680893/a979b4a7-ddd8-42bc-bbf3-fb5cb0dc34b8)

Based on the above plots, we can see that cells are scattered regardless of their cell cycle phase and mitochondrial genes expression level. So there is no need to regress out the effect of cell cycle and mitochondrial expression in this dataset.

# SCTransform
This function is useful for normalization and regressing out sources of unwanted variation at the same time.The method constructs a generalized linear model (GLM) for each gene, using UMI counts as the response variable and sequencing depth as the explanatory variable. To handle the fact that different genes have different levels of expression, information is pooled across genes with similar abundances, resulting in more accurate parameter estimates.

This regularization process yields residuals, which represent effectively normalized data values that are no longer correlated with sequencing depth.

This method is more accurate method of normalizing, estimating the variance of the raw filtered data, and identifying the most variable genes. In practice SCTransform single command replaces ```NormalizeData()```, ```ScaleData()```, and ```FindVariableFeatures()```. Since we have two group of sample we will run SCTransform on each groups after doing "integration".

# Integration
To improve clustering and downstream analyses, it can be beneficial to integrate or align samples across groups using shared highly variable genes. If cells cluster by sample, condition, batch, dataset, or modalities(scRNA, scATAC-seq), integration can help to remove these unwanted sources of variation.

For example, if we want to integrate normal samples together and BLCA samples together, we should keep each sample as a separate object and transform them accordingly for integration. This is necessary to ensure that the samples are properly aligned and that downstream analyses are meaningful. If cell types are present in one dataset, but not the other, then the cells will still appear as a separate sample-specific cluster.


## With the Complete Data Set: (Can be modified upon completing the previous steps) 
Here all the 8 samples are used. 

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
```R
# Compute percent mito ratio
merged_seurat_data_2$mitoRatio <- PercentageFeatureSet(object = merged_seurat_data_2, pattern = "^MT-")
merged_seurat_data_2$mitoRatio <- merged_seurat_data_2@meta.data$mitoRatio / 100

# adding cell column
merged_seurat_data_2$cells <- rownames(merged_seurat_data_2@meta.data)
# re-setting the rownames
rownames(merged_seurat_data_2@meta.data) <- merged_seurat_data_2@meta.data$cells
```
## Filteration:
```R
# Filteration
filtered_seurat_data_2 <- subset(merged_seurat_data_2, 
                          subset= nCount_RNA >= 1000 &
                                  nFeature_RNA <= 6000 & 
                                  mitoRatio < 0.10)
```
## Normalization and Integration: 
```R
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
```R
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
```R
merged_seurat_data_2 <- RunTSNE(merged_seurat_data_2, assay = "SCT", npcs = 50)
```
```R
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
```R
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
```R
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
```R
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
```R
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
```R
# Cluster 0
png(filename = "harmony_blca_clsuter_markers_cluster0_data_2.png", width = 16, height = 8.135, units = "in", res = 300)
plotList[[1]]
dev.off()
```
[Cluster0](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122897224) ; 
[Cluster1](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122898703) ; 
[Cluster2](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122898746) ; 
[Cluster3](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899321) ; 
[Cluster4](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899358) ; 
[Cluster5](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899411) ; 
[Cluster6](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899517) ; 
[Cluster7](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899536) ; 
[Cluster8](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899548) ; 
[Cluster9](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899574) ; 
[Cluster10](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899589) ; 
[Cluster11](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899698) ; 
[Cluster12](https://github.com/Saindhabi17/SC_RNA_Repo_Data_2/commit/e78e7852dba825a293ad641675d6a3bd9d5a079d#r122899757)

```R
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
```R
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



