##########################################################
### ENVIRONMENT CONFIGURATION PHASE ######################
##########################################################
rm(list = objects())
gc(reset = TRUE, full = TRUE)

Sys.setenv(R_MAKEVARS_USER = "~/.R/Makevars")
options(stringsAsFactors = FALSE)
set.seed(12345)

##########################################################
### VISUALIZATION PARAMETERS INITIALIZATION ##############
##########################################################
color_palette_set_1 <- c("#F19294", "#66C2A5", "#AEDDEE", "#80B1D3", "#3477A9", 
                        "#CCEBC5", "#A4D38E", "#4A9D47", "#F5B375", "#BADA90", 
                        "#FEE08B", "#F58135", "#FACB7B", "#96C3D8", "#C6DBEF",
                        "#BDA7CB", "#B4B1AE", "#00A43C", "#FDAE61", "#E6F598",
                        "#FFFFB3", "#B3DE69", "#FCCDE5")

color_palette_set_2 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                        "#CCEBC5", "#FFED6F")

##########################################################
### PACKAGE LOADING SEQUENCE #############################
##########################################################
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(tidyverse)
if(!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
library(Seurat)
if(!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)
if(!requireNamespace("harmony", quietly = TRUE)) install.packages("harmony")
library(harmony)

##########################################################
### DATA INGESTION PROCESS ###############################
##########################################################
load("DatasetA_scRNA.Rda")
DatasetA_obj <- get(ls()[grep("scRNA", ls())])
DatasetA_obj$StudyGroup <- "DatasetA"
rm(list = ls()[!ls() %in% "DatasetA_obj"])

load("DatasetB_scRNA.Rda")
DatasetB_obj <- get(ls()[grep("scRNA", ls())])
DatasetB_obj$StudyGroup <- "DatasetB"
rm(list = ls()[!ls() %in% c("DatasetA_obj", "DatasetB_obj")])

load("InternalData_scRNA.Rda")
Internal_obj <- get(ls()[grep("scRNA", ls())])
Internal_obj$StudyGroup <- "Internal"
Internal_obj$TissueSite <- "Neoplasm"
rm(list = ls()[!ls() %in% c("DatasetA_obj", "DatasetB_obj", "Internal_obj")])

##########################################################
### DATA AGGREGATION PROCEDURE ###########################
##########################################################
merged_scRNA <- merge(x = DatasetA_obj, y = list(DatasetB_obj, Internal_obj))

##########################################################
### METADATA CURATION WORKFLOW ###########################
##########################################################
expression_matrix <- merged_scRNA@assays$RNA@counts
metadata_table <- merged_scRNA@meta.data

srt_aggregated <- CreateSeuratObject(
  counts = expression_matrix,
  meta.data = metadata_table,
  min.cells = 3,
  min.features = 200
)

##########################################################
### QUALITY CONTROL IMPLEMENTATION #######################
##########################################################
srt_aggregated$pMT <- PercentageFeatureSet(
  object = srt_aggregated,
  pattern = "^MT-"
)[,1]

srt_aggregated$pHB <- PercentageFeatureSet(
  object = srt_aggregated,
  pattern = "^HBA|^HBB"
)[,1]

srt_aggregated$pRP <- PercentageFeatureSet(
  object = srt_aggregated,
  pattern = "^RPS|^RPL"
)[,1]

##########################################################
### VISUALIZATION OUTPUT GENERATION ######################
##########################################################
create_violin_plots <- function(seurat_obj) {
  vplot_1 <- VlnPlot(
    object = seurat_obj,
    features = "nCount_RNA",
    group.by = "orig.ident",
    pt.size = 0
  ) + theme(legend.position = "none")
  
  vplot_2 <- VlnPlot(
    object = seurat_obj,
    features = "nFeature_RNA",
    group.by = "orig.ident",
    pt.size = 0
  ) + theme(legend.position = "none")
  
  vplot_3 <- VlnPlot(
    object = seurat_obj,
    features = "pMT",
    group.by = "orig.ident",
    pt.size = 0
  ) + theme(legend.position = "none")
  
  CombinePlots(plots = list(vplot_1, vplot_2, vplot_3), ncol = 3)
}

generate_ridge_plots <- function(seurat_obj) {
  RidgePlot(
    object = seurat_obj,
    features = c("nCount_RNA", "nFeature_RNA", "pMT"),
    log = TRUE
  ) + theme_minimal()
}

pdf("QC_Figures.pdf", width = 16, height = 12)
print(create_violin_plots(srt_aggregated))
print(generate_ridge_plots(srt_aggregated))
dev.off()

##########################################################
### DATA FILTERING IMPLEMENTATION ########################
##########################################################
filtered_scRNA <- subset(
  x = srt_aggregated,
  subset = nCount_RNA > 1000 &
    nCount_RNA < 50000 &
    nFeature_RNA > 800 &
    nFeature_RNA < 5000 &
    pMT < 5
)

##########################################################
### CELL CYCLE PHASE DETERMINATION #######################
##########################################################
cell_cycle_genes <- list(
  S.genes = cc.genes$s.genes,
  G2M.genes = cc.genes$g2m.genes
)

filtered_scRNA <- CellCycleScoring(
  object = filtered_scRNA,
  s.features = cell_cycle_genes$S.genes,
  g2m.features = cell_cycle_genes$G2M.genes
)

##########################################################
### DATA TRANSFORMATION PIPELINE #########################
##########################################################
sct_variables <- c("S.Score", "G2M.Score")

filtered_scRNA <- SCTransform(
  object = filtered_scRNA,
  vars.to.regress = sct_variables,
  variable.features.n = 3000,
  return.only.var.genes = FALSE
)

##########################################################
### DIMENSIONALITY REDUCTION WORKFLOW ####################
##########################################################
pca_results <- RunPCA(
  object = filtered_scRNA,
  npcs = 50,
  verbose = FALSE
)

harmony_results <- RunHarmony(
  object = pca_results,
  group.by.vars = "StudyGroup",
  reduction.save = "harmony"
)

umap_coordinates <- RunUMAP(
  object = harmony_results,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 100
)

tsne_projections <- RunTSNE(
  object = harmony_results,
  reduction = "harmony",
  dims = 1:30,
  perplexity = 150
)

##########################################################
### CLUSTERING OPTIMIZATION PROCEDURE ####################
##########################################################
resolution_parameters <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for (resolution_value in resolution_parameters) {
  umap_coordinates <- FindNeighbors(
    object = umap_coordinates,
    reduction = "harmony",
    dims = 1:30
  )
  
  umap_coordinates <- FindClusters(
    object = umap_coordinates,
    resolution = resolution_value
  )
  
  cluster_plot <- DimPlot(
    object = umap_coordinates,
    reduction = "umap",
    label = TRUE
  ) + ggtitle(paste("Resolution:", resolution_value))
  
  pdf(paste0("Cluster_Resolution_", resolution_value, ".pdf"), width = 8, height = 6)
  print(cluster_plot)
  dev.off()
}

##########################################################
### SUBSET ANALYSIS WORKFLOW #############################
##########################################################
selected_celltypes <- c("Lymphoid", "Myeloid", "Stromal")

for (celltype in selected_celltypes) {
  subset_data <- subset(
    x = umap_coordinates,
    subset = CellCategory == celltype
  )
  
  subset_data <- SCTransform(subset_data)
  subset_data <- RunPCA(subset_data)
  subset_data <- RunHarmony(subset_data, group.by.vars = "StudyGroup")
  subset_data <- RunUMAP(subset_data, dims = 1:30)
  
  marker_features <- switch(celltype,
    "Lymphoid" = c("CD3D", "CD8A", "NKG7", "GNLY"),
    "Myeloid" = c("CD14", "FCGR3A", "LYZ", "CST3"),
    "Stromal" = c("COL1A1", "DCN", "LUM", "ACTA2")
  )
  
  feature_plot <- FeaturePlot(
    object = subset_data,
    features = marker_features,
    reduction = "umap"
  )
  
  pdf(paste0(celltype, "_Marker_Expression.pdf"), width = 12, height = 9)
  print(feature_plot)
  dev.off()
}

##########################################################
### DATA EXPORTATION PROCESS #############################
##########################################################
output_objects <- list(
  FullDataset = srt_aggregated,
  ProcessedData = umap_coordinates,
  CellSubsets = lapply(selected_celltypes, function(ct) {
    subset(x = umap_coordinates, subset = CellCategory == ct)
  })
)

saveRDS(object = output_objects, file = "CompleteAnalysisResults.rds")

##########################################################
### CELL SUBSET EXTRACTION PROTOCOL ######################
##########################################################
filtered_clusters <- c(0, 1, 3, 4, 5, 6)
subset_criteria <- paste0("SCT_snn_res.0.1 %in% c(", paste(filtered_clusters, collapse = ","), ")")

TNK_subpopulation <- subset(
  x = subsc,
  subset = eval(parse(text = subset_criteria))
)

print(paste("Initial dimensions:", paste(dim(subsc), collapse = " x ")))
print(paste("Filtered dimensions:", paste(dim(TNK_subpopulation), collapse = " x ")))

##########################################################
### DATA PERSISTENCE OPERATIONS ##########################
##########################################################
save(TNK_subpopulation, file = "Processed_TNK_Subset.Rda")
remove(list = objects())
gc(verbose = FALSE, full = TRUE)

##########################################################
### VISUALIZATION PARAMETERS CONFIGURATION ###############
##########################################################
color_definitions <- list(
  blue_tone = "#96C3D8",
  red_tone = "#F19294",
  green_tone = "#66C2A5",
  gray_tone = "#D9D9D9",
  extended_palette = c("#E7E1EF", "#C994C7", "#66C2A5", "#AEDDEE", 
                      "#80B1D3", "#3477A9", "#CCEBC5", "#A4D38E",
                      "#4A9D47", "#BADA90", "#FEE08B", "#F58135",
                      "#FACB7B", "#96C3D8", "#C6DBEF", "#F5B375",
                      "#BDA7CB", "#B4B1AE", "#00A43C", "#FDAE61",
                      "#E6F598", "#FFFFB3", "#B3DE69", "#FCCDE5")
)

##########################################################
### DATA REINTEGRATION WORKFLOW ##########################
##########################################################
load("Processed_TNK_Subset.Rda")
load("Integrated_Processed_Data.Rda")

cellular_categories <- list(
  Lymphoid = c("Lymphoid_Cells"),
  Myeloid = c("Myeloid_Cells"),
  Stromal = c("Stromal_Cells"),
  Epithelial = c("Epithelial_Cells"),
  Mast = c("Mast_Cells")
)

subset_collection <- list(
  Lymphoid = subset(Integrated_Processed_Data, CellCategory %in% cellular_categories$Lymphoid),
  Myeloid = subset(Integrated_Processed_Data, CellCategory %in% cellular_categories$Myeloid),
  Stromal = subset(Integrated_Processed_Data, CellCategory %in% cellular_categories$Stromal),
  Epithelial = subset(Integrated_Processed_Data, CellCategory %in% cellular_categories$Epithelial),
  Mast = subset(Integrated_Processed_Data, CellCategory %in% cellular_categories$Mast),
  TNK = TNK_subpopulation
)

##########################################################
### DATA AGGREGATION PROCESS #############################
##########################################################
merged_dataset <- merge(
  x = subset_collection[[1]],
  y = list(
    subset_collection[[2]],
    subset_collection[[3]],
    subset_collection[[4]],
    subset_collection[[5]],
    subset_collection[[6]]
  )
)

##########################################################
### METADATA CURATION PROCEDURE ##########################
##########################################################
metadata_columns <- c("SampleID", "PatientID", "DatasetSource")
curated_metadata <- merged_dataset@meta.data[, metadata_columns, drop = FALSE]

processed_object <- CreateSeuratObject(
  counts = merged_dataset@assays$RNA@counts,
  meta.data = curated_metadata
)

##########################################################
### QUALITY METRICS CALCULATION ##########################
##########################################################
mitochondrial_genes <- "^MT-"
hemoglobin_genes <- "^HBA|^HBB"
ribosomal_genes <- "^RPS|^RPL"

processed_object$mitochondrial_percent <- PercentageFeatureSet(
  object = processed_object,
  pattern = mitochondrial_genes
)

processed_object$hemoglobin_percent <- PercentageFeatureSet(
  object = processed_object,
  pattern = hemoglobin_genes
)

processed_object$ribosomal_percent <- PercentageFeatureSet(
  object = processed_object,
  pattern = ribosomal_genes
)

##########################################################
### CELL CYCLE ASSESSMENT ################################
##########################################################
cell_cycle_markers <- list(
  S_phase = cc.genes$s.genes,
  G2M_phase = cc.genes$g2m.genes
)

processed_object <- CellCycleScoring(
  object = processed_object,
  s.features = cell_cycle_markers$S_phase,
  g2m.features = cell_cycle_markers$G2M_phase
)

##########################################################
### DATA TRANSFORMATION PIPELINE #########################
##########################################################
normalization_variables <- c("mitochondrial_percent")

processed_object <- SCTransform(
  object = processed_object,
  vars.to.regress = normalization_variables,
  variable.features.n = 3000
)

##########################################################
### DIMENSIONAL REDUCTION WORKFLOW #######################
##########################################################
dimensional_reduction_steps <- list(
  PCA = list(
    method = RunPCA,
    parameters = list(npcs = 50)
  ),
  Harmony = list(
    method = RunHarmony,
    parameters = list(group.by.vars = "PatientID")
  ),
  UMAP = list(
    method = RunUMAP,
    parameters = list(dims = 1:30)
  ),
  tSNE = list(
    method = RunTSNE,
    parameters = list(perplexity = 150)
  )
)

for (step in names(dimensional_reduction_steps)) {
  processed_object <- do.call(
    what = dimensional_reduction_steps[[step]]$method,
    args = c(list(object = processed_object), dimensional_reduction_steps[[step]]$parameters)
  )
}

##########################################################
### CLUSTERING OPTIMIZATION PROCEDURE ####################
##########################################################
resolution_parameters <- seq(from = 0.1, to = 1.0, by = 0.1)

cluster_evaluation <- function(seurat_obj, resolutions) {
  for (res in resolutions) {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    print(DimPlot(seurat_obj, reduction = "umap") + ggtitle(paste("Resolution:", res)))
  }
  return(seurat_obj)
}

processed_object <- cluster_evaluation(processed_object, resolution_parameters)

##########################################################
### GENOMIC INSTABILITY ANALYSIS #########################
##########################################################
copykat_analysis <- function(seurat_obj, reference_celltype) {
  expression_matrix <- GetAssayData(seurat_obj, slot = "counts")
  reference_cells <- colnames(subset(seurat_obj, CellCategory == reference_celltype))
  
  analysis_results <- copykat(
    rawmat = expression_matrix,
    id.type = "S",
    ngene.chr = 1,
    win.size = 25,
    KS.cut = 0.1,
    norm.cell.names = reference_cells,
    sam.name = "GenomicAnalysis",
    distance = "euclidean",
    n.cores = 1
  )
  
  return(analysis_results)
}

genomic_results <- copykat_analysis(processed_object, "Reference_CellType")
save(genomic_results, file = "GenomicInstabilityResults.Rda")

##########################################################
### DATA INTERPRETATION WORKFLOW #########################
##########################################################
malignancy_classification <- data.frame(
  CellBarcode = colnames(processed_object),
  GenomicStatus = "Unclassified"
)

malignant_cells <- genomic_results$prediction$cell.names[genomic_results$prediction$copykat.pred == "aneuploid"]
nonmalignant_cells <- genomic_results$prediction$cell.names[genomic_results$prediction$copykat.pred == "diploid"]

malignancy_classification$GenomicStatus[malignancy_classification$CellBarcode %in% malignant_cells] <- "Malignant"
malignancy_classification$GenomicStatus[malignancy_classification$CellBarcode %in% nonmalignant_cells] <- "NonMalignant"

processed_object <- AddMetaData(
  object = processed_object,
  metadata = malignancy_classification$GenomicStatus,
  col.name = "GenomicStatus"
)

##########################################################
### GENOMIC CLASSIFICATION INTEGRATION ###################
##########################################################
prediction_analysis <- list(
  malignant = pred_res$cell.names[pred_res$copykat.pred == "aneuploid"],
  nonmalignant = pred_res$cell.names[pred_res$copykat.pred == "diploid"]
)

classification_matrix <- data.frame(
  CellIdentifier = colnames(expression_matrix),
  GenomicStatus = "Unclassified"
)

classification_matrix$GenomicStatus[classification_matrix$CellIdentifier %in% prediction_analysis$malignant] <- "GenomicAbnormality"
classification_matrix$GenomicStatus[classification_matrix$CellIdentifier %in% prediction_analysis$nonmalignant] <- "GenomicStability"

##########################################################
### METADATA ENHANCEMENT PROCESS #########################
##########################################################
rownames(classification_matrix) <- classification_matrix$CellIdentifier
classification_matrix$CellIdentifier <- NULL

cellular_object <- AddMetaData(
  object = cellular_object,
  metadata = classification_matrix
)

##########################################################
### CELLULAR SUBSET ISOLATION PROTOCOL ###################
##########################################################
epithelial_subpopulation <- subset(
  x = cellular_object,
  subset = CellularCategory == "Epithelial"
)

cellular_counts <- table(epithelial_subpopulation$GenomicStatus)
print(paste("Malignant cells:", cellular_counts["GenomicAbnormality"]))
print(paste("Non-malignant cells:", cellular_counts["GenomicStability"]))

##########################################################
### DATA PERSISTENCE OPERATIONS ##########################
##########################################################
save(epithelial_subpopulation, file = "Processed_Epithelial_Subset.Rda")

##########################################################
### METADATA SYNCHRONIZATION WORKFLOW ###################
##########################################################
load("Integrated_Analysis_Results.Rda")

metadata_reference <- data.frame(
  CellIdentifier = rownames(integrated_object@meta.data),
  CellularCategory = integrated_object@meta.data$CellularCategory,
  GenomicClassification = NA
)

metadata_reference$GenomicClassification[metadata_reference$CellIdentifier %in% prediction_analysis$malignant] <- "AbnormalGenome"

##########################################################
### CLASSIFICATION HIERARCHY DEVELOPMENT ################
##########################################################
classification_system <- metadata_reference %>%
  mutate(FinalClassification = coalesce(GenomicClassification, CellularCategory)) %>%
  column_to_rownames("CellIdentifier") %>%
  select(DefinitiveClassification = FinalClassification)

integrated_object <- AddMetaData(
  object = integrated_object,
  metadata = classification_system
)

##########################################################
### CELLULAR CATEGORIZATION SYSTEM #######################
##########################################################
classification_levels <- c("GenomicStability", "GenomicAbnormality", 
                          "LymphoidCells", "MyeloidCells", "StromalCells")

integrated_object$DefinitiveClassification <- factor(
  integrated_object$DefinitiveClassification,
  levels = classification_levels
)

##########################################################
### ANALYTICAL SUBSET CREATION ###########################
##########################################################
genomic_subset <- subset(
  x = integrated_object,
  subset = DefinitiveClassification %in% c("GenomicStability", "GenomicAbnormality")
)

subset_matrix <- CreateSeuratObject(
  counts = genomic_subset@assays$GeneExpression@counts,
  meta.data = genomic_subset@meta.data,
  min.cellular_observations = 3,
  min.genetic_features = 200
)

##########################################################
### DATA TRANSFORMATION PIPELINE #########################
##########################################################
normalization_variables <- c("TissueOrigin")

subset_matrix <- SCTransform(
  object = subset_matrix,
  variables.to.adjust = normalization_variables,
  variable.genetic_features = 3000
)

##########################################################
### DIMENSIONALITY ANALYSIS WORKFLOW #####################
##########################################################
dimensional_reduction_steps <- list(
  PrincipalComponents = list(
    method = RunPCA,
    parameters = list(
      npcs = 50,
      weight.by.variance = TRUE
    )
  ),
  tSNE_Projection = list(
    method = RunTSNE,
    parameters = list(
      perplexity = 150,
      dimensionality = 1:30
    )
  ),
  UMAP_Projection = list(
    method = RunUMAP,
    parameters = list(
      neighbor_count = 100,
      dimensionality = 1:30
    )
  )
)

for (analysis_step in names(dimensional_reduction_steps)) {
  subset_matrix <- do.call(
    dimensional_reduction_steps[[analysis_step]]$method,
    args = c(list(object = subset_matrix), 
            dimensional_reduction_steps[[analysis_step]]$parameters)
  )
}

##########################################################
### VISUALIZATION OUTPUT GENERATION ######################
##########################################################
generate_dimensional_plot <- function(seurat_obj, reduction_method) {
  DimPlot(
    object = seurat_obj,
    reduction = reduction_method,
    group.by = "SampleOrigin",
    label = FALSE
  ) + theme_minimal()
}

pdf("Dimensional_Visualizations.pdf", width = 11, height = 8.5)
print(generate_dimensional_plot(subset_matrix, "pca"))
print(generate_dimensional_plot(subset_matrix, "tsne"))
print(generate_dimensional_plot(subset_matrix, "umap"))
dev.off()

##########################################################
### EXPRESSION ANALYSIS PIPELINE #########################
##########################################################
aggregated_expression <- aggregate_expression_data(
  seurat_obj = epithelial_subset,
  sample_identifier = "SampleID"
)

gene_expression_profile <- analyze_gene_expression(
  expression_matrix = aggregated_expression,
  target_gene = "DNAReplicationGene"
)

##########################################################
### STRATIFICATION ANALYSIS ##############################
##########################################################
expression_thresholds <- list(
  high_expression = median(gene_expression_profile$LogExpression + 1),
  low_expression = median(gene_expression_profile$LogExpression - 1)
)

stratification_results <- gene_expression_profile %>%
  mutate(ExpressionGroup = case_when(
    LogExpression > expression_thresholds$high_expression ~ "HighExpression",
    LogExpression < expression_thresholds$low_expression ~ "LowExpression",
    TRUE ~ "MidExpression"
  ))

##########################################################
### METADATA INTEGRATION WORKFLOW ########################
##########################################################
integrated_object$ExpressionStratum <- stratification_results$ExpressionGroup[
  match(integrated_object@meta.data$SampleID, stratification_results$SampleID)
]

##########################################################
### ANALYTICAL VISUALIZATION #############################
##########################################################
visualization_palette <- c("#E7E1EF", "#C994C7", "#66C2A5", "#AEDDEE", 
                          "#80B1D3", "#3477A9", "#CCEBC5", "#A4D38E",
                          "#4A9D47", "#BADA90", "#FEE08B", "#F58135",
                          "#FACB7B", "#96C3D8", "#C6DBEF", "#F5B375",
                          "#BDA7CB", "#B4B1AE", "#00A43C", "#FDAE61",
                          "#E6F598", "#FFFFB3", "#B3DE69", "#FCCDE5")

generate_stratified_plot <- function(metadata_table) {
  ggplot(metadata_table, 
         aes(ExpressionStratum, Proportion, 
             fill = CellularCategory, 
             stratum = CellularCategory,
             alluvium = CellularCategory)) +
    geom_col(width = 0.35, color = NA) +
    geom_flow(width = 0.35, alpha = 0.25) +
    scale_fill_manual(values = visualization_palette) +
    theme_minimal() +
    theme(axis.text.radian = element_text(angle = 45)) +
    labs(x = "Expression Stratum", y = "Cellular Proportion")
}

pdf("Stratified_Analysis.pdf", width = 14, height = 10)
print(generate_stratified_plot(integrated_object@meta.data))
dev.off()