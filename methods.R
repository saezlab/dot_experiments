library(Matrix)

run_DOT <- function(ref_data, ref_annotations, srt_data, srt_coords, ratios_weight = 0.25, max_spot_size = 20)
{
  dot_ref <- DOT::setup.ref(ref_data, ref_annotations, max_genes = 10000)
  dot_srt <- DOT::setup.srt(srt_data, srt_coords)
  dot <- DOT::create.DOT(dot_srt, dot_ref)
  dot <- DOT::run.DOT.lowresolution(dot, ratios_weight = ratios_weight, max_spot_size = max_spot_size)
  
  return(dot@weights)
}

run_Harmony <- function(ref_data, ref_annotations, srt_data, K = 10)
{
  library(Seurat)
  library(dplyr)
  library(harmony)
  
  common_genes <- intersect(rownames(ref_data), rownames(srt_data))
  ref_data <- ref_data[common_genes, ]
  srt_data <- srt_data[common_genes, ]
  
  integrated_data <- 
    CreateSeuratObject(counts = cbind(ref_data, srt_data), project = "Agg", min.cells = 5) %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  integrated_data <- integrated_data  %>% 
    ScaleData(verbose = F) %>% 
    RunPCA(pc.genes = integrated_data@var.genes, npcs = 20, verbose = F)
  
  integrated_data@meta.data$technology <- c(rep("ref", ncol(ref_data)), rep("srt", ncol(srt_data)))
  integrated_data <- integrated_data %>% 
    RunHarmony("technology", plot_convergence = F)
  
  harmony_embeddings <- Embeddings(integrated_data, 'harmony')
  
  harmony_ref <- harmony_embeddings[1:ncol(ref_data), ]
  harmony_srt <- harmony_embeddings[(ncol(ref_data)+1):nrow(harmony_embeddings), ]
  
  srt_ref_dists <- fields::rdist(harmony_srt, harmony_ref)
  celltypes <- unique(ref_annotations)
  weights <- matrix(0, ncol(srt_data), length(celltypes))
  colnames(weights) <- celltypes
  for(i in 1:ncol(srt_data))
  {
    ct_i <- DOT:::int.table(ref_annotations[order(srt_ref_dists[i, ])[1:K]], TRUE)
    weights[i, names(ct_i)] <- ct_i
  }
  
  return(weights)
}

run_Seurat <- function(ref_data, ref_annotations, srt_data)
{
  library(Seurat)
  library(dplyr)
  
  ref_data <- 
    CreateSeuratObject(ref_data) %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", 
                         nfeatures = 2000, verbose = FALSE) %>%
    SCTransform()
  
  srt_data <- CreateSeuratObject(srt_data) %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", 
                         nfeatures = 2000, verbose = FALSE) %>%
    ScaleData() %>%
    RunPCA()
  
  anchors <- FindTransferAnchors(reference = ref_data, 
                                 query = srt_data, normalization.method = "SCT",
                                 npcs = 50)
  weights <- TransferData(anchorset = anchors, refdata = ref_annotations,
                                    weight.reduction = srt_data[["pca"]],
                                    prediction.assay = TRUE, we, dims = 1:50)
  weights <- t(weights@data)
  weights <- weights[, -ncol(weights)]
  return(weights)
}

run_SingleR <- function(ref_data, ref_annotations, srt_data)
{
  library(Seurat)
  library(dplyr)
  library(SingleR)
  
  ref_data <- 
    CreateSeuratObject(ref_data) %>%
    NormalizeData(verbose = F) %>%
    GetAssayData()
  
  srt_data <- 
    CreateSeuratObject(srt_data) %>%
    NormalizeData(verbose = F) %>%
    GetAssayData()
  
  pred <- SingleR(srt_data, ref_data, labels = ref_annotations)
  # pred <- data.frame(label = as.factor(pred$labels))
  pred <- model.matrix(~labels-1, data = pred)
  colnames(pred) <- stringr::str_sub(colnames(pred), nchar("labels")+1)
  
  celltypes <- unique(ref_annotations)
  weights <- matrix(0, ncol = length(celltypes), nrow = ncol(srt_data),
                    dimnames = list(colnames(srt_data), celltypes))
  weights[rownames(pred), colnames(pred)] <- pred
  return(weights)
}

run_SPOTlight <- function(ref_data, ref_annotations, srt_data, n_cells = 100)
{
  library(SPOTlight)
  library(Seurat)
  library(dplyr)
  
  common_genes <- intersect(rownames(ref_data), rownames(srt_data))
  ref_data <- ref_data[common_genes, ]
  srt_data <- srt_data[common_genes, ]
  
  ref_data <- 
    CreateSeuratObject(ref_data) %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", 
                         nfeatures = 3000, verbose = FALSE)
  
  Idents(ref_data) <- ref_annotations
  hvg <- VariableFeatures(ref_data)
  
  idx <- split(seq(ncol(ref_data)), ref_annotations)
  cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    sample(i, n_cells)
  })
  
  ref_data <- ref_data[hvg, unlist(cs_keep)]

  srt_data <- 
    CreateSeuratObject(srt_data) %>%
    NormalizeData(verbose = F)
  
  srt_data <- srt_data[hvg, ]
  
  ref_markers <- FindAllMarkers(ref_data, only.pos = TRUE)
  ref_markers <- ref_markers[which(ref_markers$avg_log2FC > 1), ]
  
  res <- SPOTlight::SPOTlight(
    x = GetAssayData(ref_data),
    y = GetAssayData(srt_data),
    groups = Idents(ref_data),
    mgs = ref_markers,
    hvg = hvg,
    weight_id = "avg_log2FC",
    group_id = "cluster",
    gene_id = "gene")
  
  return(res$mat)
}

run_CARD <- function(ref_data, ref_annotations, srt_data, srt_coords)
{
  sc_meta <- data.frame(Indetity = ref_annotations)
  sc_meta$sampleInfo <- 'sample1'
  rownames(sc_meta) <- colnames(ref_data)
  
  CARD_obj <- CARD::createCARDObject(
    sc_count = ref_data, # gene x cell
    sc_meta = sc_meta, # data.frame, with "cellType" and "sampleInfo" columns
    spatial_count = srt_data, # gene x location
    spatial_location = srt_coords, # location x (x,y)
    ct.varname = "Indetity",
    ct.select = unique(sc_meta$Indetity),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5)
  
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)
  
  return(CARD_obj@Proportion_CARD)
}

run_RCTD <- function(ref_data, ref_annotations, srt_data, srt_coords)
{
  ref_annotations <- as.factor(ref_annotations)
  names(ref_annotations) <- colnames(ref_data)
  nUMI <- colSums(counts)
  levels(ref_annotations) <- gsub("/", "-", levels(ref_annotations))
  
  query <- spacexr::SpatialRNA(srt_coords, srt_data, colSums(srt_data))
  ref <- spacexr::Reference(ref_data, ref_annotations, nUMI)
  
  RCTD <- spacexr::create.RCTD(query, ref, max_cores = 1, UMI_min = 50)
  RCTD <- spacexr::run.RCTD(RCTD, doublet_mode = "full")
  
  return(RCTD@results$results_df)
}

run_RF <- function(ref_data, ref_annotations, srt_data)
{
  library(ranger)
  
  common_genes <- intersect(rownames(ref_data), rownames(srt_data))
  ref_data <- ref_data[common_genes, ]
  srt_data <- srt_data[common_genes, ]
  
  rf <- ranger(x = t(ref_data), y = as.factor(ref_annotations), 
               classification = TRUE, probability = TRUE)
  prd <- predict(rf, t(srt_data))
  rownames(prd$predictions) <- colnames(srt_data)
  
  return(prd$predictions)
}

# Python:
#  Cell2location
#  NovoSpaRc
#  TACCO
#  Tangram

