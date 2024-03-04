library(Matrix)

load_st_her2p <- function(st_id)
{
  data_dir_st <- "../Andersson_etal_2021"

  st_counts <- read.csv(sprintf("%s/ST-cnts/%s.tsv.gz", data_dir_st, st_id), 
                        sep = "\t", row.names = 1)
  st_counts <- as(t(as.matrix(st_counts)), "Matrix")
  
  st_meta <- read.csv(sprintf("%s/ST-spotfiles/%s_selection.tsv", data_dir_st, st_id), 
                      sep = "\t")
  rownames(st_meta) <- sprintf("%dx%d", st_meta$x, st_meta$y)
  
  pat_file <- sprintf("%s/ST-pat/lbl/%s_labeled_coordinates.tsv", data_dir_st, st_id)
  
  if(file.exists(pat_file))
  {
    pat <- read.csv(pat_file, sep = "\t")
    
    pat <- pat[which(!is.na(pat$x)), ]
    pat[, c("x2", "y2")] <- round(pat[, c("x", "y")])
    rownames(pat) <- sprintf("%dx%d", pat$x2, pat$y2)
    
    pat <- pat[rownames(st_meta), ]
    
    st_meta$pat_label <- pat$label
  }
  
  common_spots <- intersect(rownames(st_meta), colnames(st_counts))
  st_meta <- st_meta[common_spots, ]
  st_counts <- st_counts[, common_spots]
  
  return(list(counts = st_counts, meta = st_meta))
}

load_vis_tnbc <- function(vis_id)
{
  data_dir_st <- "../Wu_etal_2021_vis"
  
  vis_counts <- sprintf("%s/filtered_count_matrices/%s_filtered_count_matrix/", 
                        data_dir_st, vis_id)
  vis_counts <- Seurat::Read10X(vis_counts, gene.column = 1)
  
  vis_meta <- read.csv(sprintf("%s/metadata/%s_metadata_processed.csv", 
                               data_dir_st, vis_id), row.names = 1)
  
  vis_count <- vis_count[, rownames(vis_meta)]
  
  return(list(counts = vis_count, meta = vis_meta))
}

load_sc_wu <- function(subtype = "HER2+")
{
  data_dir_sc <- "../Wu_etal_2021_sc"
  
  sc_meta <- read.csv(sprintf("%s/metadata.csv", data_dir_sc), row.names = 1)
  
  if(subtype != 'all')
  {
    sc_meta <- sc_meta[which(sc_meta$subtype == subtype), ]
  }
  
  sc_count <- Seurat::Read10X(data_dir_sc, gene.column = 1)
  
  sc_count <- sc_count[, rownames(sc_meta)]
  
  return(list(counts = sc_count, meta = sc_meta))
}

load_libd_sample <- function(sample_name)
{
  libd_data_dir <- "../LIBD"
  
  layers <- read.csv(sprintf("%s/barcode_level_layer_map.tsv", libd_data_dir), sep = "\t", header = F)
  colnames(layers) <- c("spot_name", "sample_name", "layer")
  
  sample_counts <- Seurat::Read10X_h5(sprintf("%s/%s/filtered_feature_bc_matrix.h5", libd_data_dir, sample_name))
  
  sample_layers <- layers[which(layers$sample_name == sample_name), ]
  rownames(sample_layers) <- sample_layers$spot_name
  
  loc <- read.csv(sprintf("%s/%s/tissue_positions_list.txt", libd_data_dir, sample_name), 
                  header = F, row.names = 1)
  loc <- loc[which(loc$V2 == 1),  ]
  loc <- setNames(loc[, 2:3], c("Row", "Col"))
  
  spots <- intersect(rownames(loc), rownames(sample_layers))
  loc <- cbind(loc[spots, ], sample_layers[spots, "layer", drop = FALSE])
  loc$sample <- sample_name
  sample_counts <- sample_counts[, spots]
  rs <- rowSums(sample_counts)
  sample_counts <- sample_counts[which(rs > 0), ]
  
  colnames(sample_counts) <- rownames(loc) <- sprintf("%s_%s", sample_name, spots)
  
  return(list(counts = sample_counts, meta = loc))
}

load_libd_batch <- function(sample_names)
{
  dt <- c()
  for(sample_name in sample_names)
  {
    sample_data <- load_libd_sample(sample_name)
    
    if(is.null(dt))
    {
      dt <- sample_data
    }else
    {
      cg <- intersect(rownames(dt$counts), rownames(sample_data$counts))
      dt$counts <- cbind(dt$counts[cg, ], sample_data$counts[cg, ])
      dt$meta <- rbind(dt$meta, sample_data$meta)
    }
  }
  
  return(dt)
}

load_libd_experiment <- function(ref, srt)
{
  libd_data_dir <- "../LIBD"
  
  mtdata <- read.csv(sprintf("%s/libd_clinicaldata.csv", libd_data_dir))
  rownames(mtdata) <-  mtdata$SlideID <- as.character(mtdata$SlideID)
  
  sample_names <- mtdata$SlideID
  
  srt_samples <- c()
  if(srt %in% sample_names)
  {
    srt_samples <- srt
  }else
  {
    suppressWarnings(srt_samples <- sample_names[as.integer(srt)])
  }
  
  ref_samples <- c()
  if(ref == 'aggregated')
  {
    ref_samples <- sample_names
  }else if(ref == 'brain')
  {
    ref_samples <- sample_names[which(mtdata$Brain.Number == mtdata[srt_samples, "Brain.Number"])]
  }else if(ref == 'adjacent')
  {
    ref_samples <- sample_names[which(mtdata$Brain.Number == mtdata[srt_samples, "Brain.Number"] &
                                        mtdata$Position == mtdata[srt_samples, "Position"])]
  }else{
    stop("Invalid ref.")
  }
  
  ref_samples <- setdiff(ref_samples, srt_samples)
  return(list(ref = load_libd_batch(ref_samples), 
              srt = load_libd_batch(srt_samples)))
}

produce_mop_multicell <- function(sample_id = "mouse1_sample1", sc_data = NULL, multicell_meta = NULL, 
                                  seed = 1, diameter = 100, ann_column = "cluster_label")
{
  mop_dir <- "../MOp data"
  
  if(is.null(multicell_meta))
  {
    multicell_meta <- read.csv(sprintf("%s/multicell_meta.csv", mop_dir), 
                               colClasses = c("X" = "character"), row.names = 1)
  }
  
  sc_col <- sprintf("sc_seed_%d", seed)
  if(! (sc_col %in% colnames(multicell_meta)))
    stop("Seed not found!")
  
  d_col <- sprintf("spot_%g", diameter)
  if(! (d_col %in% colnames(multicell_meta)))
    stop("Diameter not found!")
  
  if(is.null(sc_data))
  {
    sc_data <- readRDS(sprintf("%s/scRNA_10X_v2_A/processed_data.rds", mop_dir))
  }
  
  mf_cells <- which(multicell_meta$sample_id == sample_id)
  spots <- multicell_meta[mf_cells, d_col]
  sc_cells <- multicell_meta[mf_cells, sc_col]
  
  srt_counts <- sc_data$counts[, sc_cells]
  srt_counts <- as.matrix(t(srt_counts))
  srt_counts <- t(rowsum(srt_counts, spots))
  srt_counts <- as(srt_counts, "Matrix")
  
  srt_composition <- as.data.frame.matrix(table(spots, sc_data$meta[sc_cells, ann_column]))
  
  srt_meta <- stringr::str_split(rownames(srt_composition), "X")
  srt_meta <- t(sapply(srt_meta, as.numeric))
  srt_meta <- setNames(as.data.frame(srt_meta), c("x", "y"))
  
  srt_id <- as.data.frame.matrix(table(spots, multicell_meta[mf_cells, "slice_id"]))
  srt_meta$slice_id <- colnames(srt_id)[apply(srt_id, 1, which.max)]
  
  remaining_sc <- setdiff(rownames(sc_data$meta), sc_cells)
  
  return(list(ref = list(counts = sc_data$counts[, remaining_sc], 
                         meta = sc_data$meta[remaining_sc, ]),
              srt = list(counts = srt_counts, 
                         composition = srt_composition,
                         meta = srt_meta)))
}
