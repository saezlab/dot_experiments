source("process_data.R")
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)
st_id <- args[[1]] # e.g., "G2"
method <- args[[2]] # e.g., DOT, CARD, RCTD

sc_data <- load_sc_wu('HER2+')
st_data <- load_st_her2p(st_id)
annotations <- sc_data$metaa[, "celltype_major"]
locations <- st_data$meta[, c('x', 'y')]

if(method == "DOT")
{
  weights <- run_DOT(sc_data$counts, annotations, st_data$counts, locations,
                     max_spot_size = 200, normalize = FALSE) 
}else if(method == "CARD")
{
  weights <- run_CARD(sc_data$counts, annotations, st_data$counts, locations)  
}else if(method == "RCTD")
{
  weights <- run_RCTD(sc_data$counts, annotations, st_data$counts, locations)
}