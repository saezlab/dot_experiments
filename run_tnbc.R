source("process_data.R")
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)
vis_id <- args[[1]] # e.g., "1142243F"
method <- args[[2]] # e.g., DOT, CARD, RCTD

sc_data <- load_sc_wu('TNBC')
vis_data <- load_vis_tnbc(vis_id)
annotations <- sc_data$metaa[, "celltype_major"]
locations <- vis_data$meta[, c("Col", "Row")]

if(method == "DOT")
{
  weights <- run_DOT(sc_data$counts, annotations, vis_data$counts, locations, normalize = FALSE)
}else if(method == "CARD")
{
  weights <- run_CARD(sc_data$counts, annotations, vis_data$counts, locations)  
}else if(method == "RCTD")
{
  weights <- run_RCTD(sc_data$counts, annotations, vis_data$counts, locations)
}
