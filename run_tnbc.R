source("process_data.R")
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)
vis_id <- args[[1]] # e.g., "1142243F"

sc_data <- load_sc_wu('TNBC')
vis_data <- load_vis_tnbc(vis_id)

weights <- run_DOT(sc_data$counts, sc_data$meta[, "celltype_major"],
                   vis_data$counts, vis_data$meta[, c("Col", "Row")],
                   normalize = FALSE)