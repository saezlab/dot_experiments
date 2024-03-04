source("process_data.R")
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)
st_id <- args[[1]] # e.g., "G2"

sc_data <- load_sc_wu('HER2+')
st_data <- load_st_her2p(st_id)

weights <- run_DOT(sc_data$counts, sc_data$meta[, "celltype_major"],
                   st_data$counts, st_data$meta[, c('x', 'y')], max_spot_size = 200)