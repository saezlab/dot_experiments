source("process_data.R")
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[[1]] # e.g., "mouse1_sample1"
method <- args[[1]] # e.g., DOT, RCTD, SPOTlight, CARD, ...

experiment_data <- produce_mop_multicell(sample_id)
annotations <- experiment_data$ref$meta[, "cluster_label"]
locations <- experiment_data$srt$meta[, c('x', 'y')]

if(method == "DOT"){
  weights <- run_DOT(experiment_data$ref$counts, annotations,
                     experiment_data$srt$counts, locations)
}else if(method == "RCTD"){
  weights <- run_RCTD(experiment_data$ref$counts, annotations,
                      experiment_data$srt$counts, locations)
}else if(method == "SPOTlight"){
  weights <- run_SPOTlight(experiment_data$ref$counts, annotations,
                      experiment_data$srt$counts, n_cells = 100)
}else if(method == "CARD"){
  weights <- run_CARD(experiment_data$ref$counts, annotations,
                     experiment_data$srt$counts, locations)
}else if(method == "SingleR"){
  weights <- run_SingleR(experiment_data$ref$counts, annotations,
                        experiment_data$srt$counts)
}else if(method == "Seurat"){
  weights <- run_Seurat(experiment_data$ref$counts, annotations,
                         experiment_data$srt$counts)
}else if(method == "Harmony"){
  weights <- run_Harmony(experiment_data$ref$counts, annotations,
                         experiment_data$srt$counts)
}else if(method == "RF"){
  weights <- run_RF(experiment_data$ref$counts, annotations,
                         experiment_data$srt$counts)
}