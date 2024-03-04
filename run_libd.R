source('process_data.R')
source("methods.R")

args <- commandArgs(trailingOnly = TRUE)
ref_type <- args[[1]] # e.g., "brain"
sample_index <- args[[2]] # e.g., "1"
method <- args[[3]] # e.g., DOT, RCTD, CARD, ...

experiment_data <- load_libd_experiment(ref_type, sample_index)
annotations <- experiment_data$ref$meta[, "layer"]
locations <- experiment_data$srt$meta[, c('Row', 'Col')]

if(method == "DOT")
{
  weights <- run_DOT(experiment_data$ref$counts, annotations,
                     experiment_data$srt$counts, locations, 
                     ifelse(ref_type == 'aggregated', 0.25, 1))
  
}else if(method == "RCTD")
{
  weights <- run_RCTD(experiment_data$ref$counts, annotations,
                      experiment_data$srt$counts, locations)
}else if(method == "CARD")
{
  weights <- run_CARD(experiment_data$ref$counts, annotations,
                      experiment_data$srt$counts, locations)
}