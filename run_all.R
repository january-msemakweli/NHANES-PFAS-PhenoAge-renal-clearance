# Reproduce the manuscript analysis.
# From this folder: Rscript run_all.R
#
# Writes:
#   data/nhanes_raw.rds          (cached NHANES extract)
#   figures/fig1_hex_facet_CKD.jpg
#   figures/fig2_forest_marginal_effects.jpg
#   figures/fig3_mice_convergence.jpg
#   output/tables/table2_mixture.csv
#   output/tables/table3_ckd_interaction.csv
#   output/tables/table4_components.csv
#   output/tables/single_pollutant.csv

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- if (length(file_arg) && nzchar(file_arg[1])) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
setwd(root)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

message("Working directory: ", getwd())
source("R/setup.R")
message("=== 01 download ===")
source("01_download.R")
message("=== 02 construct ===")
source("02_construct.R")
message("=== 03 impute ===")
source("03_impute.R")
message("=== 04 models ===")
source("04_models.R")
message("=== 05 figures and tables ===")
source("05_figures_tables.R")
message("=== write outputs ===")
source("R/write_outputs.R")
message("Done.")
