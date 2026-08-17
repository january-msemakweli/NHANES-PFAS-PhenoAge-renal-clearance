# Shared setup: seed, packages, config, helpers.
# Sourced by run_all.R before the numbered scripts.

set.seed(20260227)

config <- list(
  cycles = c("2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"),
  exposure_family = "pfas",
  q = 4,
  min_age = 20,
  ckd_def = "egfr_or_acr",
  use_residualized_phenoage = TRUE,
  require_complete_exposure_set = FALSE,
  min_nonmissing_exposure_fraction = 0.6,
  min_nonmissing_exposures_floor = 2,
  include_clinical_covariates = TRUE,
  weight_type = "auto",
  use_mice = TRUE,
  mice_m = 24,
  mice_maxit = 20,
  mice_n_core = 10,
  mice_seed = 20260227,
  verbose = TRUE
)

required_pkgs <- c(
  "nhanesA", "dplyr", "tidyr", "purrr", "stringr", "survey", "splines",
  "ggplot2", "broom", "knitr", "qgcomp", "forcats", "tibble",
  "mice", "mitools", "future", "future.apply", "patchwork", "hexbin"
)
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package not installed: ", p)
}))

library(nhanesA)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(survey)
options(survey.lonely.psu = "adjust")
library(splines)
library(ggplot2)
library(broom)
library(knitr)
library(qgcomp)
library(forcats)
library(tibble)
library(mice)
library(mitools)
library(future)
library(future.apply)
library(patchwork)

source("R/helpers.R")
