# Renal Clearance and Albumin Binding May Distort Cross-Sectional Associations Between PFAS and Phenotypic Age Acceleration

Reproducible analysis of a four-chemical PFAS mixture (PFOS, PFOA, PFNA, PFHxS) and PhenoAge acceleration in NHANES 2005-2018, with CKD as an effect modifier.

The public entry point is `run_all.R`.

## Reproduce

```r
install.packages(c(
  "nhanesA", "dplyr", "tidyr", "purrr", "stringr",
  "survey", "splines", "ggplot2", "broom", "knitr",
  "qgcomp", "forcats", "tibble", "mice", "mitools",
  "future", "future.apply", "patchwork", "hexbin"
))
```

```r
setwd("path/to/NHANES-PFA-PhenoAge")
source("run_all.R")
```

An internet connection is required the first time. Later runs reuse `data/nhanes_raw.rds`.

## Scripts

```
R/setup.R                 # seed, packages, config
R/helpers.R               # PhenoAge, eGFR, LOD, NHANES parsers
R/write_outputs.R         # CSV tables + numeric check against the paper
01_download.R             # NHANES extract
02_construct.R            # PhenoAge, CKD, PFAS, weights, mixture score
03_impute.R               # MICE (m = 24); outcome and exposures are never imputed
04_models.R               # qgcomp, Model B, sensitivities, components
05_figures_tables.R       # flow chart, figures
run_all.R
```

## Analysis choices that match the manuscript

- Adults aged 20 and older; cycles with CRP: 2005-2010 and 2015-2018
- PhenoAge from Levine et al. (2018); acceleration is the survey-weighted residual of PhenoAge on chronological age
- CKD: 2021 CKD-EPI eGFR < 60 and/or ACR >= 30
- Adjustment: chronological age, sex, race/ethnicity, PIR, BMI, smoking, alcohol, survey cycle
- MICE m = 24, seed 20260227; primary N = 7,633

## Outputs

- `figures/fig1_hex_facet_CKD.jpg`
- `figures/fig2_forest_marginal_effects.jpg`
- `figures/fig3_mice_convergence.jpg`
- `figures/participant_flow.jpg`
- `output/tables/table1_characteristics.csv`
- `output/tables/table2_mixture.csv`
- `output/tables/table3_ckd_interaction.csv`
- `output/tables/table4_components.csv`
- `output/tables/table4_sensitivities.csv`

## Data

NHANES is public and is downloaded at runtime with [nhanesA](https://cran.r-project.org/package=nhanesA). Do not commit `data/` extracts.

## Ethics

NHANES is approved by the NCHS Ethics Review Board. IRB documentation: https://www.cdc.gov/nchs/nhanes/about/erb.html
