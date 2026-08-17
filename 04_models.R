# 04: qgcomp, Model B, sensitivities, single-pollutant, and component models.
# Model B: PhenoAgeAccel plus chronological age, so the mixture coefficient is age-adjusted.
covars_int <- c(
  if (valid_num("age", qdat)) "age",
  if (valid_cat("sex", qdat)) "factor(sex)",
  if (valid_cat("race_eth", qdat)) "factor(race_eth)",
  if (valid_num("pir", qdat)) "pir",
  if (valid_cat("smoker3", qdat)) "factor(smoker3)",
  if (valid_num("bmi", qdat)) "bmi",
  if (valid_cat("alcohol3", qdat)) "factor(alcohol3)",
  if (valid_cat("cycle", qdat)) "factor(cycle)"
)
intended_covars_int <- covars_int
build_int_formula <- function(covars_vec) {
  rhs <- if (length(covars_vec) > 0) paste("+", paste(covars_vec, collapse = " + ")) else ""
  as.formula(paste("phenoage_accel ~ mixture_qscore * ckd", rhs))
}
f_int <- build_int_formula(covars_int)

repeat {
  mf_int <- tryCatch(model.frame(f_int, data = qdat, na.action = na.omit), error = function(e) NULL)
  if (is.null(mf_int) || nrow(mf_int) == 0) {
    stop("Model B preprocessing failed: no complete-case rows available for interaction formula.")
  }
  factor_terms <- names(mf_int)[vapply(mf_int, is.factor, logical(1))]
  bad_factor_terms <- factor_terms[vapply(mf_int[factor_terms], function(x) nlevels(droplevels(x)) < 2, logical(1))]
  removable_terms <- intersect(covars_int, bad_factor_terms)
  if (length(removable_terms) == 0) break
  message("Model B: dropping collapsed factor term(s) after complete-case filtering: ",
          paste(removable_terms, collapse = ", "))
  covars_int <- setdiff(covars_int, removable_terms)
  f_int <- build_int_formula(covars_int)
}
n_modelb <- nrow(mf_int)
realized_covars_int <- covars_int
covar_report <- tibble(
  intended_term = intended_covars_int,
  realized_in_model_b = ifelse(intended_covars_int %in% realized_covars_int, "Yes", "No")
)
kable(covar_report, caption = "Model B covariate transparency: intended vs realized adjustment terms.")

# ---- Model B fit: complete-case (always) AND MI-pooled (when use_mice = TRUE) ----
fit_int_cc_error <- NULL
fit_int_cc_raw <- tryCatch(
  svyglm(f_int, design = svy_q, family = gaussian()),
  error = function(e) { fit_int_cc_error <<- conditionMessage(e); NULL }
)
if (is.null(fit_int_cc_raw)) stop("Model B complete-case (svyglm interaction) failed: ", fit_int_cc_error)
fit_int_cc <- list(fit = fit_int_cc_raw, mi = FALSE)
message("Model B complete-case fit OK. N = ", n_modelb)
print(coef(fit_int_cc_raw))

fit_int_mi_error <- NULL
fit_int_mi <- if (mice_active) {
  tryCatch(
    fit_svyglm_pool(f_int, svy_q_mi),
    error = function(e) { fit_int_mi_error <<- conditionMessage(e); NULL }
  )
} else NULL
if (mice_active && is.null(fit_int_mi)) {
  message("Model B MI pooling failed: ", fit_int_mi_error, ". Falling back to complete-case.")
}

# Primary `fit_int` exposed downstream: MI-pooled if available, else complete-case.
fit_int <- if (!is.null(fit_int_mi)) fit_int_mi else fit_int_cc
n_modelb_primary <- if (!is.null(fit_int_mi)) n_mi_base else n_modelb
modelb_primary_label <- if (!is.null(fit_int_mi)) "MI-pooled (Rubin)" else "Complete-case"
message("Model B primary inference: ", modelb_primary_label, "; N = ", n_modelb_primary,
        " (complete-case N = ", n_modelb, ")")
# Sanity check: weighted lm and svyglm (weights-only design) should match (uses complete-case data).
fit_int_lm <- lm(f_int, data = qdat_svy, weights = weight_combined)
message("Weighted lm mixture_qscore coef: ", round(coef(fit_int_lm)["mixture_qscore"], 6),
        "; svyglm CC: ", round(coef(fit_int_cc_raw)["mixture_qscore"], 6))

aqg0 <- analytic_qg %>% filter(ckd == 0)
aqg1 <- analytic_qg %>% filter(ckd == 1)
if (nrow(aqg0) > 0) aqg0$weight_normalized <- aqg0$weight_combined / mean(aqg0$weight_combined, na.rm = TRUE)
if (nrow(aqg1) > 0) aqg1$weight_normalized <- aqg1$weight_combined / mean(aqg1$weight_combined, na.rm = TRUE)
fit_qg_ckd0 <- tryCatch(
  qgcomp::qgcomp.noboot(
    f_qg, expnms = qg_exposure_vars,
    data = aqg0,
    family = gaussian(), q = q_qg,
    weights = aqg0$weight_normalized
  ),
  error = function(e) NULL
)
fit_qg_ckd1 <- tryCatch(
  qgcomp::qgcomp.noboot(
    f_qg, expnms = qg_exposure_vars,
    data = aqg1,
    family = gaussian(), q = q_qg,
    weights = aqg1$weight_normalized
  ),
  error = function(e) NULL
)


## ----tables-------------------------------------------------------------------
fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x != 0 & abs(x) < 10^-digits,
      format(x, scientific = TRUE, digits = 2),
      format(round(x, digits), nsmall = digits)
    )
  )
}
fmt_p <- function(x) ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", format(round(x, 3), nsmall = 3)))

# Table 1: characteristics by CKD -- continuous: mean (SD); categorical: n (weighted %)
t1_mean_sd <- function(formula, design) {
  est <- svymean(formula, design, na.rm = TRUE)
  v_obj <- svyvar(formula, design, na.rm = TRUE)
  sd_val <- sqrt(as.numeric(coef(v_obj)[1]))
  paste0(format(round(coef(est)[1], 2), nsmall = 2), " (", format(round(sd_val, 2), nsmall = 2), ")")
}
t1_mean_obj <- function(formula, design) {
  m <- svymean(formula, design, na.rm = TRUE)
  v <- svyvar(formula, design, na.rm = TRUE)
  list(est = unname(coef(m)[1]), se = unname(SE(m)[1]), var = unname(coef(v)[1]))
}
t1_pct_se <- function(formula, design) {
  est <- svymean(formula, design, na.rm = TRUE)
  cf <- coef(est)
  se <- SE(est)
  idx <- 1
  if (length(cf) > 1) {
    nms <- names(cf)
    idx_true <- which(str_detect(nms, "TRUE$"))
    if (length(idx_true) == 0) idx_true <- which(str_detect(nms, "1$"))
    if (length(idx_true) > 0) idx <- idx_true[1]
  }
  paste0(
    format(round(unname(cf[idx]) * 100, 1), nsmall = 1),
    " (",
    format(round(unname(se[idx]) * 100, 1), nsmall = 1),
    ")"
  )
}
t1_n_wpct <- function(cond_vec, formula, design) {
  n_unwt <- sum(cond_vec, na.rm = TRUE)
  est <- svymean(formula, design, na.rm = TRUE)
  cf <- coef(est)
  idx <- 1
  if (length(cf) > 1) {
    nms <- names(cf)
    idx_true <- which(str_detect(nms, "TRUE$"))
    if (length(idx_true) == 0) idx_true <- which(str_detect(nms, "1$"))
    if (length(idx_true) > 0) idx <- idx_true[1]
  }
  paste0(
    format(n_unwt, big.mark = ","),
    " (",
    format(round(unname(cf[idx]) * 100, 1), nsmall = 1),
    ")"
  )
}
t1_pct_obj <- function(formula, design) {
  est <- svymean(formula, design, na.rm = TRUE)
  cf <- coef(est)
  se <- SE(est)
  idx <- 1
  if (length(cf) > 1) {
    nms <- names(cf)
    idx_true <- which(str_detect(nms, "TRUE$"))
    if (length(idx_true) == 0) idx_true <- which(str_detect(nms, "1$"))
    if (length(idx_true) > 0) idx <- idx_true[1]
  }
  list(est = unname(cf[idx]), se = unname(se[idx]))
}
t1_factor_level_pct_se <- function(var, level, design) {
  est <- svymean(as.formula(paste0("~factor(", var, ")")), design, na.rm = TRUE)
  key <- paste0("factor(", var, ")", level)
  if (!(key %in% names(coef(est)))) return("NA")
  paste0(
    format(round(unname(coef(est)[key]) * 100, 1), nsmall = 1),
    " (",
    format(round(unname(SE(est)[key]) * 100, 1), nsmall = 1),
    ")"
  )
}
t1_factor_level_n_wpct <- function(var, level, design, data) {
  n_unwt <- sum(data[[var]] == level, na.rm = TRUE)
  est <- svymean(as.formula(paste0("~factor(", var, ")")), design, na.rm = TRUE)
  key <- paste0("factor(", var, ")", level)
  if (!(key %in% names(coef(est)))) {
    return(paste0(format(n_unwt, big.mark = ","), " (NA)"))
  }
  paste0(
    format(n_unwt, big.mark = ","),
    " (",
    format(round(unname(coef(est)[key]) * 100, 1), nsmall = 1),
    ")"
  )
}
t1_factor_level_obj <- function(var, level, design) {
  est <- svymean(as.formula(paste0("~factor(", var, ")")), design, na.rm = TRUE)
  key <- paste0("factor(", var, ")", level)
  if (!(key %in% names(coef(est)))) return(list(est = NA_real_, se = NA_real_))
  list(est = unname(coef(est)[key]), se = unname(SE(est)[key]))
}
# Missing category: weighted % of rows with NA on `var` (for MICE / missingness diagnostics).
t1_miss_pct_obj <- function(var, design) {
  t1_pct_obj(as.formula(paste0("~I(is.na(", var, "))")), design)
}
t1_miss_n_wpct <- function(var, design, data_sub) {
  n_unwt <- sum(is.na(data_sub[[var]]))
  o <- t1_miss_pct_obj(var, design)
  paste0(
    format(n_unwt, big.mark = ","),
    " (",
    format(round(o$est * 100, 1), nsmall = 1),
    ")"
  )
}
smd_cont <- function(o0, o1) {
  if (any(!is.finite(c(o0$est, o1$est, o0$var, o1$var)))) return(NA_real_)
  denom <- sqrt((o0$var + o1$var) / 2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  abs(o1$est - o0$est) / denom
}
smd_prop <- function(o0, o1) {
  if (any(!is.finite(c(o0$est, o1$est)))) return(NA_real_)
  denom <- sqrt((o0$est * (1 - o0$est) + o1$est * (1 - o1$est)) / 2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  abs(o1$est - o0$est) / denom
}
fmt_smd <- function(x) ifelse(is.na(x), "", format(round(x, 3), nsmall = 3))

# Build Part A (Model B PRIMARY estimation sample = MICE base; PFAS+outcome+design complete) and
# Part B (analytic-sample participants excluded from primary inference: missing PFAS / outcome / design).
# Within Part A we display Missing rows so the reader sees exactly which covariate cells are imputed by MICE.
if ("SEQN" %in% names(analytic) && exists("qdat_svy") && "SEQN" %in% names(qdat_svy)) {
  inc_idx <- analytic$SEQN %in% qdat_svy$SEQN
} else {
  inc_idx <- !is.na(analytic$ckd) & !is.na(analytic$phenoage_accel)
}
analytic_inc <- analytic[inc_idx, , drop = FALSE]
analytic_exc <- analytic[!inc_idx, , drop = FALSE]

svy_inc <- svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined,
                     nest = TRUE, data = analytic_inc)
svy_exc <- svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined,
                     nest = TRUE, data = analytic_exc)

# Block builder: returns a tibble with columns characteristic, No CKD, CKD, SMD.
# Continuous: mean (SD).  Categorical: n (weighted %).  No SE displayed anywhere.
build_t1_block <- function(design, data) {
  d0 <- subset(design, ckd == 0)
  d1 <- subset(design, ckd == 1)
  data0 <- data[!is.na(data$ckd) & data$ckd == 0, , drop = FALSE]
  data1 <- data[!is.na(data$ckd) & data$ckd == 1, , drop = FALSE]

  t1_show_miss <- function(var) var %in% names(data) && any(is.na(data[[var]]))

  age0 <- t1_mean_obj(~age, d0); age1 <- t1_mean_obj(~age, d1)
  age_miss0 <- t1_miss_pct_obj("age", d0); age_miss1 <- t1_miss_pct_obj("age", d1)
  pir0 <- t1_mean_obj(~pir, d0); pir1 <- t1_mean_obj(~pir, d1)
  pir_miss0 <- t1_miss_pct_obj("pir", d0); pir_miss1 <- t1_miss_pct_obj("pir", d1)
  bmi0 <- t1_mean_obj(~bmi, d0); bmi1 <- t1_mean_obj(~bmi, d1)
  bmi_miss0 <- t1_miss_pct_obj("bmi", d0); bmi_miss1 <- t1_miss_pct_obj("bmi", d1)

  fem0 <- t1_pct_obj(~I(sex == "Female"), d0); fem1 <- t1_pct_obj(~I(sex == "Female"), d1)
  sex_miss0 <- t1_miss_pct_obj("sex", d0); sex_miss1 <- t1_miss_pct_obj("sex", d1)
  his0 <- t1_factor_level_obj("race_eth", "Hispanic", d0); his1 <- t1_factor_level_obj("race_eth", "Hispanic", d1)
  whi0 <- t1_factor_level_obj("race_eth", "Non-Hispanic White", d0); whi1 <- t1_factor_level_obj("race_eth", "Non-Hispanic White", d1)
  blk0 <- t1_factor_level_obj("race_eth", "Non-Hispanic Black", d0); blk1 <- t1_factor_level_obj("race_eth", "Non-Hispanic Black", d1)
  asn0 <- t1_factor_level_obj("race_eth", "Non-Hispanic Asian", d0); asn1 <- t1_factor_level_obj("race_eth", "Non-Hispanic Asian", d1)
  oth0 <- t1_factor_level_obj("race_eth", "Other/Multi", d0); oth1 <- t1_factor_level_obj("race_eth", "Other/Multi", d1)
  race_miss0 <- t1_miss_pct_obj("race_eth", d0); race_miss1 <- t1_miss_pct_obj("race_eth", d1)
  edu_lths0 <- t1_factor_level_obj("education", "<HS", d0); edu_lths1 <- t1_factor_level_obj("education", "<HS", d1)
  edu_hs0   <- t1_factor_level_obj("education", "HS/GED", d0); edu_hs1   <- t1_factor_level_obj("education", "HS/GED", d1)
  edu_gt0   <- t1_factor_level_obj("education", ">HS", d0); edu_gt1   <- t1_factor_level_obj("education", ">HS", d1)
  edu_miss0 <- t1_miss_pct_obj("education", d0); edu_miss1 <- t1_miss_pct_obj("education", d1)
  smk_never0 <- t1_factor_level_obj("smoker3", "Never",   d0); smk_never1 <- t1_factor_level_obj("smoker3", "Never",   d1)
  smk_form0  <- t1_factor_level_obj("smoker3", "Former",  d0); smk_form1  <- t1_factor_level_obj("smoker3", "Former",  d1)
  smk_curr0  <- t1_factor_level_obj("smoker3", "Current", d0); smk_curr1  <- t1_factor_level_obj("smoker3", "Current", d1)
  smk_miss0 <- t1_miss_pct_obj("smoker3", d0); smk_miss1 <- t1_miss_pct_obj("smoker3", d1)
  alc_never0 <- t1_factor_level_obj("alcohol3", "Never",   d0); alc_never1 <- t1_factor_level_obj("alcohol3", "Never",   d1)
  alc_form0  <- t1_factor_level_obj("alcohol3", "Former",  d0); alc_form1  <- t1_factor_level_obj("alcohol3", "Former",  d1)
  alc_curr0  <- t1_factor_level_obj("alcohol3", "Current", d0); alc_curr1  <- t1_factor_level_obj("alcohol3", "Current", d1)
  alc_miss0 <- t1_miss_pct_obj("alcohol3", d0); alc_miss1 <- t1_miss_pct_obj("alcohol3", d1)
  pa0  <- t1_pct_obj(~I(phys_active == 1), d0); pa1  <- t1_pct_obj(~I(phys_active == 1), d1)
  pa_miss0 <- t1_miss_pct_obj("phys_active", d0); pa_miss1 <- t1_miss_pct_obj("phys_active", d1)
  dm0  <- t1_pct_obj(~I(diabetes == 1), d0); dm1  <- t1_pct_obj(~I(diabetes == 1), d1)
  dm_miss0 <- t1_miss_pct_obj("diabetes", d0); dm_miss1 <- t1_miss_pct_obj("diabetes", d1)
  htn0 <- t1_pct_obj(~I(hypertension == 1), d0); htn1 <- t1_pct_obj(~I(hypertension == 1), d1)
  htn_miss0 <- t1_miss_pct_obj("hypertension", d0); htn_miss1 <- t1_miss_pct_obj("hypertension", d1)
  cvd0 <- t1_pct_obj(~I(cvd == 1), d0); cvd1 <- t1_pct_obj(~I(cvd == 1), d1)
  cvd_miss0 <- t1_miss_pct_obj("cvd", d0); cvd_miss1 <- t1_miss_pct_obj("cvd", d1)

  show_age_miss <- t1_show_miss("age")
  show_pir_miss <- t1_show_miss("pir")
  show_bmi_miss <- t1_show_miss("bmi")
  show_sex_miss <- t1_show_miss("sex")
  show_race_miss <- t1_show_miss("race_eth")
  show_edu_miss <- t1_show_miss("education")
  show_smk_miss <- t1_show_miss("smoker3")
  show_alc_miss <- t1_show_miss("alcohol3")
  show_pa_miss <- t1_show_miss("phys_active")
  show_dm_miss <- t1_show_miss("diabetes")
  show_htn_miss <- t1_show_miss("hypertension")
  show_cvd_miss <- t1_show_miss("cvd")

  tibble(
    characteristic = c(
      "N (unweighted)", "Age, mean (SD)",
      if (show_age_miss) "Age missing, n (weighted %)" else character(0),
      "Female, n (weighted %)",
      if (show_sex_miss) "Sex missing, n (weighted %)" else character(0),
      "Hispanic, n (weighted %)", "Non-Hispanic White, n (weighted %)", "Non-Hispanic Black, n (weighted %)",
      "Non-Hispanic Asian, n (weighted %)", "Other/Multi, n (weighted %)",
      if (show_race_miss) "Race/ethnicity missing, n (weighted %)" else character(0),
      "Education <HS, n (weighted %)", "Education HS/GED, n (weighted %)", "Education >HS, n (weighted %)",
      if (show_edu_miss) "Education missing, n (weighted %)" else character(0),
      "PIR, mean (SD)",
      if (show_pir_miss) "PIR missing, n (weighted %)" else character(0),
      "BMI, mean (SD)",
      if (show_bmi_miss) "BMI missing, n (weighted %)" else character(0),
      "Never smoker, n (weighted %)", "Former smoker, n (weighted %)", "Current smoker, n (weighted %)",
      if (show_smk_miss) "Smoking status missing, n (weighted %)" else character(0),
      "Never drinker, n (weighted %)", "Former drinker, n (weighted %)", "Current drinker, n (weighted %)",
      if (show_alc_miss) "Alcohol status missing, n (weighted %)" else character(0),
      "Physically active, n (weighted %)",
      if (show_pa_miss) "Physical activity missing, n (weighted %)" else character(0),
      "Diabetes, n (weighted %)",
      if (show_dm_miss) "Diabetes missing, n (weighted %)" else character(0),
      "Hypertension, n (weighted %)",
      if (show_htn_miss) "Hypertension missing, n (weighted %)" else character(0),
      "CVD history, n (weighted %)",
      if (show_cvd_miss) "CVD history missing, n (weighted %)" else character(0)
    ),
    `No CKD` = c(
      format(nrow(data0), big.mark = ","),
      t1_mean_sd(~age, d0),
      if (show_age_miss) t1_miss_n_wpct("age", d0, data0) else character(0),
      t1_n_wpct(data0$sex == "Female", ~I(sex == "Female"), d0),
      if (show_sex_miss) t1_miss_n_wpct("sex", d0, data0) else character(0),
      t1_factor_level_n_wpct("race_eth", "Hispanic", d0, data0),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic White", d0, data0),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic Black", d0, data0),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic Asian", d0, data0),
      t1_factor_level_n_wpct("race_eth", "Other/Multi", d0, data0),
      if (show_race_miss) t1_miss_n_wpct("race_eth", d0, data0) else character(0),
      t1_factor_level_n_wpct("education", "<HS", d0, data0),
      t1_factor_level_n_wpct("education", "HS/GED", d0, data0),
      t1_factor_level_n_wpct("education", ">HS", d0, data0),
      if (show_edu_miss) t1_miss_n_wpct("education", d0, data0) else character(0),
      t1_mean_sd(~pir, d0),
      if (show_pir_miss) t1_miss_n_wpct("pir", d0, data0) else character(0),
      t1_mean_sd(~bmi, d0),
      if (show_bmi_miss) t1_miss_n_wpct("bmi", d0, data0) else character(0),
      t1_factor_level_n_wpct("smoker3", "Never",   d0, data0),
      t1_factor_level_n_wpct("smoker3", "Former",  d0, data0),
      t1_factor_level_n_wpct("smoker3", "Current", d0, data0),
      if (show_smk_miss) t1_miss_n_wpct("smoker3", d0, data0) else character(0),
      t1_factor_level_n_wpct("alcohol3", "Never",   d0, data0),
      t1_factor_level_n_wpct("alcohol3", "Former",  d0, data0),
      t1_factor_level_n_wpct("alcohol3", "Current", d0, data0),
      if (show_alc_miss) t1_miss_n_wpct("alcohol3", d0, data0) else character(0),
      t1_n_wpct(data0$phys_active == 1, ~I(phys_active == 1), d0),
      if (show_pa_miss) t1_miss_n_wpct("phys_active", d0, data0) else character(0),
      t1_n_wpct(data0$diabetes == 1, ~I(diabetes == 1), d0),
      if (show_dm_miss) t1_miss_n_wpct("diabetes", d0, data0) else character(0),
      t1_n_wpct(data0$hypertension == 1, ~I(hypertension == 1), d0),
      if (show_htn_miss) t1_miss_n_wpct("hypertension", d0, data0) else character(0),
      t1_n_wpct(data0$cvd == 1, ~I(cvd == 1), d0),
      if (show_cvd_miss) t1_miss_n_wpct("cvd", d0, data0) else character(0)
    ),
    CKD = c(
      format(nrow(data1), big.mark = ","),
      t1_mean_sd(~age, d1),
      if (show_age_miss) t1_miss_n_wpct("age", d1, data1) else character(0),
      t1_n_wpct(data1$sex == "Female", ~I(sex == "Female"), d1),
      if (show_sex_miss) t1_miss_n_wpct("sex", d1, data1) else character(0),
      t1_factor_level_n_wpct("race_eth", "Hispanic", d1, data1),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic White", d1, data1),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic Black", d1, data1),
      t1_factor_level_n_wpct("race_eth", "Non-Hispanic Asian", d1, data1),
      t1_factor_level_n_wpct("race_eth", "Other/Multi", d1, data1),
      if (show_race_miss) t1_miss_n_wpct("race_eth", d1, data1) else character(0),
      t1_factor_level_n_wpct("education", "<HS", d1, data1),
      t1_factor_level_n_wpct("education", "HS/GED", d1, data1),
      t1_factor_level_n_wpct("education", ">HS", d1, data1),
      if (show_edu_miss) t1_miss_n_wpct("education", d1, data1) else character(0),
      t1_mean_sd(~pir, d1),
      if (show_pir_miss) t1_miss_n_wpct("pir", d1, data1) else character(0),
      t1_mean_sd(~bmi, d1),
      if (show_bmi_miss) t1_miss_n_wpct("bmi", d1, data1) else character(0),
      t1_factor_level_n_wpct("smoker3", "Never",   d1, data1),
      t1_factor_level_n_wpct("smoker3", "Former",  d1, data1),
      t1_factor_level_n_wpct("smoker3", "Current", d1, data1),
      if (show_smk_miss) t1_miss_n_wpct("smoker3", d1, data1) else character(0),
      t1_factor_level_n_wpct("alcohol3", "Never",   d1, data1),
      t1_factor_level_n_wpct("alcohol3", "Former",  d1, data1),
      t1_factor_level_n_wpct("alcohol3", "Current", d1, data1),
      if (show_alc_miss) t1_miss_n_wpct("alcohol3", d1, data1) else character(0),
      t1_n_wpct(data1$phys_active == 1, ~I(phys_active == 1), d1),
      if (show_pa_miss) t1_miss_n_wpct("phys_active", d1, data1) else character(0),
      t1_n_wpct(data1$diabetes == 1, ~I(diabetes == 1), d1),
      if (show_dm_miss) t1_miss_n_wpct("diabetes", d1, data1) else character(0),
      t1_n_wpct(data1$hypertension == 1, ~I(hypertension == 1), d1),
      if (show_htn_miss) t1_miss_n_wpct("hypertension", d1, data1) else character(0),
      t1_n_wpct(data1$cvd == 1, ~I(cvd == 1), d1),
      if (show_cvd_miss) t1_miss_n_wpct("cvd", d1, data1) else character(0)
    ),
    SMD = c(
      "",
      fmt_smd(smd_cont(age0, age1)),
      if (show_age_miss) fmt_smd(smd_prop(age_miss0, age_miss1)) else character(0),
      fmt_smd(smd_prop(fem0, fem1)),
      if (show_sex_miss) fmt_smd(smd_prop(sex_miss0, sex_miss1)) else character(0),
      fmt_smd(smd_prop(his0, his1)),
      fmt_smd(smd_prop(whi0, whi1)),
      fmt_smd(smd_prop(blk0, blk1)),
      fmt_smd(smd_prop(asn0, asn1)),
      fmt_smd(smd_prop(oth0, oth1)),
      if (show_race_miss) fmt_smd(smd_prop(race_miss0, race_miss1)) else character(0),
      fmt_smd(smd_prop(edu_lths0, edu_lths1)),
      fmt_smd(smd_prop(edu_hs0, edu_hs1)),
      fmt_smd(smd_prop(edu_gt0, edu_gt1)),
      if (show_edu_miss) fmt_smd(smd_prop(edu_miss0, edu_miss1)) else character(0),
      fmt_smd(smd_cont(pir0, pir1)),
      if (show_pir_miss) fmt_smd(smd_prop(pir_miss0, pir_miss1)) else character(0),
      fmt_smd(smd_cont(bmi0, bmi1)),
      if (show_bmi_miss) fmt_smd(smd_prop(bmi_miss0, bmi_miss1)) else character(0),
      fmt_smd(smd_prop(smk_never0, smk_never1)),
      fmt_smd(smd_prop(smk_form0,  smk_form1)),
      fmt_smd(smd_prop(smk_curr0,  smk_curr1)),
      if (show_smk_miss) fmt_smd(smd_prop(smk_miss0, smk_miss1)) else character(0),
      fmt_smd(smd_prop(alc_never0, alc_never1)),
      fmt_smd(smd_prop(alc_form0,  alc_form1)),
      fmt_smd(smd_prop(alc_curr0,  alc_curr1)),
      if (show_alc_miss) fmt_smd(smd_prop(alc_miss0, alc_miss1)) else character(0),
      fmt_smd(smd_prop(pa0, pa1)),
      if (show_pa_miss) fmt_smd(smd_prop(pa_miss0, pa_miss1)) else character(0),
      fmt_smd(smd_prop(dm0, dm1)),
      if (show_dm_miss) fmt_smd(smd_prop(dm_miss0, dm_miss1)) else character(0),
      fmt_smd(smd_prop(htn0, htn1)),
      if (show_htn_miss) fmt_smd(smd_prop(htn_miss0, htn_miss1)) else character(0),
      fmt_smd(smd_prop(cvd0, cvd1)),
      if (show_cvd_miss) fmt_smd(smd_prop(cvd_miss0, cvd_miss1)) else character(0)
    )
  )
}

block_a <- build_t1_block(svy_inc, analytic_inc)
block_b <- build_t1_block(svy_exc, analytic_exc)

n_inc <- nrow(analytic_inc)
n_exc <- nrow(analytic_exc)
n_inc_0 <- sum(analytic_inc$ckd == 0, na.rm = TRUE)
n_inc_1 <- sum(analytic_inc$ckd == 1, na.rm = TRUE)
n_exc_0 <- sum(analytic_exc$ckd == 0, na.rm = TRUE)
n_exc_1 <- sum(analytic_exc$ckd == 1, na.rm = TRUE)

header_a <- tibble(
  characteristic = "PART A -- INCLUDED IN STUDY",
  `No CKD` = "", CKD = "", SMD = ""
)
header_b <- tibble(
  characteristic = "PART B -- EXCLUDED FROM STUDY",
  `No CKD` = "", CKD = "", SMD = ""
)

table1 <- bind_rows(header_a, block_a, header_b, block_b)

kable(
  table1,
  caption = "Table 1. Sample characteristics by CKD status. Part A = MICE base sample (Model B PRIMARY estimation sample; PFAS + outcome + design complete). Missing rows show pre-imputation missingness on each variable shown (continuous and categorical alike); those cells are then imputed under MICE for primary inference. Continuous: mean (SD); categorical: n (weighted %); Missing rows: n (weighted %) of NA. SMD is the absolute standardized mean difference between CKD vs No CKD within each part."
)
cat("\n*Table 1 note:* Weighted percentages use the NHANES complex survey design (strata, PSU, combined weights). ",
    "Part A = MICE base sample = participants entering Model B primary inference (non-missing PhenoAge, CKD, weight/PSU/strata, and complete primary 4-PFAS set). ",
    "Within Part A, the Missing rows quantify pre-imputation covariate missingness; those cells are imputed via MICE (m = ", config$mice_m,
    ") and pooled with Rubin's rules. ",
    "Part B = analytic-sample participants excluded from primary inference because they are missing one or more PFAS / outcome / design components (these are NOT imputed). ",
    "SMD in Part B is computed between excluded No CKD vs excluded CKD; values >0.10 are typically considered meaningful imbalance.\n", sep = "")

# Table 2: overall mixture effect (MI-pooled qgcomp via Rubin when MICE is active; CC fallback otherwise).
qg_mi_overall <- if (mice_active) {
  qgcomp_mi_pool(
    formula = f_qg, expnms = qg_exposure_vars,
    imp_list_obj = imp_list, q = q_qg,
    weight_col = "weight_combined"
  )
} else NULL

tab2_mi <- if (!is.null(qg_mi_overall)) {
  tibble(
    model = paste0("Model A: Overall mixture effect [MI-pooled, m = ", qg_mi_overall$m, "]"),
    beta = fmt(qg_mi_overall$psi),
    conf.low = fmt(qg_mi_overall$conf.low),
    conf.high = fmt(qg_mi_overall$conf.high),
    p.value = fmt_p(qg_mi_overall$p.value),
    N = round(mean(qg_mi_overall$n_per_fit)),
    source = "MI-pooled"
  )
} else NULL

tab2_cc <- safe_qgcomp_extract(
  fit_qg_overall,
  if (!is.null(qg_mi_overall)) "Model A: Overall mixture effect [Complete-case sensitivity]"
  else "Model A: Overall mixture effect [Complete-case]",
  nrow(analytic_qg)
) %>%
  mutate(
    beta = fmt(beta),
    conf.low = fmt(conf.low),
    conf.high = fmt(conf.high),
    p.value = fmt_p(p.value),
    source = "Complete-case"
  )

tab2_full <- if (!is.null(tab2_mi)) bind_rows(tab2_mi, tab2_cc) else tab2_cc

kable(
  tab2_full,
  caption = paste0(
    "Table 2. Overall association between PFAS mixture score and PhenoAge acceleration (years). ",
    "Mixture score ranges from 1 to 4; a 1-unit increase represents a 1-quantile increase averaged ",
    "across the fixed primary PFAS set. Estimates from qgcomp with normalized pooled NHANES sampling weights",
    if (!is.null(qg_mi_overall)) "; primary row is MI-pooled across imputations via Rubin's rules; CC row shown for transparency." else " (complete-case)."
  )
)
cat("\n*Table 2 note:* qgcomp model adjusted for the same covariates as Model B: age, sex, race/ethnicity, PIR, BMI, smoking (Never/Former/Current), alcohol (Never/Former/Current), cycle. ",
    "qgcomp intended covariates: ", paste(qg_covars_intended, collapse = ", "), ". ",
    "qgcomp realized covariates (after complete-case support checks): ", paste(qg_covars_realized, collapse = ", "), ". ",
    if (!is.null(qg_mi_overall)) {
      paste0("Primary row uses ", qg_mi_overall$m, " MICE imputations of model covariates; psi (overall mixture coef) and its variance are pooled via Rubin's rules. ",
             "qgcomp variance does NOT incorporate PSU/strata; statistical inference for the interaction question relies on the survey-weighted Table 3 (svyglm). ")
    } else if (!is.null(fit_qg_overall)) {
      "Inference here uses qgcomp's analytic variance on the complete-case sample; see Table 3 (svyglm) for primary survey-weighted inference. "
    } else {
      paste0("qgcomp failed: ",
             ifelse(is.null(qgcomp_overall_error), "unknown error", qgcomp_overall_error), ". ")
    },
    "PhenoAge acceleration defined as survey-weighted residual of PhenoAge on chronological age.\n", sep = "")

# Table 3: interaction results + design-based marginal effects by CKD group
if (!is.null(fit_int)) {
  int_tidy <- tidy_pool(fit_int, conf.level = 0.95)
  cf_names <- names(coef_pool(fit_int))

  int_terms_wanted <- c("mixture_qscore", cf_names[str_detect(cf_names, "^ckd$|mixture_qscore:ckd")])
  tab3_raw <- int_tidy %>%
    filter(term %in% int_terms_wanted) %>%
    transmute(
      term,
      beta = fmt(estimate),
      conf.low = fmt(conf.low),
      conf.high = fmt(conf.high),
      p.value = fmt_p(p.value),
      N = n_modelb_primary,
      model = paste0("Model B: Mixture*CKD interaction [", modelb_primary_label, "]")
    )

  # Design-based marginal effects via svycontrast with named contrast vectors (pooled across imputations if MI).

  # Robustly identify the interaction term regardless of factor coding.
  int_name <- cf_names[str_detect(cf_names, "mixture_qscore:ckd")]
  if (length(int_name) != 1) {
    warning("Could not uniquely identify mixture*ckd interaction term among: ",
            paste(cf_names, collapse = ", "), ". Skipping marginal effects.")
    tab3 <- tab3_raw
  } else {
    c_ckd0 <- setNames(rep(0, length(cf_names)), cf_names); c_ckd0["mixture_qscore"] <- 1
    c_ckd1 <- setNames(rep(0, length(cf_names)), cf_names); c_ckd1["mixture_qscore"] <- 1; c_ckd1[int_name] <- 1

    ctr0 <- contrast_pool(fit_int, c_ckd0)
    ctr1 <- contrast_pool(fit_int, c_ckd1)

    make_contrast_row <- function(ctr, label) {
      tibble(term = label, beta = fmt(ctr$estimate),
             conf.low = fmt(ctr$conf.low),
             conf.high = fmt(ctr$conf.high),
             p.value = fmt_p(ctr$p.value),
             N = n_modelb_primary, model = "Marginal effects (design-based)")
    }

    p_int_tidy <- int_tidy %>% filter(term == int_name) %>% pull(p.value)
    p_interaction <- if (length(p_int_tidy) > 0) p_int_tidy[1] else NA_real_

    tab3_marginal <- bind_rows(
      make_contrast_row(ctr0, "Mixture effect among CKD=0"),
      make_contrast_row(ctr1, "Mixture effect among CKD=1"),
      tibble(term = "P-interaction", beta = "", conf.low = "", conf.high = "",
             p.value = fmt_p(p_interaction), N = n_modelb_primary, model = "Marginal effects (design-based)")
    )
    tab3 <- bind_rows(tab3_raw, tab3_marginal)
  }
} else {
  tab3 <- tibble(
    term = c("mixture_qscore", "ckd", "mixture_qscore:ckd"),
    beta = "NA", conf.low = "NA", conf.high = "NA", p.value = "NA",
    N = n_modelb_primary, model = "Model B: Mixture*CKD interaction (not estimable)"
  )
}
kable(tab3, caption = "Table 3. Interaction model with design-based marginal effects by CKD group")
cat("\n*Table 3 note:* Model B is a survey-weighted linear model (`svyglm`) with interaction term `mixture_qscore * ckd`. ",
    "Inference: ", modelb_primary_label, " (N = ", format(n_modelb_primary, big.mark = ","),
    "); complete-case sensitivity N = ", format(n_modelb, big.mark = ","),
    if (mice_active) paste0("; MICE m = ", config$mice_m, " imputations of missing covariates (outcome/exposures/design not imputed). ") else ". ",
    "Intended covariates: ", paste(intended_covars_int, collapse = ", "), ". ",
    "Realized covariates (after complete-case support checks): ", paste(realized_covars_int, collapse = ", "), ". ",
    "Marginal effects for CKD=0 and CKD=1 are estimated using design-based contrasts (`svycontrast`); ",
    "for MI runs, contrasts are computed per imputation and pooled with Rubin's rules (mitools::MIcombine). ",
    "95% CI and p-values use t-distribution on Rubin (or survey residual) degrees of freedom.\n", sep = "")

# Table 4: sensitivity analyses -- using same covariate set as primary models.
sens_covars_rhs <- paste(realized_covars_int, collapse = " + ")
make_rhs_formula <- function(outcome, main_terms, covar_rhs) {
  rhs <- paste(main_terms, collapse = " + ")
  if (!is.null(covar_rhs) && nzchar(covar_rhs)) rhs <- paste(rhs, covar_rhs, sep = " + ")
  as.formula(paste(outcome, "~", rhs))
}
prep_svy_frame <- function(data_df) {
  data_svy <- as.data.frame(data_df)
  for (v in c("phenoage_accel", "phenoage_nokid", "phenoage_nokid_both", "phenoage_accel_unweighted",
              "mixture_qscore", "age", "pir", "bmi", "weight_combined", "psu", "strata")) {
    if (v %in% names(data_svy)) data_svy[[v]] <- suppressWarnings(as.numeric(as.character(data_svy[[v]])))
  }
  if ("ckd" %in% names(data_svy)) data_svy$ckd <- as.integer(data_svy$ckd)
  if ("ckd_sens" %in% names(data_svy)) data_svy$ckd_sens <- as.integer(data_svy$ckd_sens)
  for (v in c("sex", "race_eth", "cycle")) {
    if (v %in% names(data_svy)) data_svy[[v]] <- factor(as.character(data_svy[[v]]))
  }
  if ("smoker3"  %in% names(data_svy)) data_svy$smoker3  <- factor(as.character(data_svy$smoker3),  levels = c("Never", "Former", "Current"))
  if ("alcohol3" %in% names(data_svy)) data_svy$alcohol3 <- factor(as.character(data_svy$alcohol3), levels = c("Never", "Former", "Current"))
  data_svy$sex      <- droplevels(factor(data_svy$sex))
  data_svy$race_eth <- droplevels(factor(data_svy$race_eth))
  data_svy$cycle    <- droplevels(factor(data_svy$cycle))
  if ("smoker3"  %in% names(data_svy)) data_svy$smoker3  <- droplevels(data_svy$smoker3)
  if ("alcohol3" %in% names(data_svy)) data_svy$alcohol3 <- droplevels(data_svy$alcohol3)
  data_svy
}

# MI-aware sensitivity fit. `data_df` = complete-case subset for CC fit.
# `mi_filter_fn` (optional, used when mice_active) = function(df) -> filtered imputed df.
# Pooled MI inference is preferred over CC when available.
fit_svy_with_covar_guard <- function(data_df, outcome, main_terms, covars_terms, label,
                                     mi_filter_fn = NULL) {
  build_formula <- function(covs) {
    covar_rhs <- if (length(covs) > 0) paste(covs, collapse = " + ") else ""
    make_rhs_formula(outcome, main_terms, covar_rhs)
  }
  covs <- covars_terms
  intended <- covs
  f <- build_formula(covs)
  mf <- NULL
  repeat {
    mf <- tryCatch(model.frame(f, data = data_df, na.action = na.omit), error = function(e) NULL)
    if (is.null(mf) || nrow(mf) == 0) break
    factor_terms <- names(mf)[vapply(mf, is.factor, logical(1))]
    bad_factor_terms <- factor_terms[vapply(mf[factor_terms], function(x) nlevels(droplevels(x)) < 2, logical(1))]
    removable <- intersect(covs, bad_factor_terms)
    if (length(removable) == 0) break
    message("Sensitivity [", label, "]: dropping collapsed factor term(s): ", paste(removable, collapse = ", "))
    covs <- setdiff(covs, removable)
    f <- build_formula(covs)
  }
  if (is.null(mf) || nrow(mf) == 0) {
    return(list(
      fit = NULL, pooled = NULL, fits_mi = NULL,
      n_model = 0L, n_model_mi = NA_integer_, mi = FALSE,
      intended = intended, realized = covs,
      dropped = paste(setdiff(intended, covs), collapse = ", "),
      error = "No complete-case rows for model frame."
    ))
  }
  data_svy <- prep_svy_frame(data_df)
  des <- tryCatch(
    svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE, data = data_svy),
    error = function(e) svydesign(id = ~1, strata = NULL, weights = ~weight_combined, data = data_svy)
  )
  err <- NULL
  fit <- tryCatch(
    svyglm(f, design = des, family = gaussian()),
    error = function(e) { err <<- conditionMessage(e); NULL }
  )

  pooled <- NULL; fits_mi <- NULL; n_model_mi <- NA_integer_; err_mi <- NULL
  if (mice_active && !is.null(mi_filter_fn)) {
    mi_dfs <- tryCatch(
      lapply(seq_len(config$mice_m), function(i) {
        df <- imp_list$imputations[[i]]
        df <- mi_filter_fn(df)
        prep_svy_frame(df)
      }),
      error = function(e) { err_mi <<- conditionMessage(e); NULL }
    )
    if (!is.null(mi_dfs) && length(mi_dfs) > 0) {
      n_model_mi <- nrow(mi_dfs[[1]])
      mi_des <- tryCatch(
        survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE,
                          data = mitools::imputationList(mi_dfs)),
        error = function(e) survey::svydesign(id = ~1, strata = NULL, weights = ~weight_combined,
                                              data = mitools::imputationList(mi_dfs))
      )
      pool_obj <- tryCatch(
        fit_svyglm_pool(f, mi_des),
        error = function(e) { err_mi <<- conditionMessage(e); NULL }
      )
      if (!is.null(pool_obj) && isTRUE(pool_obj$mi)) {
        pooled <- pool_obj$pooled
        fits_mi <- pool_obj$fits
      }
    }
  }
  list(
    fit = fit, pooled = pooled, fits_mi = fits_mi,
    n_model = nrow(mf), n_model_mi = n_model_mi,
    mi = !is.null(pooled),
    intended = intended, realized = covs,
    dropped = paste(setdiff(intended, covs), collapse = ", "),
    error = if (!is.null(pooled) && !is.null(err_mi)) paste(err, "; MI:", err_mi) else err
  )
}

# Wrap a sensitivity object as a `fit_svyglm_pool`-style object for downstream tidy/contrast helpers.
sens_as_fit <- function(obj) {
  if (isTRUE(obj$mi) && !is.null(obj$pooled) && !is.null(obj$fits_mi)) {
    list(fits = obj$fits_mi, pooled = obj$pooled, mi = TRUE)
  } else if (!is.null(obj$fit)) {
    list(fit = obj$fit, mi = FALSE)
  } else {
    NULL
  }
}

sens_ckd_egfr <- qdat %>% filter(!is.na(ckd_sens))
fit_sens1_obj <- fit_svy_with_covar_guard(
  sens_ckd_egfr, "phenoage_accel", c("mixture_qscore", "ckd_sens"), realized_covars_int,
  "CKD redefined as eGFR<60",
  mi_filter_fn = function(df) df[!is.na(df$ckd_sens), , drop = FALSE]
)
fit_sens1 <- fit_sens1_obj$fit

# PhenoAge-recomputed sensitivities use the precomputed columns on qdat / qdat_svy.
fit_sens2_obj <- fit_svy_with_covar_guard(
  qdat, "phenoage_nokid", c("mixture_qscore", "ckd"), realized_covars_int,
  "PhenoAge excluding creatinine only",
  mi_filter_fn = identity
)
fit_sens2 <- fit_sens2_obj$fit

fit_sens3_obj <- fit_svy_with_covar_guard(
  qdat, "phenoage_nokid_both", c("mixture_qscore", "ckd"), realized_covars_int,
  "PhenoAge excluding creatinine and albumin",
  mi_filter_fn = identity
)
fit_sens3 <- fit_sens3_obj$fit

get_mix_row <- function(obj, label, n_model_cc = NA_integer_, term_name = "mixture_qscore") {
  fit_for_tidy <- sens_as_fit(obj)
  if (is.null(fit_for_tidy)) {
    return(tibble(sensitivity = label, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                  p.value = NA_real_, N = n_model_cc, source = "Failed"))
  }
  tb <- tryCatch(tidy_pool(fit_for_tidy, conf.level = 0.95), error = function(e) NULL)
  if (is.null(tb) || !any(tb$term == term_name)) {
    return(tibble(sensitivity = label, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                  p.value = NA_real_, N = n_model_cc, source = "Failed"))
  }
  row <- tb %>% filter(term == term_name) %>% slice(1)
  N_used <- if (isTRUE(fit_for_tidy$mi)) (obj$n_model_mi %||% n_model_cc) else n_model_cc
  src    <- if (isTRUE(fit_for_tidy$mi)) "MI-pooled" else "Complete-case"
  tibble(sensitivity = label, estimate = row$estimate, conf.low = row$conf.low,
         conf.high = row$conf.high, p.value = row$p.value, N = N_used, source = src)
}
# Null-coalesce helper (R < 4.4 doesn't have %||%)
`%||%` <- function(a, b) if (is.null(a) || (is.atomic(a) && length(a) == 1 && is.na(a))) b else a

f_sens_uw <- make_rhs_formula("phenoage_accel_unweighted", c("mixture_qscore", "ckd"), sens_covars_rhs)
f_sens_both <- make_rhs_formula("phenoage_accel", c("mixture_qscore", "ckd"), sens_covars_rhs)
fit_sens_uw_obj <- fit_svy_with_covar_guard(
  qdat, "phenoage_accel_unweighted", c("mixture_qscore", "ckd"), realized_covars_int,
  "Unweighted residualized PhenoAgeAccel (sensitivity)",
  mi_filter_fn = identity
)
fit_sens_both_markers_obj <- fit_svy_with_covar_guard(
  qdat %>% filter(!is.na(egfr), !is.na(acr)), "phenoage_accel", c("mixture_qscore", "ckd"),
  realized_covars_int, "Restricted to participants with both eGFR and ACR",
  mi_filter_fn = function(df) df[!is.na(df$egfr) & !is.na(df$acr), , drop = FALSE]
)
fit_sens_no_acute_obj <- fit_svy_with_covar_guard(
  qdat %>% filter(is.na(crp_acute) | crp_acute == 0L), "phenoage_accel", c("mixture_qscore", "ckd"),
  realized_covars_int, "Exclude participants with CRP > 1 mg/dL",
  mi_filter_fn = function(df) df[is.na(df$crp_acute) | df$crp_acute == 0L, , drop = FALSE]
)
# Strict complete-case sensitivity stays CC by definition.
fit_sens_strict_obj <- fit_svy_with_covar_guard(
  qdat_strict, "phenoage_accel", c("mixture_qscore", "ckd"), realized_covars_int,
  "Strict complete-case (qgcomp-compatible PFAS subset), same primary adjustment as Model B",
  mi_filter_fn = NULL
)

# Sensitivity: additional adjustment for cardiometabolic conditions potentially on the causal pathway.
covars_cardiometabolic <- c(
  realized_covars_int,
  if (valid_num("diabetes", qdat)) "diabetes",
  if (valid_num("hypertension", qdat)) "hypertension",
  if (valid_num("cvd", qdat)) "cvd"
)
fit_sens_cardiometabolic_obj <- fit_svy_with_covar_guard(
  qdat, "phenoage_accel", c("mixture_qscore", "ckd"), covars_cardiometabolic,
  "Additional adjustment for diabetes, hypertension, CVD (potentially on causal pathway)",
  mi_filter_fn = identity
)
fit_sens_cardiometabolic <- fit_sens_cardiometabolic_obj$fit

tab4 <- bind_rows(
  get_mix_row(fit_sens1_obj, "CKD redefined as eGFR<60", fit_sens1_obj$n_model),
  get_mix_row(fit_sens2_obj, "PhenoAge excluding creatinine only", fit_sens2_obj$n_model),
  get_mix_row(fit_sens3_obj, "PhenoAge excluding creatinine and albumin", fit_sens3_obj$n_model),
  get_mix_row(fit_sens_uw_obj, "Unweighted residualized PhenoAgeAccel (sensitivity)", fit_sens_uw_obj$n_model),
  get_mix_row(fit_sens_both_markers_obj, "Restricted to participants with both eGFR and ACR", fit_sens_both_markers_obj$n_model),
  get_mix_row(fit_sens_no_acute_obj, "Exclude participants with CRP > 1 mg/dL", fit_sens_no_acute_obj$n_model),
  get_mix_row(fit_sens_cardiometabolic_obj, "Additional adjustment for diabetes, hypertension, CVD (potentially on causal pathway)", fit_sens_cardiometabolic_obj$n_model),
  get_mix_row(fit_sens_strict_obj, "Strict complete-case (qgcomp-compatible PFAS subset), same primary adjustment as Model B", fit_sens_strict_obj$n_model)
) %>%
  transmute(
    sensitivity,
    beta = fmt(estimate),
    conf.low = fmt(conf.low),
    conf.high = fmt(conf.high),
    p.value = fmt_p(p.value),
    N,
    source
  )
kable(tab4, caption = "Table 4. Sensitivity analyses (primary adjustment: age, sex, race, PIR, BMI, smoking [3-level], alcohol [3-level], cycle).")
cat("\n*Table 4 note:* Primary and most sensitivity models use the lean adjustment set (age, sex, race/ethnicity, PIR, BMI, smoking [Never/Former/Current], alcohol [Never/Former/Current], cycle). ",
    "Each row reports the `mixture_qscore` coefficient from a survey-weighted model (`svyglm`) with NHANES design variables. ",
    "The `source` column distinguishes MI-pooled (Rubin's rules across m=", config$mice_m,
    " imputations) versus complete-case rows; the strict-CC sensitivity remains complete-case by definition.\n", sep = "")

sens_diagnostics <- bind_rows(
  tibble(sensitivity = "CKD redefined as eGFR<60", N_model = fit_sens1_obj$n_model, dropped_covariates = ifelse(fit_sens1_obj$dropped == "", "None", fit_sens1_obj$dropped), fit_error = ifelse(is.null(fit_sens1_obj$error) || fit_sens1_obj$error == "", "None", fit_sens1_obj$error)),
  tibble(sensitivity = "PhenoAge excluding creatinine only", N_model = fit_sens2_obj$n_model, dropped_covariates = ifelse(fit_sens2_obj$dropped == "", "None", fit_sens2_obj$dropped), fit_error = ifelse(is.null(fit_sens2_obj$error) || fit_sens2_obj$error == "", "None", fit_sens2_obj$error)),
  tibble(sensitivity = "PhenoAge excluding creatinine and albumin", N_model = fit_sens3_obj$n_model, dropped_covariates = ifelse(fit_sens3_obj$dropped == "", "None", fit_sens3_obj$dropped), fit_error = ifelse(is.null(fit_sens3_obj$error) || fit_sens3_obj$error == "", "None", fit_sens3_obj$error)),
  tibble(sensitivity = "Unweighted residualized PhenoAgeAccel (sensitivity)", N_model = fit_sens_uw_obj$n_model, dropped_covariates = ifelse(fit_sens_uw_obj$dropped == "", "None", fit_sens_uw_obj$dropped), fit_error = ifelse(is.null(fit_sens_uw_obj$error) || fit_sens_uw_obj$error == "", "None", fit_sens_uw_obj$error)),
  tibble(sensitivity = "Restricted to participants with both eGFR and ACR", N_model = fit_sens_both_markers_obj$n_model, dropped_covariates = ifelse(fit_sens_both_markers_obj$dropped == "", "None", fit_sens_both_markers_obj$dropped), fit_error = ifelse(is.null(fit_sens_both_markers_obj$error) || fit_sens_both_markers_obj$error == "", "None", fit_sens_both_markers_obj$error)),
  tibble(sensitivity = "Exclude participants with CRP > 1 mg/dL", N_model = fit_sens_no_acute_obj$n_model, dropped_covariates = ifelse(fit_sens_no_acute_obj$dropped == "", "None", fit_sens_no_acute_obj$dropped), fit_error = ifelse(is.null(fit_sens_no_acute_obj$error) || fit_sens_no_acute_obj$error == "", "None", fit_sens_no_acute_obj$error)),
  tibble(sensitivity = "Additional adjustment for diabetes, hypertension, CVD (potentially on causal pathway)", N_model = fit_sens_cardiometabolic_obj$n_model, dropped_covariates = ifelse(fit_sens_cardiometabolic_obj$dropped == "", "None", fit_sens_cardiometabolic_obj$dropped), fit_error = ifelse(is.null(fit_sens_cardiometabolic_obj$error) || fit_sens_cardiometabolic_obj$error == "", "None", fit_sens_cardiometabolic_obj$error)),
  tibble(sensitivity = "Strict complete-case (qgcomp-compatible PFAS subset), same primary adjustment as Model B", N_model = fit_sens_strict_obj$n_model, dropped_covariates = ifelse(fit_sens_strict_obj$dropped == "", "None", fit_sens_strict_obj$dropped), fit_error = ifelse(is.null(fit_sens_strict_obj$error) || fit_sens_strict_obj$error == "", "None", fit_sens_strict_obj$error))
)
kable(sens_diagnostics, caption = "Sensitivity model diagnostics (dropped covariates and fit errors).")

# Supplementary: single-pollutant models with same primary adjustment as Model B.
single_pollutant_out <- lapply(exposure_vars, function(v) {
  sp_obj <- fit_svy_with_covar_guard(
    qdat, "phenoage_accel", c(v, "ckd"), realized_covars_int,
    paste0("Single pollutant: ", v),
    mi_filter_fn = function(df) df[!is.na(df[[v]]), , drop = FALSE]
  )
  diag_row <- tibble(
    exposure = v,
    N_model = if (isTRUE(sp_obj$mi)) sp_obj$n_model_mi else sp_obj$n_model,
    source = if (isTRUE(sp_obj$mi)) "MI-pooled" else "Complete-case",
    dropped_covariates = ifelse(sp_obj$dropped == "", "None", sp_obj$dropped),
    fit_error = ifelse(is.null(sp_obj$error) || sp_obj$error == "", "None", sp_obj$error)
  )
  row <- get_mix_row(sp_obj, paste0("Single pollutant: ", v),
                     n_model_cc = sp_obj$n_model, term_name = v)
  list(
    row  = tibble(
      exposure = v,
      beta = row$estimate, conf.low = row$conf.low, conf.high = row$conf.high,
      p.value = row$p.value, N = row$N, source = row$source
    ),
    diag = diag_row
  )
})
single_pollutant <- bind_rows(lapply(single_pollutant_out, `[[`, "row")) %>%
  mutate(
    beta = fmt(beta),
    conf.low = fmt(conf.low),
    conf.high = fmt(conf.high),
    p.value = fmt_p(p.value)
  )
single_pollutant_diagnostics <- bind_rows(lapply(single_pollutant_out, `[[`, "diag"))

kable(
  single_pollutant,
  caption = "Supplementary Table. Single-pollutant associations with PhenoAge acceleration (survey-weighted, same primary adjustment as Model B); MI-pooled when MICE is active."
)
kable(single_pollutant_diagnostics, caption = "Supplementary model diagnostics (single-pollutant N, source, dropped covariates, fit errors).")


## ----component-driver-models--------------------------------------------------
# Map each PhenoAge biomarker to its raw column, display label/units, and any transform.
# C-reactive protein is log-transformed because that is how it enters the Levine (2018) clock.
component_specs <- tibble::tribble(
  ~outcome,           ~raw,  ~transform, ~label,
  "comp_albumin",     "alb", "identity", "Albumin (g/dL)",
  "comp_creatinine",  "scr", "identity", "Creatinine (mg/dL)",
  "comp_glucose",     "glu", "identity", "Glucose (mg/dL)",
  "comp_logcrp",      "crp", "log",      "C-reactive protein (log mg/dL)",
  "comp_lymph",       "lym", "identity", "Lymphocyte (%)",
  "comp_mcv",         "mcv", "identity", "Mean cell volume (fL)",
  "comp_rdw",         "rdw", "identity", "Red cell distribution width (%)",
  "comp_alp",         "alp", "identity", "Alkaline phosphatase (U/L)",
  "comp_wbc",         "wbc", "identity", "White blood cell count (1000/uL)"
)

# Add component outcomes (raw units; log for CRP) to any data frame (CC subset or imputed df).
# No standardization: beta is the change in the component per one-bin increase in the mixture
# quantile score, matching the interpretation of the main-model mixture coefficient.
add_comp_outcomes <- function(df) {
  for (i in seq_len(nrow(component_specs))) {
    raw <- suppressWarnings(as.numeric(as.character(df[[component_specs$raw[i]]])))
    df[[component_specs$outcome[i]]] <- if (component_specs$transform[i] == "log") log(pmax(raw, 1e-6)) else raw
  }
  df
}

qdat_comp <- add_comp_outcomes(qdat)
# Exact same covariate set as Model B (includes chronological age).
comp_covars <- realized_covars_int

component_out <- lapply(seq_len(nrow(component_specs)), function(i) {
  out <- component_specs$outcome[i]
  lab <- component_specs$label[i]
  obj <- fit_svy_with_covar_guard(
    qdat_comp, out, "mixture_qscore", comp_covars,
    paste0("Mixture vs component: ", lab),
    mi_filter_fn = function(df) add_comp_outcomes(df)
  )
  row <- get_mix_row(obj, lab, n_model_cc = obj$n_model, term_name = "mixture_qscore")
  tibble(
    component = lab,
    beta = row$estimate, conf.low = row$conf.low, conf.high = row$conf.high,
    p.value = row$p.value, N = row$N, source = row$source
  )
})
component_table <- bind_rows(component_out)

component_table_fmt <- component_table %>%
  transmute(
    `PhenoAge component` = component,
    beta = fmt(beta),
    conf.low = fmt(conf.low),
    conf.high = fmt(conf.high),
    p.value = fmt_p(p.value),
    N = N,
    source = source
  )

kable(
  component_table_fmt,
  caption = paste0(
    "Supplementary Table. PFAS mixture (quantile score) associations with each individual PhenoAge ",
    "biomarker component (survey-weighted, same primary adjustment as Model B; MI-pooled when MICE is active). ",
    "Beta is the change in the component per one-quantile increase in the mixture score."
  )
)
cat("\n*Component note:* Each row is a separate survey-weighted linear model (`svyglm`) of one ",
    "PhenoAge biomarker on the PFAS `mixture_qscore`, adjusted for the exact Model B primary ",
    "covariate set (age, sex, race/ethnicity, PIR, BMI, smoking [3-level], alcohol [3-level], cycle). ",
    "Beta is interpreted exactly as in the main models: the change in the component (raw units; log mg/dL for CRP) ",
    "per one-quantile increase in the mixture score. This does not recompute PhenoAge; it decomposes the ",
    "outcome to localize which components the mixture is associated with.\n", sep = "")


