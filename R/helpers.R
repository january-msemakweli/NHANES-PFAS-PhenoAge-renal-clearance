# Helper functions extracted from the analysis Rmd.
## ----helpers------------------------------------------------------------------
cycle_to_suffix <- function(cycle) {
  map <- c(
    "1999-2000" = "A", "2001-2002" = "B", "2003-2004" = "C", "2005-2006" = "D",
    "2007-2008" = "E", "2009-2010" = "F", "2011-2012" = "G", "2013-2014" = "H",
    "2015-2016" = "I", "2017-2018" = "J", "2019-2020" = "K"
  )
  out <- unname(map[cycle])
  if (is.na(out)) stop("Cycle not mapped: ", cycle)
  out
}

safe_nhanes <- function(table_name, verbose = TRUE) {
  out <- tryCatch(nhanesA::nhanes(table_name), error = function(e) NULL)
  if (is.null(out) && verbose) message("Unavailable table: ", table_name)
  out
}

resolve_table <- function(cycle, candidates, domain = NULL, verbose = TRUE) {
  suf <- cycle_to_suffix(cycle)
  to_try <- paste0(candidates, "_", suf)
  for (tb in to_try) {
    dat <- safe_nhanes(tb, verbose = FALSE)
    if (!is.null(dat)) {
      if (verbose) message("Using table ", tb, " for cycle ", cycle)
      return(list(code = candidates[1], table = tb, data = dat))
    }
  }
  NULL
}

coalesce_num <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  v <- lapply(df[vars], function(x) suppressWarnings(as.numeric(as.character(x))))
  dplyr::coalesce(!!!v)
}

parse_yes_no <- function(x) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  x_chr <- toupper(trimws(as.character(x)))
  case_when(
    x_num == 1 ~ 1,
    x_num == 2 ~ 0,
    str_detect(x_chr, "^YES$") ~ 1,
    str_detect(x_chr, "^NO$") ~ 0,
    TRUE ~ NA_real_
  )
}

parse_smq040 <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  case_when(
    x_num %in% c(1, 2) ~ 1,
    x_num == 3 ~ 0,
    str_detect(x_chr, "EVERY DAY|SOME DAYS") ~ 1,
    str_detect(x_chr, "^NOT AT ALL$") ~ 0,
    TRUE ~ NA_real_
  )
}

# SMQ020 = "Smoked at least 100 cigarettes in entire life?" (the gate question for SMQ040).
# 1 = Yes (lifetime smoker), 2 = No (never smoker), 7/9 = Refused/Don't know -> NA.
parse_smq020 <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  case_when(
    x_num == 1 ~ 1,
    x_num == 2 ~ 0,
    str_detect(x_chr, "^YES$") ~ 1,
    str_detect(x_chr, "^NO$")  ~ 0,
    TRUE ~ NA_real_
  )
}

# Lifetime "ever drinker" indicator harmonized across NHANES cycles.
#   ALQ101 (1999-2016) = ">=12 alcoholic drinks in any one year ever?"  (1=Yes / 2=No)
#   ALQ110 (1999-2016) = ">=12 drinks in entire life?"  asked when ALQ101=No (1=Yes / 2=No)
#   ALQ111 (2017+)     = "Ever had a drink of any kind of alcohol?"     (1=Yes / 2=No)
# Returns 1 if ANY indicator is Yes, 0 if all available indicators are No, NA otherwise
# (refused / don't know / not asked).
parse_alq_ever <- function(alq101, alq110, alq111) {
  e101 <- parse_yes_no(alq101)
  e110 <- parse_yes_no(alq110)
  e111 <- parse_yes_no(alq111)
  case_when(
    !is.na(e111) & e111 == 1 ~ 1,
    !is.na(e101) & e101 == 1 ~ 1,
    !is.na(e110) & e110 == 1 ~ 1,
    !is.na(e111) & e111 == 0 ~ 0,
    !is.na(e101) & e101 == 0 & (is.na(e110) | e110 == 0) ~ 0,
    TRUE ~ NA_real_
  )
}

# "Current drinker" (drank in past 12 months) harmonized across cycles.
#   ALQ121 (2013-2018) = "Past 12 mo, how often drink alcoholic beverages?"
#                         0 = Never in last year, 1-10 = various frequencies, 77/99 = Refused/DK.
#   ALQ120Q (1999-2016) = "How often drink alcoholic beverages?" (numeric # of times)
#                         in the unit specified by ALQ120U (week/month/year).
#                         0 = Never; 1-365 = numeric answer; 777/999 = Refused/DK.
# Returns 1 if any-frequency drinker, 0 if explicit "never in past year", NA otherwise.
parse_alq_current_yr <- function(alq121, alq120q) {
  c121 <- suppressWarnings(as.numeric(as.character(alq121)))
  c120 <- suppressWarnings(as.numeric(as.character(alq120q)))
  case_when(
    !is.na(c121) & c121 %in% 1:10 ~ 1,
    !is.na(c121) & c121 == 0 ~ 0,
    !is.na(c120) & c120 >= 1 & c120 <= 365 ~ 1,
    !is.na(c120) & c120 == 0 ~ 0,
    TRUE ~ NA_real_
  )
}

parse_phq_item <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  case_when(
    x_num == 1 ~ 0, x_num == 2 ~ 1, x_num == 3 ~ 2, x_num == 4 ~ 3,
    x_num %in% 0:3 ~ x_num,
    str_detect(x_chr, "^NOT AT ALL$") ~ 0,
    str_detect(x_chr, "^SEVERAL DAYS$") ~ 1,
    str_detect(x_chr, "^MORE THAN HALF THE DAYS$") ~ 2,
    str_detect(x_chr, "^NEARLY EVERY DAY$") ~ 3,
    TRUE ~ NA_real_
  )
}

yesno01 <- function(x) {
  xc <- toupper(trimws(as.character(x)))
  xn <- suppressWarnings(as.numeric(as.character(x)))
  case_when(
    xn == 1 ~ 1,
    xn == 2 ~ 0,
    xc %in% c("YES", "1") ~ 1,
    xc %in% c("NO", "2") ~ 0,
    TRUE ~ NA_real_
  )
}

compute_egfr_2021 <- function(scr_mgdl, age, sex_raw) {
  sex_chr <- toupper(trimws(as.character(sex_raw)))
  sex_num <- suppressWarnings(as.numeric(as.character(sex_raw)))
  male <- case_when(
    sex_num == 1 ~ TRUE,
    sex_num == 2 ~ FALSE,
    sex_chr == "MALE" ~ TRUE,
    sex_chr == "FEMALE" ~ FALSE,
    TRUE ~ NA
  )
  kappa <- ifelse(male, 0.9, 0.7)
  alpha <- ifelse(male, -0.302, -0.241)
  sex_mult <- ifelse(male, 1, 1.012)
  142 * pmin(scr_mgdl / kappa, 1)^alpha * pmax(scr_mgdl / kappa, 1)^(-1.200) * (0.9938^age) * sex_mult
}

compute_acr <- function(ua, ucr, acr_direct = NULL) {
  calc <- ifelse(!is.na(ua) & !is.na(ucr) & ucr > 0, (ua / ucr) * 100, NA_real_)
  dplyr::coalesce(acr_direct, calc)
}

compute_phenoage_components <- function(df, include_scr = TRUE, include_alb = TRUE) {
  # Levine et al. (2018) coefficients expect: albumin g/L, creatinine µmol/L, glucose mmol/L.
  # NHANES provides: albumin g/dL, creatinine mg/dL, glucose often mg/dL (we standardize to mg/dL in pipeline).
  alb_gL   <- df$alb * 10
  scr_umol <- df$scr * 88.4
  glu_mmol <- df$glu / 18.018
  alb_term <- if (include_alb) (-0.0336 * alb_gL) else 0
  scr_term <- if (include_scr) ( 0.0095 * scr_umol) else 0
  glu_term <- 0.1953 * glu_mmol

  xb <- -19.907 +
    alb_term +
    scr_term +
    glu_term +
    0.0954 * log(pmax(df$crp, 1e-6)) -
    0.0120 * df$lym +
    0.0268 * df$mcv +
    0.3306 * df$rdw +
    0.00188 * df$alp +
    0.0554 * df$wbc +
    0.0804 * df$age

  xb <- pmin(pmax(xb, -50), 50)
  gamma_const <- (exp(120 * 0.0076927) - 1) / 0.0076927
  hazard <- exp(xb) * gamma_const
  mort <- -expm1(-hazard)
  mort <- pmin(pmax(mort, 1e-12), 1 - 1e-12)
  pheno <- 141.50 + log(-0.00553 * log1p(-mort)) / 0.090165

  list(xb = xb, mort = mort, phenoage = pheno)
}

safe_phenoage <- function(df) {
  compute_phenoage_components(df, include_scr = TRUE, include_alb = TRUE)$phenoage
}

safe_qgcomp_extract <- function(fit_obj, model_label, n_obs) {
  if (is.null(fit_obj)) {
    return(tibble(model = model_label, beta = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_, N = n_obs))
  }
  cf <- coef(fit_obj)
  psi_name <- grep("^psi", names(cf), value = TRUE)[1]
  if (is.na(psi_name)) psi_name <- grep("psi", names(cf), value = TRUE)[1]
  if (is.na(psi_name)) {
    message("qgcomp terms for [", model_label, "]: no psi found in ", paste(names(cf), collapse = ", "))
    return(tibble(model = model_label, beta = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_, N = n_obs))
  }
  beta <- unname(cf[psi_name])
  vc <- tryCatch(vcov(fit_obj), error = function(e) NULL)
  if (!is.null(vc) && psi_name %in% rownames(vc) && psi_name %in% colnames(vc)) {
    se <- sqrt(vc[psi_name, psi_name])
    lo <- beta - 1.96 * se
    hi <- beta + 1.96 * se
    p_val <- 2 * pnorm(-abs(beta / se))
  } else {
    lo <- hi <- p_val <- NA_real_
  }
  tibble(model = model_label, beta = beta, conf.low = lo, conf.high = hi, p.value = p_val, N = n_obs)
}

resolve_weight_var <- function(df, weight_type = "auto", exposure_family = "pfas") {
  cols <- names(df)
  if (weight_type != "auto" && weight_type != "MEC") {
    if (!(weight_type %in% cols)) stop("Requested weight variable not found: ", weight_type)
    return(list(var = weight_type, note = paste("User-specified weight:", weight_type)))
  }
  if ("WTMEC2YR" %in% cols) return(list(var = "WTMEC2YR", note = "Generic fallback: using MEC weight WTMEC2YR because no exposure-specific subsample weight has been resolved. In PFAS analyses, final design should normally use the PFAS-specific subsample weight identified by the audit/auto-selection logic."))
  stop("No supported NHANES weight variable found.")
}

apply_lod_rule <- function(df, varname) {
  x <- suppressWarnings(as.numeric(as.character(df[[varname]])))

  base <- sub("^(LBX|URX)", "", varname)
  lc_candidates <- c(
    paste0(varname, "LC"),
    paste0("LBD", base, "LC"),
    paste0("URD", base, "LC"),
    # NHANES often uses *L (not *LC) as the detection-limit flag.
    paste0(varname, "L"),
    paste0("LBD", base, "L"),
    paste0("URD", base, "L")
  )
  lc_candidates <- lc_candidates[lc_candidates %in% names(df)]

  if (length(lc_candidates) > 0) {
    lc <- df[[lc_candidates[1]]]
    lc_num <- suppressWarnings(as.numeric(as.character(lc)))
    lc_chr <- toupper(trimws(as.character(lc)))
    below <- (!is.na(lc_num) & lc_num == 1) | str_detect(lc_chr, "BELOW|<")
    below[is.na(below)] <- FALSE
    # NHANES convention: when LC==1, the concentration column typically holds the LOD value.
    # Standard replacement: value / sqrt(2).
    x[below] <- x[below] / sqrt(2)
    # If value is still NA/non-positive after that, use smallest observed LOD / sqrt(2).
    lod_vals <- x[below & x > 0]
    fallback <- if (length(lod_vals) > 0) min(lod_vals, na.rm = TRUE) else NA_real_
    x[below & (is.na(x) | x <= 0)] <- fallback
    return(x)
  }
  # If no LOD indicator exists, do not impute missing values.
  x
}


