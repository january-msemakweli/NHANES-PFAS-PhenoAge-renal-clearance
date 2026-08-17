# 03: MICE for model covariates (outcome, exposures, and design are never imputed).
# ---- Multiple imputation (MICE) for missing model covariates ----
# We NEVER impute the outcome (phenoage_accel), exposures (PFAS / mixture_qscore), CKD, or design vars.
# Imputation base = qdat_svy (already requires non-missing PhenoAge, CKD, weights/PSU/strata, and full PFAS for mixture_qscore).
mice_target_vars <- intersect(
  c("smoker3", "alcohol3", "bmi", "pir", "education",
    "phys_active", "diabetes", "hypertension", "cvd"),
  names(qdat_svy)
)
mice_predictor_vars <- intersect(
  c("phenoage_accel", "ckd", "mixture_qscore", "age", "sex", "race_eth", "cycle"),
  names(qdat_svy)
)

mids_obj <- NULL
imp_list <- NULL
svy_q_mi <- NULL
mice_pre_n_missing <- if (length(mice_target_vars) > 0) {
  vapply(mice_target_vars, function(v) sum(is.na(qdat_svy[[v]])), integer(1))
} else integer(0)

run_mice <- isTRUE(config$use_mice) && length(mice_target_vars) > 0 && any(mice_pre_n_missing > 0)
if (run_mice) {
  mi_data <- qdat_svy[, unique(c(mice_target_vars, mice_predictor_vars)), drop = FALSE]
  # Force factor types correctly for mice methods.
  for (v in c("smoker3", "alcohol3")) {
    if (v %in% names(mi_data)) mi_data[[v]] <- factor(as.character(mi_data[[v]]),
                                                       levels = c("Never", "Former", "Current"))
  }
  for (v in intersect(c("sex", "race_eth", "cycle", "education"), names(mi_data))) {
    mi_data[[v]] <- factor(as.character(mi_data[[v]]))
  }
  for (v in intersect(c("phys_active", "diabetes", "hypertension", "cvd"), names(mi_data))) {
    mi_data[[v]] <- factor(as.character(mi_data[[v]]))
  }
  for (v in intersect(c("phenoage_accel", "mixture_qscore", "age", "bmi", "pir"), names(mi_data))) {
    mi_data[[v]] <- as.numeric(mi_data[[v]])
  }
  if ("ckd" %in% names(mi_data)) mi_data$ckd <- factor(as.character(mi_data$ckd))

  pred_mat <- matrix(0, nrow = ncol(mi_data), ncol = ncol(mi_data),
                     dimnames = list(names(mi_data), names(mi_data)))
  for (v in mice_target_vars) {
    others <- setdiff(names(mi_data), v)
    pred_mat[v, others] <- 1L
  }

  meth <- character(ncol(mi_data)); names(meth) <- names(mi_data)
  for (v in names(mi_data)) {
    if (!(v %in% mice_target_vars)) { meth[v] <- ""; next }
    if (is.factor(mi_data[[v]])) {
      meth[v] <- if (nlevels(droplevels(mi_data[[v]])) == 2) "logreg" else "polyreg"
    } else {
      meth[v] <- "pmm"
    }
  }
  mice_workers <- max(1L, min(as.integer(config$mice_n_core), parallel::detectCores()))
  message(">>> MICE: imputing ", paste(mice_target_vars, collapse = ", "),
          " with m=", config$mice_m, ", maxit=", config$mice_maxit,
          ", workers=", mice_workers, " (cap ", config$mice_n_core,
          "; machine cores=", parallel::detectCores(), ")")
  mice_err <- NULL
  mids_obj <- tryCatch(
    mice::futuremice(
      mi_data, m = config$mice_m, maxit = config$mice_maxit,
      method = meth, predictorMatrix = pred_mat,
      seed = NA,
      parallelseed = config$mice_seed,
      n.core = mice_workers,
      future.plan = "multisession",
      packages = "mice"
    ),
    error = function(e) { mice_err <<- conditionMessage(e); NULL }
  )
  if (is.null(mids_obj)) {
    message(">>> futuremice() failed (", mice_err, "); falling back to sequential mice::mice()")
    mids_obj <- tryCatch(
      mice::mice(mi_data, m = config$mice_m, maxit = config$mice_maxit,
                 method = meth, predictorMatrix = pred_mat,
                 printFlag = FALSE, seed = config$mice_seed),
      error = function(e) { message("mice() failed: ", conditionMessage(e)); NULL }
    )
  }
  if (is.null(mids_obj)) {
    stop("MICE failed. Primary inference requires the MI sample used in the manuscript. Not falling back to complete-case.")
  }
  if (!is.null(mids_obj)) {
    imp_list <- mitools::imputationList(lapply(seq_len(config$mice_m), function(i) {
      df <- qdat_svy
      completed <- mice::complete(mids_obj, i)
      for (v in mice_target_vars) df[[v]] <- completed[[v]]
      # Re-pin factor levels in the merged frame so model.frame is consistent across imputations.
      if ("smoker3"  %in% names(df)) df$smoker3  <- factor(as.character(df$smoker3),  levels = c("Never", "Former", "Current"))
      if ("alcohol3" %in% names(df)) df$alcohol3 <- factor(as.character(df$alcohol3), levels = c("Never", "Former", "Current"))
      df
    }))
    svy_q_mi <- tryCatch(
      survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined,
                        nest = TRUE, data = imp_list),
      error = function(e) survey::svydesign(id = ~1, strata = NULL,
                                            weights = ~weight_combined, data = imp_list)
    )
  }
}
n_mi_base <- nrow(qdat_svy)
mice_diag <- tibble::tibble(
  variable = mice_target_vars,
  n_missing_pre = as.integer(mice_pre_n_missing),
  pct_missing_pre = round(100 * mice_pre_n_missing / max(1, n_mi_base), 2)
)
if (length(mice_target_vars) > 0) {
  kable(mice_diag, caption = paste0(
    "MICE pre-imputation missingness (within Model B base sample, N = ",
    format(n_mi_base, big.mark = ","), ")."
  ))
}
mice_active <- isTRUE(config$use_mice) && !is.null(svy_q_mi)
message(">>> MICE active for inference: ", mice_active,
        if (mice_active) paste0(" (m=", config$mice_m, ", N=", n_mi_base, ")") else "")
if (isTRUE(config$use_mice) && !isTRUE(mice_active)) {
  stop("MICE did not become active. Refusing complete-case primary inference so the sample stays the manuscript MI sample.")
}
if (isTRUE(mice_active) && n_mi_base != 7633L) {
  warning("MICE base N is ", n_mi_base, ", not the manuscript N=7633. Check the analytic sample construction.")
}

# ---- Helpers: tidy + svycontrast pooling that handle single-fit OR MI list ----
fit_svyglm_pool <- function(formula, design, family = stats::gaussian()) {
  if (inherits(design, "svyimputationList") || inherits(design, "svy.designs.imputationList")) {
    fits <- with(design, svyglm(formula, family = family))
    pooled <- mitools::MIcombine(fits)
    list(fits = fits, pooled = pooled, mi = TRUE)
  } else {
    fit <- svyglm(formula, design = design, family = family)
    list(fit = fit, mi = FALSE)
  }
}

tidy_pool <- function(x, conf.level = 0.95) {
  alpha <- 1 - conf.level
  if (isTRUE(x$mi)) {
    p <- x$pooled
    coefs <- coef(p); ses <- sqrt(diag(vcov(p)))
    df_v <- if (!is.null(p$df) && length(p$df) == length(coefs)) as.numeric(p$df) else rep(Inf, length(coefs))
    df_v[!is.finite(df_v) | df_v <= 0] <- Inf
    t_val <- coefs / ses
    p_val <- 2 * pt(-abs(t_val), df = df_v)
    crit <- qt(1 - alpha / 2, df = df_v)
    tibble::tibble(
      term = names(coefs),
      estimate = unname(coefs),
      std.error = unname(ses),
      statistic = unname(t_val),
      p.value = unname(p_val),
      conf.low = unname(coefs - crit * ses),
      conf.high = unname(coefs + crit * ses)
    )
  } else {
    broom::tidy(x$fit, conf.int = TRUE, conf.level = conf.level)
  }
}

coef_pool <- function(x) {
  if (isTRUE(x$mi)) coef(x$pooled) else coef(x$fit)
}

contrast_pool <- function(x, c_vec, conf.level = 0.95) {
  alpha <- 1 - conf.level
  if (isTRUE(x$mi)) {
    contrasts <- lapply(x$fits, function(f) survey::svycontrast(f, c_vec))
    pooled <- mitools::MIcombine(contrasts)
    est <- as.numeric(coef(pooled)); se <- as.numeric(sqrt(diag(vcov(pooled))))
    df_v <- if (!is.null(pooled$df) && length(pooled$df) == length(est)) as.numeric(pooled$df) else Inf
    df_v[!is.finite(df_v) | df_v <= 0] <- Inf
  } else {
    ctr <- survey::svycontrast(x$fit, c_vec)
    est <- as.numeric(coef(ctr)); se <- as.numeric(SE(ctr))
    df_v <- x$fit$df.residual
    if (is.null(df_v) || !is.finite(df_v)) df_v <- Inf
  }
  t_val <- est / se
  list(
    estimate = est, se = se, df = df_v,
    p.value = 2 * pt(-abs(t_val), df = df_v),
    conf.low = est - qt(1 - alpha / 2, df = df_v) * se,
    conf.high = est + qt(1 - alpha / 2, df = df_v) * se
  )
}

# ---- qgcomp Rubin pooling: run qgcomp.noboot on each imputation, pool psi (overall mixture coef) ----
# Always writes a small diagnostic log next to the Rmd so failures are visible even when
# knitr suppresses message()/warning() output.
qgcomp_mi_pool <- function(formula, expnms, imp_list_obj, q, weight_col = "weight_combined",
                           log_path = "qgcomp_mi_pool_log.txt") {
  log_lines <- character(0)
  log <- function(...) log_lines <<- c(log_lines, paste0(...))
  on.exit({
    try(writeLines(log_lines, log_path), silent = TRUE)
  }, add = TRUE)
  log("[", format(Sys.time()), "] qgcomp_mi_pool called.")
  log("  formula: ", deparse(formula))
  log("  expnms : ", paste(expnms, collapse = ", "))
  log("  q      : ", q)
  log("  weight : ", weight_col)
  if (is.null(imp_list_obj)) { log("ABORT: imp_list_obj is NULL."); return(NULL) }
  imps <- imp_list_obj$imputations
  m <- length(imps)
  log("  m (#imputations): ", m)
  if (m == 0) { log("ABORT: empty imputation list."); return(NULL) }
  log("  imp[1] cols (first 12): ", paste(head(names(imps[[1]]), 12), collapse = ", "))
  results <- lapply(seq_len(m), function(i) {
    df <- as.data.frame(imps[[i]])
    if (!(weight_col %in% names(df))) {
      log("  imp ", i, ": missing weight column '", weight_col, "'.")
      return(NULL)
    }
    df$weight_normalized <- df[[weight_col]] / mean(df[[weight_col]], na.rm = TRUE)
    needed <- intersect(c("phenoage_accel", expnms, weight_col), names(df))
    df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
    n_unwt <- nrow(df)
    if (n_unwt == 0) {
      log("  imp ", i, ": 0 complete-case rows after filtering on outcome/exposure/weight.")
      return(NULL)
    }
    # Bind data into a temporary global so qgcomp's match.call()-based weights handling
    # resolves the weights vector cleanly outside the lapply closure.
    assign("qgcomp_data_local", df, envir = .GlobalEnv)
    on.exit(if (exists("qgcomp_data_local", envir = .GlobalEnv))
      rm("qgcomp_data_local", envir = .GlobalEnv), add = TRUE)
    fit <- tryCatch(
      do.call(
        qgcomp::qgcomp.noboot,
        list(
          f = formula, expnms = expnms, data = qgcomp_data_local,
          family = stats::gaussian(), q = q,
          weights = qgcomp_data_local$weight_normalized
        )
      ),
      error = function(e) {
        log("  imp ", i, ": qgcomp error: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    log("  imp ", i, ": OK (n_unwt=", n_unwt,
        ", psi=", tryCatch(round(fit$psi[1], 4), error = function(e) NA), ")")
    list(fit = fit, n_unwt = n_unwt)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    log("ABORT: all per-imputation fits failed.")
    return(NULL)
  }
  fits <- lapply(results, `[[`, "fit")
  m_eff <- length(fits)
  psis <- vapply(fits, function(f) unname(f$psi[1]), numeric(1))
  vars <- vapply(fits, function(f) unname(f$var.psi[1]), numeric(1))
  ns   <- vapply(results, function(r) as.integer(r$n_unwt), integer(1))
  qbar <- mean(psis)
  ubar <- mean(vars)
  Bvar <- if (m_eff > 1) stats::var(psis) else 0
  Tvar <- ubar + (1 + 1 / m_eff) * Bvar
  se   <- sqrt(Tvar)
  rubin_df <- if (Bvar > 0) (m_eff - 1) * (1 + ubar / ((1 + 1 / m_eff) * Bvar))^2 else Inf
  if (!is.finite(rubin_df) || rubin_df <= 0) rubin_df <- Inf
  t_val <- qbar / se
  log("OK: pooled m_eff=", m_eff, "/", m,
      "; psi=", round(qbar, 4), "; se=", round(se, 4),
      "; df=", round(rubin_df, 1))
  list(
    psi = qbar, se = se, df = rubin_df, m = m_eff, m_attempted = m,
    n_per_fit = ns,
    p.value = 2 * pt(-abs(t_val), df = rubin_df),
    conf.low = qbar - qt(0.975, df = rubin_df) * se,
    conf.high = qbar + qt(0.975, df = rubin_df) * se,
    fits = fits
  )
}
svy_q_strict <- if (nrow(qdat_strict) > 0) {
  qdat_strict_svy <- as.data.frame(qdat_strict)
  for (v in c("phenoage_accel", "mixture_qscore", "ckd", "age", "pir", "bmi", "weight_combined", "psu", "strata")) {
    if (v %in% names(qdat_strict_svy)) qdat_strict_svy[[v]] <- suppressWarnings(as.numeric(as.character(qdat_strict_svy[[v]])))
  }
  qdat_strict_svy$ckd <- as.integer(qdat_strict_svy$ckd)
  qdat_strict_svy$sex <- droplevels(factor(qdat_strict_svy$sex))
  qdat_strict_svy$race_eth <- droplevels(factor(qdat_strict_svy$race_eth))
  qdat_strict_svy$cycle <- droplevels(factor(qdat_strict_svy$cycle))
  if ("smoker3"  %in% names(qdat_strict_svy)) qdat_strict_svy$smoker3  <- droplevels(factor(as.character(qdat_strict_svy$smoker3),  levels = c("Never", "Former", "Current")))
  if ("alcohol3" %in% names(qdat_strict_svy)) qdat_strict_svy$alcohol3 <- droplevels(factor(as.character(qdat_strict_svy$alcohol3), levels = c("Never", "Former", "Current")))
  tryCatch(
    svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE, data = qdat_strict_svy),
    error = function(e) svydesign(id = ~1, strata = NULL, weights = ~weight_combined, data = qdat_strict_svy)
  )
} else {
  NULL
}

