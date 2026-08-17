# 05: participant flow, PFAS descriptives, manuscript figures, sessionInfo.
## ----sample-flow--------------------------------------------------------------
n0 <- nrow(nhanes_raw)
n1 <- sum(dat$age >= config$min_age, na.rm = TRUE)
# Split the old "PhenoAge AND CKD non-missing" step into two consecutive exclusions:
# (a) missing PhenoAge components -> phenoage is NA when any required biomarker is missing;
# (b) missing CKD modelling components -> ckd is NA when both eGFR and ACR cannot be computed.
n2a <- sum(dat$age >= config$min_age & !is.na(dat$phenoage), na.rm = TRUE)
n2b <- sum(dat$age >= config$min_age & !is.na(dat$phenoage) & !is.na(dat$ckd), na.rm = TRUE)
n3 <- nrow(analytic)
n4 <- if (exists("qdat")) nrow(qdat) else NA_integer_
n5 <- if (exists("analytic_qg")) nrow(analytic_qg) else NA_integer_
n6 <- if (exists("qdat_strict")) nrow(qdat_strict) else NA_integer_
n7 <- if (exists("n_modelb")) n_modelb else NA_integer_
n8 <- if (exists("n_modelb_primary")) n_modelb_primary else NA_integer_
mi_label_step <- if (exists("mice_active") && isTRUE(mice_active)) {
  paste0("Model B PRIMARY estimation sample (MICE m=", config$mice_m, "; covariates imputed)")
} else {
  "Model B PRIMARY estimation sample (complete-case; MICE off)"
}

flow <- tibble(
  step = c(
    "Downloaded NHANES participants (selected cycles)",
    paste0("Age >= ", config$min_age),
    "Non-missing PhenoAge components (all 9 biomarkers required)",
    "Non-missing CKD modelling components (eGFR or ACR available)",
    "Non-missing selected weight, PSU, and strata (analytic sample)",
    "Primary interaction model sample (fixed 4-PFAS complete-case score + design vars)",
    "Complete-case qgcomp sample",
    "Strict complete-case (qgcomp-compatible PFAS subset)",
    "Model B complete-case estimation sample (covariate complete-case)",
    mi_label_step
  ),
  included = c(n0, n1, n2a, n2b, n3, n4, n5, n6, n7, n8)
) %>%
  mutate(
    excluded_at_step = pmax(lag(included, default = first(included)) - included, 0),
    excluded_at_step = if_else(row_number() == 1, NA_integer_, as.integer(excluded_at_step)),
    excluded_label = ifelse(is.na(excluded_at_step), "", paste0("Excluded: ", excluded_at_step))
  )

kable(flow, caption = "Sample flow counts (PRISMA-style): inclusion and exclusion at each analytic step")

flow_plot <- flow %>%
  mutate(
    y = rev(seq_len(n())),
    box_label = paste0(step, "\nN = ", included)
  )

flow_fig <- ggplot(flow_plot, aes(x = 1, y = y)) +
  geom_segment(
    data = flow_plot %>% filter(y > min(y)),
    aes(x = 1, xend = 1, y = y - 0.42, yend = y - 0.95),
    arrow = arrow(length = unit(0.12, "inches")),
    linewidth = 0.5
  ) +
  geom_label(aes(label = box_label), size = 3.2, label.size = 0.2, fill = "white") +
  geom_text(aes(x = 1.62, label = excluded_label), hjust = 0, size = 3.1, color = "gray30") +
  xlim(0.5, 2.2) +
  ylim(0.5, max(flow_plot$y) + 0.8) +
  labs(
    title = "Participant flow from downloaded sample to analytic samples",
    subtitle = "PRISMA-style flow for NHANES analytic cohort development",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 11)

ggsave("figures/participant_flow.jpg", plot = flow_fig,
       width = 9, height = 7, dpi = 300, device = "jpeg", quality = 95, bg = "white")


## ----pfas-descriptives--------------------------------------------------------
# Raw serum PFAS (ng/mL) for the four primary chemicals used in the mixture.
pfas_raw_map <- c(
  "PFOS"  = "LBXPFOS",
  "PFOA"  = "LBXPFOA",
  "PFNA"  = "LBXPFNA",
  "PFHxS" = "LBXPFHS"
)
pfas_present <- pfas_raw_map[pfas_raw_map %in% names(analytic)]

# Wide -> long for plotting (one panel per chemical) using the analytic sample.
# Coerce all PFAS columns to numeric first (some NHANES vars arrive as character).
pfas_wide <- analytic %>%
  dplyr::select(dplyr::any_of(unname(pfas_present))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), function(x) {
    suppressWarnings(as.numeric(trimws(as.character(x))))
  }))

pfas_long <- pfas_wide %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "var", values_to = "value") %>%
  dplyr::filter(!is.na(value), value > 0) %>%
  dplyr::mutate(chemical = factor(names(pfas_present)[match(var, pfas_present)],
                                  levels = names(pfas_present)))

# Geometric mean / SD helpers (applied on positive values only)
gm_mean <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) NA_real_ else exp(mean(log(x)))
}
gm_sd <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) NA_real_ else exp(stats::sd(log(x)))
}

pfas_desc <- pfas_long %>%
  dplyr::group_by(Chemical = chemical) %>%
  dplyr::summarise(
    N        = sum(!is.na(value)),
    Mean     = mean(value, na.rm = TRUE),
    SD       = stats::sd(value, na.rm = TRUE),
    Median   = stats::median(value, na.rm = TRUE),
    Q1       = stats::quantile(value, 0.25, na.rm = TRUE),
    Q3       = stats::quantile(value, 0.75, na.rm = TRUE),
    IQR      = Q3 - Q1,
    Min      = min(value, na.rm = TRUE),
    Max      = max(value, na.rm = TRUE),
    GeoMean  = gm_mean(value),
    GeoSD    = gm_sd(value),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(dplyr::across(c(Mean, SD, Median, Q1, Q3, IQR, Min, Max, GeoMean, GeoSD),
                              ~ round(.x, 3)))

knitr::kable(
  pfas_desc,
  caption = "PFAS descriptive statistics in the analytic sample (serum, ng/mL). Geometric mean / SD computed on positive values only."
)

# Histogram panels (2x2): export both linear axis (true raw scale) and log10-axis (skewness readability).
hist_base <- function() {
  ggplot(pfas_long, aes(x = value)) +
    geom_histogram(bins = 40, fill = "grey80", color = "black", linewidth = 0.2) +
    facet_wrap(~ chemical, ncol = 2, scales = "free") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}
pfas_hist_linear <- hist_base() +
  labs(x = "Serum concentration (ng/mL)", y = "Count",
       title = "PFAS distributions (linear x-axis)")
pfas_hist_log10 <- hist_base() +
  scale_x_log10() +
  labs(x = "Serum concentration (ng/mL, log10 axis)", y = "Count",
       title = "PFAS distributions (log10-scaled axis)")

dir.create("figures", showWarnings = FALSE)

pfas_hist_linear / pfas_hist_log10 +
  plot_layout(ncol = 1) &
  theme(plot.title = element_text(face = "bold", size = 11))


## ----figure-------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
n_chems <- length(primary_exposure_vars)

# Figure 1: survey-weighted predicted lines per CKD stratum.
# Mirror the proven standalone logic: filter NAs, fall back design->weights-only->lm,
# normalize predict() output (svyglm returns svystat; lm returns atomic).
fit_one_stratum <- function(sub_df) {
  needed <- c("phenoage_accel", "mixture_qscore", "age", "weight_combined", "psu", "strata")
  sub_df <- sub_df[stats::complete.cases(sub_df[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(sub_df) == 0) return(NULL)
  fit <- tryCatch(
    survey::svyglm(phenoage_accel ~ mixture_qscore + age,
                   design = survey::svydesign(id = ~psu, strata = ~strata,
                                              weights = ~weight_combined, nest = TRUE,
                                              data = sub_df),
                   family = gaussian()),
    error = function(e) NULL)
  if (is.null(fit)) {
    fit <- tryCatch(
      survey::svyglm(phenoage_accel ~ mixture_qscore + age,
                     design = survey::svydesign(id = ~1, weights = ~weight_combined,
                                                data = sub_df),
                     family = gaussian()),
      error = function(e) NULL)
  }
  if (is.null(fit)) {
    fit <- tryCatch(stats::lm(phenoage_accel ~ mixture_qscore + age,
                              data = sub_df, weights = sub_df$weight_combined),
                    error = function(e) NULL)
  }
  fit
}

pred_lines <- tryCatch({
  full_score_range <- range(qdat$mixture_qscore, na.rm = TRUE)
  pdat_list <- lapply(c(0L, 1L), function(ck) {
    sub_df <- qdat %>% filter(ckd == ck)
    fit_pred <- fit_one_stratum(sub_df)
    if (is.null(fit_pred)) return(NULL)
    score_seq <- seq(full_score_range[1], full_score_range[2], length.out = 100)
    newdat <- data.frame(
      mixture_qscore = score_seq,
      age = mean(sub_df$age, na.rm = TRUE)
    )
    pr <- tryCatch(
      predict(fit_pred, newdata = newdat, se.fit = TRUE, type = "response"),
      error = function(e) NULL)
    if (is.null(pr)) return(NULL)
    fit_v <- as.numeric(if (is.list(pr)) pr$fit else pr)
    se_v  <- if (is.list(pr) && !is.null(pr$se.fit)) as.numeric(pr$se.fit)
             else as.numeric(survey::SE(pr))
    newdat$fit <- fit_v
    newdat$se  <- se_v
    newdat$ckd <- ifelse(ck == 1L, "CKD", "No CKD")
    cf <- coef(fit_pred)
    vc <- tryCatch(stats::vcov(fit_pred), error = function(e) NULL)
    slope <- as.numeric(cf["mixture_qscore"])
    se_b  <- if (!is.null(vc)) sqrt(vc["mixture_qscore", "mixture_qscore"]) else NA_real_
    newdat$slope    <- slope
    newdat$slope_lo <- if (is.finite(se_b)) slope - 1.96 * se_b else NA_real_
    newdat$slope_hi <- if (is.finite(se_b)) slope + 1.96 * se_b else NA_real_
    newdat
  })
  pdat_list <- Filter(Negate(is.null), pdat_list)
  if (length(pdat_list) == 0) NULL else bind_rows(pdat_list)
}, error = function(e) NULL)

# Overall marginal direction: as mixture increases, what happens to phenoage_accel (pooled model).
overall_direction <- tryCatch({
  svy_full <- svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE, data = qdat)
  fit_overall <- svyglm(phenoage_accel ~ mixture_qscore + age, design = svy_full, family = gaussian())
  score_seq <- seq(min(qdat$mixture_qscore, na.rm = TRUE), max(qdat$mixture_qscore, na.rm = TRUE), length.out = 100)
  newdat <- data.frame(
    mixture_qscore = score_seq,
    age = mean(qdat$age, na.rm = TRUE)
  )
  preds <- predict(fit_overall, newdata = newdat, se.fit = TRUE, type = "response")
  newdat$fit <- as.numeric(preds$fit)
  newdat
}, error = function(e) NULL)

# Faceted hex density of (mixture quantile score, PhenoAge acceleration) by CKD,
# with a thick survey-weighted trend line + 95% CI ribbon and per-facet slope label.
qdat_plot <- qdat %>%
  mutate(ckd_label = factor(ifelse(ckd == 0L, "No CKD", "CKD"),
                            levels = c("No CKD", "CKD")))
pred_plot <- if (!is.null(pred_lines)) {
  pred_lines %>%
    mutate(ckd_label = factor(ckd, levels = c("No CKD", "CKD")))
} else NULL

# Per-facet slope annotation: prefer MI-pooled marginal effects from fit_int
# (design-based contrasts already used in Table 3); fall back to per-stratum fits.
mi_ann <- tryCatch({
  if (exists("fit_int") && !is.null(fit_int)) {
    cf_names_a <- names(coef_pool(fit_int))
    int_name_a <- cf_names_a[stringr::str_detect(cf_names_a, "mixture_qscore:ckd")]
    if (length(int_name_a) == 1 && "mixture_qscore" %in% cf_names_a) {
      c0 <- setNames(rep(0, length(cf_names_a)), cf_names_a); c0["mixture_qscore"] <- 1
      c1 <- setNames(rep(0, length(cf_names_a)), cf_names_a); c1["mixture_qscore"] <- 1; c1[int_name_a] <- 1
      ctr0 <- contrast_pool(fit_int, c0)
      ctr1 <- contrast_pool(fit_int, c1)
      tibble::tibble(
        ckd_label = factor(c("No CKD", "CKD"), levels = c("No CKD", "CKD")),
        slope    = c(ctr0$estimate, ctr1$estimate),
        slope_lo = c(ctr0$conf.low, ctr1$conf.low),
        slope_hi = c(ctr0$conf.high, ctr1$conf.high)
      )
    } else NULL
  } else NULL
}, error = function(e) NULL)
ann_df <- if (!is.null(mi_ann)) mi_ann else if (!is.null(pred_plot)) {
  pred_plot %>%
    group_by(ckd_label) %>%
    summarise(slope = unique(slope), slope_lo = unique(slope_lo),
              slope_hi = unique(slope_hi), .groups = "drop")
} else NULL
if (!is.null(ann_df)) {
  ann_df <- ann_df %>%
    mutate(label = ifelse(
      is.finite(slope_lo) & is.finite(slope_hi),
      sprintf("beta = %+0.2f years per quantile increase (95%% CI %+0.2f, %+0.2f)",
              slope, slope_lo, slope_hi),
      sprintf("beta = %+0.2f years per quantile increase", slope)
    ))
}

theme_facet_box <- theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.key.height = unit(0.6, "lines"))

trend_color <- "#000000" # bright red, prints as a strong dark gray

fig1_hex <- ggplot(qdat_plot, aes(x = mixture_qscore, y = phenoage_accel)) +
  geom_hex(bins = 35) +
  scale_fill_gradient(low = "grey92", high = "black", name = "Count") +
  facet_wrap(~ ckd_label, ncol = 2) +
  scale_x_continuous(limits = range(qdat_plot$mixture_qscore, na.rm = TRUE),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Mixture quantile score",
       y = "PhenoAge acceleration (years)") +
  theme_facet_box
if (!is.null(pred_plot)) {
  fig1_hex <- fig1_hex +
    geom_ribbon(data = pred_plot,
                aes(x = mixture_qscore,
                    ymin = fit - 1.96 * se, ymax = fit + 1.96 * se,
                    group = ckd_label),
                fill = trend_color, alpha = 0.20, inherit.aes = FALSE) +
    geom_line(data = pred_plot,
              aes(x = mixture_qscore, y = fit, group = ckd_label),
              color = trend_color, linewidth = 1.3, inherit.aes = FALSE)
}
if (!is.null(ann_df)) {
  y_top <- max(qdat_plot$phenoage_accel, na.rm = TRUE)
  x_left <- min(qdat_plot$mixture_qscore, na.rm = TRUE)
  fig1_hex <- fig1_hex +
    geom_text(data = ann_df,
              aes(x = x_left, y = y_top, label = label),
              hjust = 0, vjust = 1.2, size = 3.2, color = trend_color,
              fontface = "bold", inherit.aes = FALSE)
}
fig1_hex
ggsave("figures/fig1_hex_facet_CKD.jpg", plot = fig1_hex,
       width = 9.5, height = 5, dpi = 600, device = "jpeg",
       quality = 95, bg = "white")

# Figure 2: forest plot of marginal mixture effects (design-based).
fig2 <- NULL
if (!is.null(fit_int)) {
  cf_names <- names(coef_pool(fit_int))
  int_name <- cf_names[str_detect(cf_names, "mixture_qscore:ckd")]
  if (length(int_name) == 1 && "mixture_qscore" %in% cf_names) {
    c_ckd0 <- setNames(rep(0, length(cf_names)), cf_names); c_ckd0["mixture_qscore"] <- 1
    c_ckd1 <- setNames(rep(0, length(cf_names)), cf_names); c_ckd1["mixture_qscore"] <- 1; c_ckd1[int_name] <- 1
    ctr0 <- contrast_pool(fit_int, c_ckd0)
    ctr1 <- contrast_pool(fit_int, c_ckd1)

    mk_row <- function(ctr, grp) {
      tibble(
        group = grp,
        estimate = ctr$estimate,
        conf.low = ctr$conf.low,
        conf.high = ctr$conf.high,
        p.value = ctr$p.value
      )
    }
    forest_df <- bind_rows(
      mk_row(ctr0, "No CKD"),
      mk_row(ctr1, "CKD")
    ) %>%
      mutate(group = factor(group, levels = c("No CKD", "CKD")))

    fig2 <- ggplot(forest_df, aes(x = estimate, y = group, color = group)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.8) +
      geom_point(size = 2.8) +
      geom_text(aes(label = ifelse(is.na(p.value), "",
                                   ifelse(p.value < 0.001, "p<0.001",
                                          paste0("p=", format(round(p.value, 3), nsmall = 3))))),
                nudge_y = 0.2, size = 3, color = "black", show.legend = FALSE) +
      scale_color_manual(values = c("No CKD" = "#1b9e77", "CKD" = "#d95f02")) +
      labs(
        title = "Marginal mixture effects by CKD status",
        subtitle = "Design-based contrasts from interaction model",
        x = "Beta (95% CI) for mixture score",
        y = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "none")
    fig2
    ggsave("figures/fig2_forest_marginal_effects.jpg", plot = fig2, width = 7, height = 4, dpi = 300, device = "jpeg", quality = 95, bg = "white")
  }
}

# Figure 3: MICE convergence diagnostics in MCMC-style plate
# (Trace | Autocorrelation | Kernel density) per variable, boxed panels, no overall title.
fig3 <- NULL
if (exists("mids_obj") && !is.null(mids_obj)) {
  miss_n <- vapply(mice_target_vars, function(v) sum(is.na(qdat_svy[[v]])), integer(1))
  top_vars <- names(sort(miss_n, decreasing = TRUE))[seq_len(min(4, length(miss_n)))]
  top_vars <- top_vars[miss_n[top_vars] > 0]

  if (length(top_vars) >= 1) {
    # ---- Trace: chainMean is [variable, iteration, chain] ----
    cm <- mids_obj$chainMean
    trace_df <- as.data.frame.table(cm, responseName = "value", stringsAsFactors = FALSE) %>%
      dplyr::rename(variable = Var1, iter_lab = Var2, chain_lab = Var3) %>%
      dplyr::mutate(
        iter  = as.integer(sub("[^0-9]+", "", iter_lab)),
        chain = as.integer(sub("[^0-9]+", "", chain_lab))
      ) %>%
      dplyr::filter(variable %in% top_vars, !is.na(value))

    # ---- ACF: per chain, averaged across chains ----
    max_iter <- max(trace_df$iter, na.rm = TRUE)
    lag_max  <- max(1L, min(8L, max_iter - 1L))
    acf_ci   <- 1.96 / sqrt(max_iter) # white-noise 95% CI for ACF (Bartlett)
    acf_df <- trace_df %>%
      dplyr::group_by(variable, chain) %>%
      dplyr::arrange(iter, .by_group = TRUE) %>%
      dplyr::summarise(
        lag = list(0:lag_max),
        acf = list(as.numeric(stats::acf(value, plot = FALSE, lag.max = lag_max)$acf)),
        .groups = "drop"
      ) %>%
      tidyr::unnest(c(lag, acf)) %>%
      dplyr::group_by(variable, lag) %>%
      dplyr::summarise(acf = mean(acf, na.rm = TRUE), .groups = "drop")

    # ---- Density: observed (.imp == 0) vs imputed ----
    long_full <- mice::complete(mids_obj, "long", include = TRUE) %>%
      dplyr::mutate(source = ifelse(.imp == 0, "Observed", "Imputed"))
    var_pieces <- lapply(top_vars, function(v) {
      x <- long_full[[v]]
      val <- if (is.factor(x)) as.numeric(x) else suppressWarnings(as.numeric(x))
      data.frame(source = long_full$source, variable = v, value_num = val,
                 stringsAsFactors = FALSE)
    })
    long_imp <- do.call(rbind, var_pieces) %>%
      dplyr::filter(!is.na(value_num)) %>%
      dplyr::mutate(source = factor(source, levels = c("Observed", "Imputed")))

    # ---- ResearchGate-style theme: boxed, no grid, B&W ----
    bw_boxed <- function(base_size = 9) {
      theme_bw(base_size = base_size) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5, face = "plain", size = base_size + 1),
          legend.position = "none"
        )
    }

    build_row <- function(var_name, col_titles) {
      tdf <- dplyr::filter(trace_df, variable == var_name)
      adf <- dplyr::filter(acf_df,   variable == var_name)
      ddf <- dplyr::filter(long_imp, variable == var_name)

      p_trace <- ggplot(tdf, aes(x = iter, y = value, group = chain)) +
        geom_line(alpha = 0.35, linewidth = 0.3, color = "black") +
        scale_x_continuous(expand = c(0.02, 0)) +
        labs(x = "Iteration", y = var_name,
             title = if (col_titles) "Traceplot" else NULL) +
        bw_boxed()

      p_acf <- ggplot(adf, aes(x = lag, y = acf)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.35) +
        geom_hline(yintercept = c(-acf_ci, acf_ci), color = "gray55", linewidth = 0.25) +
        geom_segment(aes(xend = lag, yend = 0), linewidth = 0.35, color = "black") +
        geom_point(size = 1.6, color = "gray20") +
        scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
        labs(x = "Lag", y = "ACF",
             title = if (col_titles) "Autocorrelation" else NULL) +
        bw_boxed()

      # KDE: observed (solid) vs imputed (dashed); vertical lines at observed 2.5/97.5%
      q_obs <- stats::quantile(ddf$value_num[ddf$source == "Observed"],
                               probs = c(0.025, 0.975), na.rm = TRUE)
      p_kde <- ggplot(ddf, aes(x = value_num, linetype = source)) +
        geom_density(color = "black", linewidth = 0.45, na.rm = TRUE) +
        geom_vline(xintercept = q_obs, color = "gray45",
                   linewidth = 0.35, linetype = "solid") +
        scale_linetype_manual(values = c("Observed" = "solid", "Imputed" = "dashed")) +
        labs(x = "Value", y = "Density",
             title = if (col_titles) "Kernel density estimation" else NULL) +
        bw_boxed()

      p_trace + p_acf + p_kde + plot_layout(ncol = 3)
    }

    rows <- lapply(seq_along(top_vars), function(i) build_row(top_vars[i], col_titles = (i == 1L)))
    fig3 <- wrap_plots(rows, ncol = 1)

    fig3
    ggsave("figures/fig3_mice_convergence.jpg", plot = fig3,
           width = 8.2, height = 1.55 * length(top_vars) + 0.6,
           dpi = 300, device = "jpeg", quality = 95, bg = "white")
  }
}


## ----session-info-------------------------------------------------------------
sessionInfo()

