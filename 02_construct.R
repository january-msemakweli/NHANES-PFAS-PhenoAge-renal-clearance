# 02: PhenoAge, CKD, PFAS, weights, analytic sample, mixture score.
## ----construct-vars-----------------------------------------------------------
# Defensive guard: ensure cycle-dependent NHANES columns exist before referencing them
# in the mutate below. cycle_data harmonization (Data Acquisition section) already pads
# cycles to the union of column names, so these will normally be present -- but if the
# user reconfigures `config$cycles` to drop certain ranges, some columns can vanish from
# the union (e.g., ALQ142 only exists in 2017+; PAD200/PAD320 only exist in 1999-2006).
# This block guarantees they exist as NA so the harmonized derivations below still run.
required_optional_cols <- c(
  # Smoking
  "SMQ020", "SMQ040",
  # Alcohol (cycle-dependent question wording -- see parser comments above)
  "ALQ101", "ALQ110", "ALQ111", "ALQ120Q", "ALQ121", "ALQ142", "ALQ151",
  # Physical activity (GPAQ in 2007+, legacy module in 1999-2006)
  "PAD200", "PAD320", "PAQ650", "PAQ665"
)
for (.col in required_optional_cols) {
  if (!.col %in% names(nhanes_raw)) nhanes_raw[[.col]] <- NA_character_
}
rm(.col)

dpq_vars <- names(nhanes_raw)[str_detect(names(nhanes_raw), "^DPQ0")]

dat <- nhanes_raw %>%
  mutate(across(any_of(dpq_vars), parse_phq_item)) %>%
  mutate(
    age = suppressWarnings(as.numeric(as.character(RIDAGEYR))),
    sex = case_when(
      suppressWarnings(as.numeric(as.character(RIAGENDR))) == 1 ~ "Male",
      suppressWarnings(as.numeric(as.character(RIAGENDR))) == 2 ~ "Female",
      toupper(trimws(as.character(RIAGENDR))) == "MALE" ~ "Male",
      toupper(trimws(as.character(RIAGENDR))) == "FEMALE" ~ "Female",
      TRUE ~ NA_character_
    ),
    race_eth = {
      r3 <- suppressWarnings(as.numeric(as.character(RIDRETH3)))
      r1 <- suppressWarnings(as.numeric(as.character(RIDRETH1)))
      r3_chr <- toupper(trimws(as.character(RIDRETH3)))
      r1_chr <- toupper(trimws(as.character(RIDRETH1)))
      case_when(
        !is.na(r3) & r3 %in% c(1, 2) ~ "Hispanic",
        !is.na(r3) & r3 == 3 ~ "Non-Hispanic White",
        !is.na(r3) & r3 == 4 ~ "Non-Hispanic Black",
        !is.na(r3) & r3 == 6 ~ "Non-Hispanic Asian",
        !is.na(r3) & r3 %in% c(5, 7) ~ "Other/Multi",
        is.na(r3) & r1 %in% c(1, 2) ~ "Hispanic",
        is.na(r3) & r1 == 3 ~ "Non-Hispanic White",
        is.na(r3) & r1 == 4 ~ "Non-Hispanic Black",
        is.na(r3) & r1 == 5 ~ "Other/Multi",
        str_detect(r3_chr, "NON-HISPANIC WHITE") ~ "Non-Hispanic White",
        str_detect(r3_chr, "NON-HISPANIC BLACK") ~ "Non-Hispanic Black",
        str_detect(r3_chr, "NON-HISPANIC ASIAN") ~ "Non-Hispanic Asian",
        str_detect(r3_chr, "OTHER RACE|MULTI") ~ "Other/Multi",
        str_detect(r3_chr, "^MEXICAN AMERICAN$|^OTHER HISPANIC$|^HISPANIC$") ~ "Hispanic",
        str_detect(r1_chr, "NON-HISPANIC WHITE") ~ "Non-Hispanic White",
        str_detect(r1_chr, "NON-HISPANIC BLACK") ~ "Non-Hispanic Black",
        str_detect(r1_chr, "OTHER RACE|MULTI") ~ "Other/Multi",
        str_detect(r1_chr, "^MEXICAN AMERICAN$|^OTHER HISPANIC$|^HISPANIC$") ~ "Hispanic",
        TRUE ~ NA_character_
      )
    },
    education = {
      e_num <- suppressWarnings(as.numeric(as.character(DMDEDUC2)))
      e_chr <- toupper(trimws(as.character(DMDEDUC2)))
      case_when(
        e_num %in% c(1, 2) ~ "<HS",
        e_num == 3 ~ "HS/GED",
        e_num %in% c(4, 5) ~ ">HS",
        str_detect(e_chr, "LESS THAN 9TH|9TH.*12TH") ~ "<HS",
        str_detect(e_chr, "HIGH SCHOOL|GED") ~ "HS/GED",
        str_detect(e_chr, "SOME COLLEGE|COLLEGE GRADUATE") ~ ">HS",
        TRUE ~ NA_character_
      )
    },
    pir = suppressWarnings(as.numeric(as.character(INDFMPIR))),
    bmi = suppressWarnings(as.numeric(as.character(BMXBMI))),
    # ---- Categorical-variable derivations ---------------------------------
    # Smoking: NHANES uses a 2-stage instrument. SMQ020 ("Ever smoked >=100 cigarettes
    # in life?") is the GATE question; SMQ040 ("Do you now smoke?") is asked only of
    # SMQ020==Yes respondents. Therefore an SMQ040-only derivation silently drops every
    # never-smoker (SMQ020==No) to NA and excludes them from any complete-case modeling
    # sample. Below we derive the standard 3-level smoker3 (Never/Former/Current) using
    # SMQ020 + SMQ040, and a true binary smoker (Current vs. Not-current) that retains
    # never-smokers as 0. This binary is what downstream models continue to use.
    smoker3 = {
      ever <- parse_smq020(SMQ020)
      smq040_chr <- toupper(trimws(as.character(SMQ040)))
      smq040_num <- suppressWarnings(as.numeric(as.character(SMQ040)))
      curr <- case_when(
        smq040_num %in% c(1, 2) ~ 1,
        smq040_num == 3 ~ 0,
        str_detect(smq040_chr, "EVERY DAY|SOME DAYS") ~ 1,
        str_detect(smq040_chr, "^NOT AT ALL$") ~ 0,
        TRUE ~ NA_real_
      )
      case_when(
        ever == 0 ~ "Never",
        ever == 1 & curr == 1 ~ "Current",
        ever == 1 & curr == 0 ~ "Former",
        TRUE ~ NA_character_
      )
    },
    smoker = case_when(
      smoker3 == "Current" ~ 1,
      smoker3 %in% c("Never", "Former") ~ 0,
      TRUE ~ NA_real_
    ),
    # Alcohol: harmonized 3-level Never/Former/Current drinker, parallel to smoker3.
    # The previous derivation used ALQ151 ("ever heavy drinker") only, which (a) silently
    # dropped never-drinkers (skip pattern), (b) returned NA for the entire 2017-2018 (J)
    # cycle (variable discontinued), and (c) measured a different construct than the
    # NHANES-literature standard. The new derivation uses the cycle-harmonized helpers:
    #   parse_alq_ever()        : lifetime ever-drinker from ALQ101 / ALQ110 / ALQ111
    #   parse_alq_current_yr()  : drank in past 12 mo from ALQ121 (2013+) or ALQ120Q (1999-2016)
    # Never-drinkers are now coded as "Never" (instead of NA), and J-cycle participants
    # get a defined value via ALQ111 + ALQ121.
    alcohol3 = {
      ever <- parse_alq_ever(ALQ101, ALQ110, ALQ111)
      curr <- parse_alq_current_yr(ALQ121, ALQ120Q)
      case_when(
        ever == 0 ~ "Never",
        ever == 1 & curr == 1 ~ "Current",
        ever == 1 & curr == 0 ~ "Former",
        TRUE ~ NA_character_
      )
    },
    # Binary alcohol covariate (used by downstream models): "Current drinker" vs.
    # "Not current drinker" (combines never + former). This matches the NHANES literature
    # convention and is far less sparse than the previous ALQ151-based binary.
    alcohol = case_when(
      alcohol3 == "Current" ~ 1,
      alcohol3 %in% c("Never", "Former") ~ 0,
      TRUE ~ NA_real_
    ),
    # Physical activity: harmonized across the 2006/2007 GPAQ transition.
    #   PAQ650 (vigorous recreational, 2007+) OR PAQ665 (moderate recreational, 2007+)
    #   PAD200 (vigorous past 30 days, 1999-2006) OR PAD320 (moderate past 30 days, 1999-2006)
    # Either questionnaire format counts as "physically active" for the harmonized binary;
    # otherwise the participant is asked one set or the other depending on cycle.
    phys_active = case_when(
      parse_yes_no(PAQ650) == 1 | parse_yes_no(PAQ665) == 1 ~ 1,
      parse_yes_no(PAQ650) == 0 & parse_yes_no(PAQ665) == 0 ~ 0,
      parse_yes_no(PAD200) == 1 | parse_yes_no(PAD320) == 1 ~ 1,
      parse_yes_no(PAD200) == 0 & parse_yes_no(PAD320) == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    # Diabetes: DIQ010 ("Doctor told you have diabetes?") codes 1=Yes, 2=No, 3=Borderline,
    # 7/9=Refused/DK. yesno01() maps 1->1, 2->0, and ALL others (including Borderline) to NA.
    # This drops borderline-diabetics from the analytic sample. If PI prefers, borderline
    # can be coded as 1 (any diabetes) or 0 (definite diabetes only).
    diabetes = yesno01(DIQ010),
    # Hypertension: BPQ020 ("Doctor told you have HBP?") asked of all adults; no gate.
    hypertension = yesno01(BPQ020),
    # CVD composite: ANY of CHF (MCQ160B), CHD (MCQ160C), heart attack (MCQ160E), stroke
    # (MCQ160F). Not asked behind a gate; respondent NA in all four -> overall NA.
    cvd = case_when(
      parse_yes_no(MCQ160B) == 1 | parse_yes_no(MCQ160C) == 1 | parse_yes_no(MCQ160E) == 1 | parse_yes_no(MCQ160F) == 1 ~ 1,
      parse_yes_no(MCQ160B) == 0 & parse_yes_no(MCQ160C) == 0 & parse_yes_no(MCQ160E) == 0 & parse_yes_no(MCQ160F) == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    phenoage_accel = NA_real_
  )

# Compute biomarkers, CKD markers, and design variables outside mutate (avoids cur_data_all()).
dat$alb <- coalesce_num(dat, c("LBXSAL"))
dat$scr <- coalesce_num(dat, c("LBXSCR", "LBDSCRSI"))
# Glucose: glu_mec / glu_glu / glu set in download_cycle (primary uses glu_mec by design).
if (!"glu" %in% names(dat)) dat$glu <- dat$glu_mec
dat$crp <- coalesce_num(dat, c("LBXCRP", "LBXHSCRP"))
dat$lym <- coalesce_num(dat, c("LBXLYPCT"))
dat$mcv <- coalesce_num(dat, c("LBXMCV", "LBXMCVSI"))
dat$rdw <- coalesce_num(dat, c("LBXRW", "LBXRDW"))
dat$alp <- coalesce_num(dat, c("LBXAPSI", "LBXSAPSI"))
dat$wbc <- coalesce_num(dat, c("LBXWBCSI", "LBXWBC"))

dat$ualb <- coalesce_num(dat, c("URXUMA"))
dat$ucr <- coalesce_num(dat, c("URXUCR", "URDUCR"))
dat$acr_direct <- coalesce_num(dat, c("URDACT"))
dat$acr <- compute_acr(dat$ualb, dat$ucr, dat$acr_direct)
dat$egfr <- compute_egfr_2021(dat$scr, dat$age, dat$RIAGENDR)
dat$ckd_egfr <- ifelse(!is.na(dat$egfr), as.integer(dat$egfr < 60), NA_integer_)
dat$ckd_acr <- ifelse(!is.na(dat$acr), as.integer(dat$acr >= 30), NA_integer_)
dat$ckd <- case_when(
  config$ckd_def == "egfr_or_acr" & (dat$ckd_egfr == 1 | dat$ckd_acr == 1) ~ 1L,
  config$ckd_def == "egfr_or_acr" & (!is.na(dat$ckd_egfr) | !is.na(dat$ckd_acr)) ~ 0L,
  config$ckd_def == "egfr_only" ~ dat$ckd_egfr,
  config$ckd_def == "acr_only" ~ dat$ckd_acr,
  TRUE ~ NA_integer_
)
dat$ckd_stage <- cut(
  dat$egfr, breaks = c(-Inf, 15, 30, 45, 60, 90, Inf),
  labels = c("G5", "G4", "G3b", "G3a", "G2", "G1"),
  right = FALSE  # KDIGO: left-closed intervals [lo, hi). E.g., 60 -> G2, 59.9 -> G3a.
)
dat$strata <- coalesce_num(dat, c("SDMVSTRA"))
dat$psu <- coalesce_num(dat, c("SDMVPSU"))

# Coding audits (raw -> derived) to prevent silent recode failures.
kable(
  dat %>% count(cycle, ridreth3_missing = is.na(RIDRETH3)) %>% arrange(cycle, ridreth3_missing),
  caption = "Race audit: availability of RIDRETH3 by cycle"
)
kable(as.data.frame(table(dat$RIDRETH1, useNA = "ifany")), caption = "Raw RIDRETH1 distribution")
kable(as.data.frame(table(dat$RIDRETH3, useNA = "ifany")), caption = "Raw RIDRETH3 distribution")
kable(as.data.frame(table(dat$race_eth, useNA = "ifany")), caption = "Derived race_eth distribution")

kable(as.data.frame(table(dat$SMQ040, useNA = "ifany")), caption = "Raw SMQ040 distribution")
kable(as.data.frame(table(dat$DIQ010, useNA = "ifany")), caption = "Raw DIQ010 distribution")
kable(as.data.frame(table(dat$BPQ020, useNA = "ifany")), caption = "Raw BPQ020 distribution")
kable(as.data.frame(table(dat$smoker, useNA = "ifany")), caption = "Derived smoker distribution")
kable(as.data.frame(table(dat$diabetes, useNA = "ifany")), caption = "Derived diabetes distribution")
kable(as.data.frame(table(dat$hypertension, useNA = "ifany")), caption = "Derived hypertension distribution")

# Codebook-aligned harmonization checks (raw coding support by cycle).
check_raw_levels <- function(df, var, allowed_patterns) {
  if (!(var %in% names(df))) return(tibble(variable = var, cycle = NA_character_, n_unexpected = NA_integer_, unexpected_examples = NA_character_))
  df %>%
    mutate(.raw = toupper(trimws(as.character(.data[[var]])))) %>%
    group_by(cycle) %>%
    summarise(
      n_nonmissing = sum(!is.na(.raw) & .raw != ""),
      n_unexpected = sum(!is.na(.raw) & .raw != "" & !str_detect(.raw, allowed_patterns)),
      unexpected_examples = paste(head(unique(.raw[!is.na(.raw) & .raw != "" & !str_detect(.raw, allowed_patterns)]), 6), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(variable = var) %>%
    select(variable, cycle, n_nonmissing, n_unexpected, unexpected_examples)
}

harmonization_checks <- bind_rows(
  check_raw_levels(dat, "RIAGENDR", "^(1|2|MALE|FEMALE)$"),
  check_raw_levels(dat, "DMDEDUC2", "^(1|2|3|4|5|7|9)$|LESS THAN 9TH|9-11TH|HIGH SCHOOL|GED|SOME COLLEGE|COLLEGE GRADUATE|DON'T KNOW|DONT KNOW|REFUSED"),
  check_raw_levels(dat, "SMQ040", "^(1|2|3|7|9)$|EVERY DAY|SOME DAYS|NOT AT ALL|REFUSED|DON'T KNOW|DONT KNOW"),
  check_raw_levels(dat, "DIQ010", "^(1|2|3|7|9)$|YES|NO|BORDERLINE|REFUSED|DON'T KNOW|DONT KNOW"),
  check_raw_levels(dat, "BPQ020", "^(1|2|7|9)$|YES|NO|REFUSED|DON'T KNOW|DONT KNOW")
)
kable(harmonization_checks, caption = "Codebook-aligned raw coding checks by cycle. Non-zero unexpected counts indicate coding values not covered by harmonization rules.")

# --- CRP unit detection and conversion ---
# PhenoAge (Levine 2018) uses CRP in mg/dL. NHANES LBXCRP is typically mg/dL;
# LBXHSCRP or CRP tables in some cycles may report mg/L.
# Step 1: check variable/table provenance.
crp_by_cycle <- dat %>%
  group_by(cycle) %>%
  summarise(
    source_table = if ("crp_table" %in% names(pick(everything()))) paste(unique(na.omit(crp_table)), collapse = ",") else "unknown",
    n_obs = sum(!is.na(crp)),
    median_crp = median(crp, na.rm = TRUE),
    p75 = quantile(crp, 0.75, na.rm = TRUE),
    p90 = quantile(crp, 0.90, na.rm = TRUE),
    .groups = "drop"
  )
kable(crp_by_cycle, caption = "CRP diagnostic by cycle -- verify units (PhenoAge expects mg/dL). If source table is HSCRP, values may be mg/L.")

# Step 2: provenance-based detection (strongest signal).
crp_tables_used <- if ("crp_table" %in% names(dat)) unique(na.omit(dat$crp_table)) else character(0)
provenance_suggests_mgl <- length(crp_tables_used) > 0 && any(str_detect(crp_tables_used, regex("HSCRP", ignore_case = TRUE)))
# Step 3: distribution-based backup -- require median > 1.5 AND p75 > 3 to trigger conversion.
crp_overall_med <- median(dat$crp, na.rm = TRUE)
crp_overall_p75 <- quantile(dat$crp, 0.75, na.rm = TRUE)
distribution_suggests_mgl <- !is.na(crp_overall_med) && crp_overall_med > 1.5 &&
                             !is.na(crp_overall_p75) && crp_overall_p75 > 3

if (provenance_suggests_mgl || distribution_suggests_mgl) {
  reason <- if (provenance_suggests_mgl) {
    paste0("CRP source table(s) [", paste(crp_tables_used, collapse = ", "), "] suggest mg/L.")
  } else {
    paste0("Distribution heuristic: median=", round(crp_overall_med, 3), ", p75=", round(crp_overall_p75, 3), ".")
  }
  message(">>> Converting CRP from mg/L to mg/dL (/ 10). Reason: ", reason)
  dat$crp <- dat$crp / 10
} else {
  message(">>> CRP assumed to be in mg/dL (median=", round(crp_overall_med, 4),
          ", p75=", round(crp_overall_p75, 4), "). Verify against NHANES documentation.")
}

# --- Glucose provenance diagnostic ---
# Coerce glu to numeric (nhanes_raw has character columns from bind_rows); needed for summarise and later unit conversion.
dat$glu <- suppressWarnings(as.numeric(as.character(dat$glu)))
glu_by_cycle <- dat %>%
  group_by(cycle) %>%
  summarise(
    source_table = if ("glu_table" %in% names(pick(everything()))) paste(unique(na.omit(glu_table)), collapse = ",") else "unknown",
    n_glu = sum(!is.na(glu)),
    median_glu = median(glu, na.rm = TRUE),
    min_glu = min(glu, na.rm = TRUE),
    max_glu = max(glu, na.rm = TRUE),
    .groups = "drop"
  )
kable(glu_by_cycle, caption = "Glucose provenance by cycle. If source_table is a fasting component (GLU_*), use a compatible fasting subsample weight.")

glu_med <- median(dat$glu, na.rm = TRUE)
if (!is.na(glu_med) && glu_med < 6) {
  warning("Glucose median is < 6 -- this looks like mmol/L. Converting to mg/dL ( x 18.018).")
  dat$glu <- dat$glu * 18.018
}

# Warn based on stored table provenance, not just weight variable existence.
glu_tables_used <- if ("glu_table" %in% names(dat)) unique(na.omit(dat$glu_table)) else character(0)
glu_is_fasting <- length(glu_tables_used) > 0 && any(str_detect(glu_tables_used, regex("^GLU_", ignore_case = TRUE)))
if (glu_is_fasting) {
  message(">>> IMPORTANT: Glucose sourced from fasting component table(s) [",
          paste(glu_tables_used, collapse = ", "),
          "]. A compatible fasting subsample weight (e.g., WTSAF2YR) is required. ",
          "Using a non-fasting or otherwise incompatible weight can bias estimates. ",
          "Review the weight audit; override config$weight_type only if auto-selection does not resolve the correct subsample weight.")
} else if (length(glu_tables_used) > 0) {
  message(">>> Glucose source table(s): ", paste(glu_tables_used, collapse = ", "),
          ". Verify whether this is a fasting component.")
}

# ------------------------------------------------------------
# CLEAN SPECIAL NHANES MISSING CODES FOR LAB VARIABLES
# ------------------------------------------------------------
clean_nhanes_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  special_codes <- c(
    7, 77, 777, 7777, 77777,
    9, 99, 999, 9999, 99999,
    888, 8888, 88888
  )
  x[x %in% special_codes] <- NA_real_
  x
}

dat$alb <- clean_nhanes_numeric(dat$alb)
dat$scr <- clean_nhanes_numeric(dat$scr)
dat$glu <- clean_nhanes_numeric(dat$glu)
dat$crp <- clean_nhanes_numeric(dat$crp)
dat$lym <- clean_nhanes_numeric(dat$lym)
dat$mcv <- clean_nhanes_numeric(dat$mcv)
dat$rdw <- clean_nhanes_numeric(dat$rdw)
dat$alp <- clean_nhanes_numeric(dat$alp)
dat$wbc <- clean_nhanes_numeric(dat$wbc)

# ------------------------------------------------------------
# HARD PLAUSIBILITY CHECK (retain CRP >1 mg/dL in primary; flag as acute)
# ------------------------------------------------------------
dat <- dat %>%
  mutate(
    crp_acute = ifelse(!is.na(crp) & crp > 1, 1L, 0L),
    wbc = ifelse(wbc > 50, NA_real_, wbc),
    glu = ifelse(glu > 600, NA_real_, glu),
    alp = ifelse(alp > 500, NA_real_, alp)
  )

# --- General biomarker plausibility ---
biomarker_ranges <- list(
  alb  = c(1, 7,    "g/dL"),
  scr  = c(0.1, 25, "mg/dL"),
  glu  = c(20, 600,  "mg/dL"),
  crp  = c(0.001, 3, "mg/dL (post-conversion if applied)"),
  lym  = c(1, 80,    "%"),
  mcv  = c(50, 150,  "fL"),
  rdw  = c(8, 30,    "%"),
  alp  = c(5, 500,   "U/L"),
  wbc  = c(1, 50,    "1000 cells/uL")
)
biomarker_check <- lapply(names(biomarker_ranges), function(nm) {
  vals <- dat[[nm]]
  rng  <- biomarker_ranges[[nm]]
  lo   <- as.numeric(rng[1]); hi <- as.numeric(rng[2]); unit <- rng[3]
  v    <- vals[!is.na(vals)]
  n_out <- sum(v < lo | v > hi)
  tibble(biomarker = nm, unit = unit, n_obs = length(v),
         min = if (length(v) > 0) min(v) else NA_real_,
         p25 = if (length(v) > 0) quantile(v, 0.25) else NA_real_,
         median = if (length(v) > 0) median(v) else NA_real_,
         p75 = if (length(v) > 0) quantile(v, 0.75) else NA_real_,
         max = if (length(v) > 0) max(v) else NA_real_,
         expected_lo = lo, expected_hi = hi,
         n_out_of_range = n_out)
}) %>% bind_rows()
kable(biomarker_check, caption = "PhenoAge biomarker ranges after unit conversion (verify units)")

for (nm in names(biomarker_ranges)) {
  v <- dat[[nm]][!is.na(dat[[nm]])]
  lo <- as.numeric(biomarker_ranges[[nm]][1])
  hi <- as.numeric(biomarker_ranges[[nm]][2])
  pct_out <- if (length(v) > 0) mean(v < lo | v > hi) else 0
  if (pct_out > 0.20) warning(sprintf(
    "Biomarker '%s': %.1f%% of values outside plausible range [%s, %s]. Check units.",
    nm, pct_out * 100, lo, hi
  ))
}

# Recompute PhenoAge after any unit conversions (CRP mg/L->mg/dL, glucose mmol/L->mg/dL).
pheno_core <- compute_phenoage_components(dat, include_scr = TRUE, include_alb = TRUE)
dat$phenoage <- pheno_core$phenoage

message(
  "PhenoAge diagnostics -- xb median=", round(stats::median(pheno_core$xb, na.rm = TRUE), 3),
  ", mort median=", format(stats::median(pheno_core$mort, na.rm = TRUE), scientific = TRUE, digits = 3),
  ", pct(mort>0.9999)=", round(mean(pheno_core$mort > 0.9999, na.rm = TRUE) * 100, 1), "%",
  ", PhenoAge SD=", round(stats::sd(dat$phenoage, na.rm = TRUE), 3)
)
ph_q <- stats::quantile(dat$phenoage, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
message(
  "PhenoAge five-number summary -- min=", round(ph_q[[1]], 3),
  ", p25=", round(ph_q[[2]], 3),
  ", median=", round(ph_q[[3]], 3),
  ", p75=", round(ph_q[[4]], 3),
  ", max=", round(ph_q[[5]], 3)
)
if (mean(pheno_core$mort > 0.9999, na.rm = TRUE) > 0.25) {
  warning("More than 25% of participants have mort > 0.9999 after stabilization. Re-check biomarker units, especially CRP and glucose.")
}



## ----exposure-processing------------------------------------------------------
pfas_vars <- c("LBXPFOS", "LBXPFOA", "LBXPFNA", "LBXPFHS", "LBXPFDA")

# --- ISOMER HARMONIZATION (2013-2018) ---
# NHANES split PFOS and PFOA into linear and branched isomers.
# We must apply the LOD/sqrt(2) rule to the individual isomers FIRST, then sum them.

isomer_map <- list(
  pfos_lin = c("SSNPFOS", "LBXNFOS"),
  pfos_br  = c("SSMPFOS", "LBXMFOS"),
  pfoa_lin = c("SSNPFOA", "LBXNFOA"),
  pfoa_br  = c("SSBPFOA", "LBXBFOA")
)

lc_map <- list(
  pfos_lin = c("SDNPFOSL", "LBDNFOSL"),
  pfos_br  = c("SDMPFOSL", "LBDMFOSL"),
  pfoa_lin = c("SDNPFOAL", "LBDNFOAL"),
  pfoa_br  = c("SDBPFOAL", "LBDBFOAL")
)

for (chem in c("pfos", "pfoa")) {
  # Get available variables for the cycle
  lin_vars <- isomer_map[[paste0(chem, "_lin")]]
  br_vars  <- isomer_map[[paste0(chem, "_br")]]
  lin_lc   <- lc_map[[paste0(chem, "_lin")]]
  br_lc    <- lc_map[[paste0(chem, "_br")]]
  
  lin_val <- coalesce_num(dat, lin_vars)
  br_val  <- coalesce_num(dat, br_vars)
  lin_lc_val <- coalesce_num(dat, lin_lc)
  br_lc_val  <- coalesce_num(dat, br_lc)
  
  # Apply LOD rule: if indicator == 1, divide by sqrt(2)
  lin_val <- ifelse(!is.na(lin_lc_val) & lin_lc_val == 1, lin_val / sqrt(2), lin_val)
  br_val  <- ifelse(!is.na(br_lc_val) & br_lc_val == 1, br_val / sqrt(2), br_val)
  
  # Sum them to recreate the total concentration
  summed_val <- lin_val + br_val
  
  # Inject back into the legacy total variables (LBXPFOS or LBXPFOA)
  target_var <- ifelse(chem == "pfos", "LBXPFOS", "LBXPFOA")
  
  if (target_var %in% names(dat)) {
    dat[[target_var]] <- coalesce_num(dat, c(target_var))
    dat[[target_var]] <- ifelse(is.na(dat[[target_var]]) & !is.na(summed_val), summed_val, dat[[target_var]])
  } else {
    dat[[target_var]] <- summed_val
  }
}
# --- END ISOMER HARMONIZATION ---

blood_metal_vars <- c("LBXBPB", "LBXBCD", "LBXTHG")
urine_metal_vars <- c("URXUAS", "URXUCD", "URXUPB")

chosen_exposures <- switch(
  config$exposure_family,
  pfas = pfas_vars,
  metals = c(blood_metal_vars, urine_metal_vars),
  both = c(pfas_vars, blood_metal_vars, urine_metal_vars),
  pfas_vars
)
chosen_exposures <- chosen_exposures[chosen_exposures %in% names(dat)]
if (length(chosen_exposures) < 3) stop("Too few exposure variables found; choose different exposure family/cycles.")

for (v in chosen_exposures) {
  vv <- apply_lod_rule(dat, v)
  dat[[paste0(v, "_log2")]] <- log2(vv)
}

if (any(urine_metal_vars %in% chosen_exposures)) {
  for (uv in intersect(urine_metal_vars, chosen_exposures)) {
    vv <- apply_lod_rule(dat, uv)
    dat[[paste0(uv, "_ucr_adj_log2")]] <- log2((vv / pmax(dat$ucr, 1e-6)) * 100)
  }
}

exposure_vars <- names(dat)[str_detect(names(dat), "_log2$|_ucr_adj_log2$")]
exposure_vars <- unique(exposure_vars)
# Prevent collinear double-counting: use creatinine-adjusted urine metals only, drop raw _log2 for them.
raw_urine_vars <- paste0(urine_metal_vars, "_log2")
exposure_vars <- setdiff(exposure_vars, raw_urine_vars)
# Primary analysis: fixed 4-PFAS subset common across all retained cycles. Sensitivity: available-exposure (adaptive) mixture run as secondary.
primary_exposure_vars <- c("LBXPFOS_log2", "LBXPFOA_log2", "LBXPFNA_log2", "LBXPFHS_log2")
primary_exposure_vars <- primary_exposure_vars[primary_exposure_vars %in% exposure_vars]
adaptive_exposure_vars <- exposure_vars

# LOD diagnostic table -- computed on dat (post-merge, pre-log transform) stratified by cycle.
lod_diag <- lapply(chosen_exposures, function(v) {
  x_raw <- suppressWarnings(as.numeric(as.character(dat[[v]])))
  base <- sub("^(LBX|URX)", "", v)
  lc_candidates <- c(
    paste0(v, "LC"), paste0("LBD", base, "LC"), paste0("URD", base, "LC"),
    paste0(v, "L"),  paste0("LBD", base, "L"),  paste0("URD", base, "L")
  )
  lc_var <- lc_candidates[lc_candidates %in% names(dat)]
  has_lc <- length(lc_var) > 0
  if (has_lc) {
    lc_vals <- suppressWarnings(as.numeric(as.character(dat[[lc_var[1]]])))
    below <- !is.na(lc_vals) & lc_vals == 1
  } else {
    below <- rep(FALSE, length(x_raw))
  }

  by_cycle <- dat %>%
    mutate(.x = x_raw, .below = below) %>%
    group_by(cycle) %>%
    summarise(
      n = sum(!is.na(.x)),
      n_below = sum(.below, na.rm = TRUE),
      pct_below = if (sum(!is.na(.x)) > 0) round(sum(.below, na.rm = TRUE) / sum(!is.na(.x)) * 100, 1) else NA_real_,
      median_below = if (sum(.below & !is.na(.x)) > 0) median(.x[.below], na.rm = TRUE) else NA_real_,
      median_above = if (sum(!.below & !is.na(.x)) > 0) median(.x[!.below & !is.na(.x)], na.rm = TRUE) else NA_real_,
      min_positive = if (any(.x > 0, na.rm = TRUE)) min(.x[.x > 0], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(exposure = v, lc_indicator = if (has_lc) lc_var[1] else "none") %>%
    select(exposure, cycle, lc_indicator, everything())
  by_cycle
}) %>% bind_rows()
kable(lod_diag, caption = "LOD diagnostic per exposure by cycle. Verify: median_below ~ LOD value; if it diverges from min_positive, review the LC->x/sqrt(2) assumption.")


## ----analytic-sample----------------------------------------------------------
analytic <- dat %>%
  filter(age >= config$min_age)

# PFAS workflow: harmonize rotating 2-year subsample weights into one unified PFAS weight (WTPFAS2YR).
analytic <- analytic %>%
  mutate(
    WTPFAS2YR = coalesce_num(cur_data(), c("WTSA2YR", "WTSB2YR", "WTSC2YR", "WTSAF2YR"))
  )

weight_candidates <- names(analytic)[stringr::str_detect(names(analytic), "^WT")]
kable(data.frame(weight_variable = weight_candidates), caption = "Weight audit: all WT* variables in analytic dataset.")

# Exposure-specific weight audit in the exposure-complete analytic sample (fixed primary PFAS set).
complete_expo <- analytic %>%
  mutate(.has_all_expo = if_all(all_of(primary_exposure_vars), ~ !is.na(.x)))
wt_cols <- names(analytic)[str_detect(names(analytic), "^WT")]
# PFAS auto-selection eligibility: exclude full-sample interview/MEC weights; keep PFAS-appropriate subsample candidates.
if (tolower(config$exposure_family) == "pfas") {
  wt_cols <- setdiff(wt_cols, c("WTINT2YR", "WTMEC2YR"))
}
if (length(wt_cols) > 0) {
  wt_audit <- complete_expo %>%
    filter(.has_all_expo) %>%
    summarise(across(all_of(wt_cols), ~ round(mean(!is.na(.x)) * 100, 1))) %>%
    pivot_longer(everything(), names_to = "weight", values_to = "pct_nonmissing_complete_expo") %>%
    arrange(desc(pct_nonmissing_complete_expo))
  kable(wt_audit, caption = "Weight availability (%) among participants with complete exposure set. Use the weight with highest coverage that matches your subsample.")
}

# Auto-weight: rank eligible candidate weights by coverage and apply the top PFAS-appropriate subsample weight.
if (config$weight_type == "auto" && exists("wt_audit") && nrow(wt_audit) > 0) {
  wt_audit_auto <- wt_audit
  if (tolower(config$exposure_family) == "pfas") {
    wt_audit_auto <- wt_audit_auto %>% filter(!(weight %in% c("WTINT2YR", "WTMEC2YR")))
  }
  if (nrow(wt_audit_auto) == 0) {
    stop("Auto-weight selection failed for PFAS: no eligible WTS* subsample weights were available in the audit table.")
  }
  top_row <- wt_audit_auto %>% slice(1)
  config$weight_type <- top_row$weight
  message(">>> Auto-weight resolved: Applying '", config$weight_type, "' (Coverage: ", top_row$pct_nonmissing_complete_expo, "%)")
}

weight_info <- resolve_weight_var(analytic, config$weight_type, config$exposure_family)

analytic$weight_2yr <- suppressWarnings(as.numeric(as.character(analytic[[weight_info$var]])))

analytic <- analytic %>%
  filter(!is.na(phenoage), !is.na(age), !is.na(ckd), !is.na(weight_2yr), !is.na(strata), !is.na(psu))

n_cycles <- length(unique(analytic$cycle))
# After resolving the PFAS-appropriate 2-year weight, pooled analyses divide by the number of contributing cycles in the final sample (NCHS multi-cycle guidance).
n_cycles_exposure_available <- analytic %>%
  filter(if_any(all_of(primary_exposure_vars), ~ !is.na(.x))) %>%
  summarise(n = n_distinct(cycle)) %>%
  pull(n)
if (!is.finite(n_cycles_exposure_available) || n_cycles_exposure_available < 1) {
  n_cycles_exposure_available <- n_cycles
}
repool_weights <- function(df) {
  ncy <- dplyr::n_distinct(df$cycle[!is.na(df$weight_2yr)])
  if (!is.finite(ncy) || ncy < 1) stop("Invalid cycle count for weight pooling.")
  df$weight_combined <- df$weight_2yr / ncy
  df
}
analytic <- repool_weights(analytic)
message("Weight pooling: repool_weights() applied to analytic (n_cycles in sample = ", n_cycles_exposure_available, ")")

# Audit B: weight availability in the final analytic sample (after all exclusions).
if (length(wt_cols) > 0) {
  wt_audit_final <- analytic %>%
    summarise(across(any_of(wt_cols), ~ round(mean(!is.na(.x)) * 100, 1))) %>%
    pivot_longer(everything(), names_to = "weight", values_to = "pct_nonmissing_final_sample") %>%
    arrange(desc(pct_nonmissing_final_sample))
  kable(wt_audit_final, caption = "Weight availability (%) in final analytic sample (post-exclusion). This is the definitive audit.")
}

# PRIMARY: survey-weighted residualization of PhenoAge on age.
# This ensures PhenoAgeAccel is a population-representative construct.
analytic$phenoage_accel_unweighted <- NA_real_
fit_age_uw <- lm(phenoage ~ age, data = analytic)
analytic$phenoage_accel_unweighted <- residuals(fit_age_uw)
keep_w <- which(!is.na(analytic$phenoage) & !is.na(analytic$age) & !is.na(analytic$weight_combined))
if (length(keep_w) > 0) {
  svy_resid <- svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE,
                          data = analytic[keep_w, , drop = FALSE])
  fit_age_wt <- svyglm(phenoage ~ age, design = svy_resid, family = gaussian())
  pred_wt <- as.numeric(predict(fit_age_wt))
  analytic$phenoage_accel <- NA_real_
  analytic$phenoage_accel[keep_w] <- analytic$phenoage[keep_w] - pred_wt
}
analytic <- analytic %>% filter(!is.na(phenoage_accel))
message(">>> PhenoAgeAccel: primary = survey-weighted residualization; unweighted version retained as sensitivity.")

ckd_marker_availability <- tibble::tibble(
  group = c("Both eGFR and ACR", "eGFR only", "ACR only", "Neither"),
  N = c(
    sum(!is.na(analytic$egfr) & !is.na(analytic$acr)),
    sum(!is.na(analytic$egfr) & is.na(analytic$acr)),
    sum(is.na(analytic$egfr) & !is.na(analytic$acr)),
    sum(is.na(analytic$egfr) & is.na(analytic$acr))
  )
)

if (!all(c(0, 1) %in% unique(analytic$ckd))) {
  stop("Both CKD groups are not represented after exclusions.")
}

svy <- svydesign(
  id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE, data = analytic
)

kable(data.frame(weight_selected = weight_info$var, note = weight_info$note), caption = "Selected NHANES weight variable")
kable(ckd_marker_availability, caption = "CKD marker availability in analytic sample")


## ----mixture-models-----------------------------------------------------------
valid_num <- function(v, d) v %in% names(d) && length(unique(stats::na.omit(d[[v]]))) >= 2
valid_cat <- function(v, d) v %in% names(d) && length(unique(stats::na.omit(as.character(d[[v]])))) >= 2

# Primary adjustment: lean set to avoid overadjustment/mediator bias (diabetes, HTN, CVD may be on pathway).
# Aligned with EHP/Env Int mixture papers: sex, race, PIR, BMI, smoking, cycle; plus alcohol for confounding control.
# Smoking and alcohol are 3-level factors (Never/Former/Current) for consistency with Table 1.
clinical_covars <- c("diabetes", "hypertension", "cvd")  # used only in sensitivity "cardiometabolic"
# Age is retained because PFAS is age-structured; residualizing only the outcome does not age-adjust the exposure.
covars <- c(
  if (valid_num("age", analytic)) "age",
  if (valid_cat("sex", analytic)) "factor(sex)",
  if (valid_cat("race_eth", analytic)) "factor(race_eth)",
  if (valid_num("pir", analytic)) "pir",
  if (valid_cat("smoker3", analytic)) "factor(smoker3)",
  if (valid_num("bmi", analytic)) "bmi",
  if (valid_cat("alcohol3", analytic)) "factor(alcohol3)",
  if (valid_cat("cycle", analytic)) "factor(cycle)"
)

# Sanity checks before modeling.
message("--- Pre-modeling sanity checks ---")
message("Analytic N: ", nrow(analytic))
message("PhenoAgeAccel: mean=", round(mean(analytic$phenoage_accel, na.rm = TRUE), 3),
        ", SD=", round(sd(analytic$phenoage_accel, na.rm = TRUE), 3))
message("CKD prevalence: ", round(mean(analytic$ckd == 1, na.rm = TRUE) * 100, 1), "%")
message("Weight selected: ", weight_info$var)

# qgcomp complete-case: use primary (fixed) exposure vars with sufficient variability and non-missing in >= 50% of analytic.
qg_min_frac <- 0.5
qg_exposure_vars_for_cc <- primary_exposure_vars[vapply(primary_exposure_vars, function(v) {
  x <- suppressWarnings(as.numeric(analytic[[v]]))
  finite <- is.finite(x)
  n_ok <- sum(finite)
  pct_ok <- n_ok / max(1, nrow(analytic))
  n_ok > 0 && dplyr::n_distinct(x[finite]) >= 2 && pct_ok >= qg_min_frac
}, logical(1))]
if (length(qg_exposure_vars_for_cc) < 2) {
  qg_exposure_vars_for_cc <- primary_exposure_vars[vapply(primary_exposure_vars, function(v) {
    x <- suppressWarnings(as.numeric(analytic[[v]]))
    x <- x[is.finite(x)]
    length(x) > 0 && dplyr::n_distinct(x) >= 2
  }, logical(1))]
  if (length(qg_exposure_vars_for_cc) < 2) qg_exposure_vars_for_cc <- primary_exposure_vars
  message("qgcomp: few primary exposures met variance + ", round(qg_min_frac * 100), "% non-missing; using variance-only subset or all primary_exposure_vars.")
}
message("qgcomp exposures for complete-case (primary fixed set, variance + >= ", round(qg_min_frac * 100), "% non-missing): ", paste(qg_exposure_vars_for_cc, collapse = ", "))
# Qgcomp uses same adjustment as Model B (age, sex, race, PIR, BMI, smoking [3-level], alcohol [3-level], cycle).
rhs_qg_cc <- paste(c(qg_exposure_vars_for_cc, covars), collapse = " + ")
f_qg_cc <- as.formula(paste("phenoage_accel ~", rhs_qg_cc))
qg_cc <- complete.cases(model.frame(f_qg_cc, data = analytic, na.action = na.pass))
analytic_qg <- analytic[qg_cc, , drop = FALSE]
analytic_qg <- repool_weights(analytic_qg)
# Normalize precision weights for qgcomp so glm does not treat sum(weights) as N (avoids N = 300M, infinitesimal SE).
analytic_qg$weight_normalized <- analytic_qg$weight_combined / mean(analytic_qg$weight_combined, na.rm = TRUE)
# Ensure covariates are proper type for qgcomp (character -> wrong model matrix, psi ~ 0).
for (v in c("phenoage_accel", "age", "pir", "bmi")) {
  if (v %in% names(analytic_qg)) analytic_qg[[v]] <- suppressWarnings(as.numeric(as.character(analytic_qg[[v]])))
}
analytic_qg$sex      <- factor(as.character(analytic_qg$sex))
analytic_qg$race_eth <- factor(as.character(analytic_qg$race_eth))
analytic_qg$cycle    <- factor(as.character(analytic_qg$cycle))
if ("smoker3"  %in% names(analytic_qg)) analytic_qg$smoker3  <- factor(as.character(analytic_qg$smoker3),  levels = c("Never", "Former", "Current"))
if ("alcohol3" %in% names(analytic_qg)) analytic_qg$alcohol3 <- factor(as.character(analytic_qg$alcohol3), levels = c("Never", "Former", "Current"))
for (v in qg_exposure_vars_for_cc) {
  if (v %in% names(analytic_qg)) analytic_qg[[v]] <- suppressWarnings(as.numeric(as.character(analytic_qg[[v]])))
}
message("qgcomp complete-case: same covariates as Model B; N = ", nrow(analytic_qg))

# Among those, keep only exposures that still have variance in the complete-case sample (in case one collapses).
qg_exposure_vars <- qg_exposure_vars_for_cc[vapply(qg_exposure_vars_for_cc, function(v) {
  x <- suppressWarnings(as.numeric(analytic_qg[[v]]))
  x <- x[is.finite(x)]
  length(x) > 0 && dplyr::n_distinct(x) >= 2
}, logical(1))]
if (length(qg_exposure_vars) < 2) {
  warning("Too few finite/variable exposures for qgcomp after complete-case filtering; Table 2 will be NA.")
}
if (length(qg_exposure_vars) > 0) {
  analytic_qg <- analytic_qg %>%
    filter(if_all(all_of(c("phenoage_accel", qg_exposure_vars)), ~ is.finite(.x)))
}
q_qg <- as.integer(config$q)
if (!is.finite(q_qg) || q_qg < 2) stop("config$q must be an integer >= 2.")

# Qgcomp model: same adjustment as Model B (includes chronological age).
covars_qg <- c(
  if (valid_num("age", analytic_qg)) "age",
  if (valid_cat("sex", analytic_qg)) "factor(sex)",
  if (valid_cat("race_eth", analytic_qg)) "factor(race_eth)",
  if (valid_num("pir", analytic_qg)) "pir",
  if (valid_cat("smoker3", analytic_qg)) "factor(smoker3)",
  if (valid_num("bmi", analytic_qg)) "bmi",
  if (valid_cat("alcohol3", analytic_qg)) "factor(alcohol3)",
  if (valid_cat("cycle", analytic_qg)) "factor(cycle)"
)
qg_covars_intended <- covars_qg
build_qg_formula <- function(exps, covs) {
  rhs <- paste(c(exps, covs), collapse = " + ")
  as.formula(paste("phenoage_accel ~", rhs))
}
f_qg <- build_qg_formula(qg_exposure_vars, covars_qg)

repeat {
  mf_qg <- tryCatch(model.frame(f_qg, data = analytic_qg, na.action = na.omit), error = function(e) NULL)
  if (is.null(mf_qg) || nrow(mf_qg) == 0) break
  factor_terms_qg <- names(mf_qg)[vapply(mf_qg, is.factor, logical(1))]
  bad_factor_terms_qg <- factor_terms_qg[vapply(mf_qg[factor_terms_qg], function(x) nlevels(droplevels(x)) < 2, logical(1))]
  removable_qg <- intersect(covars_qg, bad_factor_terms_qg)
  if (length(removable_qg) == 0) break
  message("qgcomp overall: dropping collapsed factor term(s): ", paste(removable_qg, collapse = ", "))
  covars_qg <- setdiff(covars_qg, removable_qg)
  f_qg <- build_qg_formula(qg_exposure_vars, covars_qg)
}
qg_covars_realized <- covars_qg
message("qgcomp exposures used: ", paste(qg_exposure_vars, collapse = ", "))
message("qgcomp finite analytic N: ", nrow(analytic_qg))

qgcomp_overall_error <- NULL
fit_qg_overall <- tryCatch(
  qgcomp::qgcomp.noboot(
    f = f_qg, expnms = qg_exposure_vars, data = analytic_qg,
    family = gaussian(), q = q_qg, weights = analytic_qg$weight_normalized
  ),
  error = function(e) {
    qgcomp_overall_error <<- conditionMessage(e)
    NULL
  }
)
if (!is.null(fit_qg_overall)) {
  message("qgcomp N used: ", fit_qg_overall$fit$n[1])
} else {
  warning("qgcomp overall model failed: ", qgcomp_overall_error)
}

# Survey-weighted interaction model using quantized summed score.
# Quantile cutpoints are computed within each NHANES cycle using cycle-specific sampling weights.
qdat <- analytic
qdat$ckd <- as.integer(qdat$ckd)
# Strict sample = complete-case on qgcomp-compatible PFAS subset (same as analytic_qg). Label reflects this.
qdat_strict <- analytic_qg
qdat_strict$ckd <- as.integer(qdat_strict$ckd)

# Primary model: fixed PFAS set (primary_exposure_vars). Restrict to complete case on primary set.
qdat <- qdat %>% filter(if_all(all_of(primary_exposure_vars), ~ !is.na(.x)))
if (config$require_complete_exposure_set) {
  # Optional extra: could also require full adaptive set here; we use primary only for main.
}
message("Primary interaction: fixed PFAS set = ", paste(primary_exposure_vars, collapse = ", "))
message("N primary pre-quantile sample: ", nrow(qdat))
message("N strict sample (Strict complete-case, qgcomp-compatible PFAS subset): ", nrow(qdat_strict))

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]; w <- w[keep]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    if (is.na(idx)) x[length(x)] else x[idx]
  })
}

# Primary mixture score: fixed set, cycle-specific quantiles to avoid temporal confounding (PFAS trends across retained cycles).
for (v in primary_exposure_vars) {
  qdat[[paste0(v, "_q")]] <- NA_real_
  for (cy in unique(qdat$cycle)) {
    idx <- qdat$cycle == cy
    x_cy <- suppressWarnings(as.numeric(qdat[[v]][idx]))
    w_cy <- suppressWarnings(as.numeric(qdat$weight_combined[idx]))
    probs <- weighted_quantile(x_cy, w_cy, probs = seq(0, 1, length.out = config$q + 1))
    br <- sort(unique(as.numeric(probs)))
    if (length(br) > 1) {
      br[1] <- -Inf
      br[length(br)] <- Inf
    }
    if (length(br) < 3 || any(diff(br) <= 0)) {
      r <- rank(qdat[[v]][idx], na.last = "keep", ties.method = "average")
      qdat[[paste0(v, "_q")]][idx] <- as.numeric(cut(r, breaks = config$q, labels = FALSE, include.lowest = TRUE))
    } else {
      qdat[[paste0(v, "_q")]][idx] <- as.numeric(cut(qdat[[v]][idx], breaks = br, include.lowest = TRUE, labels = FALSE))
    }
  }
}
qvars <- paste0(primary_exposure_vars, "_q")
qdat$nonmissing_q_exposures <- rowSums(!is.na(qdat[, qvars, drop = FALSE]))
qdat$mixture_qscore <- rowMeans(qdat[, qvars, drop = FALSE], na.rm = TRUE)
# Primary: require all primary PFAS non-missing (fixed set).
required_primary <- length(primary_exposure_vars)
qdat$mixture_qscore <- ifelse(qdat$nonmissing_q_exposures < required_primary, NA_real_, qdat$mixture_qscore)
# By-cycle table for reporting (primary set is fixed; this documents it).
cycle_capacity <- tibble(
  cycle = unique(qdat$cycle),
  measured_exposures = length(primary_exposure_vars),
  PFAS_measured = paste(primary_exposure_vars, collapse = ", "),
  required_nonmissing_exposures = required_primary
)

# Build strict mixture score: same exposure set as qgcomp, cycle-specific quantiles.
for (v in qg_exposure_vars) {
  qdat_strict[[paste0(v, "_q")]] <- NA_real_
  for (cy in unique(qdat_strict$cycle)) {
    idx <- qdat_strict$cycle == cy
    x_cy <- suppressWarnings(as.numeric(qdat_strict[[v]][idx]))
    w_cy <- suppressWarnings(as.numeric(qdat_strict$weight_combined[idx]))
    probs_s <- weighted_quantile(x_cy, w_cy, probs = seq(0, 1, length.out = config$q + 1))
    br_s <- sort(unique(as.numeric(probs_s)))
    if (length(br_s) > 1) {
      br_s[1] <- -Inf
      br_s[length(br_s)] <- Inf
    }
    if (length(br_s) < 3 || any(diff(br_s) <= 0)) {
      r_s <- rank(qdat_strict[[v]][idx], na.last = "keep", ties.method = "average")
      qdat_strict[[paste0(v, "_q")]][idx] <- as.numeric(cut(r_s, breaks = config$q, labels = FALSE, include.lowest = TRUE))
    } else {
      qdat_strict[[paste0(v, "_q")]][idx] <- as.numeric(cut(qdat_strict[[v]][idx], breaks = br_s, include.lowest = TRUE, labels = FALSE))
    }
  }
}
qvars_strict <- paste0(qg_exposure_vars, "_q")
qdat_strict$nonmissing_q_exposures <- rowSums(!is.na(qdat_strict[, qvars_strict, drop = FALSE]))
qdat_strict$mixture_qscore <- rowMeans(qdat_strict[, qvars_strict, drop = FALSE], na.rm = TRUE)

# Hard type enforcement to prevent silent coercion issues in svyglm.
qdat <- qdat %>%
  mutate(
    phenoage_accel = suppressWarnings(as.numeric(phenoage_accel)),
    mixture_qscore = suppressWarnings(as.numeric(mixture_qscore)),
    weight_combined = suppressWarnings(as.numeric(weight_combined)),
    psu = suppressWarnings(as.numeric(psu)),
    strata = suppressWarnings(as.numeric(strata)),
    ckd = as.integer(ckd)
  )
qdat_strict <- qdat_strict %>%
  mutate(
    phenoage_accel = suppressWarnings(as.numeric(phenoage_accel)),
    mixture_qscore = suppressWarnings(as.numeric(mixture_qscore)),
    weight_combined = suppressWarnings(as.numeric(weight_combined)),
    psu = suppressWarnings(as.numeric(psu)),
    strata = suppressWarnings(as.numeric(strata)),
    ckd = as.integer(ckd)
  )

# Primary model: fixed PFAS set (same in every cycle).
cycle_capacity_with_min <- cycle_capacity
kable(
  cycle_capacity_with_min,
  caption = "Primary analysis: fixed PFAS set by cycle. Same chemicals and required non-missing in all cycles."
)

kable(
  tibble(
    statistic = c("N", "Mean", "SD", "Min", "P25", "Median", "P75", "Max"),
    value = c(
      sum(!is.na(qdat$mixture_qscore)),
      round(mean(qdat$mixture_qscore, na.rm = TRUE), 3),
      round(sd(qdat$mixture_qscore, na.rm = TRUE), 3),
      round(min(qdat$mixture_qscore, na.rm = TRUE), 3),
      round(quantile(qdat$mixture_qscore, 0.25, na.rm = TRUE), 3),
      round(median(qdat$mixture_qscore, na.rm = TRUE), 3),
      round(quantile(qdat$mixture_qscore, 0.75, na.rm = TRUE), 3),
      round(max(qdat$mixture_qscore, na.rm = TRUE), 3)
    )
  ),
  caption = paste0("Primary mixture quantile score (fixed set, N chemicals = ", length(primary_exposure_vars),
                   "). Score = mean quantile bin across ", paste(primary_exposure_vars, collapse = ", "),
                   ". Complete case on primary set required.")
)

stopifnot(all(c("psu", "strata", "weight_combined") %in% names(qdat)))
qdat <- qdat %>% filter(!is.na(psu), !is.na(strata), !is.na(weight_combined), !is.na(phenoage_accel), !is.na(ckd), !is.na(mixture_qscore))
qdat_strict <- qdat_strict %>% filter(!is.na(psu), !is.na(strata), !is.na(weight_combined), !is.na(phenoage_accel), !is.na(ckd), !is.na(mixture_qscore))
qdat         <- repool_weights(qdat)
qdat_strict  <- repool_weights(qdat_strict)

# Ensure all model variables are proper type so svyglm does not zero out coefficients (character -> wrong design matrix).
for (v in exposure_vars) {
  if (v %in% names(qdat)) qdat[[v]] <- suppressWarnings(as.numeric(as.character(qdat[[v]])))
}
qdat$phenoage_accel <- suppressWarnings(as.numeric(as.character(qdat$phenoage_accel)))
qdat$age            <- suppressWarnings(as.numeric(as.character(qdat$age)))
qdat$pir            <- suppressWarnings(as.numeric(as.character(qdat$pir)))
qdat$bmi            <- suppressWarnings(as.numeric(as.character(qdat$bmi)))
qdat$sex            <- factor(as.character(qdat$sex))
qdat$race_eth       <- factor(as.character(qdat$race_eth))
qdat$cycle          <- factor(as.character(qdat$cycle))
if ("smoker3"  %in% names(qdat)) qdat$smoker3  <- factor(as.character(qdat$smoker3),  levels = c("Never", "Former", "Current"))
if ("alcohol3" %in% names(qdat)) qdat$alcohol3 <- factor(as.character(qdat$alcohol3), levels = c("Never", "Former", "Current"))
qdat$ckd            <- as.integer(qdat$ckd)

for (v in exposure_vars) {
  if (v %in% names(qdat_strict)) qdat_strict[[v]] <- suppressWarnings(as.numeric(as.character(qdat_strict[[v]])))
}
qdat_strict$phenoage_accel <- suppressWarnings(as.numeric(as.character(qdat_strict$phenoage_accel)))
qdat_strict$age            <- suppressWarnings(as.numeric(as.character(qdat_strict$age)))
qdat_strict$pir            <- suppressWarnings(as.numeric(as.character(qdat_strict$pir)))
qdat_strict$bmi            <- suppressWarnings(as.numeric(as.character(qdat_strict$bmi)))
qdat_strict$sex            <- factor(as.character(qdat_strict$sex))
qdat_strict$race_eth       <- factor(as.character(qdat_strict$race_eth))
qdat_strict$cycle          <- factor(as.character(qdat_strict$cycle))
if ("smoker3"  %in% names(qdat_strict)) qdat_strict$smoker3  <- factor(as.character(qdat_strict$smoker3),  levels = c("Never", "Former", "Current"))
if ("alcohol3" %in% names(qdat_strict)) qdat_strict$alcohol3 <- factor(as.character(qdat_strict$alcohol3), levels = c("Never", "Former", "Current"))
qdat_strict$ckd            <- as.integer(qdat_strict$ckd)

message("---- Model B diagnostics ----")
message("N qdat: ", nrow(qdat))
print(table(qdat$ckd, useNA = "ifany"))
message("N qdat_strict: ", nrow(qdat_strict))
print(table(qdat_strict$ckd, useNA = "ifany"))
sd_all <- sd(qdat$mixture_qscore, na.rm = TRUE)
sd_by <- tapply(qdat$mixture_qscore, qdat$ckd, sd, na.rm = TRUE)
message("SD(mixture_qscore) overall: ", sd_all)
print(sd_by)

if (length(unique(stats::na.omit(qdat$ckd))) < 2) {
  stop("After exposure/design filtering, only one CKD group remains in qdat. ",
       "Interaction mixture_qscore*ckd cannot be estimated. ",
       "Set config$require_complete_exposure_set = FALSE or reduce exposure set/cycles.")
}
if (any(is.finite(sd_by) & sd_by < 1e-8)) {
  stop("mixture_qscore has near-zero variance within at least one CKD group. ",
       "Interaction cannot be estimated; check quantile binning/ties and exposure variability.")
}

miss_compare_ckd <- qdat %>%
  mutate(
    retained_primary = complete.cases(phenoage_accel, mixture_qscore, ckd, age, sex, race_eth, pir, smoker3, bmi, alcohol3, cycle)
  ) %>%
  group_by(ckd, retained_primary) %>%
  summarise(
    N = n(),
    mean_age = round(mean(age, na.rm = TRUE), 1),
    female_pct = round(mean(sex == "Female", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  group_by(ckd) %>%
  mutate(pct_within_ckd = round(100 * N / sum(N), 1)) %>%
  ungroup() %>%
  mutate(
    ckd = ifelse(ckd == 1, "CKD", "No CKD"),
    retained_primary = ifelse(retained_primary, "Retained in Model B", "Dropped from Model B")
  )
kable(miss_compare_ckd, caption = "Primary interaction model retention by CKD status (retained vs dropped participants).")

# Use explicit data-frame copy so survey design uses correct types (avoids ~0 coefficients from reference/lazy issues).
qdat_svy <- as.data.frame(qdat)
qdat_svy$phenoage_accel <- as.numeric(qdat_svy$phenoage_accel)
qdat_svy$mixture_qscore <- as.numeric(qdat_svy$mixture_qscore)
qdat_svy$ckd            <- as.integer(qdat_svy$ckd)
qdat_svy$age            <- as.numeric(qdat_svy$age)
qdat_svy$pir            <- as.numeric(qdat_svy$pir)
qdat_svy$bmi            <- as.numeric(qdat_svy$bmi)
qdat_svy$sex            <- factor(qdat_svy$sex)
qdat_svy$race_eth       <- factor(qdat_svy$race_eth)
qdat_svy$cycle          <- factor(qdat_svy$cycle)
if ("smoker3"  %in% names(qdat_svy)) qdat_svy$smoker3  <- factor(as.character(qdat_svy$smoker3),  levels = c("Never", "Former", "Current"))
if ("alcohol3" %in% names(qdat_svy)) qdat_svy$alcohol3 <- factor(as.character(qdat_svy$alcohol3), levels = c("Never", "Former", "Current"))
qdat_svy$weight_combined <- as.numeric(qdat_svy$weight_combined)
qdat_svy$psu            <- as.numeric(qdat_svy$psu)
qdat_svy$strata         <- as.numeric(qdat_svy$strata)
# Drop unused factor levels to avoid singularity with full survey design.
qdat_svy$sex      <- droplevels(factor(qdat_svy$sex))
qdat_svy$race_eth <- droplevels(factor(qdat_svy$race_eth))
qdat_svy$cycle    <- droplevels(factor(qdat_svy$cycle))
if ("smoker3"  %in% names(qdat_svy)) qdat_svy$smoker3  <- droplevels(qdat_svy$smoker3)
if ("alcohol3" %in% names(qdat_svy)) qdat_svy$alcohol3 <- droplevels(qdat_svy$alcohol3)

# ---- Pre-compute sensitivity-derived columns on qdat_svy so they ride into MI imputations ----
# These outcomes/redefinitions depend only on biomarkers + age (not imputed), so they can be computed once.
ph2_comp_qsvy <- compute_phenoage_components(qdat_svy, include_scr = FALSE, include_alb = TRUE)
ph2_qsvy <- ph2_comp_qsvy$phenoage
qdat_svy$phenoage_nokid <- NA_real_
keep2_qsvy <- which(!is.na(ph2_qsvy) & !is.na(qdat_svy$age))
if (length(keep2_qsvy) > 0) {
  fit_nokid_qsvy <- lm(ph2_qsvy[keep2_qsvy] ~ qdat_svy$age[keep2_qsvy])
  qdat_svy$phenoage_nokid[keep2_qsvy] <- residuals(fit_nokid_qsvy)
}
ph3_comp_qsvy <- compute_phenoage_components(qdat_svy, include_scr = FALSE, include_alb = FALSE)
ph3_qsvy <- ph3_comp_qsvy$phenoage
qdat_svy$phenoage_nokid_both <- NA_real_
keep3_qsvy <- which(!is.na(ph3_qsvy) & !is.na(qdat_svy$age))
if (length(keep3_qsvy) > 0) {
  fit_nokid_both_qsvy <- lm(ph3_qsvy[keep3_qsvy] ~ qdat_svy$age[keep3_qsvy])
  qdat_svy$phenoage_nokid_both[keep3_qsvy] <- residuals(fit_nokid_both_qsvy)
}
qdat_svy$ckd_sens <- ifelse(!is.na(qdat_svy$egfr), as.integer(qdat_svy$egfr < 60), NA_integer_)

# Mirror the same derivations on qdat for complete-case sensitivities.
qdat$phenoage_nokid       <- qdat_svy$phenoage_nokid[match(qdat$SEQN, qdat_svy$SEQN)]
qdat$phenoage_nokid_both  <- qdat_svy$phenoage_nokid_both[match(qdat$SEQN, qdat_svy$SEQN)]
qdat$ckd_sens             <- qdat_svy$ckd_sens[match(qdat$SEQN, qdat_svy$SEQN)]

# Full complex survey design (PSU/strata) for valid inference; droplevels above avoids empty cells.
svy_q <- tryCatch(
  svydesign(id = ~psu, strata = ~strata, weights = ~weight_combined, nest = TRUE, data = qdat_svy),
  error = function(e) {
    message("Full design failed (", conditionMessage(e), "); falling back to weights-only design.")
    svydesign(id = ~1, strata = NULL, weights = ~weight_combined, data = qdat_svy)
  }
)

