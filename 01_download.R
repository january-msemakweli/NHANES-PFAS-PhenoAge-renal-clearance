# 01: download and bind NHANES cycles.
## ----cycle-download-----------------------------------------------------------
dir.create("data", showWarnings = FALSE, recursive = TRUE)
if (file.exists("data/nhanes_raw.rds")) {
  message("Loading cached NHANES extract: data/nhanes_raw.rds")
  nhanes_raw <- readRDS("data/nhanes_raw.rds")
} else {

# Cycle-aware: 2003-2004 (C) uses Lxx_ names (L40, L25, L06CRP, L06ALB, L24PFC, L06BMT, L06UHM); 2005+ use BIOPRO/CBC/CRP/PFC; 2013+ PFAS; I/J use HSCRP.
table_candidates <- list(
  DEMO = c("DEMO"),
  BMX = c("BMX"),
  SMQ = c("SMQ"),
  ALQ = c("ALQ"),
  PAQ = c("PAQ"),
  DIQ = c("DIQ"),
  BPQ = c("BPQ"),
  MCQ = c("MCQ"),
  BIOPRO = c("L40", "BIOPRO"),
  CBC = c("L25", "CBC"),
  CRP = c("L06CRP", "CRP", "HSCRP"),
  ALB_CR = c("L06ALB", "ALB_CR"),
  GLU = c("L10", "GLU"),
  PFAS = c("L24PFC", "PFAS", "PFC", "P_PFAS"),
  BLOOD_METALS = c("L06BMT", "PBCD"),
  URINE_METALS = c("L06UHM", "UHM", "UM")
)

download_cycle <- function(cycle, exposure_family = "pfas", verbose = TRUE) {
  core_domains <- c("DEMO", "BMX", "SMQ", "ALQ", "PAQ", "DIQ", "BPQ", "MCQ", "BIOPRO", "CBC", "CRP", "ALB_CR", "GLU")
  exp_domains <- switch(
    exposure_family,
    pfas = c("PFAS"),
    metals = c("BLOOD_METALS", "URINE_METALS"),
    both = c("PFAS", "BLOOD_METALS", "URINE_METALS"),
    c("PFAS")
  )
  all_domains <- c(core_domains, exp_domains)

  resolved <- lapply(all_domains, function(dom) resolve_table(cycle, table_candidates[[dom]], domain = dom, verbose = verbose))
  names(resolved) <- all_domains
  resolved <- resolved[!vapply(resolved, is.null, logical(1))]
  if (!"DEMO" %in% names(resolved)) return(NULL)

  dat <- resolved$DEMO$data
  if (!"SEQN" %in% names(dat)) return(NULL)

  for (nm in setdiff(names(resolved), "DEMO")) {
    d <- resolved[[nm]]$data
    # Rename glucose vars from GLU table before joining so they do not collide with BIOPRO/L40.
    if (nm == "GLU") {
      if ("LBXGLU" %in% names(d)) names(d)[names(d) == "LBXGLU"] <- "LBXGLU_GLU"
      if ("LBXSGL" %in% names(d)) names(d)[names(d) == "LBXSGL"] <- "LBXSGL_GLU"
    }
    if ("SEQN" %in% names(d)) dat <- left_join(dat, d, by = "SEQN")
  }
  # After merge: define explicit glucose sources. Primary analysis uses MEC chemistry glucose for measurement consistency in PhenoAge construction.
  dat$glu_mec  <- coalesce_num(dat, c("LBXGLU", "LBXSGL"))
  dat$glu_glu  <- coalesce_num(dat, c("LBXGLU_GLU", "LBXSGL_GLU"))
  dat$glu     <- dat$glu_mec
  dat$cycle <- cycle
  # Store provenance for key tables so we can audit weight/source compatibility.
  if ("CRP" %in% names(resolved)) dat$crp_table <- resolved$CRP$table
  if ("GLU" %in% names(resolved)) dat$glu_table <- resolved$GLU$table
  # Warn if joins created .x/.y duplicate columns (would break coalesce_num).
  dup_suffix_cols <- names(dat)[grepl("\\.x$|\\.y$", names(dat))]
  if (length(dup_suffix_cols) > 0) {
    warning("Cycle ", cycle, ": duplicate columns from join: ", paste(dup_suffix_cols, collapse = ", "))
  }
  dat
}

cycle_data <- lapply(config$cycles, download_cycle, exposure_family = config$exposure_family, verbose = config$verbose)
cycle_data <- cycle_data[!vapply(cycle_data, is.null, logical(1))]
if (length(cycle_data) == 0) stop("No cycles were successfully downloaded.")

all_cols <- sort(unique(unlist(lapply(cycle_data, names))))
cycle_data <- lapply(cycle_data, function(df) {
  miss <- setdiff(all_cols, names(df))
  if (length(miss) > 0) for (m in miss) df[[m]] <- NA
  df <- df[, all_cols, drop = FALSE]
  df %>% mutate(across(everything(), as.character))
})
nhanes_raw <- bind_rows(cycle_data)
  saveRDS(nhanes_raw, "data/nhanes_raw.rds")
  message("Wrote data/nhanes_raw.rds")
}

kable(
  data.frame(
    cycle = names(table(config$cycles))[names(table(config$cycles)) %in% unique(nhanes_raw$cycle)],
    n = as.integer(table(nhanes_raw$cycle))
  ),
  caption = "Cycle availability and downloaded participant counts"
)


