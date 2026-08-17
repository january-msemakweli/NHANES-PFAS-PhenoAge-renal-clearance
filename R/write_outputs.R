# Export manuscript tables and check them against the age-adjusted MI results.

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

if (exists("table1")) {
  write.csv(table1, "output/tables/table1_characteristics.csv", row.names = FALSE)
}
if (exists("flow")) {
  write.csv(flow, "output/tables/sample_flow.csv", row.names = FALSE)
}
if (exists("tab2_full")) {
  write.csv(tab2_full, "output/tables/table2_mixture.csv", row.names = FALSE)
}
if (exists("tab3")) {
  write.csv(tab3, "output/tables/table3_ckd_interaction.csv", row.names = FALSE)
}
if (exists("tab4")) {
  write.csv(tab4, "output/tables/table4_sensitivities.csv", row.names = FALSE)
}
if (exists("component_table_fmt")) {
  write.csv(component_table_fmt, "output/tables/table4_components.csv", row.names = FALSE)
} else if (exists("component_table")) {
  write.csv(component_table, "output/tables/table4_components.csv", row.names = FALSE)
}
if (exists("single_pollutant")) {
  write.csv(single_pollutant, "output/tables/single_pollutant.csv", row.names = FALSE)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

chk <- function(label, got, expected, tol = 0.002) {
  g <- to_num(got)[1]
  ok <- is.finite(g) && abs(g - expected) <= tol
  list(label = label, expected = expected, got = g, ok = isTRUE(ok))
}

rows <- list()
if (exists("tab2_full") && nrow(tab2_full) >= 1) {
  rows[[length(rows) + 1]] <- chk("Table 2 Model A beta", tab2_full$beta[1], -1.547)
  rows[[length(rows) + 1]] <- chk("Table 2 Model A lo", tab2_full$conf.low[1], -1.773)
  rows[[length(rows) + 1]] <- chk("Table 2 Model A hi", tab2_full$conf.high[1], -1.322)
}
if (exists("tab3")) {
  mix <- tab3[tab3$term == "mixture_qscore", , drop = FALSE]
  ckd <- tab3[tab3$term == "ckd", , drop = FALSE]
  int <- tab3[grepl("mixture_qscore:ckd", tab3$term), , drop = FALSE]
  m0 <- tab3[grepl("CKD=0", tab3$term), , drop = FALSE]
  m1 <- tab3[grepl("CKD=1", tab3$term), , drop = FALSE]
  if (nrow(mix)) rows[[length(rows) + 1]] <- chk("Table 3 mixture", mix$beta[1], -0.971)
  if (nrow(int)) rows[[length(rows) + 1]] <- chk("Table 3 interaction", int$beta[1], -1.253)
  if (nrow(m0)) rows[[length(rows) + 1]] <- chk("Table 3 no CKD", m0$beta[1], -0.971)
  if (nrow(m1)) rows[[length(rows) + 1]] <- chk("Table 3 CKD", m1$beta[1], -2.224)
}
if (exists("component_table_fmt")) {
  ct <- component_table_fmt
  pick <- function(nm) ct[ct[[1]] == nm, , drop = FALSE]
  alb <- pick("Albumin (g/dL)")
  scr <- pick("Creatinine (mg/dL)")
  glu <- pick("Glucose (mg/dL)")
  if (nrow(alb)) rows[[length(rows) + 1]] <- chk("Table 4 albumin", alb$beta[1], 0.066)
  if (nrow(scr)) rows[[length(rows) + 1]] <- chk("Table 4 creatinine", scr$beta[1], -0.014)
  if (nrow(glu)) rows[[length(rows) + 1]] <- chk("Table 4 glucose", glu$beta[1], -0.952)
}

ver <- do.call(rbind, lapply(rows, function(r) {
  data.frame(check = r$label, expected = r$expected, got = r$got, ok = r$ok)
}))

n_ok <- if (exists("n_mi_base")) n_mi_base else if (exists("n_modelb_primary")) n_modelb_primary else NA
message("MICE / Model B N: ", n_ok)

needed_figs <- c(
  "figures/fig1_hex_facet_CKD.jpg",
  "figures/fig2_forest_marginal_effects.jpg",
  "figures/fig3_mice_convergence.jpg",
  "figures/participant_flow.jpg"
)
fig_ok <- file.exists(needed_figs)
ver <- rbind(
  ver,
  data.frame(check = needed_figs, expected = 1, got = as.numeric(fig_ok), ok = fig_ok)
)

needed_tabs <- c(
  "output/tables/table1_characteristics.csv",
  "output/tables/table2_mixture.csv",
  "output/tables/table3_ckd_interaction.csv",
  "output/tables/table4_components.csv",
  "output/tables/table4_sensitivities.csv",
  "output/tables/single_pollutant.csv",
  "output/tables/sample_flow.csv"
)
tab_ok <- file.exists(needed_tabs)
ver <- rbind(
  ver,
  data.frame(check = needed_tabs, expected = 1, got = as.numeric(tab_ok), ok = tab_ok)
)
print(ver)

if (!isTRUE(all(ver$ok)) || !identical(as.integer(n_ok), 7633L)) {
  stop("Verification failed: estimates, figures, or N did not match the manuscript targets.")
}
message("Verification passed: manuscript Table 2/3/4 targets, figures, and N=7633 matched.")
