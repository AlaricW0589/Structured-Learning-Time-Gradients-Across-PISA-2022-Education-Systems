suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(broom)
  library(readr); library(ggplot2)
})

.find_helper <- function() {
  cands <- c(file.path(getwd(), "helpers", "proj_root.R"),
             file.path(getwd(), "..", "helpers", "proj_root.R"),
             "C:/Users/Moyih/auvshk/helpers/proj_root.R")
  for (p in cands) if (file.exists(p)) return(normalizePath(p))
  stop("Cannot find helpers/proj_root.R; cwd=", getwd())
}
source(.find_helper())
ROOT      <- proj_root()
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")
LOG_DIR   <- file.path(ROOT, "logs")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(LOG_DIR, "B2_imputed_composites.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[B2] === Block B v2 — composite imputation + PCA ===\n")
cat("[B2] Started at", format(Sys.time()), "\n\n")

panel <- readr::read_csv(file.path(TAB_DIR, "B1_country_panel.csv"),
                         show_col_types = FALSE) %>%
  dplyr::mutate(region = factor(region,
            levels = c("OECD (other)", "East Asia", "Nordic", "Non-OECD")))

raw_vars <- c("abgmath_mean", "stubeha_mean", "mactiv_mean",
              "edushort_mean", "stratio_mean", "private_share",
              "escs_btw_var")

cat("[B2] Country counts before imputation:\n")
print(panel %>%
        dplyr::summarise(dplyr::across(all_of(raw_vars),
                         ~ sum(!is.na(.x)))))

panel_imp <- panel %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(dplyr::across(all_of(raw_vars),
            ~ ifelse(is.na(.x), stats::median(.x, na.rm = TRUE), .x))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(all_of(raw_vars),
            ~ ifelse(is.na(.x), stats::median(.x, na.rm = TRUE), .x)))

cat("\n[B2] Country counts after regional-median imputation:\n")
print(panel_imp %>%
        dplyr::summarise(dplyr::across(all_of(raw_vars),
                         ~ sum(!is.na(.x)))))

zsc <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

panel_imp <- panel_imp %>%
  dplyr::mutate(
    exam_pressure_z  = zsc(abgmath_mean) + zsc(mactiv_mean) - zsc(stubeha_mean),
    stratification_z = zsc(escs_btw_var) + zsc(private_share),
    resource_z       = -zsc(edushort_mean) - zsc(stratio_mean)
  )

cat("\n[B2] Running PCA on the 7 raw country-level inputs ...\n")
X <- panel_imp[, raw_vars] |> as.matrix()
X_scaled <- scale(X)
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
cat(sprintf("  PC1 explains %.1f%%, PC2 %.1f%%, PC3 %.1f%% of variance\n",
            100 * pca$sdev[1]^2 / sum(pca$sdev^2),
            100 * pca$sdev[2]^2 / sum(pca$sdev^2),
            100 * pca$sdev[3]^2 / sum(pca$sdev^2)))

panel_imp$pc1_org <- pca$x[, "PC1"]
panel_imp$pc2_org <- pca$x[, "PC2"]

if (cor(panel_imp$pc1_org, panel_imp$exam_pressure_z, use = "complete.obs") < 0) {
  panel_imp$pc1_org <- -panel_imp$pc1_org
  pca$rotation[, "PC1"] <- -pca$rotation[, "PC1"]
}
if (cor(panel_imp$pc2_org, panel_imp$stratification_z, use = "complete.obs") < 0) {
  panel_imp$pc2_org <- -panel_imp$pc2_org
  pca$rotation[, "PC2"] <- -pca$rotation[, "PC2"]
}

loadings <- as.data.frame(pca$rotation[, 1:3])
loadings$variable <- rownames(loadings)
loadings <- loadings[, c("variable", "PC1", "PC2", "PC3")]
print(loadings, row.names = FALSE)
readr::write_csv(loadings, file.path(TAB_DIR, "B2_pca_loadings.csv"))

readr::write_csv(panel_imp, file.path(TAB_DIR, "B2_country_panel_imputed.csv"))
cat(sprintf("\n[B2] Saved imputed panel: %d countries x %d cols\n",
            nrow(panel_imp), ncol(panel_imp)))

fit_dat <- panel_imp %>%
  dplyr::filter(!is.na(slope_p50), !is.na(hwk_mean))
cat(sprintf("\n[B2] Fitting on %d countries (was 60 in B1)\n", nrow(fit_dat)))

models <- list(
  M0 = lm(slope_p50 ~ 1, data = fit_dat),
  M1 = lm(slope_p50 ~ hwk_mean, data = fit_dat),
  M2 = lm(slope_p50 ~ hwk_mean + region, data = fit_dat),
  M3z = lm(slope_p50 ~ hwk_mean + region + exam_pressure_z + stratification_z + resource_z,
           data = fit_dat),
  M4z = lm(slope_p50 ~ hwk_mean + exam_pressure_z + stratification_z + resource_z,
           data = fit_dat),
  M3pc = lm(slope_p50 ~ hwk_mean + region + pc1_org + pc2_org, data = fit_dat),
  M4pc = lm(slope_p50 ~ hwk_mean + pc1_org + pc2_org,           data = fit_dat),
  M5  = lm(slope_p50 ~ region, data = fit_dat),
  M6z = lm(slope_p50 ~ exam_pressure_z + stratification_z + resource_z, data = fit_dat),
  M6pc = lm(slope_p50 ~ pc1_org + pc2_org, data = fit_dat)
)

fit_summary <- purrr::imap_dfr(models, function(m, name) {
  data.frame(
    model   = name,
    n       = length(m$residuals),
    df_used = length(m$coefficients),
    r2      = summary(m)$r.squared,
    adj_r2  = summary(m)$adj.r.squared,
    aic     = stats::AIC(m),
    bic     = stats::BIC(m),
    rmse    = sqrt(mean(m$residuals^2))
  )
})
print(fit_summary, row.names = FALSE)
readr::write_csv(fit_summary, file.path(TAB_DIR, "B2_nested_models.csv"))

attenuate <- function(m_unadj, m_adj, terms) {
  tu <- broom::tidy(m_unadj); ta <- broom::tidy(m_adj)
  purrr::map_dfr(terms, function(tt) {
    a <- tu$estimate[tu$term == tt]
    b <- ta$estimate[ta$term == tt]
    if (length(a) == 0 || length(b) == 0)
      return(NULL)
    data.frame(
      adj_model       = deparse(m_adj$call$formula),
      term            = tt,
      M2_estimate     = a,
      M2_se           = tu$std.error[tu$term == tt],
      M3_estimate     = b,
      M3_se           = ta$std.error[ta$term == tt],
      attenuation_pct = 100 * (a - b) / a
    )
  })
}

region_terms <- grep("^region", names(coef(models$M2)), value = TRUE)
attn_z  <- attenuate(models$M2, models$M3z,  region_terms)
attn_pc <- attenuate(models$M2, models$M3pc, region_terms)
attn <- dplyr::bind_rows(
  attn_z  %>% dplyr::mutate(spec = "z-score composites (M3z)"),
  attn_pc %>% dplyr::mutate(spec = "PCA components (M3pc)")
)
cat("\n[B2] Region attenuation (M2 vs M3z, vs M3pc):\n\n")
print(attn[, c("spec","term","M2_estimate","M3_estimate","attenuation_pct")],
      row.names = FALSE)
readr::write_csv(attn, file.path(TAB_DIR, "B2_attenuation.csv"))

coef_tab <- purrr::imap_dfr(
  models[c("M1","M2","M3z","M3pc","M4z","M4pc","M5","M6z","M6pc")],
  function(m, nm) broom::tidy(m, conf.int = TRUE) %>% dplyr::mutate(model = nm)
) %>% dplyr::select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high)
readr::write_csv(coef_tab, file.path(TAB_DIR, "B2_coefficients.csv"))

attn_plot_df <- bind_rows(
  broom::tidy(models$M2,    conf.int = TRUE) %>% dplyr::mutate(spec = "M2: + region"),
  broom::tidy(models$M3z,   conf.int = TRUE) %>% dplyr::mutate(spec = "M3z: + region + z-composites"),
  broom::tidy(models$M3pc,  conf.int = TRUE) %>% dplyr::mutate(spec = "M3pc: + region + PC1+PC2")
) %>%
  dplyr::filter(grepl("^region", term))

p1 <- ggplot2::ggplot(attn_plot_df,
        ggplot2::aes(x = term, y = estimate, color = spec)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                         position = ggplot2::position_dodge(0.5),
                         width = 0.25, linewidth = 0.7) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Region coefficients on country slope: M2 vs composite-adjusted (n=77)",
    subtitle = "Mediation = M3 estimates shrink toward zero relative to M2",
    x = NULL, y = "Coefficient on country slope (p50)", color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "B2_region_attenuation.png"),
                p1, width = 9, height = 5.5, dpi = 150)

panel_long <- fit_dat %>%
  tidyr::pivot_longer(c(pc1_org, pc2_org),
                      names_to = "component", values_to = "score") %>%
  dplyr::mutate(component = dplyr::recode(component,
                  pc1_org = "PC1 (organizational regime)",
                  pc2_org = "PC2 (residual)"))

p2 <- ggplot2::ggplot(panel_long,
        ggplot2::aes(x = score, y = slope_p50, color = region, label = country)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm",
                       color = "black", se = TRUE, linewidth = 0.6) +
  ggplot2::geom_point(size = 2, alpha = 0.85) +
  ggplot2::geom_text(size = 2.4, vjust = -0.7, alpha = 0.7, show.legend = FALSE) +
  ggplot2::facet_wrap(~ component, scales = "free_x") +
  ggplot2::labs(
    title    = "Country slope vs PCA components of school-organizational variables",
    x = "PC score", y = "Marginal slope at country median (PISA points / hwk hour)"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "B2_pca_vs_slope.png"),
                p2, width = 11, height = 5, dpi = 150)

r2 <- function(m) summary(m)$r.squared
gain <- function(big, small) (r2(big) - r2(small)) * 100

cat("\n[B2] === DECISION SUMMARY (n=77, regional-median imputation) ===\n")
cat(sprintf("  R^2(M0)                         = %.3f\n", r2(models$M0)))
cat(sprintf("  R^2(M1: hwk_mean)                = %.3f\n", r2(models$M1)))
cat(sprintf("  R^2(M2: + region)                = %.3f\n", r2(models$M2)))
cat(sprintf("  R^2(M3z: + z-composites)         = %.3f\n", r2(models$M3z)))
cat(sprintf("  R^2(M3pc: + PC1+PC2)             = %.3f\n", r2(models$M3pc)))
cat(sprintf("  R^2(M4z: hwk + z-composites)     = %.3f\n", r2(models$M4z)))
cat(sprintf("  R^2(M4pc: hwk + PCs)             = %.3f\n", r2(models$M4pc)))
cat(sprintf("  R^2(M5: region only)             = %.3f\n", r2(models$M5)))
cat(sprintf("  R^2(M6z: z-composites only)      = %.3f\n", r2(models$M6z)))
cat(sprintf("  R^2(M6pc: PCs only)              = %.3f\n", r2(models$M6pc)))
cat(sprintf("\n  Region adds beyond hwk_mean      : +%.2f pp\n",
            gain(models$M2, models$M1)))
cat(sprintf("  z-composites add beyond M2       : +%.2f pp\n",
            gain(models$M3z, models$M2)))
cat(sprintf("  PCs add beyond M2                 : +%.2f pp\n",
            gain(models$M3pc, models$M2)))

cat(sprintf("\n  Mean region attenuation (z)      : %.1f%%\n",
            mean(attn_z$attenuation_pct, na.rm = TRUE)))
cat(sprintf("  Mean region attenuation (PC)     : %.1f%%\n",
            mean(attn_pc$attenuation_pct, na.rm = TRUE)))

if (gain(models$M3z, models$M2) > 5 ||
    mean(attn_z$attenuation_pct, na.rm = TRUE) > 20 ||
    gain(models$M3pc, models$M2) > 5 ||
    mean(attn_pc$attenuation_pct, na.rm = TRUE) > 20) {
  cat("\n[B2] v2 H4 (or some form of mediation) SUPPORTED on n=77.\n")
} else {
  cat("\n[B2] v2 H4 (mediation) STILL REJECTED even after recovering the dropped countries.\n")
}

cat("\n[B2] Saved:\n",
    "  ", file.path(TAB_DIR, "B2_country_panel_imputed.csv"), "\n",
    "  ", file.path(TAB_DIR, "B2_nested_models.csv"),         "\n",
    "  ", file.path(TAB_DIR, "B2_attenuation.csv"),            "\n",
    "  ", file.path(TAB_DIR, "B2_coefficients.csv"),           "\n",
    "  ", file.path(TAB_DIR, "B2_pca_loadings.csv"),           "\n",
    "  ", file.path(FIG_DIR, "B2_region_attenuation.png"),     "\n",
    "  ", file.path(FIG_DIR, "B2_pca_vs_slope.png"),           "\n", sep = "")

cat("\n[B2] Finished at", format(Sys.time()), "\n")
