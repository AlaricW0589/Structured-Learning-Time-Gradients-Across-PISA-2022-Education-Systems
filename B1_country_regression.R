suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  library(readr)
  library(ggplot2)
})

.find_helper <- function() {
  cands <- c(
    file.path(getwd(), "helpers", "proj_root.R"),
    file.path(getwd(), "..", "helpers", "proj_root.R"),
    "C:/Users/Moyih/auvshk/helpers/proj_root.R"
  )
  for (p in cands) if (file.exists(p)) return(normalizePath(p))
  stop("Cannot find helpers/proj_root.R; cwd=", getwd())
}
source(.find_helper())
ROOT      <- proj_root()
CLEAN_DIR <- file.path(ROOT, "data", "clean")
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")
LOG_DIR   <- file.path(ROOT, "logs")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "B1_country_regression.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[B1] === Block B — country-level decomposition of slope heterogeneity ===\n")
cat("[B1] Started at", format(Sys.time()), "\n")
cat("[B1] ROOT =", ROOT, "\n\n")

a2_path <- file.path(TAB_DIR, "A2_position_slope.csv")
if (!file.exists(a2_path)) {
  stop("A2_position_slope.csv not found. Run scripts/A2_position_slope.R first.")
}
a2 <- readr::read_csv(a2_path, show_col_types = FALSE)
cat(sprintf("[B1] Loaded A2: %d countries\n", nrow(a2)))

dat <- as.data.frame(arrow::read_parquet(
  file.path(CLEAN_DIR, "pisa2022_core.parquet")
))
cat(sprintf("[B1] Loaded parquet: %d students x %d cols\n",
            nrow(dat), ncol(dat)))

cat("\n[B1] Building country-level school-organizational aggregates ...\n")

school_country <- dat %>%
  dplyr::group_by(CNT, CNTSCHID) %>%
  dplyr::summarise(
    sch_w        = sum(W_FSTUWT, na.rm = TRUE),
    SCHSIZE      = mean(SCHSIZE,    na.rm = TRUE),
    STRATIO      = mean(STRATIO,    na.rm = TRUE),
    SC013Q01TA   = mean(SC013Q01TA, na.rm = TRUE),
    ABGMATH      = mean(ABGMATH,    na.rm = TRUE),
    STUBEHA      = mean(STUBEHA,    na.rm = TRUE),
    MACTIV       = mean(MACTIV,     na.rm = TRUE),
    EDUSHORT     = mean(EDUSHORT,   na.rm = TRUE),
    CLSIZE       = mean(CLSIZE,     na.rm = TRUE),
    sch_escs     = stats::weighted.mean(ESCS,  W_FSTUWT, na.rm = TRUE),
    sch_pvm      = stats::weighted.mean(pv_mean_math, W_FSTUWT, na.rm = TRUE),
    .groups      = "drop"
  )

country_org <- school_country %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(
    n_schools     = dplyr::n(),
    schsize_mean  = stats::weighted.mean(SCHSIZE,    sch_w, na.rm = TRUE),
    stratio_mean  = stats::weighted.mean(STRATIO,    sch_w, na.rm = TRUE),
    private_share = stats::weighted.mean(SC013Q01TA == 2, sch_w, na.rm = TRUE),
    abgmath_mean  = stats::weighted.mean(ABGMATH,    sch_w, na.rm = TRUE),
    stubeha_mean  = stats::weighted.mean(STUBEHA,    sch_w, na.rm = TRUE),
    mactiv_mean   = stats::weighted.mean(MACTIV,     sch_w, na.rm = TRUE),
    edushort_mean = stats::weighted.mean(EDUSHORT,   sch_w, na.rm = TRUE),
    clsize_mean   = stats::weighted.mean(CLSIZE,     sch_w, na.rm = TRUE),
    escs_btw_var  = stats::weighted.mean(
                      (sch_escs - stats::weighted.mean(sch_escs, sch_w, na.rm=TRUE))^2,
                      sch_w, na.rm = TRUE),
    .groups       = "drop"
  )

cat(sprintf("[B1]   Country aggregates built for %d countries\n",
            nrow(country_org)))

zsc <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

country_org <- country_org %>%
  dplyr::mutate(
    exam_pressure_z = zsc(abgmath_mean) + zsc(mactiv_mean) - zsc(stubeha_mean),
    stratification_z = zsc(escs_btw_var) + zsc(private_share),
    resource_z = -zsc(edushort_mean) - zsc(stratio_mean)
  )

panel <- a2 %>%
  dplyr::inner_join(country_org, by = c("country" = "CNT")) %>%
  dplyr::mutate(
    region = factor(region,
                    levels = c("OECD (other)", "East Asia", "Nordic", "Non-OECD"))
  )

readr::write_csv(panel, file.path(TAB_DIR, "B1_country_panel.csv"))
cat(sprintf("[B1]   Saved panel: %d countries x %d cols -> %s\n",
            nrow(panel), ncol(panel),
            file.path(TAB_DIR, "B1_country_panel.csv")))

fit_dat <- panel %>%
  dplyr::filter(
    !is.na(slope_p50), !is.na(hwk_mean),
    !is.na(exam_pressure_z), !is.na(stratification_z), !is.na(resource_z)
  )
cat(sprintf("[B1]   %d countries enter the regressions (out of %d)\n",
            nrow(fit_dat), nrow(panel)))

cat("\n[B1] Fitting nested OLS models (outcome = slope_p50) ...\n")

models <- list(
  M0 = lm(slope_p50 ~ 1, data = fit_dat),
  M1 = lm(slope_p50 ~ hwk_mean, data = fit_dat),
  M2 = lm(slope_p50 ~ hwk_mean + region, data = fit_dat),
  M3 = lm(slope_p50 ~ hwk_mean + region + exam_pressure_z + stratification_z + resource_z,
          data = fit_dat),
  M4 = lm(slope_p50 ~ hwk_mean + exam_pressure_z + stratification_z + resource_z,
          data = fit_dat),
  M5 = lm(slope_p50 ~ region, data = fit_dat),
  M6 = lm(slope_p50 ~ exam_pressure_z + stratification_z + resource_z, data = fit_dat)
)

fit_summary <- purrr::imap_dfr(models, function(m, name) {
  data.frame(
    model      = name,
    n          = length(m$residuals),
    df_used    = length(m$coefficients),
    r2         = summary(m)$r.squared,
    adj_r2     = summary(m)$adj.r.squared,
    aic        = stats::AIC(m),
    bic        = stats::BIC(m),
    rmse       = sqrt(mean(m$residuals^2))
  )
})
print(fit_summary, row.names = FALSE)
readr::write_csv(fit_summary,
                 file.path(TAB_DIR, "B1_nested_models.csv"))

cat("\n[B1] Coefficient tables for headline models:\n\n")
coef_tab <- purrr::imap_dfr(models[c("M1", "M2", "M3", "M4")], function(m, nm) {
  broom::tidy(m, conf.int = TRUE) %>%
    dplyr::mutate(model = nm)
}) %>%
  dplyr::select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high)
print(coef_tab, row.names = FALSE)
readr::write_csv(coef_tab,
                 file.path(TAB_DIR, "B1_coefficients.csv"))

cat("\n[B1] Region attenuation (M2 -> M3) -- the formal v2 H4 test:\n\n")

attenuate <- function(m_unadj, m_adj, terms) {
  tu <- broom::tidy(m_unadj)
  ta <- broom::tidy(m_adj)
  out <- purrr::map_dfr(terms, function(tt) {
    a <- tu$estimate[tu$term == tt]
    b <- ta$estimate[ta$term == tt]
    sa <- tu$std.error[tu$term == tt]
    sb <- ta$std.error[ta$term == tt]
    if (length(a) == 0 || length(b) == 0)
      return(data.frame(term = tt, M2_estimate = NA, M3_estimate = NA,
                        attenuation_pp = NA))
    data.frame(
      term            = tt,
      M2_estimate     = a,
      M2_se           = sa,
      M3_estimate     = b,
      M3_se           = sb,
      attenuation_pct = 100 * (a - b) / a
    )
  })
  out
}

region_terms <- grep("^region", names(coef(models$M2)), value = TRUE)
attn <- attenuate(models$M2, models$M3, region_terms)
print(attn, row.names = FALSE)
readr::write_csv(attn, file.path(TAB_DIR, "B1_attenuation.csv"))

cat("\n[B1] Country leave-one-out for M2 ...\n")

loo <- purrr::map_dfr(unique(fit_dat$country), function(cc) {
  m <- lm(slope_p50 ~ hwk_mean + region,
          data = fit_dat %>% dplyr::filter(country != cc))
  td <- broom::tidy(m) %>% dplyr::mutate(dropped = cc)
  td
})
loo_region <- loo %>% dplyr::filter(grepl("^region", term))

readr::write_csv(loo_region, file.path(TAB_DIR, "B1_LOO.csv"))

cat(sprintf("[B1]   LOO over %d countries; %d region-coefficient rows\n",
            length(unique(fit_dat$country)), nrow(loo_region)))

attn_plot <- bind_rows(
  broom::tidy(models$M2, conf.int = TRUE) %>%
    dplyr::mutate(spec = "M2: + region"),
  broom::tidy(models$M3, conf.int = TRUE) %>%
    dplyr::mutate(spec = "M3: + region + composites")
) %>%
  dplyr::filter(grepl("^region", term)) %>%
  dplyr::mutate(term = factor(term))

p1 <- ggplot2::ggplot(attn_plot,
        ggplot2::aes(x = term, y = estimate, color = spec)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.4), size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                         position = ggplot2::position_dodge(0.4), width = 0.2,
                         linewidth = 0.7) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Region coefficients on country slope: M2 vs M3 (composite-adjusted)",
    subtitle = "If composites mediate region, M3 estimates should shrink toward zero",
    x = NULL, y = "Coefficient on country slope (p50)",
    color = NULL,
    caption = "Outcome: country marginal slope of math achievement w.r.t. homework hours at country median."
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "B1_region_attenuation.png"),
                p1, width = 8, height = 5, dpi = 150)

p2 <- ggplot2::ggplot(loo_region, ggplot2::aes(x = estimate, y = term)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_jitter(height = 0.15, alpha = 0.4, size = 0.9) +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
                        color = "firebrick") +
  ggplot2::labs(
    title    = "Leave-one-country-out distribution of region coefficients (M2)",
    subtitle = "Each dot = one country dropped; red diamond = mean across LOO fits",
    x = "Coefficient on country slope (p50)", y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "B1_LOO.png"),
                p2, width = 8, height = 4, dpi = 150)

panel_long <- fit_dat %>%
  tidyr::pivot_longer(c(exam_pressure_z, stratification_z, resource_z),
                      names_to = "composite", values_to = "z") %>%
  dplyr::mutate(composite = dplyr::recode(composite,
    exam_pressure_z   = "Exam-pressure regime",
    stratification_z  = "School stratification",
    resource_z        = "Resource adequacy"
  ))

p3 <- ggplot2::ggplot(panel_long,
        ggplot2::aes(x = z, y = slope_p50, color = region, label = country)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm",
                       color = "black", se = TRUE, linewidth = 0.6) +
  ggplot2::geom_point(size = 2, alpha = 0.85) +
  ggplot2::geom_text(size = 2.4, vjust = -0.7, alpha = 0.7, show.legend = FALSE) +
  ggplot2::facet_wrap(~ composite, scales = "free_x") +
  ggplot2::labs(
    title    = "Country slope vs school-organizational composites",
    x = "Composite z-score", y = "Marginal slope at country median (PISA points / hwk hour)"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "B1_composites_vs_slope.png"),
                p3, width = 12, height = 4.5, dpi = 150)

r2 <- function(m) summary(m)$r.squared
gain <- function(big, small) (r2(big) - r2(small))

cat("\n[B1] === DECISION SUMMARY ===\n")
cat(sprintf("  R^2(M0)                                = %.3f\n", r2(models$M0)))
cat(sprintf("  R^2(M1: hwk_mean)                       = %.3f\n", r2(models$M1)))
cat(sprintf("  R^2(M2: hwk_mean + region)              = %.3f\n", r2(models$M2)))
cat(sprintf("  R^2(M3: + composites)                   = %.3f\n", r2(models$M3)))
cat(sprintf("  R^2(M4: hwk_mean + composites only)     = %.3f\n", r2(models$M4)))
cat(sprintf("  R^2(M5: region only)                    = %.3f\n", r2(models$M5)))
cat(sprintf("  R^2(M6: composites only)                = %.3f\n", r2(models$M6)))
cat(sprintf("  Region adds beyond hwk_mean             : +%.3f pp\n",
            100 * gain(models$M2, models$M1)))
cat(sprintf("  Composites add beyond hwk_mean+region   : +%.3f pp\n",
            100 * gain(models$M3, models$M2)))
cat(sprintf("  Composites add beyond hwk_mean alone    : +%.3f pp\n",
            100 * gain(models$M4, models$M1)))

mean_attn <- mean(attn$attenuation_pct, na.rm = TRUE)
cat(sprintf("\n  Mean region attenuation M2 -> M3        : %.1f%%\n", mean_attn))

if (gain(models$M3, models$M2) > 0.05 && mean_attn > 20) {
  cat("\n[B1] v2 H4 SUPPORTED: composites add explanatory power AND attenuate region.\n")
} else if (gain(models$M3, models$M2) > 0.05) {
  cat("\n[B1] v2 H4 PARTIAL: composites add explanatory power but do NOT attenuate region.\n")
} else {
  cat("\n[B1] v2 H4 REJECTED: composites do not add explanatory power beyond region.\n")
}

cat("\n[B1] Saved:\n",
    "  ", file.path(TAB_DIR, "B1_country_panel.csv"), "\n",
    "  ", file.path(TAB_DIR, "B1_nested_models.csv"), "\n",
    "  ", file.path(TAB_DIR, "B1_coefficients.csv"), "\n",
    "  ", file.path(TAB_DIR, "B1_attenuation.csv"), "\n",
    "  ", file.path(TAB_DIR, "B1_LOO.csv"), "\n",
    "  ", file.path(FIG_DIR, "B1_region_attenuation.png"), "\n",
    "  ", file.path(FIG_DIR, "B1_LOO.png"), "\n",
    "  ", file.path(FIG_DIR, "B1_composites_vs_slope.png"), "\n",
    sep = "")

cat("\n[B1] Finished at", format(Sys.time()), "\n")
