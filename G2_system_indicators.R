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
RAW_DIR   <- file.path(ROOT, "data", "raw")
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")
LOG_DIR   <- file.path(ROOT, "logs")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "G2_system_indicators.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[G2] === System-indicators audit ===\n")
cat("[G2] Started at", format(Sys.time()), "\n\n")

a2 <- readr::read_csv(file.path(TAB_DIR, "A2_position_slope.csv"),
                      show_col_types = FALSE) %>%
  dplyr::select(country, region, hwk_mean, slope_p50)

b1 <- readr::read_csv(file.path(TAB_DIR, "B1_country_panel.csv"),
                      show_col_types = FALSE)

sys <- readr::read_csv(file.path(RAW_DIR, "system_indicators.csv"),
                       show_col_types = FALSE) %>%
  dplyr::select(country, tracking_age, early_track, central_exit_exam,
                confucian_heritage, vocational_share_pct, instr_hrs_per_wk)

cat(sprintf("[G2] system_indicators.csv loaded with %d countries\n",
            nrow(sys)))

panel <- a2 %>%
  dplyr::inner_join(sys, by = "country") %>%
  dplyr::mutate(
    region = factor(region,
                    levels = c("OECD (other)", "East Asia", "Nordic", "Non-OECD"))
  )

readr::write_csv(panel, file.path(TAB_DIR, "G2_system_panel.csv"))

cat(sprintf("[G2] Merged panel: %d countries\n", nrow(panel)))
print(table(panel$region))

fit_dat <- panel %>%
  dplyr::filter(!is.na(slope_p50), !is.na(hwk_mean), !is.na(tracking_age),
                !is.na(central_exit_exam), !is.na(confucian_heritage),
                !is.na(vocational_share_pct), !is.na(instr_hrs_per_wk))
cat(sprintf("\n[G2] %d countries enter regressions\n", nrow(fit_dat)))

models <- list(
  M0 = lm(slope_p50 ~ 1, data = fit_dat),
  M1 = lm(slope_p50 ~ hwk_mean, data = fit_dat),
  M2 = lm(slope_p50 ~ hwk_mean + region, data = fit_dat),
  M3a = lm(slope_p50 ~ hwk_mean + tracking_age + central_exit_exam +
                       confucian_heritage + vocational_share_pct +
                       instr_hrs_per_wk, data = fit_dat),
  M3b = lm(slope_p50 ~ hwk_mean + region + tracking_age + central_exit_exam +
                       confucian_heritage + vocational_share_pct +
                       instr_hrs_per_wk, data = fit_dat),
  M4 = lm(slope_p50 ~ hwk_mean + confucian_heritage, data = fit_dat),
  M5 = lm(slope_p50 ~ hwk_mean + region + confucian_heritage, data = fit_dat),
  M6 = lm(slope_p50 ~ tracking_age + central_exit_exam + confucian_heritage +
                       vocational_share_pct + instr_hrs_per_wk, data = fit_dat)
)

fit_summary <- purrr::imap_dfr(models, function(m, name) {
  data.frame(
    model = name,
    n = length(m$residuals),
    df_used = length(m$coefficients),
    r2 = summary(m)$r.squared,
    adj_r2 = summary(m)$adj.r.squared,
    aic = stats::AIC(m),
    rmse = sqrt(mean(m$residuals^2))
  )
})

cat("\n[G2] Nested model fit summary:\n")
print(fit_summary, row.names = FALSE)
readr::write_csv(fit_summary, file.path(TAB_DIR, "G2_nested_models.csv"))

cat("\n[G2] Headline model coefficients:\n\n")
coef_tab <- purrr::imap_dfr(models[c("M2", "M3a", "M3b", "M4", "M5")],
  function(m, nm) {
    broom::tidy(m, conf.int = TRUE) %>%
      dplyr::mutate(model = nm)
  }) %>%
  dplyr::select(model, term, estimate, std.error, statistic, p.value,
                conf.low, conf.high)
print(as.data.frame(coef_tab), row.names = FALSE, digits = 3)
readr::write_csv(coef_tab, file.path(TAB_DIR, "G2_coefficients.csv"))

cat("\n[G2] Region attenuation M2 -> M3b (full system indicator set):\n")
attenuate <- function(m_unadj, m_adj, terms) {
  tu <- broom::tidy(m_unadj); ta <- broom::tidy(m_adj)
  purrr::map_dfr(terms, function(tt) {
    a <- tu$estimate[tu$term == tt]; b <- ta$estimate[ta$term == tt]
    if (length(a) == 0 || length(b) == 0)
      return(data.frame(term = tt, M2_estimate = NA, M3b_estimate = NA,
                        attenuation_pct = NA))
    data.frame(term = tt, M2_estimate = a, M3b_estimate = b,
               attenuation_pct = 100 * (a - b) / a)
  })
}
region_terms <- grep("^region", names(coef(models$M2)), value = TRUE)
attn <- attenuate(models$M2, models$M3b, region_terms)
print(attn, row.names = FALSE)
readr::write_csv(attn, file.path(TAB_DIR, "G2_attenuation.csv"))

cat("\n[G2] Region attenuation M2 -> M5 (region + Confucian heritage only):\n")
attn_conf <- attenuate(models$M2, models$M5, region_terms)
print(attn_conf, row.names = FALSE)

attn_plot <- bind_rows(
  broom::tidy(models$M2, conf.int = TRUE)  %>% dplyr::mutate(spec = "M2: region only"),
  broom::tidy(models$M5, conf.int = TRUE)  %>% dplyr::mutate(spec = "M5: + Confucian heritage"),
  broom::tidy(models$M3b, conf.int = TRUE) %>% dplyr::mutate(spec = "M3b: + 5 system indicators")
) %>%
  dplyr::filter(grepl("^region", term)) %>%
  dplyr::mutate(term = dplyr::recode(term,
    "regionEast Asia" = "East Asia",
    "regionNordic"    = "Nordic",
    "regionNon-OECD"  = "Non-OECD",
    "regionOECD (other)" = "OECD other"
  ))

p1 <- ggplot2::ggplot(attn_plot,
        ggplot2::aes(x = estimate, y = term, color = spec)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), size = 2.5) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                          position = ggplot2::position_dodge(0.5),
                          height = 0.2, linewidth = 0.7) +
  ggplot2::labs(
    title = "Region coefficients on country homework-achievement gradient: M2 vs M5 vs M3b",
    x = "Coefficient on country gradient at the country median (PISA points)",
    y = NULL, color = NULL,
    caption = "M2 = mean homework + region. M5 = M2 + Confucian heritage (collinear with East Asia by construction). M3b = M2 + five system indicators (tracking age, central exit exam, Confucian heritage, vocational share, instruction hours)."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(FIG_DIR, "G2_system_attenuation.png"),
                p1, width = 9, height = 5, dpi = 150)

panel_long <- fit_dat %>%
  tidyr::pivot_longer(c(tracking_age, central_exit_exam, confucian_heritage,
                        vocational_share_pct, instr_hrs_per_wk),
                      names_to = "indicator", values_to = "value") %>%
  dplyr::mutate(indicator = dplyr::recode(indicator,
    tracking_age          = "Tracking age",
    central_exit_exam     = "Central exit exam (0/1)",
    confucian_heritage    = "Confucian heritage (0/1)",
    vocational_share_pct  = "Vocational share (%)",
    instr_hrs_per_wk      = "Instruction hours/week"
  ))

p2 <- ggplot2::ggplot(panel_long,
        ggplot2::aes(x = value, y = slope_p50, color = region, label = country)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm",
                       color = "black", se = TRUE, linewidth = 0.6) +
  ggplot2::geom_point(size = 2, alpha = 0.85) +
  ggplot2::geom_text(size = 2.4, vjust = -0.7, alpha = 0.7, show.legend = FALSE) +
  ggplot2::facet_wrap(~ indicator, scales = "free_x", nrow = 2) +
  ggplot2::labs(
    title = "Country homework-achievement gradient by institutional indicator",
    x = NULL, y = "Country gradient at the country median (PISA points per homework hour)"
  ) +
  ggplot2::theme_minimal(base_size = 10)

ggplot2::ggsave(file.path(FIG_DIR, "G2_indicators_vs_slope.png"),
                p2, width = 13, height = 7, dpi = 150)

r2 <- function(m) summary(m)$r.squared

cat("\n[G2] === DECISION SUMMARY ===\n")
cat(sprintf("  R^2(M0)                            = %.3f\n", r2(models$M0)))
cat(sprintf("  R^2(M1: hwk_mean)                  = %.3f\n", r2(models$M1)))
cat(sprintf("  R^2(M2: hwk_mean + region)         = %.3f\n", r2(models$M2)))
cat(sprintf("  R^2(M3a: hwk_mean + 5 indicators)  = %.3f\n", r2(models$M3a)))
cat(sprintf("  R^2(M3b: + region + 5 indicators)  = %.3f\n", r2(models$M3b)))
cat(sprintf("  R^2(M4: hwk_mean + Confucian)      = %.3f\n", r2(models$M4)))
cat(sprintf("  R^2(M5: + region + Confucian)      = %.3f\n", r2(models$M5)))
cat(sprintf("  R^2(M6: 5 indicators only)         = %.3f\n", r2(models$M6)))

system_r2_gain <- r2(models$M3a) - r2(models$M1)
indicators_attn_pct <- mean(attn$attenuation_pct, na.rm = TRUE)
confucian_attn_pct  <- mean(attn_conf$attenuation_pct, na.rm = TRUE)

cat(sprintf("\n  System indicators add beyond hwk_mean alone: +%.3f pp\n",
            100 * system_r2_gain))
cat(sprintf("  System indicators add beyond hwk + region:    +%.3f pp\n",
            100 * (r2(models$M3b) - r2(models$M2))))
cat(sprintf("  Mean region attenuation M2 -> M3b (5 indic):  %.1f%%\n",
            indicators_attn_pct))
cat(sprintf("  Mean region attenuation M2 -> M5 (Confucian): %.1f%%\n",
            confucian_attn_pct))

if (system_r2_gain > 0.05 && indicators_attn_pct > 20) {
  cat("\n[G2] System indicators DO mediate the regional pattern (>5pp R^2 + >20% attenuation).\n")
  cat("     This is the publishable contrast vs PISA-school-questionnaire-only B-block result.\n")
} else if (confucian_attn_pct > 30) {
  cat("\n[G2] Confucian-heritage indicator alone substantially attenuates region.\n")
  cat("     System-level Confucian-heritage cultural-institutional features carry\n")
  cat("     a large share of the East-Asia / Nordic pattern that PISA school\n")
  cat("     questionnaire variables miss.\n")
} else {
  cat("\n[G2] System indicators do not meaningfully attenuate region either.\n")
  cat("     The descriptive null on PISA school-questionnaire variables (B-block)\n")
  cat("     extends to these crude system-level proxies.\n")
}

cat("\n[G2] Saved:\n",
    "  ", file.path(TAB_DIR, "G2_system_panel.csv"),     "\n",
    "  ", file.path(TAB_DIR, "G2_nested_models.csv"),    "\n",
    "  ", file.path(TAB_DIR, "G2_coefficients.csv"),     "\n",
    "  ", file.path(TAB_DIR, "G2_attenuation.csv"),      "\n",
    "  ", file.path(FIG_DIR, "G2_system_attenuation.png"), "\n",
    "  ", file.path(FIG_DIR, "G2_indicators_vs_slope.png"), "\n", sep = "")

cat("\n[G2] Finished at", format(Sys.time()), "\n")
