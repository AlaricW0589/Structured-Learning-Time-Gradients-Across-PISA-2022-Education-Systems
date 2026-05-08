suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(purrr)
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
CLEAN_DIR <- file.path(ROOT, "data", "clean")
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")
LOG_DIR   <- file.path(ROOT, "logs")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "D1_selection_bounds.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[D1] === Block D — selection-bound robustness ===\n")
cat("[D1] Started at", format(Sys.time()), "\n")
cat("[D1] ROOT =", ROOT, "\n\n")

cat("[D1] Loading parquet ...\n")
dat <- as.data.frame(arrow::read_parquet(
  file.path(CLEAN_DIR, "pisa2022_core.parquet")
))

dat <- dat %>%
  dplyr::filter(
    !is.na(hwk_h), !is.na(ESCS), !is.na(W_FSTUWT),
    !is.na(male), !is.na(repeated), !is.na(pv_mean_math)
  )

countries <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1000) %>%
  dplyr::pull(CNT)

dat <- dat %>% dplyr::filter(CNT %in% countries)
cat(sprintf("[D1] %d countries x %d students\n", length(countries), nrow(dat)))

a2 <- readr::read_csv(file.path(TAB_DIR, "A2_position_slope.csv"),
                      show_col_types = FALSE) %>%
  dplyr::select(country, region, hwk_mean, slope_p50)

TRIM_LEVELS <- c(0, 5, 10, 15, 20, 25, 30)

run_country_trim <- function(cc, trim_pct) {
  d <- dat %>% dplyr::filter(CNT == cc)
  if (trim_pct > 0) {
    cutoff <- stats::quantile(d$ESCS, probs = trim_pct / 100, na.rm = TRUE)
    d <- d %>% dplyr::filter(ESCS >= cutoff)
  }
  if (nrow(d) < 500) return(NULL)
  m <- stats::lm(pv_mean_math ~ hwk_h + ESCS + male + repeated,
                 data = d, weights = W_FSTUWT)
  cf <- summary(m)$coefficients
  data.frame(
    country  = cc,
    trim_pct = trim_pct,
    n        = nrow(d),
    slope    = cf["hwk_h", "Estimate"],
    se       = cf["hwk_h", "Std. Error"],
    t        = cf["hwk_h", "t value"]
  )
}

cat("\n[D1] Running ", length(countries) * length(TRIM_LEVELS),
    " country x trim regressions ...\n", sep = "")

t0 <- Sys.time()
grid <- tidyr::expand_grid(country = countries, trim_pct = TRIM_LEVELS) %>%
  purrr::pmap_dfr(function(country, trim_pct)
    run_country_trim(country, trim_pct))
cat(sprintf("[D1]   completed in %.1fs\n",
            as.numeric(Sys.time() - t0, units = "secs")))

grid <- grid %>%
  dplyr::inner_join(a2, by = "country")

readr::write_csv(grid, file.path(TAB_DIR, "D1_trim_grid.csv"))
cat(sprintf("[D1]   saved %d-row grid -> %s\n",
            nrow(grid), file.path(TAB_DIR, "D1_trim_grid.csv")))

baseline <- grid %>%
  dplyr::filter(trim_pct == 0) %>%
  dplyr::transmute(country, base_slope = slope, base_t = t)

frag <- grid %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    n_levels      = dplyr::n(),
    sign_changes  = sum(diff(sign(slope[order(trim_pct)])) != 0),
    sign_flipped  = any(sign(slope) != sign(slope[trim_pct == 0])),
    max_abs_chg   = max(abs(slope[trim_pct == 30] - slope[trim_pct == 0])),
    .groups       = "drop"
  ) %>%
  dplyr::inner_join(baseline, by = "country") %>%
  dplyr::mutate(
    pct_change = ifelse(abs(base_slope) < 1e-6, NA,
                        100 * max_abs_chg / abs(base_slope)),
    fragility = dplyr::case_when(
      sign_changes >= 2                          ~ "noisy",
      sign_flipped                               ~ "flip",
      !is.na(pct_change) & pct_change >= 50      ~ "shrink",
      TRUE                                       ~ "robust"
    )
  ) %>%
  dplyr::inner_join(a2, by = "country") %>%
  dplyr::arrange(fragility, country)

readr::write_csv(frag, file.path(TAB_DIR, "D1_robustness.csv"))

cat("\n[D1] Fragility classification:\n")
print(table(frag$fragility))
cat("\n  By region:\n")
print(table(frag$region, frag$fragility))

region_summary <- grid %>%
  dplyr::group_by(region, trim_pct) %>%
  dplyr::summarise(
    n_countries = dplyr::n_distinct(country),
    slope_mean  = mean(slope, na.rm = TRUE),
    slope_med   = stats::median(slope, na.rm = TRUE),
    slope_p25   = stats::quantile(slope, 0.25, na.rm = TRUE),
    slope_p75   = stats::quantile(slope, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(region_summary,
                 file.path(TAB_DIR, "D1_region_summary.csv"))
cat("\n[D1] Region slope means at each trim level:\n")
print(region_summary %>%
        dplyr::select(region, trim_pct, n_countries, slope_mean, slope_med),
      row.names = FALSE)

anchor <- grid %>% dplyr::filter(country %in% ANCHOR_8_PAPER)
p1 <- ggplot2::ggplot(anchor,
        ggplot2::aes(x = trim_pct, y = slope, color = country)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = slope - 1.96 * se,
                                     ymax = slope + 1.96 * se,
                                     fill = country),
                       alpha = 0.12, color = NA) +
  ggplot2::facet_wrap(~ country, ncol = 4, scales = "free_y") +
  ggplot2::labs(
    title    = "Block D — country slope of math hwk under ESCS-trimming sensitivity",
    subtitle = "Drop the bottom k% of students by ESCS within each country, refit (anchor-8)",
    x = "Trim level (% of lowest-ESCS students dropped)",
    y = "Marginal slope: pv_mean_math ~ hwk_h",
    caption = "Source: PISA 2022. Controls: ESCS, gender, repetition. Weighted by W_FSTUWT.\n95% CIs use single-PV SE; full PV+BRR pooling left to Block C."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(file.path(FIG_DIR, "D1_anchor8_trim_curves.png"),
                p1, width = 11, height = 6, dpi = 150)

p2 <- ggplot2::ggplot(region_summary,
        ggplot2::aes(x = trim_pct, y = slope_mean, color = region, group = region)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = slope_p25, ymax = slope_p75, fill = region),
                       alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    title    = "Block D — Region average slope under ESCS-trimming sensitivity",
    subtitle = "Solid line = region mean of country slopes; ribbon = country-level p25-p75 within region",
    x = "Trim level (% of lowest-ESCS students dropped)",
    y = "Country slope: pv_mean_math ~ hwk_h"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "D1_region_trim_summary.png"),
                p2, width = 9, height = 5.5, dpi = 150)

p3 <- ggplot2::ggplot(frag,
        ggplot2::aes(x = region, fill = fragility)) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(
    title    = "Slope-sign fragility under ESCS-trimming, by region",
    subtitle = "Stacked share of countries by fragility class",
    x = NULL, y = "Share of region", fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "D1_fragility_distribution.png"),
                p3, width = 8, height = 4.5, dpi = 150)

cat("\n[D1] === DECISION SUMMARY ===\n")

n_total      <- nrow(frag)
n_robust     <- sum(frag$fragility == "robust")
n_shrink     <- sum(frag$fragility == "shrink")
n_flip       <- sum(frag$fragility == "flip")
n_noisy      <- sum(frag$fragility == "noisy")

cat(sprintf("  Robust      (sign stable, |Δ|<50%%) : %2d / %2d (%.0f%%)\n",
            n_robust, n_total, 100 * n_robust / n_total))
cat(sprintf("  Shrink      (sign stable, |Δ|>=50%%): %2d / %2d (%.0f%%)\n",
            n_shrink, n_total, 100 * n_shrink / n_total))
cat(sprintf("  Flip        (sign flips somewhere) : %2d / %2d (%.0f%%)\n",
            n_flip,   n_total, 100 * n_flip   / n_total))
cat(sprintf("  Noisy       (>=2 sign changes)     : %2d / %2d (%.0f%%)\n",
            n_noisy,  n_total, 100 * n_noisy  / n_total))

fragile_share <- frag %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    n_fragile = sum(fragility %in% c("flip", "shrink", "noisy")),
    pct_fragile = 100 * n_fragile / n,
    .groups   = "drop"
  )
cat("\n  Fragile-share by region (proxy for selection sensitivity):\n")
print(fragile_share, row.names = FALSE)

ne <- fragile_share %>% dplyr::filter(region == "Non-OECD")
oa <- fragile_share %>%
  dplyr::filter(region %in% c("OECD (other)", "Nordic", "East Asia"))
ne_pct <- if (nrow(ne) == 1) ne$pct_fragile else NA_real_
oa_pct <- sum(oa$n_fragile) / sum(oa$n) * 100
if (!is.na(ne_pct) && ne_pct > oa_pct + 10) {
  cat(sprintf("\n[D1] H6 SUPPORTED: Non-OECD %.0f%% fragile vs OECD/EastAsia/Nordic %.0f%% fragile.\n",
              ne_pct, oa_pct))
} else {
  cat(sprintf("\n[D1] H6 NOT clearly supported: Non-OECD %.0f%% fragile vs other %.0f%%.\n",
              ne_pct, oa_pct))
}

cat("\n[D1] Saved:\n",
    "  ", file.path(TAB_DIR, "D1_trim_grid.csv"),       "\n",
    "  ", file.path(TAB_DIR, "D1_robustness.csv"),       "\n",
    "  ", file.path(TAB_DIR, "D1_region_summary.csv"),   "\n",
    "  ", file.path(FIG_DIR, "D1_anchor8_trim_curves.png"),    "\n",
    "  ", file.path(FIG_DIR, "D1_region_trim_summary.png"),    "\n",
    "  ", file.path(FIG_DIR, "D1_fragility_distribution.png"), "\n", sep = "")

cat("\n[D1] Finished at", format(Sys.time()), "\n")
