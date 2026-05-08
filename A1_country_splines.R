suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(mgcv)
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
cat("[A1] ROOT =", ROOT, "\n")

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "A1_country_splines.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[A1] === Block A — country-by-country GAM splines ===\n")
cat("[A1] Started at", format(Sys.time()), "\n\n")

dat_path <- file.path(CLEAN_DIR, "pisa2022_core.parquet")
dat <- as.data.frame(arrow::read_parquet(dat_path))
PV_MATH_UC <- toupper(paste0("pv", 1:10, "math"))

dat <- dat %>%
  dplyr::filter(
    !is.na(hwk_h), !is.na(ESCS), !is.na(W_FSTUWT),
    !is.na(male), !is.na(repeated)
  )

keep_countries <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1000) %>%
  dplyr::pull(CNT)
cat(sprintf("[A1] %d countries with >=1000 valid students\n", length(keep_countries)))

dat <- dat %>% dplyr::filter(CNT %in% keep_countries)

fit_one <- function(d, pv_var) {
  fmla <- as.formula(sprintf(
    "%s ~ s(hwk_h, k=4, bs='cr') + ESCS + male + repeated", pv_var
  ))
  m <- mgcv::gam(fmla, data = d, weights = W_FSTUWT)
  m
}

predict_slope <- function(m, x_grid) {
  newd <- data.frame(
    hwk_h    = x_grid,
    ESCS     = mean(m$model$ESCS, na.rm = TRUE),
    male     = mean(m$model$male, na.rm = TRUE),
    repeated = mean(m$model$repeated, na.rm = TRUE)
  )
  p_lo <- predict(m, newdata = transform(newd, hwk_h = hwk_h - 0.05))
  p_hi <- predict(m, newdata = transform(newd, hwk_h = hwk_h + 0.05))
  (p_hi - p_lo) / 0.1
}

country_block <- function(cc) {
  d <- dat %>% dplyr::filter(CNT == cc)
  if (nrow(d) < 1000) return(NULL)
  qs <- stats::quantile(d$hwk_h, c(0.25, 0.5, 0.75), na.rm = TRUE)
  cat(sprintf("[A1]   %s  n=%d  hwk[p25,p50,p75]=(%.2f, %.2f, %.2f)\n",
              cc, nrow(d), qs[1], qs[2], qs[3]))
  slopes <- purrr::map_dfr(PV_MATH_UC, function(pv) {
    m  <- fit_one(d, pv)
    s  <- predict_slope(m, qs)
    data.frame(quantile = c("p25","p50","p75"), slope = s, hwk_at = qs)
  })
  pooled <- slopes %>%
    dplyr::group_by(quantile, hwk_at) %>%
    dplyr::summarise(
      slope_mean = mean(slope),
      between    = var(slope),
      M          = dplyr::n(),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(country = cc)
  pooled
}

cat("\n[A1] Fitting GAMs (this takes a while — ~10s per country) ...\n")
res <- purrr::map_dfr(keep_countries, country_block)

readr::write_csv(res, file.path(TAB_DIR, "A1_country_slopes.csv"))
cat(sprintf("\n[A1] Saved table -> %s  (%d rows)\n",
            file.path(TAB_DIR, "A1_country_slopes.csv"), nrow(res)))

cat("\n[A1] Plotting anchor-8 spline panel ...\n")

anchor <- dat %>% dplyr::filter(CNT %in% ANCHOR_8_PAPER)

build_curve <- function(cc) {
  d <- anchor %>% dplyr::filter(CNT == cc)
  if (nrow(d) < 500) return(NULL)
  m <- mgcv::gam(PV1MATH ~ s(hwk_h, k = 4, bs = "cr") + ESCS + male + repeated,
                 data = d, weights = W_FSTUWT)
  grid <- seq(stats::quantile(d$hwk_h, 0.02, na.rm = TRUE),
              stats::quantile(d$hwk_h, 0.98, na.rm = TRUE),
              length.out = 60)
  newd <- data.frame(
    hwk_h    = grid,
    ESCS     = mean(d$ESCS, na.rm = TRUE),
    male     = mean(d$male, na.rm = TRUE),
    repeated = mean(d$repeated, na.rm = TRUE)
  )
  pred <- predict(m, newdata = newd, se.fit = TRUE)
  data.frame(
    CNT = cc, hwk = grid,
    fit = pred$fit, lo = pred$fit - 1.96*pred$se.fit,
    hi  = pred$fit + 1.96*pred$se.fit
  )
}

curves <- purrr::map_dfr(ANCHOR_8_PAPER, build_curve)

p <- ggplot2::ggplot(curves, ggplot2::aes(x = hwk, y = fit)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~ CNT, scales = "free_y", ncol = 4) +
  ggplot2::labs(
    title = "Block A — Country-specific math homework -> math achievement curves",
    subtitle = "GAM (cubic-spline) on PV1MATH, anchor-8 countries, PISA 2022",
    x = "Math homework hours/day (winsorized)",
    y = "Predicted PV1 math score",
    caption = "Source: PISA 2022 microdata. ESCS, gender, grade-repetition held at sample means.\nFor inferential SE bands use Block C output."
  ) +
  ggplot2::theme_minimal(base_size = 11)

fig_path <- file.path(FIG_DIR, "A1_country_splines.png")
ggplot2::ggsave(fig_path, p, width = 11, height = 6, dpi = 150)
cat(sprintf("[A1] Saved figure -> %s\n", fig_path))

cat("\n[A1] Finished at", format(Sys.time()), "\n")
