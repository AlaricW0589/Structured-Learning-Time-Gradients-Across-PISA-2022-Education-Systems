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

log_path <- file.path(LOG_DIR, "E3_behavioural_time_audit.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[E3] === Behavioural-time audit (response to CER Round 1 attack #1) ===\n")
cat("[E3] Started at", format(Sys.time()), "\n\n")

dat <- as.data.frame(arrow::read_parquet(
  file.path(CLEAN_DIR, "pisa2022_core.parquet")
))

dat <- dat %>%
  dplyr::filter(!is.na(W_FSTUWT), !is.na(ESCS), !is.na(male),
                !is.na(repeated), !is.na(pv_mean_math))

math_mean <- stats::weighted.mean(dat$pv_mean_math, dat$W_FSTUWT, na.rm = TRUE)
math_sd   <- sqrt(stats::weighted.mean(
  (dat$pv_mean_math - math_mean)^2, dat$W_FSTUWT, na.rm = TRUE))
dat$math_z <- (dat$pv_mean_math - math_mean) / math_sd
cat(sprintf("[E3] math global: mean=%.1f, sd=%.1f\n", math_mean, math_sd))

countries <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1000) %>%
  dplyr::pull(CNT)
dat <- dat %>% dplyr::filter(CNT %in% countries)

a2 <- readr::read_csv(file.path(TAB_DIR, "A2_position_slope.csv"),
                      show_col_types = FALSE) %>%
  dplyr::select(country, region)

INPUTS <- list(
  list(var = "hwk_h",       label = "Math homework hours (reference)"),
  list(var = "study_h",     label = "Before/after-school study (STUDYHMW)"),
  list(var = "EXERPRAC",    label = "Exercise/practice frequency"),
  list(var = "ST297Q01JA",  label = "One-on-one tutoring (frequency)"),
  list(var = "ST297Q03JA",  label = "Internet/computer tutoring"),
  list(var = "ST297Q05JA",  label = "Asynchronous video instruction"),
  list(var = "ST297Q06JA",  label = "Small-group study (2-7)"),
  list(var = "ST297Q07JA",  label = "Large-group study (>=8)")
)
INPUTS <- INPUTS[sapply(INPUTS, function(x) x$var %in% names(dat))]
cat(sprintf("[E3] Available behavioural-time inputs: %d\n", length(INPUTS)))
for (it in INPUTS) cat(sprintf("  - %s : %s\n", it$var, it$label))

run_country_input <- function(cc, input) {
  d <- dat %>% dplyr::filter(CNT == cc)
  d <- d[!is.na(d[[input]]), , drop = FALSE]
  if (nrow(d) < 500) return(NULL)
  fmla <- as.formula(paste("math_z ~", input, "+ ESCS + male + repeated"))
  m <- tryCatch(stats::lm(fmla, data = d, weights = W_FSTUWT),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients
  if (!(input %in% rownames(cf))) return(NULL)
  data.frame(
    country = cc, input = input,
    slope = cf[input, "Estimate"],
    se    = cf[input, "Std. Error"],
    n     = nrow(d)
  )
}

cat("\n[E3] Running ", length(countries) * length(INPUTS),
    " country x input regressions ...\n", sep = "")
t0 <- Sys.time()
grid <- tidyr::expand_grid(
    country = countries,
    input   = sapply(INPUTS, function(x) x$var)
  ) %>%
  purrr::pmap_dfr(function(country, input)
    run_country_input(country, input))
cat(sprintf("[E3]   completed in %.1fs (%d rows)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(grid)))

grid <- grid %>%
  dplyr::inner_join(a2, by = "country")

readr::write_csv(grid, file.path(TAB_DIR, "E3_behavioural_slopes.csv"))

region_med <- grid %>%
  dplyr::group_by(input, region) %>%
  dplyr::summarise(
    n_countries = dplyr::n_distinct(country),
    slope_med   = stats::median(slope, na.rm = TRUE),
    slope_p25   = stats::quantile(slope, 0.25, na.rm = TRUE),
    slope_p75   = stats::quantile(slope, 0.75, na.rm = TRUE),
    .groups     = "drop"
  )

ea_nordic_gap <- region_med %>%
  dplyr::filter(region %in% c("East Asia", "Nordic")) %>%
  dplyr::select(input, region, slope_med) %>%
  tidyr::pivot_wider(names_from = region, values_from = slope_med) %>%
  dplyr::mutate(
    EA_minus_Nordic = `East Asia` - Nordic,
    EA_minus_Nordic_pisa_pts = EA_minus_Nordic * math_sd
  ) %>%
  dplyr::arrange(dplyr::desc(abs(EA_minus_Nordic)))

readr::write_csv(ea_nordic_gap,
                 file.path(TAB_DIR, "E3_ea_nordic_gaps.csv"))

cat("\n[E3] East-Asia minus Nordic gap per behavioural-time input:\n")
print(as.data.frame(ea_nordic_gap), row.names = FALSE, digits = 3)

plt <- ea_nordic_gap %>%
  dplyr::mutate(
    input = factor(input, levels = ea_nordic_gap$input),
    bar_col = ifelse(abs(EA_minus_Nordic) > 0.21 * 0.5,
                     "substantial", "minor")
  )

p1 <- ggplot2::ggplot(plt,
        ggplot2::aes(x = EA_minus_Nordic_pisa_pts, y = input, fill = bar_col)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_vline(xintercept = c(-20, 20), linetype = "dotted",
                      color = "grey50") +
  ggplot2::geom_col(alpha = 0.85) +
  ggplot2::scale_fill_manual(values = c("substantial" = "#D55E00",
                                         "minor"       = "#999999")) +
  ggplot2::labs(
    title    = "Block E3 — Behavioural-time inputs: East-Asia minus Nordic median country slope",
    subtitle = "Dotted lines = ±50% of the homework reference (~±10 PISA points). Bars > dotted = substantive regional gap.",
    x = "East-Asia − Nordic median slope (PISA points per unit of input)",
    y = NULL, fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(FIG_DIR, "E3_behavioural_layer_stability.png"),
                p1, width = 9, height = 5.2, dpi = 150)

cat("\n[E3] === DECISION SUMMARY ===\n")
ref <- ea_nordic_gap$EA_minus_Nordic[ea_nordic_gap$input == "hwk_h"]
cat(sprintf("[E3] Reference (homework): EA-Nordic gap = %+.3f SD = %+.0f PISA points\n",
            ref, ref * math_sd))
cat("\n[E3] Per behavioural-time input, gap as a share of homework reference:\n")
for (i in seq_len(nrow(ea_nordic_gap))) {
  s <- ea_nordic_gap[i, ]
  pct <- abs(s$EA_minus_Nordic) / abs(ref) * 100
  flag <- if (pct > 50) "**" else if (pct > 25) " *" else "  "
  cat(sprintf("  %s %-12s  %+.3f SD (%+.0f pts)  = %3.0f%% of homework reference\n",
              flag, s$input, s$EA_minus_Nordic,
              s$EA_minus_Nordic_pisa_pts, pct))
}

n_substantial <- sum(abs(ea_nordic_gap$EA_minus_Nordic) >= 0.5 * abs(ref))
cat(sprintf("\n[E3] %d / %d behavioural-time inputs have EA-Nordic gap >= 50%% of homework reference.\n",
            n_substantial, nrow(ea_nordic_gap)))

if (n_substantial == 1) {
  cat("[E3] Only homework shows the substantial regional gap. Claim 'uniquely homework-time' is supported.\n")
} else if (n_substantial <= 3) {
  cat(sprintf("[E3] %d behavioural inputs replicate the gap. Claim 'concentrated in structured study time' is supported.\n",
              n_substantial))
} else {
  cat("[E3] Many behavioural inputs replicate the gap. Claim 'time-use broadly' is supported but 'unique' is not.\n")
}

cat("\n[E3] Saved:\n",
    "  ", file.path(TAB_DIR, "E3_behavioural_slopes.csv"),     "\n",
    "  ", file.path(TAB_DIR, "E3_ea_nordic_gaps.csv"),         "\n",
    "  ", file.path(FIG_DIR, "E3_behavioural_layer_stability.png"), "\n", sep = "")
cat("\n[E3] Finished at", format(Sys.time()), "\n")
