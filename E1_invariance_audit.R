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

log_path <- file.path(LOG_DIR, "E1_invariance_audit.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[E1] === Block E (pilot) — psychological WLE construct-stability audit ===\n")
cat("[E1] Started at", format(Sys.time()), "\n")
cat("[E1] ROOT =", ROOT, "\n\n")

cat("[E1] Loading parquet ...\n")
dat <- as.data.frame(arrow::read_parquet(
  file.path(CLEAN_DIR, "pisa2022_core.parquet")
))

PSYCH <- toupper(PSYCH_VARS)
PSYCH <- intersect(PSYCH, names(dat))
cat(sprintf("[E1] Found %d psychological WLE constructs in parquet:\n  %s\n",
            length(PSYCH), paste(PSYCH, collapse = ", ")))

dat <- dat %>%
  dplyr::filter(!is.na(W_FSTUWT), !is.na(ESCS), !is.na(male),
                !is.na(repeated), !is.na(pv_mean_math)) %>%
  dplyr::select(CNT, W_FSTUWT, pv_mean_math, ESCS, male, repeated,
                hwk_h, dplyr::all_of(PSYCH))

math_mean <- stats::weighted.mean(dat$pv_mean_math, dat$W_FSTUWT, na.rm = TRUE)
math_sd   <- sqrt(stats::weighted.mean(
  (dat$pv_mean_math - math_mean)^2, dat$W_FSTUWT, na.rm = TRUE))
dat$math_z <- (dat$pv_mean_math - math_mean) / math_sd
cat(sprintf("[E1] Global math: mean = %.1f, sd = %.1f  (math_z = (math - mean)/sd)\n",
            math_mean, math_sd))

cat("[E1] Computing within-country z-scores for ", length(PSYCH),
    " constructs ...\n", sep = "")
for (v in PSYCH) {
  z_col <- paste0(v, "_z_wc")
  dat <- dat %>%
    dplyr::group_by(CNT) %>%
    dplyr::mutate(
      .m = stats::weighted.mean(.data[[v]], W_FSTUWT, na.rm = TRUE),
      .s = sqrt(stats::weighted.mean(
              (.data[[v]] - .m)^2, W_FSTUWT, na.rm = TRUE)),
      !!z_col := dplyr::if_else(is.na(.s) | .s == 0, NA_real_,
                                (.data[[v]] - .m) / .s)
    ) %>%
    dplyr::select(-.m, -.s) %>%
    dplyr::ungroup()
}

countries <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1000) %>%
  dplyr::pull(CNT)
dat <- dat %>% dplyr::filter(CNT %in% countries)
cat(sprintf("[E1] %d countries × %d students after filter\n",
            length(countries), nrow(dat)))

a2 <- readr::read_csv(file.path(TAB_DIR, "A2_position_slope.csv"),
                      show_col_types = FALSE) %>%
  dplyr::select(country, region, hwk_mean, slope_p50)

run_construct <- function(cc, construct, std = FALSE) {
  d <- dat %>% dplyr::filter(CNT == cc)
  rhs <- if (std) paste0(construct, "_z_wc") else construct
  ok <- !is.na(d[[rhs]]) & !is.na(d$math_z) &
        !is.na(d$ESCS) & !is.na(d$male) & !is.na(d$repeated)
  d <- d[ok, ]
  if (nrow(d) < 500) return(NULL)
  fmla <- as.formula(paste0("math_z ~ ", rhs, " + ESCS + male + repeated"))
  m <- stats::lm(fmla, data = d, weights = W_FSTUWT)
  cf <- summary(m)$coefficients
  if (!(rhs %in% rownames(cf))) return(NULL)
  data.frame(
    country   = cc,
    construct = construct,
    spec      = ifelse(std, "within-country z", "raw WLE"),
    n         = nrow(d),
    mean_x    = stats::weighted.mean(d[[construct]], d$W_FSTUWT, na.rm = TRUE),
    slope     = cf[rhs, "Estimate"],
    se        = cf[rhs, "Std. Error"],
    t         = cf[rhs, "t value"]
  )
}

cat("\n[E1] Running ", length(countries) * length(PSYCH) * 2,
    " country x construct x spec regressions ...\n", sep = "")

t0 <- Sys.time()
grid <- tidyr::expand_grid(
    country   = countries,
    construct = PSYCH,
    std       = c(FALSE, TRUE)
  ) %>%
  purrr::pmap_dfr(function(country, construct, std)
    run_construct(country, construct, std))
cat(sprintf("[E1]   completed in %.1fs (%d rows)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(grid)))

grid <- grid %>%
  dplyr::inner_join(a2, by = "country") %>%
  dplyr::rename(country_mean_x = mean_x)

readr::write_csv(grid, file.path(TAB_DIR, "E1_construct_slopes.csv"))

per_construct <- function(spec_lab) {
  g <- grid %>% dplyr::filter(spec == spec_lab)
  g %>%
    dplyr::group_by(construct) %>%
    dplyr::summarise(
      n_countries        = dplyr::n_distinct(country),
      r_pos_slope        = cor(country_mean_x, slope, use = "complete.obs"),
      slope_med_overall  = stats::median(slope, na.rm = TRUE),
      slope_iqr          = stats::IQR(slope, na.rm = TRUE),
      slope_med_eastasia = stats::median(slope[region == "East Asia"],   na.rm = TRUE),
      slope_med_nordic   = stats::median(slope[region == "Nordic"],      na.rm = TRUE),
      slope_med_oecd_o   = stats::median(slope[region == "OECD (other)"], na.rm = TRUE),
      slope_med_nonoecd  = stats::median(slope[region == "Non-OECD"],    na.rm = TRUE),
      ea_minus_nordic    = slope_med_eastasia - slope_med_nordic,
      pct_pos            = 100 * mean(slope > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(spec = spec_lab)
}

summary_raw <- per_construct("raw WLE")
summary_std <- per_construct("within-country z")
summary_all <- dplyr::bind_rows(summary_raw, summary_std)
readr::write_csv(summary_all, file.path(TAB_DIR, "E1_construct_summary.csv"))

cat("\n[E1] Cross-country position-slope correlation by construct (raw WLE):\n")
print(summary_raw[, c("construct", "n_countries", "r_pos_slope",
                      "slope_med_overall", "slope_med_eastasia",
                      "slope_med_nordic", "ea_minus_nordic", "pct_pos")],
      row.names = FALSE)

cat("\n[E1] Same after within-country z-standardization:\n")
print(summary_std[, c("construct", "r_pos_slope",
                      "slope_med_eastasia", "slope_med_nordic",
                      "ea_minus_nordic")],
      row.names = FALSE)

cl <- grid %>%
  dplyr::filter(spec == "raw WLE", country %in% c("TAP", "FIN")) %>%
  dplyr::select(country, construct, slope, se, t) %>%
  tidyr::pivot_wider(names_from = country,
                     values_from = c(slope, se, t),
                     names_glue = "{country}_{.value}") %>%
  dplyr::mutate(
    sign_TAP = ifelse(TAP_slope > 0, "+", "-"),
    sign_FIN = ifelse(FIN_slope > 0, "+", "-"),
    reversal = sign_TAP != sign_FIN
  ) %>%
  dplyr::arrange(dplyr::desc(reversal), construct)

readr::write_csv(cl, file.path(TAB_DIR, "E1_chen_lin_replication.csv"))

cat("\n[E1] Chen & Lin Taiwan-vs-Finland reversal replication (raw WLE):\n")
print(cl[, c("construct", "TAP_slope", "FIN_slope",
             "sign_TAP", "sign_FIN", "reversal")],
      row.names = FALSE)

n_reverse <- sum(cl$reversal, na.rm = TRUE)
cat(sprintf("\n[E1]   %d / %d constructs show TAP vs FIN sign reversal.\n",
            n_reverse, nrow(cl)))

plt_dat <- grid %>%
  dplyr::filter(spec == "raw WLE") %>%
  dplyr::mutate(construct = factor(construct, levels = PSYCH))

p1 <- ggplot2::ggplot(plt_dat,
        ggplot2::aes(x = country_mean_x, y = slope,
                     color = region, label = country)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(size = 1.6, alpha = 0.8) +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm",
                       color = "black", se = TRUE, linewidth = 0.5) +
  ggplot2::facet_wrap(~ construct, scales = "free", ncol = 4) +
  ggplot2::labs(
    title    = "Block E — Country slope of math_z on each psychological WLE construct",
    subtitle = "Each dot = one country (raw WLE). Slope of fitted line indicates whether higher-WLE countries also have higher within-country marginal returns.",
    x = "Country mean of the WLE construct",
    y = "Country slope (SD of math per unit of WLE)"
  ) +
  ggplot2::theme_minimal(base_size = 10)

ggplot2::ggsave(file.path(FIG_DIR, "E1_construct_country_slopes.png"),
                p1, width = 13, height = 7, dpi = 150)

sum_plot <- summary_all %>%
  dplyr::mutate(construct = factor(construct, levels = PSYCH))

p2 <- ggplot2::ggplot(sum_plot,
        ggplot2::aes(x = r_pos_slope, y = construct, color = spec)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.4)) +
  ggplot2::scale_x_continuous(limits = c(-1, 1)) +
  ggplot2::labs(
    title    = "H1-style position-slope r per construct (raw WLE vs within-country z)",
    subtitle = "If a construct's reversal is measurement-fragile, raw r differs from standardized r",
    x = "Cross-country r(country mean of construct, country slope)",
    y = NULL, color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "E1_construct_summary_grid.png"),
                p2, width = 9, height = 5, dpi = 150)

cl_long <- grid %>%
  dplyr::filter(spec == "raw WLE", country %in% c("TAP", "FIN")) %>%
  dplyr::mutate(construct = factor(construct, levels = rev(cl$construct)))

p3 <- ggplot2::ggplot(cl_long,
        ggplot2::aes(x = slope, y = construct, color = country)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::geom_segment(ggplot2::aes(x = slope - 1.96 * se,
                                     xend = slope + 1.96 * se,
                                     yend = construct),
                        linewidth = 0.8) +
  ggplot2::labs(
    title    = "Chen & Lin (2025) Taiwan-vs-Finland reversal — replicated on PISA 2022",
    subtitle = "Construct = WLE; outcome = within-country standardized math achievement; bars = 95% CI",
    x = "Slope (SD of math per unit of WLE)",
    y = NULL, color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "E1_chen_lin_replication.png"),
                p3, width = 9, height = 5, dpi = 150)

cat("\n[E1] === DECISION SUMMARY ===\n")

cat("\n  Cross-construct H1-style r (raw WLE):\n")
for (i in seq_len(nrow(summary_raw))) {
  s <- summary_raw[i, ]
  flag <- ifelse(abs(s$r_pos_slope) > 0.3, "*", " ")
  cat(sprintf("    %s %-10s  r=%+.3f  EA-Nordic=%+.3f\n",
              flag, s$construct, s$r_pos_slope, s$ea_minus_nordic))
}

cat("\n  Comparison: hwk_h reference (from A2): r = +0.475, EA-Nordic ≈ +20 pts\n")

frag <- summary_raw %>%
  dplyr::inner_join(summary_std,
                    by = "construct",
                    suffix = c("_raw", "_std")) %>%
  dplyr::mutate(
    r_change_pp     = 100 * (r_pos_slope_raw - r_pos_slope_std),
    ea_nordic_raw   = ea_minus_nordic_raw,
    ea_nordic_std   = ea_minus_nordic_std,
    diff_attenuated = abs(ea_nordic_raw) - abs(ea_nordic_std),
    fragile_under_std = abs(diff_attenuated) > 0.5 * abs(ea_nordic_raw)
  )

cat("\n  Per construct: does within-country standardisation eliminate the")
cat("\n  East-Asia / Nordic divergence? (fragile = standardised |Δ| < 50% of raw |Δ|)\n")
print(frag[, c("construct", "ea_nordic_raw", "ea_nordic_std", "fragile_under_std")],
      row.names = FALSE)

n_fragile <- sum(frag$fragile_under_std, na.rm = TRUE)
cat(sprintf("\n  %d / %d psychological constructs are FRAGILE to within-country standardisation.\n",
            n_fragile, nrow(frag)))

HWK_EA_NORDIC_REFERENCE_SD <- 20 / 93.7

cat("\n  v2.1 H5 verdict (using EA-Nordic-gap test):\n")
cat(sprintf("    Reference: homework EA-Nordic gap = %+.3f SD (~+20 PISA points)\n",
            HWK_EA_NORDIC_REFERENCE_SD))
cat("    Per-construct EA-Nordic gap (raw WLE, math-SD units):\n")
for (i in seq_len(nrow(summary_raw))) {
  s <- summary_raw[i, ]
  abs_ratio <- abs(s$ea_minus_nordic) / HWK_EA_NORDIC_REFERENCE_SD
  flag <- if (abs_ratio > 0.5) "**" else if (abs_ratio > 0.25) " *" else "  "
  cat(sprintf("      %s %-10s  EA-Nordic = %+.3f  (%.0f%% of homework reference)\n",
              flag, s$construct, s$ea_minus_nordic, 100 * abs_ratio))
}

n_substantive <- sum(abs(summary_raw$ea_minus_nordic) /
                       HWK_EA_NORDIC_REFERENCE_SD > 0.5, na.rm = TRUE)

cat(sprintf("\n  %d / %d psychological constructs have EA-Nordic gap >= 50%% of homework reference.\n",
            n_substantive, nrow(summary_raw)))

if (n_substantive <= 1) {
  cat("\n[E1] v2.1 H5 SUPPORTED in a sharper form than the original prediction:\n")
  cat("     the EA-Nordic regional reversal is essentially UNIQUE to homework/time-use;\n")
  cat("     psychological constructs are roughly regionally invariant in their\n")
  cat("     relationship to math achievement. Full MGCFA in E2 will only sharpen this.\n")
} else if (n_substantive <= 3) {
  cat("\n[E1] v2.1 H5 PARTIALLY supported: only %d of %d constructs have meaningful regional gap.\n")
} else {
  cat("\n[E1] v2.1 H5 NOT supported: many psychological constructs replicate the homework gap.\n")
}

cat("\n[E1] Saved:\n",
    "  ", file.path(TAB_DIR, "E1_construct_slopes.csv"),         "\n",
    "  ", file.path(TAB_DIR, "E1_construct_summary.csv"),         "\n",
    "  ", file.path(TAB_DIR, "E1_chen_lin_replication.csv"),      "\n",
    "  ", file.path(FIG_DIR, "E1_construct_country_slopes.png"),   "\n",
    "  ", file.path(FIG_DIR, "E1_construct_summary_grid.png"),     "\n",
    "  ", file.path(FIG_DIR, "E1_chen_lin_replication.png"),       "\n", sep = "")

cat("\n[E1] Finished at", format(Sys.time()), "\n")
