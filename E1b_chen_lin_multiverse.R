suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(ggplot2); library(sandwich); library(lmtest)
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

log_path <- file.path(LOG_DIR, "E1b_chen_lin_multiverse.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[E1b] === Chen & Lin (2025) replication multiverse ===\n")
cat("[E1b] Started at", format(Sys.time()), "\n\n")

dat <- as.data.frame(arrow::read_parquet(
  file.path(CLEAN_DIR, "pisa2022_core.parquet")
))

dat <- dat %>%
  dplyr::filter(!is.na(W_FSTUWT), !is.na(ESCS), !is.na(male),
                !is.na(repeated), !is.na(SDLEFF))

math_mean <- stats::weighted.mean(dat$pv_mean_math, dat$W_FSTUWT, na.rm = TRUE)
math_sd   <- sqrt(stats::weighted.mean(
  (dat$pv_mean_math - math_mean)^2, dat$W_FSTUWT, na.rm = TRUE))
dat$math_z <- (dat$pv_mean_math - math_mean) / math_sd

dat <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::mutate(
    .m = stats::weighted.mean(SDLEFF, W_FSTUWT, na.rm = TRUE),
    .s = sqrt(stats::weighted.mean(
            (SDLEFF - .m)^2, W_FSTUWT, na.rm = TRUE)),
    SDLEFF_z_wc = dplyr::if_else(is.na(.s) | .s == 0, NA_real_,
                                 (SDLEFF - .m) / .s)
  ) %>%
  dplyr::select(-.m, -.s) %>%
  dplyr::ungroup()

m1_mean <- stats::weighted.mean(dat$PV1MATH, dat$W_FSTUWT, na.rm = TRUE)
m1_sd   <- sqrt(stats::weighted.mean(
  (dat$PV1MATH - m1_mean)^2, dat$W_FSTUWT, na.rm = TRUE))
dat$pv1_z <- (dat$PV1MATH - m1_mean) / m1_sd

cat(sprintf("[E1b] Total students (with non-missing SDLEFF): %d\n", nrow(dat)))
cat(sprintf("[E1b] Taiwan students: %d\n", sum(dat$CNT == "TAP")))
cat(sprintf("[E1b] Finland students: %d\n", sum(dat$CNT == "FIN")))

PV_MATH_FULL <- paste0("PV", 1:10, "MATH")

fit_one <- function(d, outcome, regressor, controls,
                    use_pv = FALSE, cluster_var = NULL) {
  controls <- intersect(controls, names(d))
  rhs_terms <- c(regressor, controls)
  fmla <- as.formula(paste(outcome, "~", paste(rhs_terms, collapse = " + ")))

  if (use_pv) {
    fits <- lapply(PV_MATH_FULL, function(pv) {
      dd <- d
      pv_g_mean <- stats::weighted.mean(dd[[pv]], dd$W_FSTUWT, na.rm = TRUE)
      pv_g_sd   <- sqrt(stats::weighted.mean(
        (dd[[pv]] - pv_g_mean)^2, dd$W_FSTUWT, na.rm = TRUE))
      dd$.outcome_z <- (dd[[pv]] - pv_g_mean) / pv_g_sd
      f2 <- as.formula(paste(".outcome_z ~", paste(rhs_terms, collapse = " + ")))
      stats::lm(f2, data = dd, weights = W_FSTUWT)
    })
    coefs <- sapply(fits, function(m) coef(m)[regressor])
    ses   <- sapply(fits, function(m) {
      if (!is.null(cluster_var) && cluster_var %in% names(d)) {
        sqrt(diag(sandwich::vcovCL(m, cluster = d[[cluster_var]])))[regressor]
      } else {
        sqrt(diag(vcov(m)))[regressor]
      }
    })
    theta_bar <- mean(coefs)
    U_bar     <- mean(ses^2)
    B         <- if (length(coefs) > 1) var(coefs) else 0
    se_total  <- sqrt(U_bar + (1 + 1 / length(coefs)) * B)
    return(data.frame(estimate = theta_bar, se = se_total,
                      n = length(fits[[1]]$residuals)))
  } else {
    m <- stats::lm(fmla, data = d, weights = W_FSTUWT)
    cf <- summary(m)$coefficients
    if (!(regressor %in% rownames(cf)))
      return(data.frame(estimate = NA, se = NA, n = nrow(d)))
    if (!is.null(cluster_var) && cluster_var %in% names(d)) {
      vc <- sandwich::vcovCL(m, cluster = d[[cluster_var]])
      se_cl <- sqrt(diag(vc))[regressor]
      return(data.frame(estimate = cf[regressor, "Estimate"],
                        se = se_cl, n = nrow(m$model)))
    }
    data.frame(estimate = cf[regressor, "Estimate"],
               se       = cf[regressor, "Std. Error"],
               n        = nrow(m$model))
  }
}

run_spec <- function(spec) {
  d <- dat %>% dplyr::filter(CNT %in% spec$countries)

  ctrl <- if (spec$controls == "reduced")  c("ESCS", "male") else
          if (spec$controls == "standard") c("ESCS", "male", "repeated") else
          character()

  reg <- if (spec$wle == "raw") "SDLEFF" else "SDLEFF_z_wc"

  out <- spec$countries %>% purrr::map_dfr(function(cc) {
    dd <- d %>% dplyr::filter(CNT == cc)
    if (nrow(dd) < 200) return(NULL)
    if (spec$pv == "averaged") {
      r <- fit_one(dd, outcome = "math_z", regressor = reg,
                   controls = ctrl, use_pv = FALSE,
                   cluster_var = spec$cluster_var)
    } else if (spec$pv == "pv1") {
      r <- fit_one(dd, outcome = "pv1_z", regressor = reg,
                   controls = ctrl, use_pv = FALSE,
                   cluster_var = spec$cluster_var)
    } else {
      r <- fit_one(dd, outcome = "math_z", regressor = reg,
                   controls = ctrl, use_pv = TRUE,
                   cluster_var = spec$cluster_var)
    }
    r$country <- cc
    r
  })

  out
}

specs <- list(
  list(label = "S1: Chen-Lin original (TAP, FIN, averaged PV, raw, reduced controls)",
       countries = c("TAP", "FIN"), pv = "pv1", wle = "raw",
       controls = "reduced", cluster_var = NULL),
  list(label = "S2: Chen-Lin + cluster-robust SE on school",
       countries = c("TAP", "FIN"), pv = "pv1", wle = "raw",
       controls = "reduced", cluster_var = "CNTSCHID"),
  list(label = "S3: Chen-Lin + Rubin PV pooling",
       countries = c("TAP", "FIN"), pv = "rubin", wle = "raw",
       controls = "reduced", cluster_var = NULL),
  list(label = "S4: Chen-Lin + Rubin PV + cluster-robust",
       countries = c("TAP", "FIN"), pv = "rubin", wle = "raw",
       controls = "reduced", cluster_var = "CNTSCHID"),
  list(label = "S5: Within-country z of SDLEFF",
       countries = c("TAP", "FIN"), pv = "rubin", wle = "z_wc",
       controls = "reduced", cluster_var = "CNTSCHID"),
  list(label = "S6: Standard controls (+ repeated) instead of reduced",
       countries = c("TAP", "FIN"), pv = "rubin", wle = "raw",
       controls = "standard", cluster_var = "CNTSCHID")
)

cat("\n[E1b] Running ", length(specs), " specifications x 2 countries ...\n",
    sep = "")

mv <- purrr::imap_dfr(specs, function(s, idx) {
  res <- run_spec(s)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$spec_id    <- paste0("S", idx)
  res$spec_label <- s$label
  res
})

mv <- mv %>%
  dplyr::transmute(spec_id, spec_label, country,
                   estimate = round(estimate, 4),
                   se = round(se, 4),
                   t = round(estimate / se, 2),
                   n)
print(mv, row.names = FALSE)

readr::write_csv(mv, file.path(TAB_DIR, "E1b_chen_lin_multiverse.csv"))

reversal_check <- mv %>%
  tidyr::pivot_wider(id_cols = c(spec_id, spec_label),
                     names_from = country,
                     values_from = c(estimate, se)) %>%
  dplyr::mutate(
    sign_TAP = ifelse(estimate_TAP > 0, "+", "-"),
    sign_FIN = ifelse(estimate_FIN > 0, "+", "-"),
    chen_lin_pattern = sign_TAP == "+" & sign_FIN == "-"
  )

cat("\n[E1b] Per specification: does Chen-Lin (TAP+, FIN-) sign pattern appear?\n\n")
print(reversal_check %>%
        dplyr::select(spec_id, estimate_TAP, estimate_FIN,
                      sign_TAP, sign_FIN, chen_lin_pattern),
      row.names = FALSE)

n_replicates <- sum(reversal_check$chen_lin_pattern, na.rm = TRUE)
cat(sprintf("\n[E1b] Chen-Lin sign pattern (TAP+, FIN-) appears in %d / %d specs.\n",
            n_replicates, nrow(reversal_check)))

plt <- mv %>%
  dplyr::mutate(spec_id = factor(spec_id, levels = paste0("S", 1:length(specs))))

p1 <- ggplot2::ggplot(plt,
        ggplot2::aes(x = estimate, y = spec_id, color = country)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_point(size = 3.5,
                      position = ggplot2::position_dodge(0.4)) +
  ggplot2::geom_segment(ggplot2::aes(x = estimate - 1.96 * se,
                                     xend = estimate + 1.96 * se,
                                     yend = spec_id),
                        linewidth = 0.8,
                        position = ggplot2::position_dodge(0.4)) +
  ggplot2::labs(
    title    = "Chen & Lin (2025) replication multiverse — TAP vs FIN, 6 specifications",
    subtitle = "If the sign-reversal claim is methodological, only S1 should reproduce it (Chen-Lin estimation choices).",
    x = "Slope of math_z on SDLEFF (raw or within-country z)",
    y = NULL, color = NULL,
    caption = "S1 = Chen-Lin original; S2 + cluster-robust; S3 + Rubin PV pooling; S4 = S2+S3 combined; S5 = within-country z SDLEFF; S6 = + repetition control."
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(file.path(FIG_DIR, "E1b_chen_lin_multiverse.png"),
                p1, width = 11, height = 5, dpi = 150)

cat("\n[E1b] Saved:\n",
    "  ", file.path(TAB_DIR, "E1b_chen_lin_multiverse.csv"), "\n",
    "  ", file.path(FIG_DIR, "E1b_chen_lin_multiverse.png"), "\n", sep = "")
cat("\n[E1b] Finished at", format(Sys.time()), "\n")
