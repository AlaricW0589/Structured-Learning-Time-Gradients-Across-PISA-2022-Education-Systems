suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(intsvy)
  library(ggplot2)
  library(broom)
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
cat("[smoke] ROOT =", ROOT, "\n")

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "02_smoke_test.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[smoke] === Block A smoke test — Niu et al. 2025 replication ===\n")
cat("[smoke] Started at", format(Sys.time()), "\n\n")

dat_path <- file.path(CLEAN_DIR, "pisa2022_core.parquet")
if (!file.exists(dat_path)) {
  stop("Cleaned data not found: ", dat_path,
       "\nRun `Rscript scripts/01_ingest.R` first.")
}

cat("[smoke] Loading", dat_path, "...\n")
dat <- as.data.frame(arrow::read_parquet(dat_path))
cat(sprintf("[smoke]   %d students, %d cols\n", nrow(dat), ncol(dat)))

countries_4 <- ANCHOR_4_NIU
cat(sprintf("[smoke] Restricting to %s\n", paste(countries_4, collapse = ", ")))
sub <- dat %>% dplyr::filter(CNT %in% countries_4)
cat(sprintf("[smoke]   %d students retained\n", nrow(sub)))

cat("\n[smoke] Running per-country PV+BRR regressions...\n")
PV_MATH_UC <- toupper(paste0("pv", 1:10, "math"))

PV_MATH_FULL <- paste0("PV", 1:10, "MATH")

run_country <- function(cc) {
  cat(sprintf("[smoke]   %s ", cc))
  d <- sub %>% dplyr::filter(CNT == cc)
  fit <- intsvy::pisa.reg.pv(
    pvlabel = PV_MATH_FULL,
    x       = c("hwk_h", "ESCS", "male", "esl_native", "repeated"),
    data    = d
  )
  out <- as.data.frame(fit$reg)
  out$country  <- cc
  out$variable <- rownames(out)
  rownames(out) <- NULL
  cat("ok\n")
  out
}

raw <- purrr::map_dfr(countries_4, run_country)
hwk <- raw %>% dplyr::filter(variable == "hwk_h")
print(hwk)

readr::write_csv(hwk, file.path(TAB_DIR, "smoke_homework_4country.csv"))
cat(sprintf("\n[smoke] Saved table -> %s\n",
            file.path(TAB_DIR, "smoke_homework_4country.csv")))

decision <- hwk %>%
  dplyr::transmute(
    country = country,
    estimate = `Estimate`,
    expected_sign = dplyr::case_when(
      country %in% c("SGP", "KOR") ~ "+",
      country %in% c("FIN", "DNK") ~ "-",
      TRUE ~ "?"
    ),
    observed_sign = dplyr::case_when(
      estimate > 0 ~ "+",
      estimate < 0 ~ "-",
      TRUE ~ "0"
    ),
    matches_niu = expected_sign == observed_sign
  )

cat("\n[smoke] Sign comparison vs Niu et al. (2025):\n")
print(decision)

n_match <- sum(decision$matches_niu)
cat(sprintf("\n[smoke] %d/%d countries match the Niu sign pattern.\n",
            n_match, nrow(decision)))

if (n_match >= 3) {
  cat("[smoke] PASS — pipeline reproduces the published reversal.\n")
} else {
  cat("[smoke] FAIL — fewer than 3 countries match. Stop and debug.\n")
}

cat("\n[smoke] Plotting homework -> math binned curves...\n")

plt <- sub %>%
  dplyr::mutate(hwk_bin = round(hwk_h * 2) / 2) %>%
  dplyr::filter(!is.na(hwk_bin), !is.na(pv_mean_math)) %>%
  dplyr::group_by(CNT, hwk_bin) %>%
  dplyr::summarise(
    math = stats::weighted.mean(pv_mean_math, W_FSTUWT, na.rm = TRUE),
    n    = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n >= 30) %>%
  ggplot2::ggplot(ggplot2::aes(x = hwk_bin, y = math, color = CNT)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.6) +
  ggplot2::facet_wrap(~ CNT, ncol = 4) +
  ggplot2::labs(
    title = "Smoke test — Math homework hours vs PISA 2022 math (Niu et al. anchor 4)",
    subtitle = "Weighted mean of PV1-10 math per 0.5h homework bin (n>=30 per bin)",
    x = "Math homework hours/day (winsorized)",
    y = "Mean math plausible-value score",
    caption = "Source: PISA 2022 microdata. Author analysis."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "none")

fig_path <- file.path(FIG_DIR, "smoke_homework_curves.png")
ggplot2::ggsave(fig_path, plt, width = 10, height = 3.5, dpi = 150)
cat(sprintf("[smoke] Saved figure -> %s\n", fig_path))

cat("\n[smoke] Finished at", format(Sys.time()), "\n")
