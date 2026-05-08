suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

.find_helper <- function() {
  cands <- c(file.path(getwd(), "helpers", "proj_root.R"),
             file.path(getwd(), "..", "helpers", "proj_root.R"),
             "C:/Users/Moyih/auvshk/helpers/proj_root.R")
  for (p in cands) if (file.exists(p)) return(normalizePath(p))
  stop("Cannot find helpers/proj_root.R; cwd=", getwd())
}
source(.find_helper()); ROOT <- proj_root()
CLEAN_DIR <- file.path(ROOT, "data", "clean")
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")

source(file.path(ROOT, "helpers", "pisa_io.R"))

slopes  <- readr::read_csv(file.path(TAB_DIR, "A1_country_slopes.csv"),
                           show_col_types = FALSE)
summary <- readr::read_csv(file.path(CLEAN_DIR, "country_summary.csv"),
                           show_col_types = FALSE)

slope_summary <- slopes %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    slope_p25  = slope_mean[quantile == "p25"][1],
    slope_p50  = slope_mean[quantile == "p50"][1],
    slope_p75  = slope_mean[quantile == "p75"][1],
    slope_mid  = mean(slope_mean, na.rm = TRUE),
    concavity  = slope_p25 - slope_p75,
    .groups    = "drop"
  )

merged <- slope_summary %>%
  dplyr::inner_join(summary, by = c("country" = "CNT")) %>%
  dplyr::mutate(
    region = dplyr::case_when(
      east_asia == 1 ~ "East Asia",
      nordic    == 1 ~ "Nordic",
      oecd      == 1 ~ "OECD (other)",
      TRUE           ~ "Non-OECD"
    )
  )

readr::write_csv(merged, file.path(TAB_DIR, "A2_position_slope.csv"))

cat("\n=== H1: position-slope prediction ===\n")
test_h1 <- function(slope_col, label) {
  ok <- !is.na(merged[[slope_col]]) & !is.na(merged$hwk_mean)
  r <- cor(merged$hwk_mean[ok], merged[[slope_col]][ok])
  s <- cor.test(merged$hwk_mean[ok], merged[[slope_col]][ok])
  cat(sprintf("  %-20s  r = %+.3f  (n=%d)  p=%.4g  %s\n",
              label, r, sum(ok), s$p.value,
              ifelse(r < 0, "✓ thesis-consistent",
                     ifelse(r > 0, "✗ refutes", "= flat")))
  )
  invisible(r)
}
test_h1("slope_p25", "slope at p25")
test_h1("slope_p50", "slope at median")
test_h1("slope_p75", "slope at p75")
test_h1("slope_mid", "mean across p25-p75")

cat("\n=== Concavity prevalence ===\n")
n_concave  <- sum(merged$concavity > 0,  na.rm = TRUE)
n_convex   <- sum(merged$concavity < 0,  na.rm = TRUE)
n_flat     <- sum(merged$concavity == 0, na.rm = TRUE)
cat(sprintf("  concave (slope falls): %d / %d  (%.0f%%)\n",
            n_concave, nrow(merged), 100 * n_concave / nrow(merged)))
cat(sprintf("  convex  (slope rises): %d / %d  (%.0f%%)\n",
            n_convex,  nrow(merged), 100 * n_convex  / nrow(merged)))
cat(sprintf("  flat                 : %d / %d\n", n_flat, nrow(merged)))

p1 <- ggplot2::ggplot(merged,
        ggplot2::aes(x = hwk_mean, y = slope_p50, color = region, label = country)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm",
                       color = "black", se = TRUE, linewidth = 0.6) +
  ggplot2::geom_point(size = 2.5, alpha = 0.85) +
  ggplot2::geom_text(size = 3, vjust = -0.7, alpha = 0.8, show.legend = FALSE) +
  ggplot2::labs(
    title = "Country mean homework intensity vs country homework-achievement gradient",
    subtitle = "Each dot = one PISA 2022 system; gradient = derivative of country-specific GAM at the country median of homework hours",
    x = "Country mean math homework hours per day",
    y = "Country homework-achievement gradient at the country median (PISA points per hour)"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(FIG_DIR, "A2_position_vs_slope.png"),
                p1, width = 10, height = 6.5, dpi = 150)

p2 <- ggplot2::ggplot(merged,
        ggplot2::aes(x = slope_p25, y = slope_p75, color = region, label = country)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(size = 2.2, alpha = 0.85) +
  ggplot2::geom_text(size = 2.6, vjust = -0.7, alpha = 0.8, show.legend = FALSE) +
  ggplot2::labs(
    title    = "Concavity: slope at p25 vs slope at p75 of country's homework distribution",
    subtitle = "Below dashed 45° line = concave (slope falls with higher hwk). Lower-right quadrant = within-country reversal.",
    x = "Marginal slope at country p25",
    y = "Marginal slope at country p75"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(FIG_DIR, "A2_concavity_test.png"),
                p2, width = 9, height = 8, dpi = 150)

cat("\n[A2] Saved:\n",
    "  ", file.path(TAB_DIR, "A2_position_slope.csv"), "\n",
    "  ", file.path(FIG_DIR, "A2_position_vs_slope.png"), "\n",
    "  ", file.path(FIG_DIR, "A2_concavity_test.png"), "\n", sep = "")
