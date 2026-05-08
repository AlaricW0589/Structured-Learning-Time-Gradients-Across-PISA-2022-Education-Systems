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
source(.find_helper())
ROOT      <- proj_root()
TAB_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR   <- file.path(ROOT, "results", "figures")
LOG_DIR   <- file.path(ROOT, "logs")

source(file.path(ROOT, "helpers", "pisa_io.R"))

log_path <- file.path(LOG_DIR, "G1_within_region.log")
sink(log_path, split = TRUE)
on.exit(sink())

cat("[G1] === Within-region heterogeneity + LOO EA-Nordic gap ===\n")
cat("[G1] Started at", format(Sys.time()), "\n\n")

a2 <- readr::read_csv(file.path(TAB_DIR, "A2_position_slope.csv"),
                      show_col_types = FALSE)
frag <- readr::read_csv(file.path(TAB_DIR, "D1_robustness.csv"),
                        show_col_types = FALSE)
wf <- readr::read_csv(file.path(TAB_DIR, "F1_waterfall.csv"),
                      show_col_types = FALSE)

master <- a2 %>%
  dplyr::select(country, region, hwk_mean, slope_p25, slope_p50, slope_p75,
                concavity) %>%
  dplyr::left_join(frag %>% dplyr::select(country, fragility, base_slope, pct_change),
                   by = "country") %>%
  dplyr::mutate(
    slope_p50_pisa_pts = slope_p50,
    fragility = ifelse(is.na(fragility), "robust", fragility)
  )

math_sd <- 93.7

ea <- master %>% dplyr::filter(region == "East Asia") %>%
  dplyr::arrange(dplyr::desc(slope_p50))
nordic <- master %>% dplyr::filter(region == "Nordic") %>%
  dplyr::arrange(slope_p50)

cat("\n[G1] East Asia (n=", nrow(ea), "):\n", sep = "")
print(ea %>% dplyr::select(country, hwk_mean, slope_p25, slope_p50, slope_p75,
                            concavity, fragility), row.names = FALSE)

cat("\n[G1] Nordic (n=", nrow(nordic), "):\n", sep = "")
print(nordic %>% dplyr::select(country, hwk_mean, slope_p25, slope_p50, slope_p75,
                                concavity, fragility), row.names = FALSE)

within_region <- master %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    n_countries  = dplyr::n(),
    slope_min    = min(slope_p50, na.rm = TRUE),
    slope_max    = max(slope_p50, na.rm = TRUE),
    slope_med    = stats::median(slope_p50, na.rm = TRUE),
    slope_iqr    = stats::IQR(slope_p50, na.rm = TRUE),
    range        = slope_max - slope_min,
    .groups = "drop"
  ) %>%
  dplyr::arrange(slope_med)

cat("\n[G1] Within-region heterogeneity:\n")
print(as.data.frame(within_region), row.names = FALSE)

readr::write_csv(master,
                 file.path(TAB_DIR, "G1_within_region.csv"))

ea_countries <- ea$country
nordic_countries <- nordic$country

base_gap <- median(ea$slope_p50, na.rm = TRUE) -
            median(nordic$slope_p50, na.rm = TRUE)
base_gap_pisa <- base_gap
cat(sprintf("\n[G1] Baseline EA-Nordic median gap: %+.2f PISA points\n",
            base_gap_pisa))

loo_ea <- purrr::map_dfr(ea_countries, function(cc) {
  ea_minus <- ea %>% dplyr::filter(country != cc)
  gap <- median(ea_minus$slope_p50, na.rm = TRUE) -
         median(nordic$slope_p50, na.rm = TRUE)
  data.frame(dropped_region = "East Asia", dropped_country = cc,
             new_gap = gap, change = gap - base_gap)
})

loo_nordic <- purrr::map_dfr(nordic_countries, function(cc) {
  nordic_minus <- nordic %>% dplyr::filter(country != cc)
  gap <- median(ea$slope_p50, na.rm = TRUE) -
         median(nordic_minus$slope_p50, na.rm = TRUE)
  data.frame(dropped_region = "Nordic", dropped_country = cc,
             new_gap = gap, change = gap - base_gap)
})

loo <- dplyr::bind_rows(loo_ea, loo_nordic) %>%
  dplyr::arrange(dplyr::desc(abs(change)))

cat("\n[G1] Leave-one-country-out impact on the EA-Nordic median gap:\n")
print(loo, row.names = FALSE)

readr::write_csv(loo, file.path(TAB_DIR, "G1_loo_ea_nordic_gap.csv"))

loo_min <- min(loo$new_gap)
loo_max <- max(loo$new_gap)
cat(sprintf("\n[G1] LOO range of EA-Nordic gap: [%+.2f, %+.2f] PISA points (baseline %+.2f)\n",
            loo_min, loo_max, base_gap_pisa))

if (loo_min > 0.5 * base_gap_pisa) {
  cat("[G1] No single country drives the EA-Nordic gap.\n")
  cat("     Even excluding the most extreme country, the gap remains > 50% of baseline.\n")
} else {
  cat("[G1] The EA-Nordic gap depends materially on at least one country.\n")
  most_drivers <- loo[abs(loo$change) > 0.5 * abs(base_gap_pisa), ]
  cat("     Most influential drops:\n")
  print(most_drivers, row.names = FALSE)
}

plt_dat <- master %>%
  dplyr::filter(region %in% c("East Asia", "Nordic")) %>%
  dplyr::mutate(country = factor(country,
                                 levels = country[order(region, slope_p50)]))

p1 <- ggplot2::ggplot(plt_dat,
        ggplot2::aes(x = country, y = slope_p50, fill = region)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", slope_p50)),
                     vjust = ifelse(plt_dat$slope_p50 > 0, -0.3, 1.2),
                     size = 3.2) +
  ggplot2::facet_grid(. ~ region, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_manual(values = c("East Asia" = "#D55E00",
                                         "Nordic"    = "#56B4E9")) +
  ggplot2::labs(
    title = "Within-region heterogeneity in homework-achievement gradients",
    subtitle = "Country gradient at the country median of homework hours (PISA points per hour); PISA 2022",
    x = NULL,
    y = "Country gradient at the country median (PISA points per hour)",
    caption = "East-Asia analytic group n=6 spans MAC +4.5 to TAP +31.5; Nordic n=4 spans DNK -6.2 to SWE -30.2.\nThe East-Asia minus Nordic median gap of +23.4 PISA points survives leave-one-country-out exclusion of any member (see Section 4.6)."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

ggplot2::ggsave(file.path(FIG_DIR, "G1_within_region_heterogeneity.png"),
                p1, width = 10, height = 5, dpi = 150)

loo_plot <- loo %>%
  dplyr::mutate(label = paste0("drop ", dropped_country, " (", dropped_region, ")"),
                label = factor(label, levels = label[order(new_gap)]))

p2 <- ggplot2::ggplot(loo_plot,
        ggplot2::aes(x = new_gap, y = label, fill = dropped_region)) +
  ggplot2::geom_vline(xintercept = base_gap_pisa, linetype = "dashed",
                      color = "grey50") +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", new_gap)),
                     hjust = -0.05, size = 3.2) +
  ggplot2::scale_fill_manual(values = c("East Asia" = "#D55E00",
                                         "Nordic"    = "#56B4E9")) +
  ggplot2::labs(
    title    = "Leave-one-country-out impact on the East-Asia minus Nordic median gap",
    subtitle = sprintf("Dashed line = baseline gap (%+.1f PISA points). Each bar drops one country and recomputes the regional medians.",
                       base_gap_pisa),
    x = "EA-Nordic median gap, with the named country excluded (PISA points)",
    y = NULL, fill = "Region of dropped country"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(FIG_DIR, "G1_loo_ea_nordic_gap.png"),
                p2, width = 9, height = 5.5, dpi = 150)

cat("\n[G1] Saved:\n",
    "  ", file.path(TAB_DIR, "G1_within_region.csv"),         "\n",
    "  ", file.path(TAB_DIR, "G1_loo_ea_nordic_gap.csv"),     "\n",
    "  ", file.path(FIG_DIR, "G1_within_region_heterogeneity.png"), "\n",
    "  ", file.path(FIG_DIR, "G1_loo_ea_nordic_gap.png"),     "\n", sep = "")

cat("\n[G1] Finished at", format(Sys.time()), "\n")
