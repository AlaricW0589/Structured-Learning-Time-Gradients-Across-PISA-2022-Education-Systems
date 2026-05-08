suppressPackageStartupMessages({
  library(intsvy)
  library(dplyr)
})

set.seed(42)

build <- function(cc, slope) {
  n <- 2000
  hwk  <- pmax(0, rgamma(n, 2, 1))
  escs <- rnorm(n)
  male <- rbinom(n, 1, 0.5)
  base <- list(SGP = 520, KOR = 510, FIN = 480, DNK = 470)[[cc]]
  mu   <- base + slope * hwk + 30 * escs - 5 * male
  pv   <- replicate(10, mu + rnorm(n, 0, 70))
  colnames(pv) <- paste0("PV", 1:10, "MATH")
  fw   <- runif(n, 0.5, 2)
  rw   <- replicate(80, fw * runif(n, 0.5, 1.5))
  colnames(rw) <- paste0("W_FSTURWT", 1:80)
  data.frame(
    CNT        = cc,
    CNTSCHID   = sample(100:300, n, replace = TRUE),
    W_FSTUWT   = fw,
    hwk_h      = hwk,
    ESCS       = escs,
    male       = male,
    esl_native = rbinom(n, 1, 0.9),
    repeated   = rbinom(n, 1, 0.05),
    pv,
    rw
  )
}

df <- dplyr::bind_rows(
  build("SGP", +6),
  build("KOR", +4),
  build("FIN", -3),
  build("DNK", -2)
)
cat(sprintf("synthetic dataset: %d rows, %d cols\n", nrow(df), ncol(df)))

PV_MATH_FULL <- paste0("PV", 1:10, "MATH")

run_country <- function(cc) {
  d <- df[df$CNT == cc, ]
  fit <- intsvy::pisa.reg.pv(
    pvlabel = PV_MATH_FULL,
    x       = c("hwk_h", "ESCS", "male", "esl_native", "repeated"),
    data    = d
  )
  out <- as.data.frame(fit$reg)
  out$country  <- cc
  out$variable <- rownames(out)
  rownames(out) <- NULL
  out[out$variable == "hwk_h",
      c("country", "variable", "Estimate", "Std. Error", "t value")]
}

cat("\nRunning intsvy::pisa.reg.pv on each fake country...\n\n")
hwk <- dplyr::bind_rows(lapply(c("SGP", "KOR", "FIN", "DNK"), run_country))
print(hwk)

decision <- hwk %>%
  dplyr::mutate(
    expected_sign = ifelse(country %in% c("SGP", "KOR"), "+", "-"),
    observed_sign = ifelse(Estimate > 0, "+", "-"),
    matches_niu   = expected_sign == observed_sign
  )
cat("\nSign comparison vs synthetic ground truth:\n")
print(decision)

if (sum(decision$matches_niu) == 4L) {
  cat("\nPASS: synthetic pilot reproduces the embedded reversal.\n")
  cat("      Pipeline (intsvy + helpers) works end-to-end.\n")
} else {
  cat("\nFAIL: synthetic pilot did NOT recover the reversal -- DEBUG.\n")
}
