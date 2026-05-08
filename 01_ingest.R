suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(arrow)
  library(readr)
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
RAW_DIR   <- file.path(ROOT, "data", "raw")
CLEAN_DIR <- file.path(ROOT, "data", "clean")
LOG_DIR   <- file.path(ROOT, "logs")
dir.create(CLEAN_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,   showWarnings = FALSE, recursive = TRUE)
cat("[ingest] ROOT =", ROOT, "\n")

source(file.path(ROOT, "helpers", "pisa_io.R"))

PV_MATH_UC <- toupper(PV_MATH)
DEMOG_UC   <- toupper(DEMOG_VARS)
TIME_UC    <- toupper(TIME_VARS)
TUTOR_UC   <- toupper(TUTOR_VARS)
PSYCH_UC   <- toupper(PSYCH_VARS)
SCHOOL_UC  <- toupper(SCHOOL_VARS)

WEIGHT_UC  <- c("W_FSTUWT", paste0("W_FSTURWT", 1:80))

KEEP_STU <- unique(c(
  "CNT", "CNTSCHID", "CNTSTUID",
  PV_MATH_UC, WEIGHT_UC,
  DEMOG_UC, TIME_UC, TUTOR_UC, PSYCH_UC
))
KEEP_SCH <- unique(c("CNT", "CNTSCHID", SCHOOL_UC))

stu_path <- file.path(RAW_DIR, "CY08MSP_STU_QQQ.sas7bdat")
sch_path <- file.path(RAW_DIR, "CY08MSP_SCH_QQQ.sas7bdat")

if (!file.exists(stu_path)) {
  stop("Student file not found: ", stu_path,
       "\nFollow scripts/00_download_pisa.md to acquire the raw data.")
}
if (!file.exists(sch_path)) {
  stop("School file not found: ", sch_path)
}

cat("[ingest] Reading student file ...\n")
t0 <- Sys.time()
stu <- haven::read_sas(stu_path, col_select = any_of(KEEP_STU))
stu <- haven::zap_labels(stu)
cat(sprintf("[ingest]   %d rows x %d cols loaded in %.1fs\n",
            nrow(stu), ncol(stu), as.numeric(Sys.time() - t0, units = "secs")))

cat("[ingest] Reading school file ...\n")
t0 <- Sys.time()
sch <- haven::read_sas(sch_path, col_select = any_of(KEEP_SCH))
sch <- haven::zap_labels(sch)
cat(sprintf("[ingest]   %d rows x %d cols loaded in %.1fs\n",
            nrow(sch), ncol(sch), as.numeric(Sys.time() - t0, units = "secs")))

cat("[ingest] Merging student x school on (CNT, CNTSCHID) ...\n")
dat <- stu %>%
  dplyr::left_join(sch, by = c("CNT", "CNTSCHID"), suffix = c("", "_sch"))

cat("[ingest] Harmonizing covariates ...\n")
dat <- dat %>%
  dplyr::mutate(
    male       = dplyr::if_else(ST004D01T == 2, 1L, 0L),
    esl_native = dplyr::if_else(IMMIG == 1, 1L, 0L),
    repeated   = dplyr::if_else(`REPEAT` >= 1, 1L, 0L),
    grade_ctr  = GRADE - 0L,
    hwk_h      = winsorize(ST296Q01JA, p = c(0.005, 0.995)),
    study_h    = winsorize(STUDYHMW,   p = c(0.005, 0.995)),
    east_asia  = as.integer(CNT %in% EAST_ASIA),
    nordic     = as.integer(CNT %in% NORDIC),
    oecd       = as.integer(CNT %in% OECD_38)
  )

dat <- dat %>%
  dplyr::mutate(
    pv_mean_math = rowMeans(dplyr::across(dplyr::all_of(PV_MATH_UC)),
                            na.rm = TRUE)
  )

cat("[ingest] Computing country summary table ...\n")
country_summary <- dat %>%
  dplyr::group_by(CNT) %>%
  dplyr::summarise(
    n_stud       = dplyr::n(),
    n_school     = dplyr::n_distinct(CNTSCHID),
    pv_math_mean = stats::weighted.mean(pv_mean_math, w = W_FSTUWT, na.rm = TRUE),
    hwk_mean     = stats::weighted.mean(hwk_h,        w = W_FSTUWT, na.rm = TRUE),
    study_mean   = stats::weighted.mean(study_h,      w = W_FSTUWT, na.rm = TRUE),
    escs_mean    = stats::weighted.mean(ESCS,         w = W_FSTUWT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    east_asia = as.integer(CNT %in% EAST_ASIA),
    nordic    = as.integer(CNT %in% NORDIC),
    oecd      = as.integer(CNT %in% OECD_38)
  ) %>%
  dplyr::arrange(dplyr::desc(pv_math_mean))

readr::write_csv(country_summary,
                 file.path(CLEAN_DIR, "country_summary.csv"))
cat(sprintf("[ingest]   %d countries\n", nrow(country_summary)))

out_path <- file.path(CLEAN_DIR, "pisa2022_core.parquet")
cat(sprintf("[ingest] Writing Parquet to %s ...\n", out_path))
arrow::write_parquet(dat, out_path, compression = "zstd")
cat(sprintf("[ingest]   wrote %.1f MB.\n",
            file.info(out_path)$size / 1024^2))
cat("[ingest] DONE.\n")
