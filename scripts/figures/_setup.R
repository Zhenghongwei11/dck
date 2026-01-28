source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()
source("scripts/utils/manifest.R")

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

ezhu_dir_create("plots/publication")
ezhu_dir_create("plots/publication/png")
ezhu_dir_create("results/figures")

par(family = "Helvetica")

fmt_env <- Sys.getenv("EZHU_FIGURE_FORMAT", unset = "pdf,png")
TARGET_FORMATS <- trimws(unlist(strsplit(fmt_env, "[,;\\s]+")))
TARGET_FORMATS <- TARGET_FORMATS[TARGET_FORMATS %in% c("pdf", "png")]
TARGET_FORMATS <- unique(TARGET_FORMATS)
TARGET_FORMATS <- TARGET_FORMATS[order(match(TARGET_FORMATS, c("pdf", "png")))]
if (length(TARGET_FORMATS) == 0) TARGET_FORMATS <- c("pdf")
