disable_renv_raw <- Sys.getenv("EZHU_DISABLE_RENV", unset = "")
disable_renv <- tolower(trimws(disable_renv_raw)) %in% c("1", "true", "yes", "y")

if (!disable_renv) {
  source("renv/activate.R")
}
