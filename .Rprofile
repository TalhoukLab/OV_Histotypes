if (renv:::renv_platform_windows()) {
  NETWORK_DRIVE <- "i:"
} else {
  NETWORK_DRIVE <- "/Volumes/files.ubc.ca"
}
Sys.setenv("RENV_PATHS_CELLAR" = file.path(NETWORK_DRIVE, "cellar"))
Sys.setenv("RENV_PATHS_ROOT" = file.path(NETWORK_DRIVE, "renv"))
rm(NETWORK_DRIVE)

source("renv/activate.R")
