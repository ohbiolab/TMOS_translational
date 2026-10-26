# ============================================================================
# Project path
#
# Every figure script sources this file. Set TMOS_DIR to the project folder
# that holds the processed Seurat/Xenium objects and the output folders
# (Clean/, Single_cell/, Spatial/), either by editing the default below or
# with the TMOS_DIR environment variable (e.g. in ~/.Renviron).
# ============================================================================

TMOS_DIR <- Sys.getenv("TMOS_DIR", unset = "path/to/project")

if (!dir.exists(TMOS_DIR)) {
  warning("TMOS_DIR does not exist: ", TMOS_DIR,
          "\nEdit config/paths.R or set the TMOS_DIR environment variable.")
}
