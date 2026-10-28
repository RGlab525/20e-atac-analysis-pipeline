# Shared paths for all R scripts.
#   ATAC_PROJECT : folder that holds the input files (BAMs, peak files, DESeq2 tables, ...).
#                  Relative input names in the scripts are resolved against it.
#   ATAC_OUT     : folder for outputs; each script writes into its own subfolder.
# Run every script from the repository root, e.g.
#   ATAC_PROJECT=/path/to/data Rscript 02_differential_accessibility/04_cmeans_TSS_groups.R
# and start each script with: source("config/config.R")

`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a

PROJECT_DIR <- Sys.getenv("ATAC_PROJECT", ".")
OUT_ROOT    <- Sys.getenv("ATAC_OUT", "results")

# input file in the project folder
inp <- function(...) file.path(PROJECT_DIR, ...)

# output folder for one script (created on first use)
out_dir <- function(name) {
    d <- file.path(OUT_ROOT, name)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
}

TIMES <- c("0h", "1h", "3h", "6h", "12h", "24h")

# DAR cutoffs used throughout the manuscript
FC_CUT   <- log2(1.5)
PADJ_CUT <- 0.05
