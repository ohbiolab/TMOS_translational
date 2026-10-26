# Install the R packages used across the analysis scripts.
# Analyses were run with R >= 4.3 and Seurat v5.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))     install.packages("remotes")

cran <- c(
  "Seurat", "harmony", "Matrix", "dplyr", "tidyr", "tibble", "stringr", "forcats",
  "purrr", "rlang", "tidyverse", "data.table", "readxl", "openxlsx", "writexl",
  "ggplot2", "ggrepel", "ggpubr", "ggrastr", "ggvenn", "patchwork", "cowplot",
  "scales", "viridis", "RColorBrewer", "pheatmap", "circlize", "treemapify",
  "magick", "rstatix", "broom", "broom.mixed", "coin", "lme4", "lmerTest",
  "lmtest", "sandwich", "logistf", "glmnet", "pROC", "msigdbr", "FNN", "dbscan",
  "zoo", "tidytext", "FactoMineR", "factoextra", "future", "future.apply",
  "hdf5r", "SoupX", "magrittr"
)

bioc <- c(
  "AnnotationDbi", "org.Hs.eg.db", "reactome.db", "ReactomePA", "clusterProfiler",
  "enrichplot", "fgsea", "DESeq2", "apeglm", "edgeR", "limma", "MAST",
  "ComplexHeatmap", "SingleR", "celldex", "scuttle", "SummarizedExperiment",
  "zellkonverter"
)

github <- c(
  "jinworks/CellChat",
  "chris-mcginnis-ucsf/DoubletFinder",
  "satijalab/seurat-data",
  "davidsjoberg/ggsankey",
  "JiaLiVUMC/scMRMA"
)

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
for (g in github) remotes::install_github(g, upgrade = "never")
