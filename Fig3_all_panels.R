# ============================================================================
# Fig3_all_panels.R
# Figure 3 and Supplementary Figure S3A-B - all panels in one script.
#
#   Panel(s)            Part   Source script (lines)
#   3A, 3B, 3J          -      Xenium Explorer images, no code
#   3C (UMAP + bar)     1      spatial_heart.R (536-586, 645-796)
#   3K                  2      spatial_heart.R (1574-1628)
#   3H, 3I; S3B         3      spatial_heart.R (3427-4182)
#   3E (UMAP + bar), 3G 4      spatial_skeletal.R (617-803, 811-982)
#   3D, 3F              5      Xenium marker dot-plot script
#   S3A (niche UMAP)    -      made interactively; no script found
#   S3C-S3F             6a/6b  Clean/fig4/fig4.R (SingleR and subtype sections)
# ============================================================================
if (!exists("TMOS_DIR")) source(file.path("config", "paths.R"))


##############################################################################
##  PART 1 - Fig 3C: heart Xenium UMAP and cell-type proportions per patient
##  (Spatial/Heart/spatial_heart.R)
##############################################################################

setwd(file.path(TMOS_DIR, "Spatial/Heart"))
xo_heart <- readRDS("xo_heart_processed.rds")
outdir <- file.path(TMOS_DIR, "Spatial/Heart")

library(Seurat)
library(ggplot2)
library(ggrastr)

major_celltype_colors <- c(
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "Cardiomyocyte"        = "#d73027",
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#bdbdbd",
  "Mixed"                = "#4d4d4d"
)


# ----------------------------
# Theme (same style)
# ----------------------------
umap_theme_border <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    plot.title   = element_text(hjust = 0.5)
  )

pt.size <- 0.01

# ----------------------------
# 1) Major cell type UMAP (NO labels; TRUE colors; small PDF)
# ----------------------------
xo_heart$major_celltype <- factor(xo_heart$major_celltype, levels = names(major_celltype_colors))

p_celltype <- DimPlot(
  xo_heart,
  reduction = "umap",
  group.by  = "major_celltype",
  cols      = major_celltype_colors,
  pt.size   = pt.size,
  label     = FALSE,
  raster    = FALSE        # <-- KEY: prevent Seurat fading
) +
  ggtitle("Heart Xenium — Major Cell Types") +
  umap_theme_border

p_celltype_r <- rasterise(p_celltype, dpi = 600, dev = "cairo")  # <-- KEY: one clean raster pass

pdf(file.path(outdir, "UMAP_Heart_Major_Celltypes_NoLabels.pdf"), width = 7, height = 5)
print(p_celltype_r)
dev.off()

############################################################
## HEART: stacked barplots of MAJOR celltype proportions by
##  - pt_id
##  - sample_id
## Uses SAME fixed coloring + formatting style you provided
############################################################

library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# Colors (fixed order = stacking order)
# Make sure keys match levels in xo_heart$major_celltype EXACTLY
# ------------------------------------------------------------
major_celltype_colors <- c(
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "Cardiomyocyte"        = "#d73027",
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#bdbdbd",
  "Mixed"                = "#4d4d4d"
)

# ------------------------------------------------------------
# Pull metadata + enforce fixed celltype order
# ------------------------------------------------------------
stopifnot("major_celltype" %in% colnames(xo_heart@meta.data))
stopifnot("pt_id" %in% colnames(xo_heart@meta.data))
stopifnot("sample_id" %in% colnames(xo_heart@meta.data))

meta_df <- xo_heart@meta.data %>%
  mutate(
    major_celltype = factor(major_celltype, levels = names(major_celltype_colors))
  )

# ------------------------------------------------------------
# Compute proportions (NO re-ordering)
# ------------------------------------------------------------
compute_proportions_fixed_major <- function(df, group_var) {
  df %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(major_celltype)) %>%
    dplyr::count(.data[[group_var]], major_celltype, name = "n") %>%
    dplyr::rename_with(~"group", all_of(group_var)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
}

# ------------------------------------------------------------
# Stacked bar plotting function (fixed order + strong border)
# ------------------------------------------------------------
plot_stacked_bar_major <- function(df,
                                   title,
                                   filename,
                                   outdir,
                                   width = 9,
                                   height = 4) {
  
  p <- ggplot(df, aes(x = group, y = prop, fill = major_celltype)) +
    geom_col(width = 0.8, color = NA) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = major_celltype_colors,
      drop   = FALSE,
      breaks = names(major_celltype_colors)   # legend order fixed
    ) +
    labs(
      x = NULL,
      y = "Cell proportion (%)",
      title = title,
      fill = "Major cell type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks  = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  ggsave(
    filename = filename,
    plot     = p,
    path     = outdir,
    device   = cairo_pdf,
    width    = width,
    height   = height
  )
  
  return(p)
}

# ------------------------------------------------------------
# Output directory (Heart)
# ------------------------------------------------------------
outdir_heart_bar <- file.path(outdir, "StackedBars_MajorCellType")
dir.create(outdir_heart_bar, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Plot dimensions
# (sample_id usually needs wider)
# ------------------------------------------------------------
plot_dims <- list(
  pt_id     = c(5, 4),
  sample_id = c(6, 4.5)
)

# ------------------------------------------------------------
# Generate plots for pt_id + sample_id
# ------------------------------------------------------------
for (gv in names(plot_dims)) {
  
  if (!gv %in% colnames(meta_df)) {
    message("Skipping ", gv, " (not found in metadata).")
    next
  }
  
  prop_df <- compute_proportions_fixed_major(meta_df, gv)
  
  plot_stacked_bar_major(
    df       = prop_df,
    title    = paste("Major cell-type proportions by", gv),
    filename = paste0("Heart_StackedBar_MajorCellType_by_", gv, ".pdf"),
    outdir   = outdir_heart_bar,
    width    = plot_dims[[gv]][1],
    height   = plot_dims[[gv]][2]
  )
}

message("Saved stacked bar PDFs to: ", outdir_heart_bar)

## end of Part 1
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 2 - Fig 3K: spatial niches (N1-N4)
##  (Spatial/Heart/spatial_heart.R)
##############################################################################

setwd(file.path(TMOS_DIR, "Spatial/Heart"))
xo_heart <- readRDS("xo_heart_processed.rds")
outdir <- file.path(TMOS_DIR, "Spatial/Heart")

library(Seurat)
library(ggplot2)
library(ggrastr)
library(grid)
outdir_spatial <- file.path(outdir, "figures_spatial_niches")
dir.create(outdir_spatial, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2) Manual niche colors (fixed forever)
##    (uses Set3 by default, but now saved as explicit vector)
############################################################
# If you want to hard-code specific colors, replace the line above with e.g.:


niche_colors_manual <- c(
  "N1: CS" = "#0072B2",  # strong blue (unchanged)
  "N2: VR" = "#009E73",  # bluish green (unchanged)
  "N3: MIS" = "#F4A261",  # lighter orange
  "N4: IE" = "#E76F51"   # lighter red / coral
)
# The saved object stores the niches as N1-N4, so key the same colours by
# those labels (otherwise every niche is drawn grey)
names(niche_colors_manual) <- c("N1", "N2", "N3", "N4")

############################################################
## 3) Spatial centroids — Niches (N1..Nk) using manual colors
############################################################
img_arg <- if ("fov" %in% names(formals(ImageDimPlot))) "fov" else
  if ("images" %in% names(formals(ImageDimPlot))) "images" else
    stop("ImageDimPlot() lacks 'fov'/'images' args.")
has_geom_arg   <- "geom"   %in% names(formals(ImageDimPlot))
has_raster_arg <- "raster" %in% names(formals(ImageDimPlot))

args <- list(
  object = xo_heart,
  group.by = "niches",
  size = 0.18,
  dark.background = FALSE
)
args[[img_arg]] <- "fov"
if (has_geom_arg)   args$geom   <- "point"
if (has_raster_arg) args$raster <- TRUE

p_niches <- do.call(ImageDimPlot, args) +
  scale_fill_manual(values = niche_colors_manual, drop = FALSE) +
  scale_color_manual(values = niche_colors_manual, drop = FALSE) +
  guides(fill = guide_legend(title = "Niche", override.aes = list(size = 3))) +
  ggtitle("Spatial centroids — niches (N1..)") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.6, "cm")
  )

for (i in seq_along(p_niches$layers)) {
  p_niches$layers[[i]] <- ggrastr::rasterise(p_niches$layers[[i]], dpi = 180)
}

ggsave(
  filename = file.path(outdir_spatial, "Spatial_Centroids_Niches_Nlabels.pdf"),
  plot     = p_niches,
  width    = 14, height = 6.5,
  device   = cairo_pdf
)

## end of Part 2
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 3 - Fig 3H, 3I; Fig S3B: niche composition and niche Reactome pathways
##  (Spatial/Heart/spatial_heart.R)
##############################################################################

setwd(file.path(TMOS_DIR, "Spatial/Heart"))
outdir <- file.path(TMOS_DIR, "Spatial/Heart")
xo_heart <-readRDS("xo_heart_processed.rds")


# Niche Heatmap


############################################################
## HEATMAP: niche composition by major cell type
## - builds count + proportion matrices
## - saves 2 PDFs:
##    (A) Column-normalized proportions (within-niche composition)
##    (B) Row-wise Z-score (cell-type enrichment across niches)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
})

# -----------------------------
# 0) Sanity
# -----------------------------
stopifnot(exists("xo_heart"))
stopifnot(all(c("niches","major_celltype") %in% colnames(xo_heart@meta.data)))

# -----------------------------
# 1) Fixed orders + colors
# -----------------------------
niche_order <- c("N1","N2","N3","N4")

niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

# Use your major celltype palette if it exists; otherwise auto
if (!exists("major_celltype_colors")) {
  ct_levels <- sort(unique(as.character(xo_heart$major_celltype)))
  major_celltype_colors <- setNames(grDevices::hcl.colors(length(ct_levels), "Dark 3"), ct_levels)
}

celltype_order <- names(major_celltype_colors)

# -----------------------------
# 2) Build count table
# -----------------------------
meta <- xo_heart@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  filter(!is.na(niches), !is.na(major_celltype)) %>%
  mutate(
    niches = factor(as.character(niches), levels = niche_order),
    major_celltype = factor(as.character(major_celltype), levels = celltype_order)
  )

count_df <- meta %>%
  count(major_celltype, niches, name = "n")

count_mat <- count_df %>%
  tidyr::pivot_wider(names_from = niches, values_from = n, values_fill = 0) %>%
  as.data.frame()

rownames(count_mat) <- count_mat$major_celltype
count_mat$major_celltype <- NULL
count_mat <- as.matrix(count_mat)

# enforce column order
count_mat <- count_mat[, niche_order, drop = FALSE]

# -----------------------------
# 3) Column-normalized proportions (within-niche composition)
#    each column sums to 1
# -----------------------------
prop_col_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/")
prop_col_mat[is.na(prop_col_mat)] <- 0

# -----------------------------
# 4) Row-wise Z-score (enrichment of each cell type across niches)
#    each row centered/scaled across niches
# -----------------------------
z_row_mat <- t(scale(t(prop_col_mat)))
z_row_mat[is.na(z_row_mat)] <- 0

# -----------------------------
# 5) Annotations (niche colors across columns)
# -----------------------------
annotation_col <- data.frame(
  Niche = factor(colnames(prop_col_mat), levels = niche_order)
)
rownames(annotation_col) <- colnames(prop_col_mat)

annotation_colors <- list(
  Niche = niche_colors_manual
)

# -----------------------------
# 6) Output directory + save PDFs
# -----------------------------
outdir_heat <- file.path(outdir, "Heatmaps_NicheComposition")
dir.create(outdir_heat, recursive = TRUE, showWarnings = FALSE)



# B) Cell-type enrichment across niches (row Z-score)
pdf(file.path(outdir_heat, "Heatmap_NicheComposition_RowZ_Enrichment.pdf"),
    width = 6, height = 5.5)
pheatmap(
  z_row_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Cell-type enrichment across niches\n(row-wise Z-score of col-normalized proportions)",
  angle_col = 0,
  fontsize_row = 12,
  fontsize_col = 12
)
dev.off()

# -----------------------------
# 7) Save the underlying matrices for reproducibility
# -----------------------------
write.csv(count_mat,     file.path(outdir_heat, "Matrix_Counts_CellType_by_Niche.csv"))
write.csv(prop_col_mat,  file.path(outdir_heat, "Matrix_ColNormalizedProportions_CellType_by_Niche.csv"))
write.csv(z_row_mat,     file.path(outdir_heat, "Matrix_RowZ_CellTypeEnrichment_by_Niche.csv"))

message("✅ Saved heatmaps + matrices to: ", outdir_heat)









############################################################
## Stacked bar — major cell type proportions by niche (N1–N4)
## - Fixed niche order + fixed cell-type color order
## - Same formatting style (theme_bw + bold border, no grid)
## - Saves PDF to: outdir_spatial
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

############################################################
## Define output directory (once per session)
############################################################

outdir_spatial <- file.path(TMOS_DIR, "Spatial/Heart")

dir.create(outdir_spatial, showWarnings = FALSE, recursive = TRUE)

message("Output directory set to: ", outdir_spatial)
############################################################
## Stacked bar — major cell type proportions by niche
## Using YOUR fixed color palette
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})



# ----------------------------------------------------------
# 2) Define niche colors + order
# ----------------------------------------------------------
niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)
niche_levels <- names(niche_colors_manual)

# ----------------------------------------------------------
# 3) Your major cell type palette (EXACT)
# ----------------------------------------------------------
major_celltype_colors <- c(
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "Cardiomyocyte"        = "#d73027",
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#bdbdbd",
  "Mixed"                = "#4d4d4d"
)

celltype_levels <- names(major_celltype_colors)

# ----------------------------------------------------------
# 4) Clean metadata + enforce matching labels
#    (IMPORTANT: harmonize underscores vs hyphens)
# ----------------------------------------------------------

meta_df <- xo_heart@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  mutate(
    niches = factor(as.character(niches), levels = niche_levels),
    
    # harmonize naming to match your color keys
    major_celltype = as.character(major_celltype),
    major_celltype = gsub("-", "_", major_celltype),
    major_celltype = gsub(" ", "_", major_celltype),
    major_celltype = factor(major_celltype, levels = celltype_levels)
  )

# ----------------------------------------------------------
# 5) Compute proportions
# ----------------------------------------------------------
prop_df <- meta_df %>%
  filter(!is.na(niches), !is.na(major_celltype)) %>%
  count(niches, major_celltype, name = "n") %>%
  group_by(niches) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ----------------------------------------------------------
# 6) Plot (publication style)
# ----------------------------------------------------------
p_bar <- ggplot(prop_df, aes(x = niches, y = prop, fill = major_celltype)) +
  geom_col(width = 0.8, color = NA) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = major_celltype_colors,
    drop   = FALSE,
    breaks = celltype_levels
  ) +
  labs(
    x = "Niche",
    y = "Cell proportion (%)",
    title = "Major cell-type proportions by niche",
    fill = "Major cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5, color = "black"),
    axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.ticks   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line    = element_blank(),
    plot.margin  = margin(8, 8, 8, 8)
  )

# ----------------------------------------------------------
# 7) Save PDF
# ----------------------------------------------------------
ggsave(
  filename = file.path(outdir_spatial, "StackedBar_MajorCellType_by_Niche.pdf"),
  plot     = p_bar,
  width    = max(6, 0.55 * length(niche_levels) + 3),
  height   = 5,
  device   = cairo_pdf
)

message("✅ Saved stacked bar to Figure 3 folder")





suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(writexl)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# 0) INPUTS / PATHS
# ------------------------------------------------------------
stopifnot(exists("xo_heart"))
stopifnot(all(c("niches","major_celltype") %in% colnames(xo_heart@meta.data)))

# set your output root (edit if you want a subfolder)
outdir <- file.path(TMOS_DIR, "Spatial/Heart")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outdir_niche_labels <- file.path(outdir, "Niche_Labeling")
dir.create(outdir_niche_labels, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) COLORS (your palette)
# ------------------------------------------------------------
major_celltype_colors <- c(
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "Cardiomyocyte"        = "#d73027",
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#bdbdbd",
  "Mixed"                = "#4d4d4d"
)

niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

# enforce consistent factor order
xo_heart$niches <- factor(as.character(xo_heart$niches), levels = names(niche_colors_manual))
xo_heart$major_celltype <- factor(as.character(xo_heart$major_celltype), levels = names(major_celltype_colors))

# ------------------------------------------------------------
# 2) CELL-TYPE COMPOSITION PER NICHE  (counts + proportions)
# ------------------------------------------------------------
meta_df <- xo_heart@meta.data %>%
  tibble::as_tibble(rownames = "cell_id") %>%
  dplyr::filter(!is.na(niches), !is.na(major_celltype)) %>%
  dplyr::mutate(
    niches = factor(as.character(niches), levels = names(niche_colors_manual)),
    major_celltype = factor(as.character(major_celltype), levels = names(major_celltype_colors))
  )

celltype_counts_props <- meta_df %>%
  dplyr::count(niches, major_celltype, name = "n_cells") %>%
  dplyr::group_by(niches) %>%
  dplyr::mutate(prop = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

# Wide versions (counts + props)
celltype_counts_wide <- celltype_counts_props %>%
  dplyr::select(niches, major_celltype, n_cells) %>%
  tidyr::pivot_wider(names_from = major_celltype, values_from = n_cells, values_fill = 0)

celltype_props_wide <- celltype_counts_props %>%
  dplyr::select(niches, major_celltype, prop) %>%
  tidyr::pivot_wider(names_from = major_celltype, values_from = prop, values_fill = 0)

# Save to Excel
xlsx_celltype <- file.path(outdir_niche_labels, "Niche_Celltype_Composition_Counts_and_Props.xlsx")
writexl::write_xlsx(
  list(
    celltype_counts_long = celltype_counts_props,
    celltype_counts_wide = celltype_counts_wide,
    celltype_props_wide  = celltype_props_wide
  ),
  path = xlsx_celltype
)
message("✅ Saved cell-type composition Excel: ", xlsx_celltype)

# ------------------------------------------------------------
# 3) HEATMAP (cell-type proportions by niche) -> PDF
# ------------------------------------------------------------
mat_prop <- celltype_counts_props %>%
  dplyr::select(niches, major_celltype, prop) %>%
  tidyr::pivot_wider(names_from = niches, values_from = prop, values_fill = 0) %>%
  as.data.frame()

rownames(mat_prop) <- mat_prop$major_celltype
mat_prop$major_celltype <- NULL
mat_prop <- as.matrix(mat_prop)

# order rows in your palette order
mat_prop <- mat_prop[names(major_celltype_colors), names(niche_colors_manual), drop = FALSE]

annotation_col <- data.frame(Niche = colnames(mat_prop))
rownames(annotation_col) <- colnames(mat_prop)

annotation_colors <- list(Niche = niche_colors_manual)

pdf(file.path(outdir_niche_labels, "Heatmap_CelltypeProportions_by_Niche.pdf"), width = 6.5, height = 4.8)
pheatmap::pheatmap(
  mat_prop,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Major cell-type proportions by niche",
  angle_col = 0
)
dev.off()
message("✅ Saved heatmap PDF: Heatmap_CelltypeProportions_by_Niche.pdf")

# ------------------------------------------------------------
# 4) NICHE LABELING USING:
#    (a) top enriched cell type (highest prop)
#    (b) Reactome enrichment of UP genes from niche-vs-rest DE
# ------------------------------------------------------------

DefaultAssay(xo_heart) <- "Xenium"
Idents(xo_heart) <- "niches"

padj_cut <- 0.05
lfc_cut  <- 0.25   # use 0.25 for enrichment sensitivity; change to 0.58 if you want stricter

# Helpers
safe_as_tibble_enrich <- function(enr_obj) {
  if (is.null(enr_obj)) return(NULL)
  df <- tryCatch(as.data.frame(enr_obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)  # force base data.frame
  tibble::as_tibble(df)
}

convert_to_entrez <- function(genes) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) return(NULL)
  
  suppressMessages({
    m <- clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  })
  if (is.null(m) || nrow(m) == 0) return(NULL)
  m %>% dplyr::distinct(ENTREZID)
}

run_reactome <- function(entrez_df) {
  if (is.null(entrez_df) || nrow(entrez_df) == 0) return(NULL)
  suppressMessages({
    ReactomePA::enrichPathway(
      gene          = entrez_df$ENTREZID,
      organism      = "human",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      readable      = TRUE
    )
  })
}

niches_use <- levels(xo_heart$niches)
niches_use <- niches_use[!is.na(niches_use) & niches_use != ""]

reactome_all_sheets <- list()
top_label_rows <- list()
de_tables <- list()

for (nch in niches_use) {
  
  message("Processing niche: ", nch)
  
  # ---- (a) top cell type by proportion
  ct_df <- celltype_counts_props %>% dplyr::filter(niches == nch)
  top_ct <- ct_df %>% dplyr::arrange(dplyr::desc(prop)) %>% dplyr::slice_head(n = 1)
  
  # ---- (b) DE niche vs rest
  de_raw <- FindMarkers(
    object          = xo_heart,
    ident.1         = nch,
    ident.2         = NULL,          # rest
    assay           = "Xenium",
    slot            = "data",
    test.use        = "wilcox",
    logfc.threshold = 0,
    min.pct         = 0
  )
  
  fc_col <- if ("avg_log2FC" %in% colnames(de_raw)) "avg_log2FC" else "avg_logFC"
  
  de_df <- de_raw %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      p_val_adj = as.numeric(p_val_adj),
      fc = as.numeric(.data[[fc_col]])
    ) %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj != 1, !is.na(fc))
  
  de_tables[[paste0("DE_", nch, "_vs_rest")]] <- de_df
  
  # ---- UP genes in this niche
  genes_up <- de_df %>%
    dplyr::filter(p_val_adj < padj_cut, fc > lfc_cut) %>%
    dplyr::pull(gene) %>% unique()
  
  # ---- Reactome enrichment
  entrez <- convert_to_entrez(genes_up)
  enr <- run_reactome(entrez)
  enr_tbl <- safe_as_tibble_enrich(enr)
  
  # Save ALL pathways for this niche (sheet)
  if (is.null(enr_tbl)) {
    reactome_all_sheets[[paste0("Reactome_", nch)]] <- data.frame()
    top_rx   <- NA_character_
    top_rx_p <- NA_real_
  } else {
    reactome_all_sheets[[paste0("Reactome_", nch)]] <- as.data.frame(enr_tbl)
    
    # pick top pathway WITHOUT slice() (avoids Rle/list issues)
    enr_tbl2 <- enr_tbl %>% dplyr::arrange(p.adjust)
    top_rx   <- as.character(enr_tbl2$Description[1])
    top_rx_p <- as.numeric(enr_tbl2$p.adjust[1])
  }
  
  # ---- Top label table row
  top_label_rows[[nch]] <- data.frame(
    niches = nch,
    top_celltype = as.character(top_ct$major_celltype[1]),
    top_celltype_prop = as.numeric(top_ct$prop[1]),
    reactome_top = top_rx,
    reactome_padj = top_rx_p,
    n_cells_in_niche = sum(xo_heart$niches == nch),
    n_up_genes_used_for_reactome = length(genes_up),
    stringsAsFactors = FALSE
  )
}

top_label_table <- dplyr::bind_rows(top_label_rows) %>%
  dplyr::arrange(niches)

# ------------------------------------------------------------
# 5) SAVE: ALL REACTOME PATHWAYS + LABEL TABLE + DE TABLES
# ------------------------------------------------------------
xlsx_reactome <- file.path(outdir_niche_labels, "Niche_Reactome_AllPathways_and_TopLabels.xlsx")

writexl::write_xlsx(
  c(
    list(NICHE_TOP_LABELS = top_label_table),
    reactome_all_sheets,
    de_tables
  ),
  path = xlsx_reactome
)

message("✅ Saved Reactome+labels Excel: ", xlsx_reactome)

# Print the niche labels in console
print(top_label_table)







suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(forcats)
  library(patchwork)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

# Point this to the Excel you already generated
xlsx_reactome <- file.path(outdir_niche_labels, "Niche_Reactome_AllPathways_and_TopLabels.xlsx")

# Output
out_pdf <- file.path(outdir_niche_labels, "Reactome_Top5Pathways_by_Niche_2x2.pdf")
out_png <- file.path(outdir_niche_labels, "Reactome_Top5Pathways_by_Niche_2x2.png")

# ------------------------------------------------------------
# Helper to read one niche sheet + build plot
# ------------------------------------------------------------
make_top5_plot <- function(nch) {
  sh <- paste0("Reactome_", nch)
  
  df <- readxl::read_excel(xlsx_reactome, sheet = sh) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  # Handle empty sheets gracefully
  if (nrow(df) == 0 || !all(c("Description","p.adjust") %in% colnames(df))) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(nch, ": no Reactome terms"))
    )
  }
  
  top5 <- df %>%
    mutate(
      p.adjust = as.numeric(p.adjust),
      Description = as.character(Description)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Description), Description != "") %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    mutate(
      neglog10padj = -log10(p.adjust),
      Description = forcats::fct_reorder(Description, neglog10padj)  # for horizontal bars
    )
  
  ggplot(top5, aes(x = Description, y = neglog10padj)) +
    geom_col(fill = niche_colors_manual[[nch]], width = 0.75) +
    coord_flip() +
    labs(
      title = nch,
      x = NULL,
      y = expression(-log[10]~"(adj p)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}

# ------------------------------------------------------------
# Build the 2x2 panel plot (N1–N4)
# ------------------------------------------------------------
p1 <- make_top5_plot("N1")
p2 <- make_top5_plot("N2")
p3 <- make_top5_plot("N3")
p4 <- make_top5_plot("N4")

p_all <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Top 5 Reactome pathways per niche (lowest adjusted p-value)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(out_pdf, p_all, width = 14, height = 8, device = cairo_pdf)
ggsave(out_png, p_all, width = 14, height = 8, dpi = 300)

message("✅ Saved: ", out_pdf)
message("✅ Saved: ", out_png)






suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(forcats)
  library(patchwork)
})

niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

xlsx_reactome <- file.path(outdir_niche_labels,
                           "Niche_Reactome_AllPathways_and_TopLabels.xlsx")

out_pdf <- file.path(outdir_niche_labels,
                     "Reactome_Top5Pathways_by_Niche_2x2_LabelOnBar.pdf")

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------
make_top5_plot <- function(nch) {
  
  sh <- paste0("Reactome_", nch)
  
  df <- readxl::read_excel(xlsx_reactome, sheet = sh) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  if (nrow(df) == 0) {
    return(
      ggplot() + theme_void() +
        labs(title = paste0(nch, ": no Reactome terms"))
    )
  }
  
  top5 <- df %>%
    dplyr::mutate(
      p.adjust = as.numeric(p.adjust),
      Description = as.character(Description)
    ) %>%
    dplyr::filter(!is.na(p.adjust), Description != "") %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::mutate(
      neglog10padj = -log10(p.adjust),
      Description = forcats::fct_reorder(Description, neglog10padj)
    )
  
  max_height <- max(top5$neglog10padj)
  
  ggplot(top5, aes(x = Description, y = neglog10padj)) +
    geom_col(fill = niche_colors_manual[[nch]], width = 0.75) +
    geom_text(
      aes(
        y = 0.02 * max_height,   # start near left edge of bar
        label = Description
      ),
      hjust = 0,
      color = "black",
      size = 3.6
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = nch,
      x = NULL,
      y = expression(-log[10]~"(adj p)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}

# ------------------------------------------------------------
# Combine into 2×2 layout
# ------------------------------------------------------------
p_all <- (make_top5_plot("N1") | make_top5_plot("N2")) /
  (make_top5_plot("N3") | make_top5_plot("N4")) +
  plot_annotation(
    title = "Top 5 Reactome pathways per niche"
  ) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(out_pdf, p_all, width = 7, height = 8, device = cairo_pdf)

message("✅ Saved: ", out_pdf)

## end of Part 3
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 4 - Fig 3E, 3G: skeletal Xenium UMAP, proportions per patient, P53 heart vs skeletal
##  (Spatial/Skeletal/spatial_skeletal.R)
##############################################################################

setwd(file.path(TMOS_DIR, "Spatial/Skeletal"))
xo_skeletal <- readRDS("skeletal_combined_352026.rds")
outdir <- file.path(TMOS_DIR, "Spatial/Skeletal")

library(Seurat)
library(ggplot2)
library(ggrastr)

# ----------------------------
# 1) Setup Skeletal Specific Colors
# ----------------------------
# Mapping your predicted_celltype to the project's color scheme
skeletal_celltype_colors <- c(
  "Endothelial"         = "#1b9e77",
  "Fibroblast/Stromal"  = "#7570b3",
  "Smooth Muscle"       = "#66a61e",
  "Skeletal Muscle"     = "#800000", # Mapping to the red color used for Cardiomyocytes
  "T/NK Lymphocytes"    = "#1f78b4",
  "Myeloid"             = "#e6ab02",
  "B cells"             = "#d95f02",
  "Cycling"             = "#bdbdbd",
  "Myogenic Progenitor" = "#b2df8a"
)

# ----------------------------
# 2) Theme & Styling (Matching Heart Formatting)
# ----------------------------
umap_theme_border <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )

pt.size <- 0.01

# ----------------------------
# 3) Process and Plot
# ----------------------------
# Ensure the cell types are factored for consistent legend ordering
xo_skeletal$predicted_celltype <- factor(
  xo_skeletal$predicted_celltype, 
  levels = names(skeletal_celltype_colors)
)

p_celltype <- DimPlot(
  xo_skeletal,
  reduction = "umap",
  group.by  = "predicted_celltype",
  cols      = skeletal_celltype_colors,
  pt.size   = pt.size,
  label     = FALSE,
  raster    = FALSE        # Using ggrastr manually for better quality
) +
  ggtitle("Skeletal Xenium — Major Cell Types") +
  umap_theme_border

# ----------------------------
# 4) High-Quality Raster Export
# ----------------------------
# Rasterizing the points at 600 DPI keeps the file size small but the quality sharp
p_celltype_r <- rasterise(p_celltype, dpi = 600, dev = "cairo")

out_path <- file.path(outdir, "UMAP_Skeletal_Major_Celltypes_NoLabels.pdf")
pdf(out_path, width = 7, height = 5)
print(p_celltype_r)
dev.off()

message("✅ UMAP Saved: ", out_path)






library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Colors (Fixed Order for Stacking)
# ------------------------------------------------------------
skeletal_celltype_colors <- c(
  "Endothelial"         = "#1b9e77",
  "Fibroblast/Stromal"  = "#7570b3",
  "Smooth Muscle"       = "#66a61e",
  "Skeletal Muscle"     = "#800000", 
  "T/NK Lymphocytes"    = "#1f78b4",
  "Myeloid"             = "#e6ab02",
  "B cells"             = "#d95f02",
  "Cycling"             = "#bdbdbd",
  "Myogenic Progenitor" = "#b2df8a"
)

# ------------------------------------------------------------
# 2. Pull Metadata + Enforce Fixed Order
# ------------------------------------------------------------
# Using xo_skeletal and the SingleR-predicted labels
meta_df <- xo_skeletal@meta.data %>%
  mutate(
    predicted_celltype = factor(predicted_celltype, levels = names(skeletal_celltype_colors))
  )

# ------------------------------------------------------------
# 3. Proportion Calculation Function
# ------------------------------------------------------------
compute_proportions_skeletal <- function(df, group_var) {
  df %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(predicted_celltype)) %>%
    dplyr::count(.data[[group_var]], predicted_celltype, name = "n") %>%
    dplyr::rename(group = all_of(group_var)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
}

# ------------------------------------------------------------
# 4. Plotting Function (Matching Heart Formatting)
# ------------------------------------------------------------
plot_stacked_bar_skeletal <- function(df, title, filename, outdir, width, height) {
  
  p <- ggplot(df, aes(x = group, y = prop, fill = predicted_celltype)) +
    geom_col(width = 0.8, color = NA) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = skeletal_celltype_colors,
      drop   = FALSE,
      breaks = names(skeletal_celltype_colors)
    ) +
    labs(
      x = NULL, y = "Cell proportion (%)",
      title = title, fill = "Cell Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.ticks   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin  = margin(8, 8, 8, 8)
    )
  
  ggsave(
    filename = filename,
    plot     = p,
    path     = outdir,
    device   = cairo_pdf,
    width    = width,
    height   = height
  )
  
  return(p)
}

# ------------------------------------------------------------
# 5. Execution
# ------------------------------------------------------------
outdir_skel_bar <- file.path(outdir, "StackedBars_Skeletal_Proportions")
dir.create(outdir_skel_bar, recursive = TRUE, showWarnings = FALSE)

# Define variables to plot (ensure these exist in your xo_skeletal metadata)
vars_to_plot <- c("pt_id", "sample_id")

for (gv in vars_to_plot) {
  if (!gv %in% colnames(meta_df)) {
    warning("Column ", gv, " not found in xo_skeletal metadata. Skipping.")
    next
  }
  
  prop_df <- compute_proportions_skeletal(meta_df, gv)
  
  plot_stacked_bar_skeletal(
    df       = prop_df,
    title    = paste("Skeletal proportions by", gv),
    filename = paste0("Skeletal_StackedBar_by_", gv, ".pdf"),
    outdir   = outdir_skel_bar,
    width    = ifelse(gv == "sample_id", 7, 5), # sample_id often needs more width
    height   = 4.5
  )
}

message("✅ Proportions plots saved to: ", outdir_skel_bar)




xo_heart <- readRDS(file.path(TMOS_DIR, "Spatial/Heart/xo_heart_processed.rds"))

##################### pT 153 ###############################



library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------
# 1. Setup Shared Colors & Mapping
# ------------------------------------------------------------
major_celltype_colors <- c(
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "Cardiomyocyte"        = "#d73027", # Use this for Skeletal Muscle too
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#bdbdbd"
)

# ------------------------------------------------------------
# 2. Extract and Harmonize Metadata
# ------------------------------------------------------------

# Process Heart Metadata
heart_df <- xo_heart@meta.data %>%
  dplyr::filter(pt_id == "153") %>%
  transmute(
    tissue = "Heart",
    celltype = as.character(major_celltype_broad)
  )

# Process Skeletal Metadata & Map to Heart-like strings
skel_df <- xo_skeletal@meta.data %>%
  dplyr::filter(pt_id == "153") %>%
  transmute(
    tissue = "Skeletal",
    celltype = case_when(
      predicted_celltype == "Skeletal Muscle"    ~ "Cardiomyocyte", # Re-mapping for side-by-side
      predicted_celltype == "Fibroblast/Stromal" ~ "Fibroblast_stromal",
      predicted_celltype == "Smooth Muscle"     ~ "Smoothe_muscle_mural",
      predicted_celltype == "T/NK Lymphocytes"  ~ "T_NK_lymphocytes",
      predicted_celltype == "B cells"            ~ "B_cells",
      predicted_celltype == "Cycling"            ~ "Cycling_immune",
      TRUE                                       ~ as.character(predicted_celltype)
    )
  )

# Combine and Filter for "Matching" types (excluding Myogenic/Mixed/NAs)
combined_meta <- bind_rows(heart_df, skel_df) %>%
  dplyr::filter(celltype %in% names(major_celltype_colors)) %>%
  mutate(celltype = factor(celltype, levels = names(major_celltype_colors)))

# ------------------------------------------------------------
# 3. Compute Proportions
# ------------------------------------------------------------
# ------------------------------------------------------------
# 3. Compute Proportions (Fixed for n() error)
# ------------------------------------------------------------
plot_data <- combined_meta %>%
  dplyr::group_by(tissue, celltype) %>%
  dplyr::summarise(count_val = dplyr::n(), .groups = 'drop') %>% 
  dplyr::group_by(tissue) %>%
  dplyr::mutate(prop = count_val / sum(count_val))

# Verify the data exists before plotting
if (exists("plot_data")) {
  print(head(plot_data))
} else {
  stop("plot_data failed to create. Check if combined_meta is empty.")
}

# ------------------------------------------------------------
# 4. Generate Faceted Bar Plot
# ------------------------------------------------------------
p_compare <- ggplot(plot_data, aes(x = tissue, y = prop, fill = celltype)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~tissue, scales = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.01)) +
  scale_fill_manual(values = major_celltype_colors) +
  labs(
    title = "Cell Type Proportion Comparison (Patient 153)",
    subtitle = "Skeletal Muscle mapped to 'Cardiomyocyte' for visual alignment",
    x = NULL,
    y = "Proportion of Total Cells (%)",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(), # Tissue name is in the facet label
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# Save the comparison
ggsave(file.path(outdir, "PT153_Tissue_Comparison_StackedBar.pdf"), 
       p_compare, width = 6, height = 6, device = cairo_pdf)

print(p_compare)




library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Prepare Data & Set 50% Limit
# ------------------------------------------------------------
pyramid_data <- plot_data %>%
  mutate(
    # Flip Heart to the left (negative)
    display_prop = ifelse(tissue == "Heart", -prop, prop),
    # Ensure the Y-axis (Cell Types) is ordered correctly
    celltype = factor(celltype, levels = rev(names(major_celltype_colors)))
  )

# ------------------------------------------------------------
# 2. Generate the Plot
# ------------------------------------------------------------
p_pyramid_fixed <- ggplot(pyramid_data, aes(x = celltype, y = display_prop, fill = celltype)) +
  # Draw the bars with a thin black outline for each bar
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  # Add the vertical center line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  # Flip coordinates to make bars horizontal
  coord_flip() +
  # Set horizontal axis to 50% max each side (-0.5 to 0.5)
  scale_y_continuous(
    limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, 0.1),
    labels = c("50%", "40%", "30%", "20%", "10%", "0", "10%", "20%", "30%", "40%", "50%"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = major_celltype_colors, guide = "none") +
  theme_classic(base_size = 12) +
  theme(
    # Add the requested border around the entire figure
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), # Remove default lines since we have a border
    axis.text.y = element_text(face = "plain", color = "black", size = 8.5),
    axis.text.x = element_text(color = "black", size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid.major.x = element_line(color = "grey95", linetype = "solid"),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    title = "Patient 153: Tissue Composition (Tornado Plot)",
    subtitle = "← HEART  |  SKELETAL →",
    x = NULL,
    y = "Proportion of Total Cells (%)"
  )

# ------------------------------------------------------------
# 3. Add Tissue Header Labels
# ------------------------------------------------------------
p_pyramid_fixed <- p_pyramid_fixed +
  annotate("label", x = length(unique(pyramid_data$celltype)) + 0.5, y = -0.25, 
           label = "HEART", fill = "white", fontface = "bold", size = 4) +
  annotate("label", x = length(unique(pyramid_data$celltype)) + 0.5, y = 0.25, 
           label = "SKELETAL", fill = "white", fontface = "bold", size = 4)

# Save the plot
out_pyramid_50 <- file.path(outdir, "PT153_Tornado_50pct_Border.pdf")
ggsave(out_pyramid_50, p_pyramid_fixed, width = 4.5, height = 4.5, device = cairo_pdf)

## end of Part 4
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 5 - Fig 3D, 3F: heart and skeletal marker dot plots
##  (fig3/fig3DF_xenium_marker_dotplots.R)
##############################################################################

# ============================================================================
# fig3DF_xenium_marker_dotplots.R
# Xenium heart and skeletal marker dot plots
#
# Figure panels : Fig 3D, 3F
# Note          : Fig 3A, 3B, 3J are Xenium Explorer images (no code).
# ============================================================================

setwd(file.path(TMOS_DIR, "Clean/fig3"))
xo_heart <- readRDS("xo_heart_processed.rds")
xo_skeletal <- readRDS("skeletal_combined_352026.rds")


library(Seurat)
library(ggplot2)

# 0) Make sure we use gene expression, not the niche assay
DefaultAssay(xo_heart) <- "Xenium"

# 1) Idents and ordering
stopifnot("major_celltype" %in% colnames(xo_heart@meta.data))
Idents(xo_heart) <- "major_celltype"

celltype_order_compartment <- c(
  "B_cells",
  "T_NK_lymphocytes",
  "Myeloid",
  "Cycling_immune",
  "Endothelial",
  "Fibroblast_stromal",
  "Smoothe_muscle_mural",
  "Cardiomyocyte",
  "Mixed"
)

xo_heart$major_celltype <- factor(xo_heart$major_celltype, levels = celltype_order_compartment)
Idents(xo_heart) <- "major_celltype"

# 2) Markers
marker_list_major <- list(
  B_cells              = c("MZB1", "POU2AF1", "CD38"),
  T_NK_lymphocytes     = c("CD8A", "GZMK", "CTLA4"),
  Myeloid              = c("CD14", "CD163", "FCGR3A"),
  Cycling_immune       = c("MKI67", "TOP2A", "CDK1"),
  Endothelial          = c("PECAM1", "S1PR1", "ENG"),
  Fibroblast_stromal   = c("PDGFRA", "LTBP2", "COL5A1"),
  Smoothe_muscle_mural = c("PDGFRB", "RGS5", "KCNJ8"),
  Cardiomyocyte        = c("TNNI3", "TNNT2", "ACTN2"),
  Mixed                = c("DUX4", "SOX2-OT", "HOXB7")
)

features_ordered <- unlist(marker_list_major[celltype_order_compartment], use.names = FALSE)

# 3) Handle Xenium naming quirks (often "-" becomes "." in feature names)
feature_universe <- rownames(xo_heart[["Xenium"]])

features_mapped <- ifelse(
  features_ordered %in% feature_universe, features_ordered,
  ifelse(gsub("-", ".", features_ordered) %in% feature_universe, gsub("-", ".", features_ordered), NA)
)

features_use <- unique(na.omit(features_mapped))

# Print what didn't map (so you can fix spellings if needed)
missing <- features_ordered[is.na(features_mapped)]
if (length(missing) > 0) {
  message("Still missing in Xenium assay:\n  ", paste(missing, collapse = ", "))
}

stopifnot(length(features_use) > 0)

# 4) Build DotPlot then redraw with genes on columns, cells on rows
byr_cols <- c("#4FA3D9", "#FFF2A6", "#F15A4A")

dp <- DotPlot(xo_heart, features = features_use, dot.scale = 6)

# Ensure factors are set for the desired order
# Note: For y-axis (rows), we usually use the original order so the first item is at the top
dp$data$id <- factor(dp$data$id, levels = rev(celltype_order_compartment)) 
dp$data$features.plot <- factor(dp$data$features.plot, levels = features_use)

p_genes_on_cols <- ggplot(dp$data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = byr_cols) +
  scale_size(range = c(0, 6)) +
  theme_classic(base_size = 12) +
  theme(
    # Genes are now on X, so we rotate these labels
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y  = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    plot.margin  = ggplot2::margin(t = 10, r = 10, b = 10, l = 10) 
  ) +
  labs(
    x = "Markers",
    y = "Cell Type",
    color = "Average expression",
    size  = "Percent expressed"
  )

p_genes_on_cols

# Saving the updated orientation
ggsave(
  filename = "DotPlot_Heart_MajorCellType_GenesOnCols_CellsOnRows.pdf",
  plot     = p_genes_on_cols,
  width    = 9,   # Swapped width/height since genes (usually more numerous) are now horizontal
  height   = 3.5,
  device   = cairo_pdf
)













# Ensure we are using the correct assay and identity
DefaultAssay(xo_skeletal) <- "Xenium"
Idents(xo_skeletal) <- "predicted_celltype"

# Find markers (using 'wilcox' is standard, 'only.pos' ensures we get upregulated genes)
all_markers <- FindAllMarkers(xo_skeletal, 
                              only.pos = TRUE, 
                              min.pct = 0.1, 
                              logfc.threshold = 0.25)

# Extract top 3 markers per cluster by Fold Change
library(dplyr)
top3_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

# Create the marker list for the plot
marker_list_skeletal <- split(top3_markers$gene, top3_markers$cluster)

# View the list to confirm markers make biological sense
print(marker_list_skeletal)




# =========================================================
# 1. Setup Plotting Parameters
# =========================================================
celltype_order <- c("T/NK Lymphocytes", "B cells", "Myeloid", "Cycling", 
                    "Endothelial", "Fibroblast/Stromal", "Smooth Muscle", "Skeletal Muscle",
                    "Myogenic Progenitor")

byr_cols <- c("#4FA3D9", "#FFF2A6", "#F15A4A")

# =========================================================
# 1. Setup Skeletal-Specific Parameters
# =========================================================
celltype_order_skeletal <- c(
  "B cells", "T/NK Lymphocytes", 
  "Myeloid", 
  "Cycling", 
  "Endothelial", 
  "Fibroblast/Stromal", 
  "Smooth Muscle", 
  "Skeletal Muscle",
  "Myogenic Progenitor"
)

byr_cols <- c("#4FA3D9", "#FFF2A6", "#F15A4A")

# =========================================================
# 2. Flexible DotPlot Function
# =========================================================
generate_skeletal_dotplot <- function(obj, group_col, marker_list, ordered_levels) {
  
  # 1. Set Identity and Factor Levels for Y-axis (Cells)
  # rev() ensures the first item in your list is at the TOP of the plot
  obj[[group_col]] <- factor(obj@meta.data[[group_col]], levels = rev(ordered_levels))
  Idents(obj) <- group_col
  
  # 2. Order Genes (X-axis) based on the provided order
  ordered_marker_names <- ordered_levels[ordered_levels %in% names(marker_list)]
  
  # Extract genes in that specific order
  markers_ordered <- unname(unlist(marker_list[ordered_marker_names]))
  
  # Filter for genes actually present in the Xenium assay
  markers_ordered <- unique(markers_ordered[markers_ordered %in% rownames(obj)])
  
  # 3. Create Plot
  p <- DotPlot(
    obj,
    features  = markers_ordered,
    assay     = "Xenium",
    dot.scale = 6,
    scale     = TRUE 
  ) 
  
  # Force X-axis (genes) to respect the marker_list order
  p$data$features.plot <- factor(p$data$features.plot, levels = markers_ordered)
  
  # Update aesthetics and fix the margin namespace error
  p <- p +
    scale_color_gradientn(colors = byr_cols) + 
    theme_classic(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y     = element_text(size = 10, color = "black"),
      legend.position = "right",
      axis.line       = element_blank(),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin     = ggplot2::margin(8, 8, 8, 8), # Explicitly ggplot2::
      plot.title      = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x     = "Top Markers (Ordered by Cell Type)", 
      y     = "Cell Type",
            color = "Avg Expression",
      size  = "% Expressed"
    )
  
  return(p)
}

# =========================================================
# 3. Execution for Skeletal
# =========================================================

# Make sure marker_list_skeletal has keys that match celltype_order_skeletal EXACTLY
p_skeletal <- generate_skeletal_dotplot(
  obj           = xo_skeletal, 
  group_col     = "predicted_celltype", 
  marker_list   = marker_list_skeletal,
  ordered_levels = celltype_order_skeletal
)

# Save the plot
ggsave(
  filename = "DotPlot_Skeletal_Final_Ordered.pdf", 
  plot     = p_skeletal, 
  width    = 9, 
  height   = 3.5, 
  device   = cairo_pdf
)

print(p_skeletal)






library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# =========================================================
# 1. Define the Standardized Theme & Scales
# =========================================================
# =========================================================
# 1. Define the Standardized Theme & Scales (Updated Font)
# =========================================================
shared_setup <- list(
  scale_color_gradientn(
    colors = byr_cols, 
    limits = c(-2.5, 2.5), 
    oob = scales::squish
  ),
  scale_size(
    range = c(0, 6), 
    limits = c(0, 100)
  ),
  theme_classic(base_size = 12),
  theme(
    # Decreased font sizes here:
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y     = element_text(size = 8, color = "black"),
    # Optional: decrease axis titles (Markers / Cell Type) as well
    axis.title      = element_text(size = 8),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line       = element_blank(),
    plot.margin     = ggplot2::margin(5, 5, 5, 5)
  )
)


# =========================================================
# 2. Build Heart Plot
# =========================================================
# Use your heart-specific features
p_heart <- DotPlot(xo_heart, features = features_use, assay = "Xenium", scale = TRUE) + 
  shared_setup +
  labs(title = "Xenium Heart", x = NULL, y = "Cell Type")

# Force the specific heart order
p_heart$data$id <- factor(p_heart$data$id, levels = rev(celltype_order_compartment))
p_heart$data$features.plot <- factor(p_heart$data$features.plot, levels = features_use)

# =========================================================
# 3. Build Skeletal Plot
# =========================================================
# Use your skeletal-specific markers
# Ensure skeletal_features is defined from your marker_list_skeletal
skeletal_features <- unname(unlist(marker_list_skeletal[celltype_order_skeletal]))
skeletal_features <- unique(skeletal_features[skeletal_features %in% rownames(xo_skeletal)])

p_skeletal <- DotPlot(xo_skeletal, features = skeletal_features, assay = "Xenium", scale = TRUE) + 
  shared_setup +
  labs(title = "Xenium Skeletal", x = "Markers", y = "Cell Type")

# Force the specific skeletal order
p_skeletal$data$id <- factor(p_skeletal$data$id, levels = rev(celltype_order_skeletal))
p_skeletal$data$features.plot <- factor(p_skeletal$data$features.plot, levels = skeletal_features)




# =========================================================
# 4. Combine with Patchwork (The "Joint" Part)
# =========================================================
# '/' stacks them vertically
# 'guides = "collect"' pulls the legends into one single shared legend
p_final <- (p_heart / p_skeletal) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# =========================================================
# 5. Save as PDF
# =========================================================
ggsave(
  filename = "Stacked_DotPlots_Matched_Scales.pdf",
  plot     = p_final,
  width    = 8,   # Adjust based on your longest gene list
  height   = 6, 
  device   = cairo_pdf
)

print(p_final)

## end of Part 5
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())



##############################################################################
##  PART 6a - Fig S3C-S3F inputs: SingleR (Monaco) subtyping of T/NK and myeloid cells -> tnk_myeloid_heart_st.rds
##  (Clean/fig4/fig4.R, SingleR section; also scores the IL-17 program used in Fig 4)
##############################################################################

setwd(file.path(TMOS_DIR, "Clean/fig4"))
xo_heart <- readRDS("xo_heart_processed.rds")
# This part uses only the expression data and metadata. The spatial image (FOV)
# was saved with an older SeuratObject and fails validation in subset() with
# newer versions, so it is dropped here (results are unchanged).
xo_heart@images <- list()

library(reactome.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(openxlsx)

library(SingleR)
library(celldex)

# ── 0. IL-17 SIGNALING-RELATED PROGRAM SCORE (43 genes) ──
# Scored on the full object so the checkpoint carries IL17_Full, which Fig 4
# (fig4/Fig4_all_panels.R Parts 1, 2 and 7) uses. Same steps as Fig 4 Part 1.
cols_to_drop <- c("Reactome_MAPK_1", "Reactome_MAPK", "Reactome_IL17_1", "Reactome_IL17",
                  "Reactome_Overlap_1", "Reactome_Overlap", "IL17_Full_1", "IL17_Full")
xo_heart@meta.data <- xo_heart@meta.data[, !(colnames(xo_heart@meta.data) %in% cols_to_drop)]
reac_to_entrez <- as.list(reactomePATHID2EXTID)
sym_il17  <- AnnotationDbi::select(org.Hs.eg.db, keys = reac_to_entrez[["R-HSA-448424"]],
                                   columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
il17_full <- unique(sym_il17[sym_il17 %in% rownames(xo_heart)])
xo_heart <- AddModuleScore(xo_heart, features = list(il17_full), name = "IL17_Full_", assay = "Xenium")
xo_heart$IL17_Full <- xo_heart$IL17_Full_1
xo_heart@meta.data <- xo_heart@meta.data[, !grepl("_1$", colnames(xo_heart@meta.data))]

# ── 1. SUBSET TO TNK + MYELOID ──
refined_order <- c("CD4_Tcell", "CD8_Tcell", "Exhausted_CD8_Tcell",
                   "Cycling_Tcell", "NK_cell", "Tfh", "Treg",
                   "COLEC12_P2RY14_Macro", "CXCL10_C1QC_Macro", "LYVE1_MRC1_Macro",
                   "FCAR_IL1R2_PMN", "S100A12_FCN1_Mono", "STAT4_GPR183_Mono")

immune_sub <- subset(xo_heart, subset = major_celltype_refined %in% refined_order)

# ── 2. GET EXPRESSION MATRIX ──
expr_mat <- GetAssayData(immune_sub, assay = "Xenium", layer = "data")

# ── 3. LOAD REFERENCE ──
# Monaco is best for immune cells (T, NK, Myeloid subtypes)
ref <- celldex::MonacoImmuneData()

# ── 4. RUN SINGLER ──
pred <- SingleR(
  test      = expr_mat,
  ref       = ref,
  labels    = ref$label.fine   # fine = detailed subtypes
)

cat("SingleR label distribution:\n")
print(table(pred$labels))


# ── 1. ADD SINGLER LABELS TO METADATA ──
# First add to immune_sub
immune_sub$SingleR_label <- pred$labels

# ── 2. REMOVE B CELLS ──
b_cell_labels <- c(
  "Naive B cells",
  "Non-switched memory B cells", 
  "Switched memory B cells",
  "Exhausted B cells",
  "Plasmablasts"
)

immune_sub_filtered <- subset(immune_sub, 
                              subset = SingleR_label %in% b_cell_labels, 
                              invert = TRUE)

cat("Cells before B cell removal:", ncol(immune_sub), "\n")
cat("Cells after B cell removal: ", ncol(immune_sub_filtered), "\n")
print(table(immune_sub_filtered$SingleR_label))

# ── 3. VERIFY METADATA COLUMN IS THERE ──
cat("\nChecking metadata columns:\n")
cat("SingleR_label" %in% colnames(immune_sub_filtered@meta.data))

# ── 4. SAVE RDS CHECKPOINT ──
saveRDS(immune_sub_filtered, 
        file = "tnk_myeloid_heart_st.rds")

## end of Part 6a
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 6b - Fig S3C, S3D, S3E, S3F: subtype marker dot plots and composition per niche
##  (uses tnk_myeloid_heart_st.rds written by Part 6a)
##############################################################################

setwd(file.path(TMOS_DIR, "Clean/fig4"))
library(reactome.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(openxlsx)
library(SingleR)
library(celldex)
immune_sub_filtered <- readRDS("tnk_myeloid_heart_st.rds")
# Only expression and metadata are used here; drop the spatial image (FOV),
# which was saved with an older SeuratObject and fails validation in subset()
immune_sub_filtered@images <- list()






library(Seurat)
library(ggplot2)

# ── 1. DEFINE WHICH LABELS ARE TNK VS MYELOID ──
tnk_labels <- c(
  "Central memory CD8 T cells",
  "Effector memory CD8 T cells",
  "Naive CD8 T cells",
  "Terminal effector CD8 T cells",
  "Naive CD4 T cells",
  "Terminal effector CD4 T cells",
  "Th1 cells",
  "Th2 cells",
  "Th17 cells",
  "Th1/Th17 cells",
  "T regulatory cells",
  "Follicular helper T cells",
  "MAIT cells",
  "Non-Vd2 gd T cells",
  "Vd2 gd T cells",
  "Natural killer cells",
  "Progenitor cells"
)

myeloid_labels <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "Non classical monocytes",
  "Myeloid dendritic cells",
  "Plasmacytoid dendritic cells",
  "Low-density basophils",
  "Low-density neutrophils"
)

# ── 2. SUBSET EACH ──
tnk_only     <- subset(immune_sub_filtered, subset = SingleR_label %in% tnk_labels)
myeloid_only <- subset(immune_sub_filtered, subset = SingleR_label %in% myeloid_labels)

Idents(tnk_only)     <- "SingleR_label"
Idents(myeloid_only) <- "SingleR_label"

# ── 3. FIND TOP MARKERS ──
cat("Finding T/NK markers...\n")
tnk_markers <- FindAllMarkers(
  tnk_only,
  assay           = "Xenium",
  test.use        = "wilcox",
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.25
)

cat("Finding Myeloid markers...\n")
myeloid_markers <- FindAllMarkers(
  myeloid_only,
  assay           = "Xenium",
  test.use        = "wilcox",
  only.pos        = TRUE,
  min.pct         = 0.10,
  logfc.threshold = 0.25
)

library(Seurat)
library(ggplot2)
library(dplyr)

# ── 1. GET TOP 3 MARKERS PER SUBTYPE (DEDUPLICATED) ──

# T/NK — deduplicated, ordered by cluster
top3_tnk_dedup <- tnk_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%  # grab extra to allow for dedup
  ungroup() %>%
  # Keep gene for the cluster where it has highest FC
  arrange(desc(avg_log2FC)) %>%
  filter(!duplicated(gene)) %>%
  # Now take top 3 per cluster after dedup
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup() %>%
  # Order by cluster factor level for clean flow
  mutate(cluster = factor(cluster, levels = tnk_labels)) %>%
  arrange(cluster)

# Myeloid — same logic
top3_myeloid_dedup <- myeloid_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup() %>%
  arrange(desc(avg_log2FC)) %>%
  filter(!duplicated(gene)) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = myeloid_labels)) %>%
  arrange(cluster)

# ── 2. GENE ORDER FOLLOWS CELL TYPE ORDER ──
genes_tnk_ordered     <- unique(top3_tnk_dedup$gene)
genes_myeloid_ordered <- unique(top3_myeloid_dedup$gene)

cat("T/NK genes (ordered):\n");    print(genes_tnk_ordered)
cat("Myeloid genes (ordered):\n"); print(genes_myeloid_ordered)

# ── UPDATED THEME WITH BORDER ──
dotplot_theme <- function() {
  theme_classic() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
      axis.text.y      = element_text(size = 8),
      axis.title       = element_blank(),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position  = "right",
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line        = element_blank()
    )
}

# ── T/NK DOTPLOT ──
dp_tnk <- DotPlot(
  tnk_only,
  features  = genes_tnk_ordered,
  assay     = "Xenium",
  cols      = c("blue", "red"),   # placeholder — overridden below
  dot.scale = 8
) +
  # Override color scale with 3-point gradient
  scale_color_gradientn(
    colors = c("#4FA3D9", "#FFF2A6", "#F15A4A"),
    name   = "Average\nExpression"
  ) +
  scale_y_discrete(limits = rev(tnk_labels)) +
  labs(title = "T/NK Subtype Markers") +
  guides(
    size  = guide_legend(title = "% Expressed", override.aes = list(color = "grey50")),
    color = guide_colorbar(title = "Average\nExpression", barwidth = 1, barheight = 6)
  ) +
  dotplot_theme()

ggsave("TNK_Subtypes_DotPlot_Top3Markers.pdf",
       dp_tnk, width = 14, height = 7)
cat("✅ T/NK dot plot saved.\n")

# ── MYELOID DOTPLOT ──
dp_myeloid <- DotPlot(
  myeloid_only,
  features  = genes_myeloid_ordered,
  assay     = "Xenium",
  cols      = c("blue", "red"),   # placeholder — overridden below
  dot.scale = 8
) +
  scale_color_gradientn(
    colors = c("#4FA3D9", "#FFF2A6", "#F15A4A"),
    name   = "Average\nExpression"
  ) +
  scale_y_discrete(limits = rev(myeloid_labels)) +
  labs(title = "Myeloid Subtype Markers") +
  guides(
    size  = guide_legend(title = "% Expressed", override.aes = list(color = "grey50")),
    color = guide_colorbar(title = "Average\nExpression", barwidth = 1, barheight = 6)
  ) +
  dotplot_theme()

ggsave("Myeloid_Subtypes_DotPlot_Top3Markers.pdf",
       dp_myeloid, width = 10, height = 5)
cat("✅ Myeloid dot plot saved.\n")










library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)

# ── 1. DEFINE LABELS & COLORS ──
tnk_labels <- c(
  "Central memory CD8 T cells",
  "Effector memory CD8 T cells",
  "Naive CD8 T cells",
  "Terminal effector CD8 T cells",
  "Naive CD4 T cells",
  "Terminal effector CD4 T cells",
  "Th1 cells",
  "Th2 cells",
  "Th17 cells",
  "Th1/Th17 cells",
  "T regulatory cells",
  "Follicular helper T cells",
  "MAIT cells",
  "Non-Vd2 gd T cells",
  "Vd2 gd T cells",
  "Natural killer cells",
  "Progenitor cells"
)

myeloid_labels <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "Non classical monocytes",
  "Myeloid dendritic cells",
  "Plasmacytoid dendritic cells",
  "Low-density basophils",
  "Low-density neutrophils"
)

niche_order <- c("N1", "N2", "N3", "N4")

niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

niche_labels_named <- c(
  "N1" = "N1: CS",
  "N2" = "N2: VR",
  "N3" = "N3: MIS",
  "N4" = "N4: IE"
)

tnk_colors <- c(
  "Central memory CD8 T cells"    = "#1F77B4",
  "Effector memory CD8 T cells"   = "#AEC7E8",
  "Naive CD8 T cells"             = "#4E9AC7",
  "Terminal effector CD8 T cells" = "#0A4C82",
  "Naive CD4 T cells"             = "#2CA02C",
  "Terminal effector CD4 T cells" = "#98DF8A",
  "Th1 cells"                     = "#D62728",
  "Th2 cells"                     = "#FF9896",
  "Th17 cells"                    = "#FF7F0E",
  "Th1/Th17 cells"                = "#FFBB78",
  "T regulatory cells"            = "#9467BD",
  "Follicular helper T cells"     = "#C5B0D5",
  "MAIT cells"                    = "#8C564B",
  "Non-Vd2 gd T cells"            = "#C49C94",
  "Vd2 gd T cells"                = "#E377C2",
  "Natural killer cells"          = "#17BECF",
  "Progenitor cells"              = "#BCBD22"
)

myeloid_colors <- c(
  "Classical monocytes"           = "#1F77B4",
  "Intermediate monocytes"        = "#FF7F0E",
  "Non classical monocytes"       = "#2CA02C",
  "Myeloid dendritic cells"       = "#D62728",
  "Plasmacytoid dendritic cells"  = "#9467BD",
  "Low-density basophils"         = "#8C564B",
  "Low-density neutrophils"       = "#E377C2"
)

# ── 2. COUNTS & % WITHIN LINEAGE PER NICHE ──

# Total T/NK cells per niche
tnk_niche_totals <- immune_sub_filtered@meta.data %>%
  filter(SingleR_label %in% tnk_labels, !is.na(niches)) %>%
  group_by(niches) %>%
  summarise(lineage_total = n(), .groups = "drop")

# Total Myeloid cells per niche
myeloid_niche_totals <- immune_sub_filtered@meta.data %>%
  filter(SingleR_label %in% myeloid_labels, !is.na(niches)) %>%
  group_by(niches) %>%
  summarise(lineage_total = n(), .groups = "drop")

# T/NK subtype counts + %
tnk_niche_pct <- immune_sub_filtered@meta.data %>%
  filter(SingleR_label %in% tnk_labels, !is.na(niches)) %>%
  group_by(niches, SingleR_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(tnk_niche_totals, by = "niches") %>%
  mutate(
    pct           = (n / lineage_total) * 100,
    niches        = factor(niches, levels = niche_order),
    SingleR_label = factor(SingleR_label, levels = tnk_labels),
    niche_label   = paste0(niche_labels_named[as.character(niches)],
                           "\n(n=", lineage_total, ")")
  )

# Myeloid subtype counts + %
myeloid_niche_pct <- immune_sub_filtered@meta.data %>%
  filter(SingleR_label %in% myeloid_labels, !is.na(niches)) %>%
  group_by(niches, SingleR_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(myeloid_niche_totals, by = "niches") %>%
  mutate(
    pct           = (n / lineage_total) * 100,
    niches        = factor(niches, levels = niche_order),
    SingleR_label = factor(SingleR_label, levels = myeloid_labels),
    niche_label   = paste0(niche_labels_named[as.character(niches)],
                           "\n(n=", lineage_total, ")")
  )

# Fix niche_label factor order
tnk_niche_pct <- tnk_niche_pct %>%
  mutate(niche_label = factor(niche_label,
                              levels = unique(niche_label[order(niches)])))

myeloid_niche_pct <- myeloid_niche_pct %>%
  mutate(niche_label = factor(niche_label,
                              levels = unique(niche_label[order(niches)])))

# ── 3. PLOT FUNCTION ──
make_lineage_bar <- function(data, color_map, title, lineage) {
  
  niche_lvls <- levels(data$niches)
  x_colors   <- niche_colors_manual[niche_lvls]
  
  ggplot(data, aes(x = niche_label, y = pct, fill = SingleR_label)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = color_map) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    labs(
      title = title,
      x     = NULL,
      y     = paste0("% of ", lineage, " Cells"),
      fill  = "Subtype"
    ) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x     = element_text(size = 10, face = "bold",
                                     color = x_colors),
      axis.text.y     = element_text(size = 9),
      axis.title.y    = element_text(size = 10),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(3.5, "mm"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line       = element_blank()
    )
}

# ── 4. GENERATE PLOTS ──
p_tnk <- make_lineage_bar(
  tnk_niche_pct,
  tnk_colors,
  "T/NK Subtype Composition per Niche",
  "T/NK"
)

p_myeloid <- make_lineage_bar(
  myeloid_niche_pct,
  myeloid_colors,
  "Myeloid Subtype Composition per Niche",
  "Myeloid"
)

# ── 5. SAVE PDFs ──
ggsave("TNK_Subtype_Pct_of_TNK_perNiche.pdf",
       p_tnk, width = 7, height = 6)

ggsave("Myeloid_Subtype_Pct_of_Myeloid_perNiche.pdf",
       p_myeloid, width = 7, height = 6)

cat("✅ Lineage composition plots saved.\n")

# ── 6. SAVE EXCEL ──
wb <- createWorkbook()

addWorksheet(wb, "TNK_Pct_per_Niche")
writeData(wb, "TNK_Pct_per_Niche",
          tnk_niche_pct %>%
            dplyr::select(niches, SingleR_label, n, lineage_total, pct) %>%
            arrange(niches, SingleR_label))
setColWidths(wb, "TNK_Pct_per_Niche", cols = 1:5, widths = "auto")

addWorksheet(wb, "Myeloid_Pct_per_Niche")
writeData(wb, "Myeloid_Pct_per_Niche",
          myeloid_niche_pct %>%
            dplyr::select(niches, SingleR_label, n, lineage_total, pct) %>%
            arrange(niches, SingleR_label))
setColWidths(wb, "Myeloid_Pct_per_Niche", cols = 1:5, widths = "auto")

saveWorkbook(wb, "Lineage_Subtype_Composition.xlsx", overwrite = TRUE)
cat("✅ Excel saved.\n")




# ── COMBINED STACKED BAR ──

# Combine both with a lineage label
combined_niche_pct <- bind_rows(
  tnk_niche_pct     %>% mutate(lineage = "T/NK"),
  myeloid_niche_pct %>% mutate(lineage = "Myeloid")
) %>%
  mutate(
    # Recalculate pct as % of ALL immune cells (T/NK + Myeloid) per niche
    lineage       = factor(lineage, levels = c("T/NK", "Myeloid")),
    SingleR_label = factor(SingleR_label, levels = c(tnk_labels, myeloid_labels))
  )

# Combined colors
all_colors <- c(tnk_colors, myeloid_colors)

# Recalculate totals as combined T/NK + Myeloid per niche
combined_totals <- immune_sub_filtered@meta.data %>%
  filter(SingleR_label %in% c(tnk_labels, myeloid_labels), !is.na(niches)) %>%
  group_by(niches) %>%
  summarise(combined_total = n(), .groups = "drop")

combined_niche_pct <- combined_niche_pct %>%
  left_join(combined_totals, by = "niches") %>%
  mutate(
    pct_combined = (n / combined_total) * 100,
    niche_label  = paste0(niche_labels_named[as.character(niches)],
                          "\n(n=", combined_total, ")")
  ) %>%
  mutate(niche_label = factor(niche_label,
                              levels = unique(niche_label[order(niches)])))

# ── PLOT ──
niche_lvls <- levels(combined_niche_pct$niches)
x_colors   <- niche_colors_manual[niche_lvls]

p_combined <- ggplot(combined_niche_pct,
                     aes(x = niche_label, y = pct_combined,
                         fill = SingleR_label)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  
  # Dividing line between T/NK and Myeloid sections within each bar
  # via facet wrap on lineage — shows them as separate stacked segments
  labs(
    title = "T/NK & Myeloid Subtype Composition per Niche",
    x     = NULL,
    y     = "% of Combined Immune Cells",
    fill  = "Subtype"
  ) +
  # Split legend into two columns for readability
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x     = element_text(size = 10, face = "bold",
                                   color = x_colors),
    axis.text.y     = element_text(size = 9),
    axis.title.y    = element_text(size = 10),
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8, face = "bold"),
    legend.key.size = unit(3, "mm"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line       = element_blank()
  )

ggsave("Combined_TNK_Myeloid_StackedBar_Niche.pdf",
       p_combined, width = 8, height = 8)

cat("✅ Combined stacked bar saved.\n")




make_count_bar <- function(data, color_map, title) {
  ggplot(data, aes(x = SingleR_label, y = n, fill = niches)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(
      values = niche_colors_manual,
      labels = niche_labels_named
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = title,
      x     = NULL,
      y     = "Cell Count",
      fill  = "Niche"
    ) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 9),
      axis.title.y    = element_text(size = 10),
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(3.5, "mm"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line       = element_blank()
    )
}

# The count plots use tnk_count_plot / myeloid_count_plot, which were created
# interactively and are not defined in the original script. They are the
# per-niche subtype tables built above (counts n + niche labels with totals):
tnk_count_plot     <- tnk_niche_pct
myeloid_count_plot <- myeloid_niche_pct

p_tnk_count     <- make_count_bar(tnk_count_plot,     tnk_colors,     "T/NK Subtype Counts per Niche")
p_myeloid_count <- make_count_bar(myeloid_count_plot, myeloid_colors, "Myeloid Subtype Counts per Niche")

ggsave("TNK_Subtype_Counts_perNiche.pdf",
       p_tnk_count, width = 10, height = 6)

ggsave("Myeloid_Subtype_Counts_perNiche.pdf",
       p_myeloid_count, width = 8, height = 6)

cat("✅ Stacked count bar plots saved.\n")




make_count_bar <- function(data, color_map, title) {
  ggplot(data, aes(x = niche_label, y = n, fill = SingleR_label)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(values = color_map) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = title,
      x     = NULL,
      y     = "Cell Count",
      fill  = "Subtype"
    ) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x     = element_text(size = 10, face = "bold",
                                     color = niche_colors_manual[niche_lvls]),
      axis.text.y     = element_text(size = 9),
      axis.title.y    = element_text(size = 10),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(3.5, "mm"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line       = element_blank()
    )
}

niche_lvls <- c("N1", "N2", "N3", "N4")

p_tnk_count     <- make_count_bar(tnk_count_plot,     tnk_colors,     "T/NK Subtype Counts per Niche")
p_myeloid_count <- make_count_bar(myeloid_count_plot, myeloid_colors, "Myeloid Subtype Counts per Niche")

ggsave("TNK_Subtype_Counts_perNiche.pdf",
       p_tnk_count, width = 7, height = 6)

ggsave("Myeloid_Subtype_Counts_perNiche.pdf",
       p_myeloid_count, width = 7, height = 6)

## end of Part 6b
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())
