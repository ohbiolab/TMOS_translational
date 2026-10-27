# ============================================================================
# Fig4_all_panels.R
# Figure 4, Supplementary Figure S4 and Supplementary Figure S5 - all panels in one script.
#
#   Panel(s)                 Part   Source script
#   4A, 4G                   -      Xenium Explorer images, no code
#   4B, 4C; S4A              1      Clean/fig4/fig4.R
#   4D, 4E, 4H; S4B, S4D     2a     Clean/fig4/fig4_prox.R
#   S4H                      2b     Clean/fig4/fig4_prox.R (last section)
#   4F                       3      Clean/reviewer_comments/rev2.3_spatial.R
#   S5A, S5B                 4a     Clean/reviewer_comments/spatial_adjustment.R
#   4I, 4J; S5C-S5G          4b     Clean/reviewer_comments/spatial_adjustment.R (standalone section)
#   4K                       5      Clean/reviewer_comments/R1.1_spatial_to_sc.R
#   4L; S4E, S4G             6      Clean/fig4/niche_to_pbmc/fig4_spatial_to_pbmc_il17.R
#   S4C                      7      Clean/fig4/th17_homeo_patho.R
#   S4F                      8      Clean/fig4/il17_receptor_ligand.R
# ============================================================================
if (!exists("TMOS_DIR")) source(file.path("config", "paths.R"))

# Xenium objects saved with an older SeuratObject lack FOV slots added later;
# fix_fov() fills them with defaults so subset() works (coordinates unchanged).
fix_fov <- function(obj) {
  patch <- function(x) {
    if (!isS4(x)) return(x)
    sl <- methods::getSlots(class(x))
    proto <- methods::getClass(class(x))@prototype
    for (s in names(sl)) if (!methods::.hasSlot(x, s)) {
      # use the class default (e.g. compact = FALSE) when it has one
      v <- attr(proto, s, exact = TRUE)
      if (!is.null(v)) { attr(x, s) <- v; next }
      # class unions such as OptionalList cannot be instantiated; use a member class
      cls <- c(sl[[s]], setdiff(names(methods::getClass(sl[[s]])@subclasses), "NULL"))
      for (k in cls) {
        v <- tryCatch(methods::new(k), error = function(e) NULL)
        if (!is.null(v)) { attr(x, s) <- v; break }
      }
    }
    x
  }
  for (i in names(obj@images)) {
    f <- patch(obj@images[[i]])
    if (methods::.hasSlot(f, "boundaries")) f@boundaries <- lapply(f@boundaries, patch)
    obj@images[[i]] <- f
  }
  obj
}


##############################################################################
##  PART 1 - Fig 4B, 4C; Fig S4A: IL-17 signaling-related program (43 genes) by niche and cell type
##  (Clean/fig4/fig4.R; S3C-S3F section moved to fig3 Part 6)
##############################################################################

setwd(file.path(TMOS_DIR, "Clean/fig4"))
xo_heart <- readRDS("xo_heart_processed.rds")

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

# ── 1. CLEANUP ──
cols_to_drop <- c("Reactome_MAPK_1", "Reactome_MAPK", "Reactome_IL17_1", "Reactome_IL17",
                  "Reactome_Overlap_1", "Reactome_Overlap", "IL17_Full_1", "IL17_Full")
xo_heart@meta.data <- xo_heart@meta.data[, !(colnames(xo_heart@meta.data) %in% cols_to_drop)]

# ── 2. PULL GENE SETS ──
reac_to_entrez <- as.list(reactomePATHID2EXTID)

entrez_mapk <- reac_to_entrez[["R-HSA-450294"]]
sym_mapk    <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_mapk, columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL

entrez_il17 <- reac_to_entrez[["R-HSA-448424"]]
sym_il17    <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_il17, columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL

valid_genes      <- rownames(xo_heart)
mapk_filtered    <- unique(sym_mapk[sym_mapk %in% valid_genes])
il17_filtered    <- unique(sym_il17[sym_il17 %in% valid_genes])
overlap_filtered <- intersect(mapk_filtered, il17_filtered)

# All 43 IL-17 genes (the full IL-17 panel = 38 shared + 5 unique)
il17_full <- il17_filtered  # all 43

cat("MAPK genes in panel:        ", length(mapk_filtered), "\n")
cat("IL-17 genes in panel:       ", length(il17_filtered), "\n")
cat("Overlap (shared) genes:     ", length(overlap_filtered), "\n")
cat("IL-17 unique only genes:    ", length(setdiff(il17_filtered, mapk_filtered)), "\n")

# (MAPK / IL-17 Venn diagram not used in the figures; omitted)


# ── 4. SCORE ALL CELLS (IL-17 full 43 genes only) ──
xo_heart <- AddModuleScore(xo_heart, features = list(il17_full), name = "IL17_Full_", assay = "Xenium")
xo_heart$IL17_Full <- xo_heart$IL17_Full_1
xo_heart@meta.data <- xo_heart@meta.data[, !grepl("_1$", colnames(xo_heart@meta.data))]

# ── 5. SUMMARY: By Cell Type AND Niche ──
celltype_niche_summary <- xo_heart@meta.data %>%
  filter(!is.na(major_celltype_broad), !is.na(niches)) %>%
  group_by(major_celltype_broad, niches) %>%
  summarise(
    n             = n(),
    IL17Full_Mean = mean(IL17_Full),
    .groups       = "drop"
  )

# ── 6. HEATMAP SETUP ──
ct_order    <- c("T_NK_lymphocytes", "Myeloid", "B_cells", "Cycling_immune",
                 "Endothelial", "Fibroblast_stromal", "Smoothe_muscle_mural", "Cardiomyocyte")
niche_order <- c("N1", "N2", "N3", "N4")

niche_colors_manual <- c(
  "N1" = "#0072B2",
  "N2" = "#009E73",
  "N3" = "#F4A261",
  "N4" = "#E76F51"
)

niche_labels <- c(
  "N1" = "N1: CS",
  "N2" = "N2: VR",
  "N3" = "N3: MIS",
  "N4" = "N4: IE"
)

# ── 7. BUILD MATRIX ──
mat_il17full <- celltype_niche_summary %>%
  filter(major_celltype_broad %in% ct_order) %>%
  mutate(
    major_celltype_broad = factor(major_celltype_broad, levels = ct_order),
    niches               = factor(niches, levels = niche_order)
  ) %>%
  dplyr::select(major_celltype_broad, niches, IL17Full_Mean) %>%
  pivot_wider(names_from = niches, values_from = IL17Full_Mean) %>%
  column_to_rownames("major_celltype_broad") %>%
  as.matrix()

mat_il17full_z <- t(scale(t(mat_il17full)))

# ── 8. COLORS & ANNOTATIONS ──
col_fun <- colorRamp2(c(-2, 0, 2), c("#4682B4", "white", "#b30000"))

top_anno <- HeatmapAnnotation(
  Niche                = niche_order,
  col                  = list(Niche = niche_colors_manual),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  simple_anno_size     = unit(3, "mm")
)

# ── 9. BUILD HEATMAP ──
ht_il17full <- Heatmap(
  mat_il17full_z,
  name             = "Z-score",
  col              = col_fun,
  top_annotation   = top_anno,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  row_names_gp     = gpar(fontsize = 9),
  column_labels    = niche_labels[colnames(mat_il17full_z)],
  column_names_gp  = gpar(fontsize = 9, fontface = "bold", col = "black"),
  column_names_rot = 45,
  border           = FALSE,
  rect_gp          = gpar(col = NA),
  column_title     = "IL-17 Signaling (All 43 Genes)",
  column_title_gp  = gpar(fontsize = 11, fontface = "bold")
)

# ── 10. SAVE HEATMAP PDF ──
pdf("AllCellTypes_IL17_Full43_NicheHeatmap.pdf", width = 5, height = 5)
draw(
  ht_il17full,
  merge_legend = TRUE,
  padding      = unit(c(5, 15, 5, 5), "mm")
)
dev.off()
cat("✅ Heatmap saved.\n")

# ── 11. SAVE GENE LIST TO EXCEL ──
wb <- createWorkbook()

# Sheet 1: All 43 IL-17 genes
il17_unique   <- setdiff(il17_filtered, mapk_filtered)
overlap_genes <- intersect(il17_filtered, mapk_filtered)

gene_table <- data.frame(
  Gene     = il17_full,
  Category = ifelse(il17_full %in% overlap_genes, "Shared with MAPK", "IL-17 Unique")
)

addWorksheet(wb, "IL17_All43_Genes")
writeData(wb, "IL17_All43_Genes", gene_table)
setColWidths(wb, "IL17_All43_Genes", cols = 1:2, widths = "auto")

# Sheet 2: Overlap genes only
addWorksheet(wb, "Shared_MAPK_IL17_38genes")
writeData(wb, "Shared_MAPK_IL17_38genes", data.frame(Gene = overlap_genes))
setColWidths(wb, "Shared_MAPK_IL17_38genes", cols = 1, widths = "auto")

# Sheet 3: IL-17 unique only
addWorksheet(wb, "IL17_Unique_5genes")
writeData(wb, "IL17_Unique_5genes", data.frame(Gene = il17_unique))
setColWidths(wb, "IL17_Unique_5genes", cols = 1, widths = "auto")

# Sheet 4: MAPK full gene list
addWorksheet(wb, "MAPK_All38_Genes")
writeData(wb, "MAPK_All38_Genes", data.frame(Gene = mapk_filtered))
setColWidths(wb, "MAPK_All38_Genes", cols = 1, widths = "auto")

saveWorkbook(wb, "MAPK_IL17_GeneList.xlsx", overwrite = TRUE)

# ── T/NK + MYELOID SUBTYPES (SingleR checkpoint, fig3/Fig3_all_panels.R Part 6a) ──
immune_sub_filtered <- readRDS("tnk_myeloid_heart_st.rds")
# IL-17 program score from the full object scored above (same cells)
immune_sub_filtered$IL17_Full <- xo_heart$IL17_Full[colnames(immune_sub_filtered)]
rm(xo_heart); invisible(gc())

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

# ── N4 CELLS FOR FIG 4C ──
n4_data <- immune_sub_filtered@meta.data %>%
  filter(niches == "N4", SingleR_label %in% c(tnk_labels, myeloid_labels), !is.na(IL17_Full)) %>%
  mutate(cell_class = ifelse(SingleR_label %in% tnk_labels, "T/NK", "Myeloid"))



library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# ── 1. SUMMARY BY SINGLER SUBTYPE AND NICHE ──
singler_niche_scores <- immune_sub_filtered@meta.data %>%
  filter(!is.na(SingleR_label), !is.na(niches),
         SingleR_label %in% c(tnk_labels, myeloid_labels)) %>%
  group_by(SingleR_label, niches) %>%
  summarise(
    n          = n(),
    IL17_Mean  = mean(IL17_Full, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  filter(n >= 10)

# ── 2. BUILD MATRIX ──
mat_il17 <- singler_niche_scores %>%
  mutate(
    SingleR_label = factor(SingleR_label, levels = c(tnk_labels, myeloid_labels)),
    niches        = factor(niches, levels = niche_order)
  ) %>%
  dplyr::select(SingleR_label, niches, IL17_Mean) %>%
  pivot_wider(names_from = niches, values_from = IL17_Mean) %>%
  column_to_rownames("SingleR_label") %>%
  as.matrix()

# Keep only complete rows
mat_il17 <- mat_il17[complete.cases(mat_il17), ]

# Z-score rows
mat_il17_z <- t(scale(t(mat_il17)))

cat("Subtypes included:", nrow(mat_il17_z), "\n")
print(rownames(mat_il17_z))

# ── 3. ROW ANNOTATION ──
row_class <- ifelse(rownames(mat_il17_z) %in% tnk_labels, "T/NK", "Myeloid")
row_split <- factor(row_class, levels = c("T/NK", "Myeloid"))

row_anno <- rowAnnotation(
  Lineage              = row_class,
  col                  = list(Lineage = c("T/NK" = "#1f78b4", "Myeloid" = "#e6ab02")),
  show_annotation_name = FALSE,
  simple_anno_size     = unit(3, "mm")
)

# ── 4. TOP ANNOTATION ──
top_anno <- HeatmapAnnotation(
  Niche                = niche_order,
  col                  = list(Niche = niche_colors_manual),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  simple_anno_size     = unit(3, "mm")
)

# ── 5. COLOR SCALE ──
col_fun <- colorRamp2(c(-2, 0, 2), c("#4682B4", "white", "#b30000"))

# ── 6. HEATMAP ──
ht <- Heatmap(
  mat_il17_z,
  name             = "Z-score",
  col              = col_fun,
  top_annotation   = top_anno,
  left_annotation  = row_anno,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  row_split        = row_split,
  row_title_gp     = gpar(fontsize = 10, fontface = "bold"),
  show_row_names   = TRUE,
  row_names_gp     = gpar(fontsize = 8),
  column_labels    = niche_labels[colnames(mat_il17_z)],
  column_names_gp  = gpar(fontsize = 9, fontface = "bold", col = "black"),
  column_names_rot = 45,
  border           = FALSE,
  rect_gp          = gpar(col = NA),
  column_title     = "IL-17 Signaling — All 43 Genes",
  column_title_gp  = gpar(fontsize = 11, fontface = "bold")
)

# ── 7. SAVE ──
pdf("TNK_Myeloid_Subtypes_IL17Full_Heatmap.pdf", width = 6, height = 9)
draw(ht, merge_legend = TRUE, padding = unit(c(5, 15, 5, 5), "mm"))
dev.off()

cat("✅ IL-17 Full 43 genes subtype heatmap saved.\n")



library(ggplot2)
library(dplyr)

# ── 1. SUMMARY STATS FOR PLOTTING ──
n4_plot_data <- n4_data %>%
  group_by(SingleR_label, cell_class) %>%
  summarise(
    n      = n(),
    Mean   = mean(IL17_Full),
    SEM    = sd(IL17_Full) / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean)) %>%
  mutate(SingleR_label = factor(SingleR_label, levels = SingleR_label))
# (highest mean first, as in Fig 4C; a later reversed-order copy of this block was removed)

# ── 2. PLOT ──
p_n4_dot <- ggplot(n4_plot_data,
                   aes(y = Mean, x = SingleR_label, color = cell_class)) +
  
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.3, linewidth = 0.6) +
  
  geom_point(aes(size = n), alpha = 0.9) +
  
  geom_text(aes(y = Mean + SEM + 0.003,
                label = paste0("n=", n)),
            vjust = 0, hjust = 0.5, size = 2.8, color = "grey40") +
  
  scale_color_manual(values = c("T/NK" = "#1f78b4", "Myeloid" = "#e6ab02")) +
  
  scale_size_continuous(range  = c(2, 8),
                        breaks = c(50, 500, 2000, 5000),
                        name   = "Cell Count") +
  
  labs(
    title    = "IL-17 Signaling Score — N4: IE Niche",
    subtitle = "Mean ± SEM | Module score: 43 IL-17 pathway genes",
    y        = "IL-17 Module Score (Mean ± SEM)",
    x        = NULL,
    color    = "Lineage"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y   = element_text(size = 9),
    axis.title.y  = element_text(size = 10),
    legend.title  = element_text(size = 8, face = "bold"),
    legend.text   = element_text(size = 8),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line     = element_blank()
  )

ggsave("N4_IL17Full_Score_DotPlot.pdf",
       p_n4_dot, width = 7, height = 4.5)

## end of Part 1
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 2a - Fig 4D, 4E, 4H; Fig S4B, S4D: IL-17 hotspots, immune proximity, proximal vs distal CMs
##  (Clean/fig4/fig4_prox.R; writes xo_heart_processed_withspatialil17labels.rds)
##############################################################################


# ─── 0. Setup ─────────────────────────────────────────────────────────────────
setwd(file.path(TMOS_DIR, "Clean/fig4"))

library(Seurat)
library(dplyr)
library(tibble)
library(FNN)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(reactome.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(stringr)
library(patchwork)

xo_heart             <- readRDS("xo_heart_processed.rds")
tnk_myeloid_heart_st <- readRDS("tnk_myeloid_heart_st.rds")
xo_heart             <- fix_fov(xo_heart)
tnk_myeloid_heart_st <- fix_fov(tnk_myeloid_heart_st)




# Convert to character to avoid factor level errors
xo_heart$updated_cell_labels <- as.character(xo_heart$major_celltype)

# Extract labels as character vector
labels_to_add <- as.character(tnk_myeloid_heart_st$SingleR_label)
names(labels_to_add) <- rownames(tnk_myeloid_heart_st@meta.data)

# Transfer: overwrite using barcodes
xo_heart$updated_cell_labels[names(labels_to_add)] <- labels_to_add

# 2. Precise Removal of Contaminated Cells
# We keep only the specific SingleR subtypes AND the non-immune major celltypes
valid_subtypes <- unique(as.character(tnk_myeloid_heart_st$SingleR_label))
non_immune_types <- c("Endothelial", "Fibroblast_stromal", "Smoothe_muscle_mural", 
                      "Cardiomyocyte", "B_cells", "Cycling_immune", "Mixed")

# Logic: Keep if the label is one of the new subtypes OR one of the original structural types
keep_cells <- xo_heart$updated_cell_labels %in% c(valid_subtypes, non_immune_types)

xo_heart_cleaned <- subset(xo_heart, cells = colnames(xo_heart)[keep_cells])

# Verify
table(xo_heart_cleaned$updated_cell_labels)



library(FNN)

# Prepare coordinates
coords <- GetTissueCoordinates(xo_heart_cleaned, which = "centroids")
df_spatial <- cbind(coords, xo_heart_cleaned@meta.data[, c("updated_cell_labels", "sample_id")])

immune_subtypes <- unique(as.character(tnk_myeloid_heart_st$SingleR_label))
all_cell_types <- unique(df_spatial$updated_cell_labels)




# Add niches to the spatial dataframe
df_spatial <- cbind(coords, xo_heart_cleaned@meta.data[, c("updated_cell_labels", "sample_id", "niches")])

# ── 1. FORCE METADATA SYNC ──
# Explicitly ensure the column is in the metadata dataframe slot
xo_heart_cleaned@meta.data$updated_cell_labels <- as.character(xo_heart_cleaned$updated_cell_labels)

# ── 2. SUBSET AND COUNT (BASE R APPROACH) ──
# This avoids 'object not found' by accessing the columns directly via $
n4_rows <- xo_heart_cleaned@meta.data[xo_heart_cleaned@meta.data$niches == "N4", ]

# Create the count table and convert to a clean dataframe
actual_counts_n4 <- as.data.frame(table(n4_rows$updated_cell_labels))

# ── 3. RENAME FOR YOUR DOWNSTREAM PIPELINE ──
colnames(actual_counts_n4) <- c("query", "true_n")

# ── 4. VERIFY ──
# Sort by count so you can quickly check the top cell types
actual_counts_n4 <- actual_counts_n4[order(-actual_counts_n4$true_n), ]
print(head(actual_counts_n4, 10))

# Optional: verify total Niche 4 cell count
cat("\nTotal cells in Niche 4:", sum(actual_counts_n4$true_n), "\n")

dist_results_n4 <- list()

for (samp in unique(df_spatial$sample_id)) {
  samp_df <- df_spatial[df_spatial$sample_id == samp, ]
  
  # Queries: Strictly cells in Niche 4
  n4_immune_df <- samp_df %>% filter(niches == "N4", updated_cell_labels %in% immune_subtypes)
  
  for (immune_query in unique(n4_immune_df$updated_cell_labels)) {
    query_coords <- n4_immune_df[n4_immune_df$updated_cell_labels == immune_query, c("x", "y")]
    if (nrow(query_coords) < 5) next
    
    for (target_type in all_cell_types) {
      # Targets: All cells of this type in the sample (regardless of niche)
      target_coords <- samp_df[samp_df$updated_cell_labels == target_type, c("x", "y")]
      if (nrow(target_coords) < 5) next
      
      # FIX: Check if we are comparing a cell type to itself
      if (immune_query == target_type) {
        # Search for 2 neighbors; the 1st is the cell itself, the 2nd is the true neighbor
        nn <- get.knnx(data = target_coords, query = query_coords, k = 2)
        dist_val <- nn$nn.dist[, 2] 
      } else {
        # Search for the 1st neighbor
        nn <- get.knnx(data = target_coords, query = query_coords, k = 1)
        dist_val <- nn$nn.dist[, 1]
      }
      
      dist_results_n4[[length(dist_results_n4) + 1]] <- data.frame(
        sample = samp,
        query = immune_query,
        target = target_type,
        distance = dist_val
      )
    }
  }
}

final_dist_n4_df <- dplyr::bind_rows(dist_results_n4)




# (per-subtype proximity plots not used in the figures; removed)




library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# 1. PREPARE THE DATA
# Aggregate to median distance per pair
dist_summary <- final_dist_n4_df %>%
  group_by(query, target) %>%
  summarise(med_dist = median(distance), .groups = "drop")

# Pivot and convert to matrix
mat_dist <- dist_summary %>%
  pivot_wider(names_from = target, values_from = med_dist) %>%
  column_to_rownames("query") %>%
  as.matrix()

# 2. Z-SCORE SCALING
# We scale across targets (rows) to see which targets a specific query is "attracted" to
mat_dist_z <- t(scale(t(mat_dist)))



# ── 1. FINAL TARGET ORDERING PER SPECIFICATIONS ──

# ── 1. CALCULATE THE MISSING VECTOR ──
# Pull mean IL17 scores for the subtypes in Niche 4
il17_means_n4 <- tnk_myeloid_heart_st@meta.data %>%
  dplyr::filter(niches == "N4") %>% 
  group_by(SingleR_label) %>%
  summarise(mean_il17 = mean(IL17_Full, na.rm = TRUE)) %>%
  tibble::deframe() 

# Map the scores to the rownames of your distance matrix
# This creates the 'row_il17_vector' that the error is asking for
row_il17_vector <- il17_means_n4[rownames(mat_dist_z)]



# Highly Proximal (Group 1)
t1 <- c("Intermediate monocytes", "Classical monocytes", "Endothelial", 
        "Low-density neutrophils", "Natural killer cells", "Terminal effector CD8 T cells", 
        "Fibroblast_stromal", "Myeloid dendritic cells", "Cardiomyocyte", 
        "Effector memory CD8 T cells", "Cycling_immune", "Non classical monocytes")

# Moderately Proximal (Group 3/Cluster C - The 'bluer' slice)
t2 <- c("Mixed", "Non-Vd2 gd T cells", "MAIT cells", "Th1 cells", 
        "Th1/Th17 cells", "Terminal effector CD4 T cells", "Vd2 gd T cells")

# Distally Associated (Group 2/Cluster D - The 'whiter' slice)
t3 <- c("Th2 cells", "T regulatory cells", "B_cells", "Smoothe_muscle_mural", 
        "Th17 cells", "Central memory CD8 T cells", "Follicular helper T cells")

# Spatially Excluded (Group 4)
t4 <- c("Low-density basophils", "Plasmacytoid dendritic cells", "Progenitor cells", 
        "Naive CD4 T cells", "Naive CD8 T cells")

target_order <- c(t1, t2, t3, t4)
mat_dist_z_ordered <- mat_dist_z[, target_order]

# ── 2. SPLIT VECTOR ──
column_split_vec <- factor(rep(c("Highly Proximal", "Moderately Proximal", 
                                 "Distally Associated", "Spatially Excluded"), 
                               times = c(length(t1), length(t2), length(t3), length(t4))),
                           levels = c("Highly Proximal", "Moderately Proximal", 
                                      "Distally Associated", "Spatially Excluded"))

# (all-subtype heatmap not used in the figures; removed)


# ── 1. IDENTIFY THE TOP 5 IL-17 PRODUCERS (N > 100) ──
# We filter the metadata to identify the specific drivers of the story
top_il17_subtypes <- tnk_myeloid_heart_st@meta.data %>%
  dplyr::filter(niches == "N4") %>%
  group_by(SingleR_label) %>%
  summarise(
    mean_il17 = mean(IL17_Full, na.rm = TRUE),
    cell_count = n()
  ) %>%
  dplyr::filter(cell_count > 100) %>%
  slice_max(mean_il17, n = 5) %>%
  pull(SingleR_label)

# ── 2. FILTER DISTANCE MATRIX ──
# Subset the matrix to only include these 5 primary "query" cells
mat_dist_z_filtered <- mat_dist_z[rownames(mat_dist_z) %in% top_il17_subtypes, ]

# Ensure the columns follow your established 4-tier spatial logic
mat_dist_z_ordered <- mat_dist_z_filtered[, target_order]

# ── 3. RE-CALCULATE ANNOTATION FOR SUBSET ──
row_il17_vector_sub <- il17_means_n4[rownames(mat_dist_z_ordered)]

row_anno <- rowAnnotation(
  `IL17 Score` = row_il17_vector_sub, 
  col = list(`IL17 Score` = colorRamp2(
    c(min(row_il17_vector_sub, na.rm=TRUE), max(row_il17_vector_sub, na.rm=TRUE)), 
    c("white", "#008080") 
  )),
  show_annotation_name = TRUE,
  annotation_legend_param = list(`IL17 Score` = list(title = "Mean IL17"))
)

# ── 4. BUILD THE FOCUSED HEATMAP ──
col_fun <- colorRamp2(c(-2, 0, 2), c("#b30000", "white", "#4682B4"))

ht_focused <- Heatmap(
  mat_dist_z_ordered,
  name             = "Dist Z-score",
  col              = col_fun,
  left_annotation  = row_anno,        
  cluster_rows     = TRUE,
  cluster_columns  = FALSE, 
  column_split     = column_split_vec,
  
  # Styling for a clean, publication-ready look
  column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp     = gpar(fontsize = 10, fontface = "bold"), # Bold queries for emphasis
  column_names_gp  = gpar(fontsize = 9, fontface = "plain"), 
  column_names_rot = 45,
  border           = TRUE,
  rect_gp          = gpar(col = NA)
)

# ── 5. SAVE TO PDF ──
pdf("Niche4_Story_Focused_Heatmap.pdf", width = 11.5, height = 3.35)

draw(
  ht_focused, 
  merge_legend = TRUE, 
  column_title = "The Pathogenic Axis: Top IL-17 Producers & Spatial Targets", 
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  padding = unit(c(5, 10, 5, 15), "mm")
)

dev.off()

cat("✅ Focused 'Story' Heatmap generated with top 5 IL-17 producers.\n")



# (Set1-colour version of S4B not used; removed)



library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── 1. IDENTIFY TOP 5 IL-17 PRODUCERS (N > 100) ──
# Based on your previous criteria for highest mean IL17_Full
top_5_summary <- tnk_myeloid_heart_st@meta.data %>%
  dplyr::filter(niches == "N4") %>%
  group_by(SingleR_label) %>%
  summarise(
    mean_il17 = mean(IL17_Full, na.rm = TRUE),
    cell_count = n()
  ) %>%
  dplyr::filter(cell_count > 100) %>%
  slice_max(mean_il17, n = 5)

top_5_labels <- top_5_summary$SingleR_label

# ── 2. PREPARE AGGREGATED DATA ──
agg_data_top5 <- final_dist_n4_df %>%
  dplyr::filter(query %in% top_5_labels) %>%
  group_by(query, target) %>%
  summarise(med_dist = median(distance), .groups = "drop")

# ── 3. DEFINE COLORS & PLOT ──
nature_cols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C8DBCFF", "#7E6148FF")

p_joint <- ggplot(agg_data_top5, aes(x = med_dist, y = reorder(target, -med_dist))) +
  # Keep the previous formatting: thicker grey tracks
  geom_line(aes(group = target), color = "grey92", linewidth = 1.5) +
  # Add dots with requested Nature colors
  geom_point(aes(color = query), size = 4.5, alpha = 0.9) +
  scale_color_manual(values = nature_cols) +
  labs(
    x = "Median Distance (µm)",
    y = "Target Cell Type",
    color = "Top IL-17 Producers"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8.5),
    axis.text.x = element_text(size = 9),
    # Optimized legend for 6-inch width to prevent clipping
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  # Use 2-3 rows to keep text from running off the edge
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

# ── 4. SAVE TO PDF (6x8) ──
pdf("Niche4_Top5_Proximity_NatureColors.pdf", width = 5, height = 8)

p_final <- p_joint + 
  plot_annotation(
    title = "Spatial Convergence: Top IL-17 Producers",
    subtitle = "Aggregated median distance in Niche 4",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

print(p_final)
dev.off()

cat("✅ Joint Proximity Plot saved with Nature colors and 6x8 optimized formatting.\n")



# ─── 0. Setup ─────────────────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(tibble)
library(FNN)
library(ggplot2)
library(ggrepel)
library(ReactomePA)
library(org.Hs.eg.db)

setwd(file.path(TMOS_DIR, "Clean/fig4"))

# ─── 1. Identify IL17 Hubs (Top 25% within Niche 4) ──────────────────────────
il17_scores <- tnk_myeloid_heart_st@meta.data %>%
  rownames_to_column("cell_id") %>%
  filter(!is.na(IL17_Full) & niches == "N4") %>%
  dplyr::select(cell_id, sample_id, IL17_Full)

il17_thresholds <- il17_scores %>%
  group_by(sample_id) %>%
  summarise(threshold = quantile(IL17_Full, 0.75), .groups = "drop")

hub_cells <- il17_scores %>%
  left_join(il17_thresholds, by = "sample_id") %>%
  filter(IL17_Full >= threshold)

# ─── 2. Extract and Standardize Spatial Coordinates ──────────────────────────
all_coords <- as.data.frame(GetTissueCoordinates(xo_heart_cleaned))  # UPDATED

if ("cell" %in% colnames(all_coords)) {
  all_coords$cell_id <- all_coords[["cell"]]
} else {
  all_coords$cell_id <- rownames(all_coords)
}

x_col <- grep("x", colnames(all_coords), ignore.case = TRUE, value = TRUE)[1]
y_col <- grep("y", colnames(all_coords), ignore.case = TRUE, value = TRUE)[1]

all_coords$x <- all_coords[[x_col]]
all_coords$y <- all_coords[[y_col]]

all_coords <- all_coords[, c("cell_id", "x", "y")]

# ─── 3. Run the Per-Sample Distance Loop ─────────────────────────────────────
samples <- unique(xo_heart_cleaned$sample_id)  # UPDATED
all_cm_results <- list()

for (samp in samples) {
  cat("\n--- Processing sample:", as.character(samp), "---\n")
  
  # Filter Hubs for this sample
  samp_hub_ids <- hub_cells %>% filter(sample_id == samp) %>% pull(cell_id)
  coords_hub <- all_coords %>% filter(cell_id %in% samp_hub_ids) %>% na.omit()
  
  # Filter Cardiomyocytes for this sample
  samp_cm_ids <- rownames(xo_heart_cleaned@meta.data)[                           # UPDATED
    xo_heart_cleaned@meta.data$updated_cell_labels == "Cardiomyocyte" &  # UPDATED
      xo_heart_cleaned@meta.data$sample_id == samp]                 # UPDATED
  coords_cm <- all_coords %>% filter(cell_id %in% samp_cm_ids) %>% na.omit()
  
  if (nrow(coords_hub) == 0 | nrow(coords_cm) == 0) {
    cat("  Skipping: No Hubs or CMs found.\n")
    next
  }
  
  cat("  Calculating distances for", nrow(coords_cm), "CMs vs", nrow(coords_hub), "Hubs...\n")
  
  nn_res <- FNN::knnx.dist(data = as.matrix(coords_hub[, c("x", "y")]), 
                           query = as.matrix(coords_cm[, c("x", "y")]), k = 1)
  
  coords_cm$dist_to_hub <- nn_res[, 1]
  coords_cm$sample_id <- samp
  all_cm_results[[samp]] <- coords_cm
}

cat("\n✅ Distance loop finished for all samples.\n")

# ─── 4. Classify CMs into Close / Intermediate / Far ─────────────────────────
coords_cm_all <- bind_rows(all_cm_results) %>%
  group_by(sample_id) %>%
  mutate(
    close_thr = quantile(dist_to_hub, 0.25),
    far_thr   = quantile(dist_to_hub, 0.75),
    cm_group  = case_when(
      dist_to_hub <= close_thr ~ "Close_to_IL17_hub",
      dist_to_hub >= far_thr   ~ "Far_from_IL17_hub",
      TRUE                     ~ "Intermediate"
    )
  ) %>% ungroup()


# Map back CMs
xo_heart_cleaned$cm_IL17_group <- "Other"
xo_heart_cleaned@meta.data[coords_cm_all$cell_id, "cm_IL17_group"] <- coords_cm_all$cm_group

# ADD THIS: Label the hub cells
xo_heart_cleaned@meta.data[hub_cells$cell_id, "cm_IL17_group"] <- "IL17_Hub"


# ─── 5. Differential Expression: Close vs Far ────────────────────────────────
cm_sub <- subset(xo_heart_cleaned, subset = cm_IL17_group %in%                   # UPDATED
                   c("Close_to_IL17_hub", "Far_from_IL17_hub"))
Idents(cm_sub) <- "cm_IL17_group"

de_results <- FindMarkers(
  cm_sub,
  ident.1 = "Close_to_IL17_hub",
  ident.2 = "Far_from_IL17_hub",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.1,
  assay = "Xenium"
)

de_results$gene <- rownames(de_results)
write.csv(de_results, "DE_CM_Close_vs_Far.csv")

# ─── 6. Export Annotation CSV for Xenium Explorer ────────────────────────────
annotation_df <- data.frame(cell_id = rownames(xo_heart_cleaned@meta.data)) %>%  # UPDATED
  left_join(dplyr::select(coords_cm_all, cell_id, group = cm_group), by = "cell_id") %>%
  mutate(
    group = ifelse(cell_id %in% hub_cells$cell_id, "IL17_Hub", group),
    group = ifelse(is.na(group), "Other", group),
    color = case_when(
      group == "IL17_Hub"          ~ "#FF6B35",
      group == "Close_to_IL17_hub" ~ "#b30000",
      group == "Far_from_IL17_hub" ~ "#4682B4",
      TRUE                         ~ "#D3D3D3"
    )
  )

write.csv(annotation_df, "CM_IL17_Hub_CellAnnotation.csv", row.names = FALSE, quote = FALSE)




library(ggplot2)
library(dplyr)
library(scales)

# ── HISTOGRAM: Distance distribution colored by CM group, per sample ──────────

group_cols <- c(
  "Close_to_IL17_hub" = "#b30000",
  "Intermediate"      = "#d9d9d9",
  "Far_from_IL17_hub" = "#4682B4"
)

# ── BUILD dist_plot_df from coords_cm_all ──────────────────────────────────────
dist_plot_df <- coords_cm_all %>%
  dplyr::filter(cm_group %in% c("Close_to_IL17_hub", "Intermediate", "Far_from_IL17_hub")) %>%
  dplyr::select(cell_id, sample_id, dist_to_hub, cm_group) %>%
  dplyr::mutate(cm_group = factor(cm_group, 
                                  levels = c("Close_to_IL17_hub", 
                                             "Intermediate", 
                                             "Far_from_IL17_hub")))

# (pooled by-sample histograms not used; per-sample histograms below are S4D)





# ── BUILD supporting objects for the per-sample loop ──────────────────────────

# 1. Hub counts per sample
hub_counts <- hub_cells %>%
  group_by(sample_id) %>%
  summarise(n_hubs = n(), .groups = "drop") %>%
  as.data.frame()

# 2. Q25 / Q75 threshold lines per sample (for the vline annotations)
threshold_lines <- coords_cm_all %>%
  group_by(sample_id) %>%
  summarise(
    Q25 = quantile(dist_to_hub, 0.25),
    Q75 = quantile(dist_to_hub, 0.75),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c("Q25", "Q75"),
                      names_to  = "threshold",
                      values_to = "value")



# ── PER-SAMPLE HISTOGRAMS WITH PROPORTIONAL X AXIS ───────────────────────────
for (samp in unique(dist_plot_df$sample_id)) {
  
  samp_df         <- dist_plot_df %>% dplyr::filter(sample_id == samp)
  samp_thresholds <- threshold_lines %>% dplyr::filter(sample_id == samp)
  
  # Stats for subtitle
  n_close <- sum(samp_df$cm_group == "Close_to_IL17_hub")
  n_int   <- sum(samp_df$cm_group == "Intermediate")
  n_far   <- sum(samp_df$cm_group == "Far_from_IL17_hub")
  n_hubs  <- hub_counts$n_hubs[hub_counts$sample_id == samp]
  
  q25_val <- samp_thresholds$value[samp_thresholds$threshold == "Q25"]
  q75_val <- samp_thresholds$value[samp_thresholds$threshold == "Q75"]
  
  p_samp <- ggplot(samp_df, aes(x = dist_to_hub, fill = cm_group)) +
    geom_histogram(bins = 60, color = "white", linewidth = 0.15) +
    
    # Threshold lines
    geom_vline(
      data      = samp_thresholds,
      aes(xintercept = value, linetype = threshold),
      color     = "black",
      linewidth = 0.7
    ) +
    
    # Threshold value labels
    geom_text(
      data        = samp_thresholds,
      aes(x       = value, 
          y       = Inf, 
          label   = paste0(threshold, "\n", round(value, 1), " µm")),
      vjust       = 1.3,
      hjust       = -0.1,
      size        = 3.5,
      color       = "black",
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(
      values = group_cols,
      labels = c(
        "Close_to_IL17_hub" = paste0("Close ≤Q25 (n=", n_close, ")"),
        "Intermediate"      = paste0("Intermediate (n=", n_int, ")"),
        "Far_from_IL17_hub" = paste0("Far ≥Q75 (n=", n_far, ")")
      )
    ) +
    scale_linetype_manual(
      values = c("Q25" = "dashed", "Q75" = "dotted"),
      labels = c("Q25" = paste0("Q25 = ", round(q25_val, 1), " µm"),
                 "Q75" = paste0("Q75 = ", round(q75_val, 1), " µm"))
    ) +
    
    # Proportional x axis to actual data range of this sample
    scale_x_continuous(
      labels = comma,
      limits = c(0, max(samp_df$dist_to_hub) * 1.05)
    ) +
    scale_y_continuous(labels = comma) +
    
    labs(
      title    = paste0("CM Distance to Nearest IL-17 Hub — ", samp),
      subtitle = paste0("Total CMs: ", format(nrow(samp_df), big.mark = ","),
                        "  |  IL-17 Hubs: ", format(n_hubs, big.mark = ","),
                        "  |  Q25: ", round(q25_val, 1), " µm",
                        "  |  Q75: ", round(q75_val, 1), " µm"),
      x        = "Distance to Nearest IL-17 Hub (µm)",
      y        = "Number of Cardiomyocytes",
      fill     = "CM Group",
      linetype = "Threshold"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(color = "grey40", size = 9),
      legend.position  = "bottom",
      legend.box       = "vertical",
      panel.grid.minor = element_blank()
    )
  
  filename <- paste0("CM_Distance_Histogram_", samp, ".pdf")
  
  pdf(filename, width = 8, height = 6)
  print(p_samp)
  dev.off()
  
  cat("✅ Saved:", filename, "\n")
}



# ─── Volcano Plot: Close vs Far Cardiomyocytes ────────────────────────────────

# ── 1. PREPARE DATA ──
volcano_df <- de_results %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),
    
    direction = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up in Close",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Up in Far",
      TRUE                                  ~ "NS"
    ),
    
    label = case_when(
      direction != "NS" ~ gene,
      TRUE              ~ NA_character_
    )
  )

# ── 2. COUNTS FOR LEGEND ──
n_up_close <- sum(volcano_df$direction == "Up in Close")
n_up_far   <- sum(volcano_df$direction == "Up in Far")

cat("Up in Close:", n_up_close, "\n")
cat("Up in Far:  ", n_up_far,   "\n")
cat("Labeled:    ", sum(!is.na(volcano_df$label)), "\n")

# ── 3. COLORS ──
cols <- c(
  "Up in Close" = "#b30000",
  "Up in Far"   = "#4682B4",
  "NS"          = "grey60"
)

# ── 4. PLOT ──
p_volcano <- ggplot(volcano_df, aes(x = avg_log2FC, y = neg_log10_p, color = direction)) +
  
  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5),  linetype = "dashed", color = "grey50", linewidth = 0.4) +
  
  # NS points first so significant sit on top
  # Points — make bigger
  geom_point(data = filter(volcano_df, direction == "NS"), size = 2.0, alpha = 0.4) +   # was 1.2
  geom_point(data = filter(volcano_df, direction != "NS"), size = 3.5, alpha = 0.9) +   # was 2.0
  
  # Labels — make font bigger
  geom_text_repel(
    aes(label = label),
    size               = 5.0,    # was 2.5
    fontface           = "italic",
    color              = "black",
    box.padding        = 0.3,
    point.padding      = 0.2,
    max.overlaps       = Inf,
    force              = 2,
    force_pull         = 0.5,
    min.segment.length = 0,
    segment.color      = "grey60",
    segment.size       = 0.3,
    na.rm              = TRUE
  ) +
  
  # Colors
  scale_color_manual(
    values = cols,
    labels = c(
      "Up in Close" = paste0("Up near IL17 hub (n=", n_up_close, ")"),
      "Up in Far"   = paste0("Up far from hub (n=",  n_up_far,   ")"),
      "NS"          = "Not significant"
    )
  ) +
  
  # Axis labels
  labs(
    title    = "Cardiomyocyte DE: Proximity to IL-17 Hubs",
    subtitle = "Close vs Far cardiomyocytes | Xenium spatial transcriptomics",
    x        = expression("Average log"[2]*" Fold Change"),
    y        = expression("-log"[10]*" Adjusted P-value"),
    color    = NULL
  ) +
  
  # Theme
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey95")
  ) +
  
  scale_x_continuous(limits = c(
    min(volcano_df$avg_log2FC) * 1.1,
    max(volcano_df$avg_log2FC) * 1.1
  ))

# ── 5. SAVE ──
pdf("Volcano_CM_Close_vs_Far_IL17Hub.pdf", width = 6, height = 8)
print(p_volcano)
dev.off()

cat("✅ Volcano plot saved.\n")



library(ReactomePA)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

# ─── 1. EXTRACT GENE LISTS FROM VOLCANO ───────────────────────────────────────
genes_up_close <- volcano_df %>% filter(direction == "Up in Close") %>% pull(gene)
genes_up_far   <- volcano_df %>% filter(direction == "Up in Far")   %>% pull(gene)

cat("Genes Up in Close:", length(genes_up_close), "\n")
cat("Genes Up in Far:  ", length(genes_up_far),   "\n")

# ─── 2. CONVERT GENE SYMBOLS TO ENTREZ IDS ────────────────────────────────────
convert_to_entrez <- function(gene_symbols) {
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys    = gene_symbols,
    column  = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  entrez <- entrez[!is.na(entrez)]
  cat("  Mapped:", length(entrez), "of", length(gene_symbols), "genes to Entrez IDs\n")
  return(entrez)
}

entrez_up_close <- convert_to_entrez(genes_up_close)
entrez_up_far   <- convert_to_entrez(genes_up_far)

# ─── 3. RUN REACTOME ORA ──────────────────────────────────────────────────────
run_reactome <- function(entrez_ids, label) {
  if (length(entrez_ids) < 3) {
    cat("  Skipping", label, "— too few genes\n")
    return(NULL)
  }
  
  res <- enrichPathway(
    gene          = entrez_ids,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 5,      # was 10
    maxGSSize     = 1000,   # UPDATED — was 500, catches broad terms like "Immune System"
    readable      = TRUE
  )
  
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat("  No significant pathways for", label, "\n")
    return(NULL)
  }
  
  df <- as.data.frame(res) %>%
    mutate(FoldEnrichment = as.numeric(sub("/.*", "", GeneRatio)) /
             as.numeric(sub("/.*", "", BgRatio)))
  
  cat("  Found", nrow(df), "significant pathways for", label, "\n")
  return(df)
}

df_up_close <- run_reactome(entrez_up_close, "Up in Close")
df_up_far   <- run_reactome(entrez_up_far,   "Up in Far")

# ─── 4. PLOTTING FUNCTIONS (YOUR EXISTING STYLE) ──────────────────────────────
prep_bar_df <- function(df, wrap_width = 55) {
  df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(
      neglog10_padj = -log10(p.adjust),
      Pathway_wrap  = str_wrap(Description, width = wrap_width),
      Pathway_wrap  = factor(Pathway_wrap, levels = rev(Pathway_wrap))
    )
}


# ─── 5. MAKE PLOTS ────────────────────────────────────────────────────────────
col_close <- "#b30000"   # Red  — Up in Close (near IL17 hub)
col_far   <- "#4682B4"   # Blue — Up in Far

# (gradient-fill versions of 4H not used; solid versions below are 4H)




make_bar_solid <- function(df, title_text, base_color) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df2 <- prep_bar_df(df)
  
  ggplot(df2, aes(x = neglog10_padj, y = Pathway_wrap)) +
    geom_col(width = 0.86, fill = base_color) +
    geom_text(
      aes(label = Pathway_wrap),
      x = 0.02, hjust = 0,
      size = 4.2, color = "white"    # white text reads better on solid color
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = title_text,
      x     = expression(-log[10]("adjusted p-value")),
      y     = NULL
    )
}

# ── MAKE AND SAVE ──────────────────────────────────────────────────────────────
p_close_solid <- make_bar_solid(
  df_up_close,
  "Reactome ORA — Up near IL17 Hub (Close CMs)",
  col_close
)

p_far_solid <- make_bar_solid(
  df_up_far,
  "Reactome ORA — Up far from IL17 Hub (Far CMs)",
  col_far
)

if (!is.null(p_close_solid)) {
  ggsave(
    filename = "Reactome_ORA_Top10_Up_Close_IL17Hub_solid.pdf",
    plot     = p_close_solid,
    width    = 4.5,
    height   = 5,
    device   = cairo_pdf
  )
  cat("✅ Saved: Reactome_ORA_Top10_Up_Close_IL17Hub_solid.pdf\n")
}

if (!is.null(p_far_solid)) {
  ggsave(
    filename = "Reactome_ORA_Top10_Up_Far_IL17Hub_solid.pdf",
    plot     = p_far_solid,
    width    = 4.5,
    height   = 5,
    device   = cairo_pdf
  )
  cat("✅ Saved: Reactome_ORA_Top10_Up_Far_IL17Hub_solid.pdf\n")
}




library(openxlsx)

# ─── SAVE REACTOME RESULTS TO EXCEL ──────────────────────────────────────────
wb <- createWorkbook()

if (!is.null(df_up_close) && nrow(df_up_close) > 0) {
  addWorksheet(wb, "Up_Close_IL17Hub")
  writeData(wb, "Up_Close_IL17Hub", df_up_close)
}

if (!is.null(df_up_far) && nrow(df_up_far) > 0) {
  addWorksheet(wb, "Up_Far_IL17Hub")
  writeData(wb, "Up_Far_IL17Hub", df_up_far)
}

saveWorkbook(wb, "Reactome_ORA_CM_Close_vs_Far.xlsx", overwrite = TRUE)
cat("✅ Saved: Reactome_ORA_CM_Close_vs_Far.xlsx\n")



saveRDS(xo_heart_cleaned, "xo_heart_processed_withspatialil17labels.rds")


# ─── STACKED BAR: NICHE 4 MAJOR CELL TYPE COMPOSITION ────────────────────────

n4_major <- xo_heart_cleaned@meta.data %>%
  dplyr::filter(niches == "N4") %>%
  dplyr::count(major_celltype) %>%
  dplyr::mutate(
    major_celltype = as.character(major_celltype),
    pct            = n / sum(n)
  ) %>%
  dplyr::filter(major_celltype != "Mixed")  # REMOVE Mixed

# Recalculate pct after removing Mixed
n4_major <- n4_major %>%
  dplyr::mutate(pct = n / sum(n),
                label_txt = ifelse(pct >= 0.03,
                                   scales::percent(pct, accuracy = 0.1),
                                   ""))

cols_major <- c(
  "Cardiomyocyte"        = "#d73027",
  "Endothelial"          = "#1b9e77",
  "Fibroblast_stromal"   = "#7570b3",
  "Smoothe_muscle_mural" = "#66a61e",
  "T_NK_lymphocytes"     = "#1f78b4",
  "Myeloid"              = "#e6ab02",
  "B_cells"              = "#d95f02",
  "Cycling_immune"       = "#118ab2"
)

p_bar_n4 <- ggplot(n4_major, aes(x = "", y = pct, fill = major_celltype)) +
  geom_col(width = 0.5) +
  geom_text(
    aes(label = label_txt),
    position   = position_stack(vjust = 0.5),
    size       = 3.5,
    color      = "white",
    fontface   = "bold",
    lineheight = 0.9
  ) +
  scale_fill_manual(values = cols_major) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "Cell Type Composition — Niche 4",
    subtitle = paste0("n = ", format(sum(n4_major$n), big.mark = ",")),
    x        = NULL,
    y        = "Proportion",
    fill     = NULL
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40"),
    legend.position = "right",
    legend.text     = element_text(size = 9),
    plot.margin     = margin(10, 20, 10, 20)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave("StackedBar_N4_MajorCellType.pdf",
       plot   = p_bar_n4,
       width  = 4,
       height = 7,
       device = cairo_pdf)

cat("✅ Saved: StackedBar_N4_MajorCellType.pdf\n")


## end of Part 2a
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 2b - Fig S4H: CM signature gradient by distance to IL-17 hubs
##  (Clean/fig4/fig4_prox.R, last section; uses the object written by Part 2a)
##############################################################################

setwd(file.path(TMOS_DIR, "Clean/fig4"))
library(Seurat)
library(dplyr)
library(FNN)
library(openxlsx)
xo_heart_cleaned <- fix_fov(readRDS("xo_heart_processed_withspatialil17labels.rds"))

# CM distance to the nearest IL-17 hub cell in the same section
# (same calculation as Part 2a, from the saved hub labels)
all_coords <- as.data.frame(GetTissueCoordinates(xo_heart_cleaned))
all_coords$cell_id <- if ("cell" %in% colnames(all_coords)) all_coords[["cell"]] else rownames(all_coords)
x_col <- grep("x", colnames(all_coords), ignore.case = TRUE, value = TRUE)[1]
y_col <- grep("y", colnames(all_coords), ignore.case = TRUE, value = TRUE)[1]
all_coords$x <- all_coords[[x_col]]
all_coords$y <- all_coords[[y_col]]
all_coords <- all_coords[, c("cell_id", "x", "y")]

md <- xo_heart_cleaned@meta.data
all_cm_results <- list()
for (samp in unique(xo_heart_cleaned$sample_id)) {
  hub_ids <- rownames(md)[md$cm_IL17_group %in% "IL17_Hub" & md$sample_id == samp]
  cm_ids  <- rownames(md)[md$updated_cell_labels == "Cardiomyocyte" & md$sample_id == samp]
  coords_hub <- all_coords %>% filter(cell_id %in% hub_ids) %>% na.omit()
  coords_cm  <- all_coords %>% filter(cell_id %in% cm_ids) %>% na.omit()
  if (nrow(coords_hub) == 0 | nrow(coords_cm) == 0) next
  coords_cm$dist_to_hub <- FNN::knnx.dist(data  = as.matrix(coords_hub[, c("x", "y")]),
                                          query = as.matrix(coords_cm[, c("x", "y")]), k = 1)[, 1]
  coords_cm$sample_id <- samp
  all_cm_results[[as.character(samp)]] <- coords_cm
}
coords_cm_all <- bind_rows(all_cm_results)
rm(all_coords, md); invisible(gc())



library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(msigdbr)
library(stringr)

# ── 1. BUILD SIGNATURES FROM MSIGDB (OBJECTIVE) ───────────────────────────────
panel_genes <- rownames(xo_heart_cleaned[["Xenium"]])

# Hallmark sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

sets_of_interest <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

signatures_hallmark <- hallmark_sets %>%
  filter(gs_name %in% sets_of_interest) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol), .groups = "drop") %>%
  tibble::deframe()

# Reactome muscle/cardiac sets
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

muscle_sigs <- c2_sets %>%
  filter(gs_name %in% c("REACTOME_STRIATED_MUSCLE_CONTRACTION")) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol), .groups = "drop") %>%
  tibble::deframe()

# Combine and filter to Xenium panel
all_sigs_raw <- c(signatures_hallmark, muscle_sigs)
all_sigs <- lapply(all_sigs_raw, function(g) g[g %in% panel_genes])

# Remove sets with fewer than 3 panel genes
all_sigs <- all_sigs[sapply(all_sigs, length) >= 3]

# Clean names
names(all_sigs) <- names(all_sigs) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("REACTOME_", "", .) %>%
  gsub("_", " ", .) %>%
  str_to_title()

cat("Final signatures:\n")
print(sapply(all_sigs, length))

# ── 2. ENSURE dist_to_hub IS IN METADATA ─────────────────────────────────────
xo_heart_cleaned@meta.data$dist_to_hub <- NA_real_
xo_heart_cleaned@meta.data[coords_cm_all$cell_id, "dist_to_hub"] <- coords_cm_all$dist_to_hub

# ── 3. SUBSET TO CMs ──────────────────────────────────────────────────────────
cm_cells <- rownames(xo_heart_cleaned@meta.data)[
  xo_heart_cleaned@meta.data$updated_cell_labels == "Cardiomyocyte" &
    !is.na(xo_heart_cleaned@meta.data$dist_to_hub)
]
cat("CM cells:", length(cm_cells), "\n")
cm_obj <- subset(xo_heart_cleaned, cells = cm_cells)

# ── 4. SCORE WITH AddModuleScore ──────────────────────────────────────────────
# Remove any old score columns first
old_cols <- grep("1$", colnames(cm_obj@meta.data), value = TRUE)
if (length(old_cols) > 0) {
  cm_obj@meta.data <- cm_obj@meta.data[, !colnames(cm_obj@meta.data) %in% old_cols]
}

for (sig_name in names(all_sigs)) {
  cm_obj <- AddModuleScore(
    cm_obj,
    features = list(all_sigs[[sig_name]]),
    assay    = "Xenium",
    name     = sig_name
  )
}

# Safely get score columns
score_cols <- tail(colnames(cm_obj@meta.data), length(all_sigs))
cat("\nScore columns:\n")
print(score_cols)

# ── 5. BUILD SCORE MATRIX SORTED BY DISTANCE ─────────────────────────────────
cm_meta <- cm_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  arrange(dist_to_hub) %>%
  mutate(sample_id = factor(sample_id,
                            levels = c("Heart_153_A", "Heart_153_B", "Heart_155_A")))

score_mat <- cm_meta[, score_cols, drop = FALSE] %>%
  as.matrix() %>%
  t()

rownames(score_mat) <- names(all_sigs)
cat("\nMatrix dimensions:", nrow(score_mat), "x", ncol(score_mat), "\n")

# ── 6. Z-SCORE AND CLIP ───────────────────────────────────────────────────────
score_mat_z <- t(scale(t(score_mat)))
score_mat_z <- pmin(pmax(score_mat_z, -2), 2)

# ── 7. SPEARMAN CORRELATIONS ──────────────────────────────────────────────────
cor_with_dist <- apply(score_mat_z, 1, function(x) {
  cor(x, cm_meta$dist_to_hub, method = "spearman", use = "complete.obs")
})

cat("\nSpearman correlations with distance:\n")
print(round(sort(cor_with_dist), 3))

# ── 8. SELECT SIGNATURES WITH MEANINGFUL ASSOCIATION ─────────────────────────
# Use |rho| >= 0.10 — more lenient since these are unbiased gene sets
# which will naturally have lower correlations than hand-picked genes
keep_sigs <- names(cor_with_dist[abs(cor_with_dist) >= 0.10])
cat("\nKeeping (|rho| >= 0.10):", paste(keep_sigs, collapse = ", "), "\n")

# Order: strongest negative first, then positive
htl <- names(sort(cor_with_dist[keep_sigs][cor_with_dist[keep_sigs] <= 0]))
lth <- names(sort(cor_with_dist[keep_sigs][cor_with_dist[keep_sigs] >  0],
                  decreasing = TRUE))

sig_order    <- c(htl, lth)
score_sorted <- score_mat_z[sig_order, , drop = FALSE]

row_split <- factor(
  ifelse(cor_with_dist[sig_order] <= 0, "High to Low", "Low to High"),
  levels = c("High to Low", "Low to High")
)

cat("\nHigh to Low:", paste(htl, collapse = ", "), "\n")
cat("Low to High:", paste(lth, collapse = ", "), "\n")

# ── 9. SAMPLE POSITIONS FOR TEXT LABELS ──────────────────────────────────────
sample_positions <- cm_meta %>%
  mutate(col_idx = row_number()) %>%
  group_by(sample_id) %>%
  summarise(mid = round(mean(col_idx)), .groups = "drop")

# ── 10. ANNOTATIONS ───────────────────────────────────────────────────────────
library(viridis)

# ── 10. ANNOTATIONS ───────────────────────────────────────────────────────────
# Distance bar using viridis
dist_vals <- cm_meta$dist_to_hub
viridis_cols <- viridis(100)

dist_anno <- HeatmapAnnotation(
  `Distance` = dist_vals,
  col = list(`Distance` = colorRamp2(
    seq(min(dist_vals), max(dist_vals), length.out = 100),
    viridis_cols
  )),
  show_legend          = TRUE,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  height               = unit(8, "mm"),
  annotation_legend_param = list(
    Distance = list(
      title     = "Distance to\nIL-17 Hub",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
)

# ── 11. COLOR SCALE ───────────────────────────────────────────────────────────
col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#4682B4", "lightcyan", "white", "orange", "#b30000")
)

# ── 12. BUILD HEATMAP ─────────────────────────────────────────────────────────
# ── 12. BUILD HEATMAP ─────────────────────────────────────────────────────────
ht_final <- Heatmap(
  score_sorted,
  name             = "Signature\nZ-score",
  col              = col_fun,
  top_annotation   = dist_anno,
  
  cluster_columns  = FALSE,
  cluster_rows     = FALSE,
  
  row_split        = row_split,
  row_title_gp     = gpar(fontsize = 11, fontface = "bold"),
  row_title_rot    = 0,
  row_gap          = unit(4, "mm"),
  
  show_column_names = FALSE,
  row_names_gp     = gpar(fontsize = 10),
  row_names_side   = "right",
  
  right_annotation = rowAnnotation(
    ` ` = anno_text(
      paste0("r = ", round(cor_with_dist[sig_order], 2)),
      gp   = gpar(fontsize = 9, fontface = "italic"),
      just = "left"
    )
  ),
  
  border          = TRUE,
  use_raster      = TRUE,
  raster_quality  = 5,
  
  column_title      = "Proximal                              Distal",
  column_title_side = "bottom",
  column_title_gp   = gpar(fontsize = 10, fontface = "italic")
)

# ── 13. SAVE ──────────────────────────────────────────────────────────────────
pdf("CM_SignatureGradient_MSigDB_Final.pdf", width = 10, height = 2.5)
draw(
  ht_final,
  merge_legend = FALSE,
  padding      = unit(c(5, 45, 10, 25), "mm")
)
dev.off()

cat("✅ Saved: CM_SignatureGradient_MSigDB_Final.pdf\n")

# ── 14. SAVE SUMMARY TABLE ────────────────────────────────────────────────────
cor_df <- data.frame(
  Signature      = names(cor_with_dist),
  N_genes_panel  = sapply(all_sigs[names(cor_with_dist)], length),
  Spearman_rho   = round(cor_with_dist, 4),
  Direction      = ifelse(cor_with_dist <= 0, "High to Low", "Low to High"),
  Kept           = names(cor_with_dist) %in% keep_sigs
) %>% arrange(Spearman_rho)

wb3 <- createWorkbook()
addWorksheet(wb3, "MSigDB_Signature_Correlations")
writeData(wb3, "MSigDB_Signature_Correlations", cor_df)
setColWidths(wb3, "MSigDB_Signature_Correlations", cols = 1:5, widths = "auto")
saveWorkbook(wb3, "CM_SignatureGradient_MSigDB_Correlations.xlsx", overwrite = TRUE)

cat("✅ Saved: CM_SignatureGradient_MSigDB_Correlations.xlsx\n")



## end of Part 2b
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 3 - Fig 4F: per-gene gradient in cardiomyocytes vs distance to IL-17 hubs
##  (Clean/reviewer_comments/rev2.3_spatial.R)
##############################################################################


## =====================================================================
## R2.3 – Per-gene gradient in cardiomyocytes vs distance to IL-17 hubs
## 14-gene IL-17 signaling–related core + IL-17 receptor genes
## All CMs pooled across sections, ordered closest -> farthest
## Display: rolling mean along distance; % of CMs expressing each gene
## =====================================================================
library(Seurat)
library(dplyr)
library(tibble)
library(FNN)
library(zoo)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(viridis)
library(openxlsx)

## ── 0. LOAD XENIUM OBJECT ────────────────────────────────────────────────
base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.3_CM_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

# written by Part 2 (an identical copy is kept in Clean/reviewer_comments)
xo_file <- file.path(TMOS_DIR, "Clean/fig4", "xo_heart_processed_withspatialil17labels.rds")
if (!file.exists(xo_file)) {
  print(list.files(base_dir, pattern = "\\.rds$"))
  stop("Set xo_file to the Xenium heart .rds listed above")
}
xo <- fix_fov(readRDS(xo_file))
DefaultAssay(xo) <- "Xenium"

## ── 1. GENES ─────────────────────────────────────────────────────────────
core14 <- c("TAB2","MAPK8","UBE2V1","NFKB1","MAP2K4","BTRC","CREB1",
            "PPP2CB","RPS6KA5","MEF2A","RPS6KA3","FBXW11","MAP2K1","DUSP4")
il17r  <- c("IL17RA","IL17RB","IL17RC","IL17RD","IL17RE","TRAF3IP2")

panel  <- rownames(xo[["Xenium"]])
g_core <- intersect(core14, panel)
g_rec  <- intersect(il17r, panel)
cat("Core genes on panel:    ", paste(g_core, collapse = ", "), "\n")
cat("Receptor genes on panel:", paste(g_rec,  collapse = ", "), "\n")
cat("Not on panel:           ", paste(setdiff(c(core14, il17r), panel), collapse = ", "), "\n")
genes <- c(g_core, g_rec)

## ── 2. DISTANCE OF EACH CM TO NEAREST IL-17 HUB (hubs within same section) ──
coords <- as.data.frame(GetTissueCoordinates(xo))
coords$cell_id <- if ("cell" %in% colnames(coords)) coords[["cell"]] else rownames(coords)
x_col <- grep("x", colnames(coords), ignore.case = TRUE, value = TRUE)[1]
y_col <- grep("y", colnames(coords), ignore.case = TRUE, value = TRUE)[1]
coords$x <- coords[[x_col]]; coords$y <- coords[[y_col]]
coords <- coords[, c("cell_id", "x", "y")]

meta <- xo@meta.data %>% rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, sample_id, updated_cell_labels, cm_IL17_group) %>%
  left_join(coords, by = "cell_id")

dist_list <- list()
for (samp in unique(meta$sample_id)) {
  hubs <- meta %>% filter(sample_id == samp, cm_IL17_group == "IL17_Hub") %>%
    dplyr::select(x, y) %>% na.omit()
  cms  <- meta %>% filter(sample_id == samp, updated_cell_labels == "Cardiomyocyte") %>%
    dplyr::select(cell_id, sample_id, x, y) %>% na.omit()
  if (nrow(hubs) == 0 || nrow(cms) == 0) next
  cms$dist_to_hub <- FNN::knnx.dist(data = as.matrix(hubs), query = as.matrix(cms[, c("x", "y")]), k = 1)[, 1]
  dist_list[[samp]] <- cms
}

## all CMs pooled, ordered closest -> farthest
cm_df <- bind_rows(dist_list) %>%
  mutate(sample_id = factor(sample_id, levels = c("Heart_153_A", "Heart_153_B", "Heart_155_A"))) %>%
  arrange(dist_to_hub)
cat("CMs with distance:", nrow(cm_df), "\n")

## ── 3. EXPRESSION FOR THESE GENES ONLY ───────────────────────────────────
expr <- t(as.matrix(LayerData(xo, assay = "Xenium", layer = "data")[genes, cm_df$cell_id, drop = FALSE]))
rm(xo, coords, meta, dist_list); gc()

## ── 4. % OF CMs EXPRESSING EACH GENE ────────────────────────────────────
pct_expr <- colMeans(expr > 0) * 100

## ── 5. DISPLAY MATRIX: rolling mean along pooled distance order, z-score ──
win <- 501
sm  <- apply(expr, 2, function(x) zoo::rollmean(x, k = min(win, length(x)), fill = "extend"))
mat_z <- t(scale(sm))
mat_z[is.na(mat_z)] <- 0
mat_z <- pmin(pmax(mat_z, -2), 2)

## ── 6. ROW ORDER / SPLIT ────────────────────────────────────────────────
ord_core   <- g_core
ord_rec    <- g_rec
gene_order <- c(ord_core, ord_rec)
row_split  <- factor(c(rep("14-gene core", length(ord_core)), rep("IL-17 receptor", length(ord_rec))),
                     levels = c("14-gene core", "IL-17 receptor"))

## ── 7. ANNOTATIONS ───────────────────────────────────────────────────────
dist_vals <- cm_df$dist_to_hub

top_anno <- HeatmapAnnotation(
  Distance = dist_vals,
  col = list(Distance = colorRamp2(seq(min(dist_vals), max(dist_vals), length.out = 100), viridis(100))),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size     = unit(6, "mm"),
  annotation_legend_param = list(Distance = list(title = "Distance to\nIL-17 hub (µm)")))

lab <- sprintf("%.1f%% CMs", pct_expr[gene_order])
right_anno <- rowAnnotation(
  r = anno_text(lab, gp = gpar(fontsize = 8.5, fontface = "italic"), just = "left"),
  show_annotation_name = FALSE)

col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("#4682B4", "lightcyan", "white", "orange", "#b30000"))

## ── 8. HEATMAP ───────────────────────────────────────────────────────────
ht <- Heatmap(
  mat_z[gene_order, , drop = FALSE],
  name = "Expression\nZ-score\n(smoothed)", col = col_fun,
  top_annotation = top_anno,
  cluster_columns = FALSE, cluster_rows = FALSE,
  row_split = row_split, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 0, row_gap = unit(3, "mm"),
  show_column_names = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  border = TRUE, use_raster = TRUE, raster_quality = 5,
  column_title = paste0("Proximal  →  Distal (all ", format(nrow(cm_df), big.mark = ","),
                        " cardiomyocytes, ordered by distance to nearest IL-17 hub)"),
  column_title_side = "bottom", column_title_gp = gpar(fontsize = 10, fontface = "italic"))

h_in <- 2 + 0.3 * length(gene_order)
pdf("CM_GeneGradient_Core14_IL17R_pooled.pdf", width = 12, height = h_in)
draw(ht, merge_legend = TRUE, padding = unit(c(5, 5, 8, 5), "mm")); dev.off()

## ── 9. TABLE ─────────────────────────────────────────────────────────────
gene_df <- data.frame(gene = gene_order, set = as.character(row_split),
                      pct_CM_expressing = round(pct_expr[gene_order], 2),
                      check.names = FALSE, row.names = NULL)
print(gene_df)

wb <- createWorkbook()
addWorksheet(wb, "Gene_gradient"); writeData(wb, "Gene_gradient", gene_df)
saveWorkbook(wb, "CM_GeneGradient_Core14_IL17R_pooled.xlsx", overwrite = TRUE)
cat("Saved: CM_GeneGradient_Core14_IL17R_pooled.pdf and .xlsx\n")



## end of Part 3
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 4a - Fig S5A, S5B: Close vs Far CMs stratified by immune density and inflammation
##  (Clean/reviewer_comments/spatial_adjustment.R)
##############################################################################


# =========================================================
# R1.3: CM proximity to IL-17 hubs — adjustment approaches A, B, C
# Main figures: all sections pooled | Supplementary: per section
# Standalone
# =========================================================
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(FNN)
library(dbscan)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(openxlsx)

setwd(file.path(TMOS_DIR, "Clean/reviewer_comments"))
out_dir <- "spatial_adjust"; dir.create(out_dir, showWarnings = FALSE)

# ---- settings ----
r_um         <- 50   # radius for immune density (µm)
k_nbr        <- 15   # nearest non-CM cells for neighborhood inflammation
immune_types <- c("T_NK_lymphocytes", "Myeloid", "B_cells", "Cycling_immune")
grp_cols     <- c("Close" = "#b30000", "Far" = "#4682B4")

il17_module <- c("NOD1","TAB1","IL17F","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","TAB2",
                 "FBXW11","IL17RA","IKBKB","IL17A","IRAK1","IRAK2","MEF2A","MEF2C","ATF1","NFKB1",
                 "PPP2CA","PPP2R1A","IL17RB","MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
                 "MAP2K1","NOD2","MAP2K4","IL25","SKP1","MAP3K7","TRAF6","CUL1","IKBKG","RIPK2",
                 "BTRC","RPS6KA5","MAPKAPK2")

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "p<0.001", paste0("p=", formatC(p, format = "f", digits = 3))))

# =========================================================
# 1) Load: counts + existing normalized data + metadata, then drop the object
# =========================================================
# written by Part 2 (an identical copy is kept in Clean/reviewer_comments)
xen      <- readRDS(file.path(TMOS_DIR, "Clean/fig4", "xo_heart_processed_withspatialil17labels.rds"))
counts   <- LayerData(xen, assay = "Xenium", layer = "counts")
data_mat <- LayerData(xen, assay = "Xenium", layer = "data")
meta <- xen@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, sample_id, updated_cell_labels, major_celltype_broad,
                cm_IL17_group, niches, x = x_centroid, y = y_centroid)
rm(xen); gc()

meta <- meta %>%
  mutate(is_cm     = updated_cell_labels == "Cardiomyocyte",
         is_immune = major_celltype_broad %in% immune_types,
         is_hub    = cm_IL17_group == "IL17_Hub")

# =========================================================
# 2) General inflammation score (Hallmark IFNG + TNFA/NFKB + IL6/JAK-STAT3, minus IL-17 module)
# =========================================================
hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
inflam_genes <- hm %>%
  filter(gs_name %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_IL6_JAK_STAT3_SIGNALING")) %>%
  pull(gene_symbol) %>% unique()
inflam_genes <- setdiff(intersect(inflam_genes, rownames(data_mat)), il17_module)
cat("Inflammation genes on panel:", length(inflam_genes), "\n")
meta$inflam_score <- Matrix::colMeans(data_mat[inflam_genes, meta$cell_id])

# =========================================================
# 3) Per-CM spatial covariates (computed within each section)
# =========================================================
xy <- function(d) as.matrix(d[, c("x", "y")])

cm_df <- bind_rows(lapply(unique(meta$sample_id), function(s) {
  m   <- meta %>% filter(sample_id == s)
  cm  <- m %>% filter(is_cm)
  hub <- m %>% filter(is_hub)
  imm <- m %>% filter(is_immune)
  nc  <- m %>% filter(!is_cm)
  if (nrow(cm) == 0 || nrow(hub) == 0) return(NULL)
  q  <- xy(cm)
  nb <- frNN(xy(imm), eps = r_um, query = q)$id
  cm$dist_hub       <- get.knnx(xy(hub), q, k = 1)$nn.dist[, 1]
  cm$immune_density <- lengths(nb)
  cm$n_hub_50       <- sapply(nb, function(j) sum(imm$is_hub[j]))
  nn <- get.knnx(xy(nc), q, k = k_nbr)$nn.index
  cm$inflam_nbr     <- rowMeans(matrix(nc$inflam_score[nn], ncol = k_nbr))
  cm
}))

cm_df <- cm_df %>%
  mutate(sample_id = factor(sample_id),
         z_dist    = as.numeric(scale(log1p(dist_hub))),
         z_dens    = as.numeric(scale(log1p(immune_density))),
         z_inflam  = as.numeric(scale(inflam_nbr)))

# =========================================================
# 4) Close-CM program: original DE method
# =========================================================
cm_obj <- CreateSeuratObject(counts = counts[, cm_df$cell_id], assay = "Xenium",
                             meta.data = cm_df %>% column_to_rownames("cell_id"))
LayerData(cm_obj, assay = "Xenium", layer = "data") <- data_mat[, cm_df$cell_id]
rm(counts, data_mat); gc()

cf_obj <- subset(cm_obj, subset = cm_IL17_group %in% c("Close_to_IL17_hub", "Far_from_IL17_hub"))
Idents(cf_obj) <- "cm_IL17_group"
de_orig <- FindMarkers(cf_obj, ident.1 = "Close_to_IL17_hub", ident.2 = "Far_from_IL17_hub",
                       test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1, assay = "Xenium")
de_orig$gene <- rownames(de_orig)
rm(cf_obj); gc()

genes_up_close <- de_orig %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)
cat("Close-CM program genes:", length(genes_up_close), "\n")

cm_obj <- AddModuleScore(cm_obj, features = list(genes_up_close), assay = "Xenium", name = "Prox", seed = 1)
cm_df$score_close <- cm_obj$Prox1[cm_df$cell_id]
rm(cm_obj); gc()

# =========================================================
# Helper: pooled test within a stratum (score ~ group + section) + Wilcoxon
# =========================================================
pooled_test <- function(d, grp_col, lvl_hi) {
  d <- d %>% mutate(g = factor(.data[[grp_col]] == lvl_hi, levels = c(FALSE, TRUE)))
  if (n_distinct(d$g) < 2) return(data.frame(beta = NA, se = NA, p_adj_section = NA, p_wilcox = NA))
  f  <- if (n_distinct(d$sample_id) > 1) score_close ~ g + sample_id else score_close ~ g
  co <- summary(lm(f, data = d))$coefficients
  data.frame(beta = co["gTRUE", 1], se = co["gTRUE", 2], p_adj_section = co["gTRUE", 4],
             p_wilcox = wilcox.test(score_close ~ g, data = d)$p.value)
}

# =========================================================
# A) Close vs Far within immune-density bins and inflammation tertiles
# =========================================================
cf2 <- cm_df %>%
  filter(cm_IL17_group %in% c("Close_to_IL17_hub", "Far_from_IL17_hub")) %>%
  mutate(group    = factor(ifelse(cm_IL17_group == "Close_to_IL17_hub", "Close", "Far"), levels = c("Close", "Far")),
         dens_bin = cut(immune_density, c(-Inf, 0, 5, 15, Inf), labels = c("0", "1-5", "6-15", ">15"))) %>%
  group_by(sample_id) %>%
  mutate(inflam_bin = cut(inflam_nbr, quantile(inflam_nbr, c(0, 1/3, 2/3, 1)),
                          labels = c("Low", "Mid", "High"), include.lowest = TRUE)) %>%
  ungroup()

stratA_fun <- function(d, bin_col, by_section = FALSE) {
  grp_vars <- if (by_section) c("sample_id", bin_col) else bin_col
  d %>%
    group_by(across(all_of(grp_vars))) %>%
    filter(sum(group == "Close") >= 20, sum(group == "Far") >= 20) %>%
    group_modify(~ bind_cols(
      data.frame(n_close = sum(.x$group == "Close"), n_far = sum(.x$group == "Far"),
                 median_close = median(.x$score_close[.x$group == "Close"]),
                 median_far   = median(.x$score_close[.x$group == "Far"])),
      pooled_test(.x %>% mutate(sample_id = if (by_section) factor(.y$sample_id) else sample_id),
                  "group", "Close"))) %>%
    ungroup() %>%
    rename(bin = all_of(bin_col)) %>%
    mutate(stratifier = bin_col)
}

stratA_main <- bind_rows(stratA_fun(cf2, "dens_bin"), stratA_fun(cf2, "inflam_bin"))
print(stratA_main)

# (pooled boxplot version of panel A not used; S5A is the line plot below)

# supplementary A: per section
pA_supp1 <- ggplot(cf2, aes(dens_bin, score_close, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = grp_cols) + facet_wrap(~ sample_id, scales = "free_y") +
  labs(x = "Immune cells within 50 µm", y = "Close-CM program score", fill = NULL,
       title = "Supplementary: Close vs Far within immune-density strata, per section") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
pA_supp2 <- ggplot(cf2, aes(inflam_bin, score_close, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = grp_cols) + facet_wrap(~ sample_id, scales = "free_y") +
  labs(x = "Neighborhood inflammation (within-section tertile)", y = "Close-CM program score", fill = NULL,
       title = "Supplementary: Close vs Far within inflammation strata, per section") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# (adjusted-model panel B not used in the figures; removed)


saveRDS(cf2,         file.path(out_dir, "cf2_close_far_cells.rds"))
saveRDS(stratA_main, file.path(out_dir, "stratA_main.rds"))

inflam_tbl <- hm %>%
  dplyr::filter(gs_name %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                               "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                               "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
                gene_symbol %in% inflam_genes) %>%
  dplyr::distinct(gs_name, gene_symbol)

length(inflam_genes)
table(inflam_tbl$gs_name)
write.csv(inflam_tbl, file.path(out_dir, "inflammation_score_genes.csv"), row.names = FALSE)




# ── Fig S5A: Close vs Far CM program score across matched strata ──
fmt_p2 <- function(p) ifelse(p < 0.001, "p<0.001", paste0("p=", formatC(p, format = "f", digits = 3)))

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9),
        legend.position = "bottom")

plotA_line <- function(bin_col, xlab, title) {
  d <- cf2 %>%
    dplyr::filter(!is.na(.data[[bin_col]])) %>%
    dplyr::group_by(bin = .data[[bin_col]], group) %>%
    dplyr::summarise(n = dplyr::n(), mean = mean(score_close),
                     se = sd(score_close) / sqrt(n), .groups = "drop")
  
  rng  <- range(c(d$mean - 1.96 * d$se, d$mean + 1.96 * d$se))
  ylab <- rng[2] + 0.12 * diff(rng)
  
  lab <- stratA_main %>%
    dplyr::filter(stratifier == bin_col) %>%
    dplyr::mutate(bin = factor(bin, levels = levels(d$bin)),
                  label = fmt_p2(p_adj_section))
  
  ggplot(d, aes(bin, mean, colour = group, group = group)) +
    geom_line(linewidth = 0.8) +
    geom_pointrange(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                    size = 0.45, linewidth = 0.7) +
    geom_text(data = lab, aes(x = bin, y = ylab, label = label),
              inherit.aes = FALSE, size = 3) +
    scale_colour_manual(values = grp_cols, name = NULL,
                        labels = c(Close = "Close to IL-17 hub", Far = "Far from IL-17 hub")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = xlab, y = "Close-CM program score\n(mean ± 95% CI)", title = title) +
    bt
}

pA1 <- plotA_line("dens_bin",   "Immune cells within 50 µm",           "Matched immune density")
pA2 <- plotA_line("inflam_bin", "Neighborhood inflammation (tertile)", "Matched inflammation")

pA_lines <- (pA1 | pA2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
pA_lines

ggsave(file.path(out_dir, "Main_A_stratified_lines.pdf"), pA_lines, width = 7.5, height = 3.8, device = cairo_pdf)

# ── Fig S5B: Close vs Far per section, within immune-density and inflammation strata ──
pdf(file.path(out_dir, "Supp_ABC_per_section.pdf"), width = 10, height = 4.5)
print(pA_supp1); print(pA_supp2)
dev.off()

saveRDS(meta,  file.path(out_dir, "meta_cells.rds"))
saveRDS(cm_df, file.path(out_dir, "cm_df_scores.rds"))

## end of Part 4a
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 4b - Fig 4I, 4J; Fig S5C-S5G: IL-17 hub-centred vs matched hub-free immune neighbourhoods
##  (Clean/reviewer_comments/spatial_adjustment.R, standalone section; reads meta_cells.rds
##   and cm_df_scores.rds written by Part 4a)
##############################################################################

## =========================================================
## Spatial: IL-17 hub-centred vs matched hub-free immune neighbourhoods
## Greedy matching (best match first; each control centre used once)
## Standalone
## =========================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(dbscan)
library(FNN)

setwd(file.path(TMOS_DIR, "Clean/reviewer_comments"))
out_dir <- "spatial_adjust"

## clear any results from earlier runs
suppressWarnings(rm(list = intersect(ls(), c("cand", "pairs_list", "res_list", "fill_max"))))

meta  <- readRDS(file.path(out_dir, "meta_cells.rds"))
cm_df <- readRDS(file.path(out_dir, "cm_df_scores.rds"))
stopifnot(all(c("inflam_score", "is_cm", "is_immune", "is_hub", "x", "y") %in% names(meta)))

## ---- settings ----
R_win      <- 100     # radius for matching counts / inflammation (µm)
R_ctrl     <- 150     # control centre: no hub within this distance
R_map      <- 150     # outermost ring / drawn window radius
min_cm     <- 15
z_tol      <- 0.1
z_tol_both <- 0.25
max_use    <- 1       # each control centre used at most once
ring_b     <- c(0, 25, 50, 100, 150)
ring_l     <- c("0–25", "25–50", "50–100", "100–150")
col_hub    <- "#b30000"
col_ctrl   <- "#4682B4"
win_cols   <- c("Hub-centred" = col_hub, "Control" = col_ctrl)

cm_scores <- cm_df %>% dplyr::select(cell_id, score_close)
cm_tab <- meta %>% dplyr::filter(is_cm) %>%
  dplyr::inner_join(cm_scores, by = "cell_id") %>%
  dplyr::select(sample_id, x, y, score_close)

## =========================================================
## 0) Candidate centres
## =========================================================
cand <- dplyr::bind_rows(lapply(unique(as.character(meta$sample_id)), function(s) {
  m   <- meta %>% dplyr::filter(sample_id == s)
  hub <- m %>% dplyr::filter(is_hub)
  imm <- m %>% dplyr::filter(is_immune)
  if (nrow(hub) < 5) return(NULL)
  
  cl <- dbscan::dbscan(as.matrix(hub[, c("x", "y")]), eps = 30, minPts = 5)$cluster
  hub_c <- hub %>% dplyr::mutate(cl = cl) %>% dplyr::filter(cl > 0) %>%
    dplyr::group_by(cl) %>% dplyr::summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>%
    dplyr::transmute(cx, cy, type = "Hub-centred")
  
  d_hub  <- FNN::get.knnx(as.matrix(hub[, c("x", "y")]), as.matrix(imm[, c("x", "y")]), k = 1)$nn.dist[, 1]
  ctrl_c <- imm[d_hub > R_ctrl, ] %>% dplyr::transmute(cx = x, cy = y, type = "Control")
  
  centres <- dplyr::bind_rows(hub_c, ctrl_c)
  nb <- dbscan::frNN(as.matrix(m[, c("x", "y")]), eps = R_win,
                     query = as.matrix(centres[, c("cx", "cy")]))$id
  centres %>% dplyr::mutate(
    sample_id = s,
    n_imm  = sapply(nb, function(j) sum(m$is_immune[j])),
    n_hub  = sapply(nb, function(j) sum(m$is_hub[j])),
    n_cm   = sapply(nb, function(j) sum(m$is_cm[j])),
    inflam = sapply(nb, function(j) { j <- j[!m$is_cm[j]]; if (length(j)) mean(m$inflam_score[j]) else NA_real_ }))
})) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(z_inflam = (inflam - mean(inflam, na.rm = TRUE)) / sd(inflam, na.rm = TRUE)) %>%
  dplyr::ungroup()

## =========================================================
## 1) Greedy pairing
## =========================================================
build_pairs <- function(match) {
  hubs_c <- cand %>% dplyr::filter(type == "Hub-centred", n_cm >= min_cm)
  ctrl_c <- cand %>% dplyr::filter(type == "Control", n_cm >= min_cm, n_hub == 0) %>%
    dplyr::mutate(ctrl_id = dplyr::row_number())
  p <- dplyr::inner_join(hubs_c, ctrl_c, by = "sample_id", suffix = c("_hub", "_ctrl"),
                         relationship = "many-to-many")
  imm_ok <- function(d) dplyr::filter(d, n_imm_ctrl / n_imm_hub >= 0.8, n_imm_ctrl / n_imm_hub <= 1.25)
  p <- switch(match,
              immune       = imm_ok(p),
              inflammation = dplyr::filter(p, !is.na(z_inflam_hub), !is.na(z_inflam_ctrl),
                                           abs(z_inflam_ctrl - z_inflam_hub) <= z_tol),
              both         = imm_ok(p) %>% dplyr::filter(!is.na(z_inflam_hub), !is.na(z_inflam_ctrl),
                                                         abs(z_inflam_ctrl - z_inflam_hub) <= z_tol_both))
  
  ## match quality: closer immune count + inflammation = better; then more CMs
  p <- p %>% dplyr::mutate(hub_id = paste(sample_id, round(cx_hub), round(cy_hub)),
                           dist_q = abs(log(n_imm_ctrl / n_imm_hub)) +
                             abs(dplyr::coalesce(z_inflam_ctrl - z_inflam_hub, 0)),
                           min_cm_pair = pmin(n_cm_hub, n_cm_ctrl)) %>%
    dplyr::arrange(dist_q, dplyr::desc(min_cm_pair))
  
  ## greedy matching: one control per hub, each control centre used at most max_use times
  used_ctrl <- integer(0); used_hub <- character(0); keep <- logical(nrow(p))
  for (i in seq_len(nrow(p))) {
    if (p$hub_id[i] %in% used_hub) next
    if (sum(used_ctrl == p$ctrl_id[i]) >= max_use) next
    keep[i] <- TRUE
    used_hub  <- c(used_hub, p$hub_id[i])
    used_ctrl <- c(used_ctrl, p$ctrl_id[i])
  }
  p[keep, ] %>%
    dplyr::arrange(dplyr::desc(min_cm_pair))
}

## =========================================================
## 2) Ring profiles + stats (paired Wilcoxon)
## =========================================================
ring_profile <- function(s, cx, cy) {
  cm_tab %>% dplyr::filter(sample_id == s) %>%
    dplyr::mutate(r = sqrt((x - cx)^2 + (y - cy)^2)) %>%
    dplyr::filter(r <= max(ring_b)) %>%
    dplyr::mutate(ring = cut(r, ring_b, labels = ring_l, include.lowest = TRUE)) %>%
    dplyr::group_by(ring) %>%
    dplyr::summarise(n_cm = dplyr::n(), score = mean(score_close), .groups = "drop") %>%
    tidyr::complete(ring, fill = list(n_cm = 0, score = NA))
}

run_profile <- function(pairs) {
  prof <- dplyr::bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i, ]
    dplyr::bind_rows(
      ring_profile(p$sample_id, p$cx_hub,  p$cy_hub)  %>% dplyr::mutate(type = "Hub-centred"),
      ring_profile(p$sample_id, p$cx_ctrl, p$cy_ctrl) %>% dplyr::mutate(type = "Control")) %>%
      dplyr::mutate(pair = i, sample_id = p$sample_id)
  })) %>% dplyr::mutate(type = factor(type, levels = c("Hub-centred", "Control")))
  
  stats <- prof %>% dplyr::select(pair, ring, type, score) %>%
    tidyr::pivot_wider(names_from = type, values_from = score) %>%
    dplyr::filter(!is.na(`Hub-centred`), !is.na(Control)) %>%
    dplyr::mutate(d = `Hub-centred` - Control) %>%
    dplyr::group_by(ring) %>%
    dplyr::group_modify(~ {
      p_w <- if (nrow(.x) >= 5) wilcox.test(.x$`Hub-centred`, .x$Control, paired = TRUE)$p.value else NA
      tibble::tibble(n_pairs = nrow(.x), diff = mean(.x$d), p = p_w)
    }) %>% dplyr::ungroup()
  
  sumry <- prof %>% dplyr::filter(!is.na(score)) %>%
    dplyr::group_by(type, ring) %>%
    dplyr::summarise(n = dplyr::n(), mean = mean(score), se = sd(score) / sqrt(n), .groups = "drop")
  
  list(prof = prof, stats = stats, sum = sumry)
}

## =========================================================
## 3) Shared styling
## =========================================================
base_txt <- 9
th_txt <- theme(plot.title = element_text(face = "bold", size = base_txt + 1),
                plot.subtitle = element_text(size = base_txt - 1, colour = "grey35"),
                legend.title = element_text(size = base_txt - 1),
                legend.text = element_text(size = base_txt - 1.5))

circ_df <- function(cx, cy, r, id, n = 120) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t), id = id)
}

## =========================================================
## 4) Panel A: section map (no linking lines)
## =========================================================
make_panelA <- function(pairs) {
  sec <- names(which.max(table(as.character(pairs$sample_id))))
  sp  <- pairs %>% dplyr::filter(sample_id == sec)
  cu  <- dplyr::distinct(sp, cx_ctrl, cy_ctrl)
  
  circ_all <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(seq_len(nrow(sp)), function(i)
      circ_df(sp$cx_hub[i], sp$cy_hub[i], R_map, paste0("h", i)))) %>%
      dplyr::mutate(type = "IL-17 hub-centred window"),
    dplyr::bind_rows(lapply(seq_len(nrow(cu)), function(i)
      circ_df(cu$cx_ctrl[i], cu$cy_ctrl[i], R_map, paste0("c", i)))) %>%
      dplyr::mutate(type = "Matched hub-free immune window")) %>%
    dplyr::mutate(type = factor(type, levels = c("IL-17 hub-centred window", "Matched hub-free immune window")))
  
  m_sec <- meta %>% dplyr::filter(sample_id == sec)
  set.seed(1); m_bg <- dplyr::slice_sample(m_sec, prop = 0.4)
  xr <- range(m_sec$x); yr <- range(m_sec$y)
  sb <- data.frame(x = xr[1], xend = xr[1] + 500, y = yr[1] - 0.05 * diff(yr))
  
  ggplot() +
    geom_point(data = m_bg, aes(x, y), colour = "grey88", size = 0.02, stroke = 0) +
    geom_point(data = dplyr::filter(m_sec, is_hub), aes(x, y), colour = "grey20", size = 0.05, stroke = 0) +
    geom_polygon(data = circ_all, aes(x, y, group = id, fill = type, colour = type, linetype = type),
                 alpha = 0.07, linewidth = 0.45) +
    geom_point(data = sp, aes(cx_hub, cy_hub), colour = col_hub, size = 0.6) +
    geom_point(data = cu, aes(cx_ctrl, cy_ctrl), colour = col_ctrl, size = 0.6) +
    geom_segment(data = sb, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.7) +
    geom_text(data = sb, aes(x = (x + xend) / 2, y = y - 0.04 * diff(yr)), label = "500 µm", size = 2.5) +
    scale_fill_manual(values = c(col_hub, col_ctrl), name = NULL) +
    scale_colour_manual(values = c(col_hub, col_ctrl), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    coord_fixed(xlim = xr + c(-0.01, 0.01) * diff(xr),
                ylim = c(yr[1] - 0.1 * diff(yr), yr[2] + 0.01 * diff(yr)), expand = FALSE) +
    labs(title = paste0("Matched windows, section ", sec),
         subtitle = "Dark points: IL-17 hub cells; circles: 150 µm windows") +
    theme_void(base_size = base_txt) + th_txt +
    theme(legend.position = "bottom", legend.key.size = unit(0.45, "cm"),
          legend.margin = margin(t = -4))
}

## =========================================================
## 5) Panel B: average bullseyes with outside callouts
## =========================================================
make_panelB <- function(res, fill_max) {
  bull_lv <- c("Around IL-17 hubs", "Around hub-free\nimmune clusters")
  r_draw  <- c(0, 1, 2, 3, 4)
  avg_ring <- res$sum %>%
    dplyr::transmute(window = factor(ifelse(type == "Hub-centred", bull_lv[1], bull_lv[2]), levels = bull_lv),
                     ring = factor(ring, levels = ring_l), score = mean)
  
  t_seq <- seq(0, 2 * pi, length.out = 240)
  discs <- dplyr::bind_rows(lapply(bull_lv, function(w)
    dplyr::bind_rows(lapply(rev(seq_along(ring_l)), function(i)
      data.frame(window = w, ring = ring_l[i], r_out = r_draw[i + 1],
                 x = r_draw[i + 1] * cos(t_seq), y = r_draw[i + 1] * sin(t_seq)))))) %>%
    dplyr::mutate(window = factor(window, levels = bull_lv), ring = factor(ring, levels = ring_l)) %>%
    dplyr::left_join(avg_ring, by = c("window", "ring")) %>%
    dplyr::mutate(grp = factor(paste(window, r_out), levels = unique(paste(window, r_out))))
  
  r_mid <- c(0.5, 1.5, 2.5, 3.5)
  val_lab <- avg_ring %>%
    dplyr::mutate(x = 0, y = -r_mid[as.integer(ring)] - ifelse(ring == "0–25", 0.05, 0),
                  lab = sprintf("%.2f", score))
  
  ang  <- c(62, 42, 22, 2) * pi / 180
  call <- data.frame(window = factor(bull_lv[2], levels = bull_lv),
                     x0 = r_mid * cos(ang), y0 = r_mid * sin(ang),
                     x1 = 4.35 * cos(ang) + 0.25, y1 = 4.35 * sin(ang),
                     lab = paste(ring_l, "µm"))
  
  icon_t <- seq(0, 2 * pi, length.out = 7)[-7]
  cpts <- dplyr::bind_rows(
    data.frame(window = bull_lv[1], x = c(0, 0.22 * cos(icon_t)), y = c(0.12, 0.12 + 0.22 * sin(icon_t)), kind = "hub"),
    data.frame(window = bull_lv[2], x = c(0, 0.22 * cos(icon_t)), y = c(0.12, 0.12 + 0.22 * sin(icon_t)), kind = "imm")) %>%
    dplyr::mutate(window = factor(window, levels = bull_lv))
  
  ggplot() +
    geom_polygon(data = discs, aes(x, y, group = grp, fill = score), colour = "white", linewidth = 0.7) +
    geom_point(data = dplyr::filter(cpts, kind == "hub"), aes(x, y),
               fill = "black", colour = "white", shape = 21, size = 1.3, stroke = 0.25) +
    geom_point(data = dplyr::filter(cpts, kind == "imm"), aes(x, y),
               fill = "grey50", colour = "white", shape = 21, size = 1.3, stroke = 0.25) +
    geom_text(data = val_lab, aes(x, y, label = lab), size = 2.7, fontface = "bold",
              colour = ifelse(val_lab$score > 0.35 * fill_max, "white", "grey15")) +
    geom_segment(data = call, aes(x = x0, y = y0, xend = x1 - 0.2, yend = y1),
                 colour = "grey45", linewidth = 0.25) +
    geom_point(data = call, aes(x0, y0), colour = "grey45", size = 0.5) +
    geom_text(data = call, aes(x1, y1, label = lab), size = 2.2, colour = "grey30", hjust = 0) +
    facet_wrap(~ window) +
    scale_fill_gradient(low = "#f7f7f7", high = col_hub, limits = c(0, fill_max),
                        oob = scales::squish, name = "Mean Close-CM\nprogram score",
                        guide = guide_colourbar(barheight = unit(2.5, "cm"), barwidth = unit(0.3, "cm"))) +
    coord_fixed(xlim = c(-4.2, 6.6), ylim = c(-4.2, 4.2), expand = FALSE, clip = "off") +
    labs(title = "Average CM program by distance from window centre",
         subtitle = paste0("Mean across all matched pairs (n = ", length(unique(res$prof$pair)),
                           "); ring widths not to scale")) +
    theme_void(base_size = base_txt) + th_txt +
    theme(strip.text = element_text(face = "bold", size = base_txt, margin = margin(b = 4)),
          legend.position = "right", panel.spacing = unit(-0.6, "cm"))
}

## =========================================================
## 6) Supplementary: paired line plot
## =========================================================
make_supp <- function(res, title_txt, balance_txt) {
  s  <- res$sum %>% dplyr::mutate(xi = as.integer(factor(ring, levels = ring_l)))
  rw <- s %>% dplyr::select(type, xi, mean) %>% tidyr::pivot_wider(names_from = type, values_from = mean)
  st <- res$stats %>% dplyr::mutate(xi = as.integer(factor(ring, levels = ring_l)))
  ytop <- max(s$mean + 1.96 * s$se, na.rm = TRUE)
  ybot <- min(0, min(s$mean - 1.96 * s$se, na.rm = TRUE))
  
  ggplot(s, aes(xi, mean, colour = type)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    geom_ribbon(data = rw, aes(x = xi, ymin = Control, ymax = `Hub-centred`),
                inherit.aes = FALSE, fill = col_hub, alpha = 0.08) +
    geom_line(linewidth = 0.9) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.08, linewidth = 0.6) +
    geom_point(size = 2.4) +
    geom_text(data = st, aes(x = xi, y = ytop + 0.12 * (ytop - ybot),
                             label = ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))),
              inherit.aes = FALSE, size = 2.7, colour = "grey20") +
    scale_x_continuous(breaks = 1:4, labels = ring_l, expand = expansion(add = 0.3)) +
    scale_y_continuous(limits = c(ybot, ytop + 0.18 * (ytop - ybot)), expand = expansion(mult = c(0.02, 0))) +
    scale_colour_manual(values = win_cols, name = NULL,
                        labels = c("Hub-centred" = "Around IL-17 hubs",
                                   "Control" = "Around matched hub-free immune clusters")) +
    labs(x = "Distance from window centre (µm)", y = "Close-CM program score\n(mean ± 95% CI)",
         title = title_txt,
         subtitle = paste0("n = ", length(unique(res$prof$pair)), " pairs; paired Wilcoxon test per ring\n", balance_txt)) +
    theme_classic(base_size = base_txt) + th_txt +
    theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
          axis.text = element_text(size = base_txt - 0.5, colour = "black"),
          legend.position = "bottom", legend.key.width = unit(0.6, "cm"))
}

## =========================================================
## 7) Run all three versions
## =========================================================
versions <- c(immune       = "immune density",
              inflammation = "neighbourhood inflammation",
              both         = "immune density and inflammation")

pairs_list <- lapply(names(versions), build_pairs); names(pairs_list) <- names(versions)
res_list   <- lapply(pairs_list, run_profile)
fill_max   <- max(sapply(res_list, function(r) max(r$sum$mean, na.rm = TRUE)))

# Spatial_main_both = Fig 4I-J; _immune = S5C; _inflammation = S5D
# Supp_Spatial_profile_both = S5E; _immune = S5F; _inflammation = S5G
for (v in names(versions)) {
  pr  <- pairs_list[[v]]
  bal <- sprintf("Balance (hub vs control): immune cells %.0f vs %.0f; inflammation z %.2f vs %.2f",
                 mean(pr$n_imm_hub), mean(pr$n_imm_ctrl),
                 mean(pr$z_inflam_hub, na.rm = TRUE), mean(pr$z_inflam_ctrl, na.rm = TRUE))
  cat("\n==", v, "| pairs:", nrow(pr),
      "| distinct control windows:", nrow(dplyr::distinct(pr, sample_id, cx_ctrl, cy_ctrl)),
      "==\n")
  print(table(pr$sample_id)); cat(bal, "\n"); print(res_list[[v]]$stats)
  
  ttl <- paste0("Matched for ", versions[[v]])
  
  fig_main <- (make_panelA(pr) / make_panelB(res_list[[v]], fill_max)) +
    plot_layout(heights = c(1, 0.9)) +
    plot_annotation(tag_levels = "A",
                    title = paste0("IL-17 hub-centred vs hub-free immune neighbourhoods: ", tolower(ttl)),
                    subtitle = bal,
                    theme = theme(plot.title = element_text(face = "bold", size = base_txt + 2),
                                  plot.subtitle = element_text(size = base_txt - 0.5, colour = "grey30"))) &
    theme(plot.tag = element_text(face = "bold", size = 12))
  print(fig_main)
  ggsave(file.path(out_dir, paste0("Spatial_main_", v, "_matched.pdf")), fig_main,
         width = 8.5, height = 8, device = cairo_pdf)
  
  fig_supp <- make_supp(res_list[[v]], ttl, bal)
  ggsave(file.path(out_dir, paste0("Supp_Spatial_profile_", v, "_matched.pdf")), fig_supp,
         width = 5, height = 4.6, device = cairo_pdf)
  
  write.csv(res_list[[v]]$stats, file.path(out_dir, paste0("Spatial_ring_stats_", v, "_matched.csv")), row.names = FALSE)
  write.csv(pr, file.path(out_dir, paste0("Spatial_pairs_", v, "_matched.csv")), row.names = FALSE)
}



## end of Part 4b
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 5 - Fig 4K: IL-17 signaling-related program in N4-like T/NK and myeloid cells, TMOS vs MC (heart scRNA-seq)
##  (Clean/reviewer_comments/R1.1_spatial_to_sc.R)
##############################################################################

# =========================================================
# IL-17 signaling-related program in N4-like cells: TMOS vs MC (scRNA heart)
# T/NK and myeloid cells analysed separately
# =========================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

setwd(file.path(TMOS_DIR, "Clean/reviewer_comments"))

# ---- settings ----
grp_cols    <- c("MC" = "#4682B4", "TMOS" = "#b30000")
cts         <- c("T/NK Lymphocytes", "Myeloid")
min_cells   <- 10    # min cells of a type per donor
min_n4like  <- 5     # min N4-like cells per donor
n4_top_frac <- 0.25  # top 25% of N4 score = N4-like

# ---- gene sets ----
module_genes <- c("NOD1","TAB1","IL17F","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","TAB2",
                  "FBXW11","IL17RA","IKBKB","IL17A","IRAK1","IRAK2","MEF2A","MEF2C","ATF1","NFKB1",
                  "PPP2CA","PPP2R1A","IL17RB","MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
                  "MAP2K1","NOD2","MAP2K4","IL25","SKP1","MAP3K7","TRAF6","CUL1","IKBKG","RIPK2",
                  "BTRC","RPS6KA5","MAPKAPK2")

n4_genes <- c("CTSL","MSR1","FCGR2A","CD163","SNX10","GZMB","THBS1","WARS","CXCL9","GBP1",
              "ANXA2","CXCL10","ANXA5","GSTO1","PDIA3","GBP5","HIF1A","PKM","LAMP1","PSMA6",
              "ITGB2","LAG3","ARF1","CCL7","CCL8","MAP3K7CL","SLAMF7","AQP9","CDKN1A","CXCL16",
              "PTPN1","CD300E","VAMP5","PLEKHO1","PLAUR","TYMP","ARL8B","CD47","RTN4","FPR2",
              "MAP2K1","FCGR1A","PTPN2","PLSCR1","ZYX","LGALS9","STAT1","SOCS3","GNS","TLR2",
              "SLC39A8","PGD","GLIPR2","CSNK2B","ADGRE5","PSTPIP2","HNRNPF","CYBA","ITGAX","DAZAP2",
              "GCH1","CD38","DEK","BCL2A1")

# =========================================================
# 1) load and score
# =========================================================
sc <- readRDS("heart_3232026.rds")
sc <- DietSeurat(sc, layers = c("counts", "data"))
md <- sc@meta.data

genes_found <- data.frame(
  gene_set = c("IL-17 signaling-related program (43 genes)", "N4 signature"),
  n_total  = c(length(module_genes), length(n4_genes)),
  n_found  = c(sum(module_genes %in% rownames(sc)), sum(n4_genes %in% rownames(sc))),
  missing  = c(paste(setdiff(module_genes, rownames(sc)), collapse = ", "),
               paste(setdiff(n4_genes,     rownames(sc)), collapse = ", "))
)
print(genes_found)

# abundance of each cell type as a fraction of all cardiac cells, per donor
burden <- md %>%
  group_by(sample_id) %>%
  summarise(`T/NK Lymphocytes` = mean(cell_type_unified == "T/NK Lymphocytes"),
            Myeloid            = mean(cell_type_unified == "Myeloid"), .groups = "drop") %>%
  tidyr::pivot_longer(-sample_id, names_to = "ct", values_to = "celltype_frac")

cell_scores <- bind_rows(lapply(cts, function(ct) {
  obj <- subset(sc, subset = cell_type_unified == ct & m_group %in% c("MC", "MC+MS+MG"))
  obj <- AddModuleScore(obj,
                        features = list(intersect(module_genes, rownames(obj)),
                                        intersect(n4_genes,     rownames(obj))),
                        name = "Score", seed = 1)
  data.frame(cell = colnames(obj), ct = ct,
             sample_id = obj$sample_id, m_group = obj$m_group,
             il17 = obj$Score1, n4 = obj$Score2)
}))
rm(sc); gc()

cell_scores <- cell_scores %>%
  group_by(ct) %>%
  mutate(n4_like = n4 >= quantile(n4, 1 - n4_top_frac)) %>%
  ungroup() %>%
  mutate(group = factor(ifelse(m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS")))

# =========================================================
# 2) donor-level values
# =========================================================
donor <- cell_scores %>%
  group_by(ct, sample_id, group) %>%
  summarise(n_cells        = n(),
            n_n4like       = sum(n4_like),
            il17_in_n4like = if (sum(n4_like) >= min_n4like) mean(il17[n4_like]) else NA_real_,
            .groups = "drop") %>%
  filter(n_cells >= min_cells, !is.na(il17_in_n4like)) %>%
  left_join(burden, by = c("sample_id", "ct")) %>%
  arrange(ct, group, sample_id)

# =========================================================
# 3) stats: exact Wilcoxon + linear model adjusted for the abundance of the
#    corresponding cell type (T/NK or myeloid fraction of all cells)
# =========================================================
stats <- donor %>%
  group_by(ct) %>%
  group_modify(~ {
    fit <- summary(lm(il17_in_n4like ~ group + celltype_frac, data = .x))$coefficients
    data.frame(n_MC                 = sum(.x$group == "MC"),
               n_TMOS               = sum(.x$group == "TMOS"),
               median_MC            = median(.x$il17_in_n4like[.x$group == "MC"]),
               median_TMOS          = median(.x$il17_in_n4like[.x$group == "TMOS"]),
               p_wilcox             = wilcox.test(il17_in_n4like ~ group, data = .x)$p.value,
               TMOS_effect_adj      = fit["groupTMOS",     "Estimate"],
               p_TMOS_adj           = fit["groupTMOS",     "Pr(>|t|)"],
               celltype_frac_effect = fit["celltype_frac", "Estimate"],
               p_celltype_frac      = fit["celltype_frac", "Pr(>|t|)"])
  }) %>% ungroup()
print(stats)

# =========================================================
# 4) Fig 4K: donor-level score in N4-like T/NK (left) and myeloid (right) cells
# =========================================================
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
ylab  <- "IL-17 signaling-related program\n(N4-like cells, donor mean)"

make_box <- function(cell_type) {
  d  <- filter(donor, ct == cell_type)
  st <- filter(stats, ct == cell_type)
  yr <- range(d$il17_in_n4like); ypad <- diff(yr) * 0.28
  ggplot(d, aes(group, il17_in_n4like, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.6) +
    geom_jitter(aes(color = group), width = 0.12, size = 2.2, show.legend = FALSE) +
    annotate("text", x = 0.55, y = yr[2] + ypad * 0.85, hjust = 0, size = 3,
             label = paste0("Wilcoxon ", fmt_p(st$p_wilcox), "\nAdjusted ", fmt_p(st$p_TMOS_adj))) +
    scale_fill_manual(values = grp_cols) + scale_color_manual(values = grp_cols) +
    scale_x_discrete(labels = c(MC = "MC\nalone", TMOS = "TMOS")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
    facet_wrap(~ ct) + labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95"))
}

ggsave("Fig4K_N4like_TNK_Myeloid.pdf",
       make_box("T/NK Lymphocytes") + make_box("Myeloid"),
       width = 6.4, height = 3.4)

# =========================================================
# 5) Excel
# =========================================================
settings <- data.frame(
  parameter = c("Comparison", "Cell types", "Unit of analysis", "Scoring",
                "N4-like definition", "Min cells per donor (cell type)", "Min N4-like cells per donor",
                "Test 1", "Test 2"),
  value = c("TMOS (MC+MS+MG) vs MC alone", paste(cts, collapse = "; "), "Donor (mean of N4-like cells)",
            "Seurat AddModuleScore (seed = 1), within each cell type",
            paste0("Top ", n4_top_frac * 100, "% of 64-gene N4 signature score within each cell type (MC + TMOS cells)"),
            min_cells, min_n4like,
            "Wilcoxon rank-sum (two-sided, exact)",
            "Linear model: score ~ group + fraction of the corresponding cell type among all cells")
)

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#EDEDED", border = "bottom")
add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
  freezePane(wb, name, firstRow = TRUE)
}
add_sheet("Stats",        stats)
add_sheet("Donor_values", donor)
add_sheet("Settings",     settings)
add_sheet("Genes_found",  genes_found)
add_sheet("Gene_sets",    data.frame(IL17_program = c(module_genes, rep(NA, length(n4_genes) - length(module_genes))),
                                     N4_signature = n4_genes))
saveWorkbook(wb, "Fig4K_N4like_TNK_Myeloid_TMOS_vs_MC.xlsx", overwrite = TRUE)


## end of Part 5
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 6 - Fig 4L; Fig S4E, S4G: N4 and IL-17 signatures in longitudinal PBMC
##  (Clean/fig4/niche_to_pbmc/fig4_spatial_to_pbmc_il17.R)
##############################################################################


# ── 0. SET WORKING DIRECTORY & LOAD LIBRARIES ──
setwd(file.path(TMOS_DIR, "Clean/fig4/niche_to_pbmc"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reactome.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ── 1. LOAD OBJECTS ──
pbmc     <- readRDS("pbmc_642025.rds")
# keep only the default assay (reductions/graphs are not used) to save memory
pbmc     <- DietSeurat(pbmc, assays = DefaultAssay(pbmc), dimreducs = NULL, graphs = NULL)
# (the Xenium heart object is loaded further down, where it is first used)

# ── 2. PREP PBMC METADATA ──
pbmc@meta.data <- pbmc@meta.data %>%
  mutate(
    deg_group = case_when(
      m_group == "MC"                          ~ "MC",
      m_group %in% c("MC+MS", "MC+MS+MG")     ~ "MCO",
      TRUE                                      ~ NA_character_
    )
  )

pbmc_map <- c(
  "B and plasma" = "B Cells",
  "CD4"          = "T/NK Lymphocytes",
  "CD8 and NK"   = "T/NK Lymphocytes",
  "MNP"          = "Myeloid",
  "cDCs"         = "Myeloid",
  "pDCs"         = "Myeloid"
)

pbmc@meta.data$cell_type_unified <- pbmc_map[pbmc@meta.data$lineage]

# ── 3. SUBSET PBMC TO MYELOID + TNK ──
pbmc_immune <- subset(pbmc,
                      subset = cell_type_unified %in% c("Myeloid", "T/NK Lymphocytes"))
cat("PBMC immune cells:", ncol(pbmc_immune), "\n")
print(table(pbmc_immune@meta.data$cell_type_unified))
rm(pbmc); invisible(gc())

# ── 4. PULL IL-17 GENE LIST FROM REACTOME ──
# ── 4. PULL IL-17 GENE LIST FROM REACTOME ──
reac_to_entrez <- as.list(reactomePATHID2EXTID)

entrez_il17 <- reac_to_entrez[["R-HSA-448424"]]
sym_il17    <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = entrez_il17,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)$SYMBOL

# ALL genes from Reactome — no panel filtering
il17_full <- unique(na.omit(sym_il17))

cat("IL-17 total genes from Reactome:", length(il17_full), "\n")
print(il17_full)

# ── 5. CHECK GENES IN PBMC ──
genes_il17_pbmc <- il17_full[il17_full %in% rownames(pbmc_immune)]
cat("IL-17 genes found in PBMC:  ", length(genes_il17_pbmc), "/", length(il17_full), "\n")
cat("Missing from PBMC:          ", paste(il17_full[!il17_full %in% rownames(pbmc_immune)], collapse = ", "), "\n")
# ── 6. SCORE PBMC ──
pbmc_immune <- AddModuleScore(
  pbmc_immune,
  features = list(genes_il17_pbmc),
  name     = "IL17_sig",
  seed     = 42
)
pbmc_immune$IL17_sig  <- pbmc_immune$IL17_sig1
pbmc_immune$IL17_sig1 <- NULL

cat("Scoring done\n")
summary(pbmc_immune$IL17_sig)

# ── 7. DEFINE TIMEPOINTS & DONORS ──
required_tp <- c("pre_ici", "on_ici", "pre_steroid", "post_steroid")

donors_with_pre_ici <- pbmc_immune@meta.data %>%
  filter(timepoint_cat == "pre_ici") %>%
  pull(donor) %>%
  unique()

cat("Donors with pre_ici timepoint:", length(donors_with_pre_ici), "\n")

# (separate IL-17 plots, post-steroid comparison, fold change and LME plots not used; the final combined figure below (4L, S4E, S4G) recomputes what it needs)
# ── N4 COMBINED TNK & MYELOID SIGNATURE ──

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(openxlsx)
library(tidyr)

# ── 1. SUBSET HEART XO TO TNK AND MYELOID ──
heart_xo <- fix_fov(readRDS("xo_heart_processed.rds"))
tnk_myeloid_xo <- subset(heart_xo,
                         subset = major_celltype_broad %in% c("T_NK_lymphocytes", "Myeloid"))

cat("T/NK + Myeloid cells in Xenium:\n")
print(table(tnk_myeloid_xo@meta.data$niches,
            tnk_myeloid_xo@meta.data$major_celltype_broad))

# ── 2. DE: N4 vs ALL OTHER NICHES — TNK ──
tnk_xo <- subset(tnk_myeloid_xo,
                 subset = major_celltype_broad == "T_NK_lymphocytes")
Idents(tnk_xo)       <- "niches"
DefaultAssay(tnk_xo)  <- "Xenium"

cat("\nFinding T/NK N4 markers...\n")
tnk_n4_markers <- FindMarkers(
  tnk_xo,
  ident.1         = "N4",
  ident.2         = NULL,
  test.use        = "wilcox",
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  only.pos        = TRUE
)
tnk_n4_markers$gene <- rownames(tnk_n4_markers)

tnk_n4_sig <- tnk_n4_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC))

cat("T/NK N4 signature genes:", nrow(tnk_n4_sig), "\n")
print(tnk_n4_sig[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

# ── 3. DE: N4 vs ALL OTHER NICHES — MYELOID ──
myeloid_xo <- subset(tnk_myeloid_xo,
                     subset = major_celltype_broad == "Myeloid")
Idents(myeloid_xo)       <- "niches"
DefaultAssay(myeloid_xo)  <- "Xenium"

cat("\nFinding Myeloid N4 markers...\n")
myeloid_n4_markers <- FindMarkers(
  myeloid_xo,
  ident.1         = "N4",
  ident.2         = NULL,
  test.use        = "wilcox",
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  only.pos        = TRUE
)
myeloid_n4_markers$gene <- rownames(myeloid_n4_markers)

myeloid_n4_sig <- myeloid_n4_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC))

cat("Myeloid N4 signature genes:", nrow(myeloid_n4_sig), "\n")
print(myeloid_n4_sig[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

# ── 4. COMBINE INTO ONE SIGNATURE ──
# Union of both gene lists — deduplicated
combined_sig_genes <- union(tnk_n4_sig$gene, myeloid_n4_sig$gene)

cat("\nCombined N4 signature (T/NK + Myeloid union):", length(combined_sig_genes), "genes\n")

# Tag each gene with its origin
gene_origin <- data.frame(
  gene   = combined_sig_genes,
  origin = case_when(
    combined_sig_genes %in% tnk_n4_sig$gene     & 
      combined_sig_genes %in% myeloid_n4_sig$gene ~ "Shared",
    combined_sig_genes %in% tnk_n4_sig$gene       ~ "T/NK only",
    combined_sig_genes %in% myeloid_n4_sig$gene   ~ "Myeloid only"
  )
)
print(table(gene_origin$origin))

# ── 5. CHECK GENE COVERAGE IN PBMC ──
# pbmc_immune = combined Myeloid + T/NK from PBMC
genes_combined_pbmc <- combined_sig_genes[combined_sig_genes %in% rownames(pbmc_immune)]

cat("\nCombined genes found in PBMC: ", length(genes_combined_pbmc),
    "/", length(combined_sig_genes), "\n")
cat("Missing from PBMC: ",
    paste(combined_sig_genes[!combined_sig_genes %in% rownames(pbmc_immune)],
          collapse = ", "), "\n")

# ── 6. SCORE PBMC (combined Myeloid + TNK) ──
pbmc_immune <- AddModuleScore(
  pbmc_immune,
  features = list(genes_combined_pbmc),
  name     = "N4_combined_sig",
  seed     = 42
)
pbmc_immune$N4_combined_sig  <- pbmc_immune$N4_combined_sig1
pbmc_immune$N4_combined_sig1 <- NULL

cat("Scoring done\n")
print(summary(pbmc_immune$N4_combined_sig))

# (separate N4 plots not used; see combined figure below)
# ── 10. SAVE EXCEL ──
wb <- createWorkbook()

addWorksheet(wb, "Combined_Sig_Genes")
writeData(wb, "Combined_Sig_Genes", gene_origin)
setColWidths(wb, "Combined_Sig_Genes", cols = 1:2, widths = "auto")

addWorksheet(wb, "TNK_N4_Markers")
writeData(wb, "TNK_N4_Markers",
          tnk_n4_sig %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val_adj))
setColWidths(wb, "TNK_N4_Markers", cols = 1:5, widths = "auto")

addWorksheet(wb, "Myeloid_N4_Markers")
writeData(wb, "Myeloid_N4_Markers",
          myeloid_n4_sig %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val_adj))
setColWidths(wb, "Myeloid_N4_Markers", cols = 1:5, widths = "auto")

saveWorkbook(wb, "N4_Combined_Signature.xlsx", overwrite = TRUE)
cat("✅ Excel saved.\n")
# (N4 LME plot and variance diagnostics not used; see combined figure below)


# ══════════════════════════════════════════════════════════════════
# COMBINED FIGURE — TRAJECTORY (MC | MCO | LME) + LOG2FC
# IL-17 (row 1) and N4 (row 2)
# Shared y-axis per row, strip labels on trajectory panels
# ══════════════════════════════════════════════════════════════════

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(tidyr)
library(lme4)
library(lmerTest)

col_vals_mco <- c("MC" = "#4682B4", "MCO" = "#b30000")
shape_vals   <- c("MC" = 16, "MC+MS" = 17, "MC+MS+MG" = 15)
required_tp  <- c("pre_ici", "on_ici", "pre_steroid", "post_steroid")

pbmc_meta <- pbmc_immune@meta.data %>%
  filter(timepoint_cat %in% required_tp,
         donor %in% donors_with_pre_ici,
         m_group %in% c("MC", "MC+MS", "MC+MS+MG"))

# ══════════════════════════════════════════════════════════════════
# COMPUTE SHARED Y LIMITS PER SIGNATURE
# ══════════════════════════════════════════════════════════════════

get_ylim <- function(score_col, lme_model) {
  
  raw <- pbmc_meta %>%
    mutate(
      group_plot = case_when(
        m_group == "MC"                      ~ "MC",
        m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
        TRUE                                  ~ NA_character_
      )
    ) %>%
    filter(!is.na(group_plot)) %>%
    group_by(donor, timepoint_cat, group_plot) %>%
    summarise(mean_score = mean(.data[[score_col]]), .groups = "drop") %>%
    pull(mean_score)
  
  pred <- predict(lme_model,
                  newdata = expand.grid(
                    tp_num     = 1:4,
                    group_plot = factor(c("MC", "MCO"), levels = c("MC", "MCO"))
                  ),
                  re.form = NA)
  
  all_vals <- c(raw, pred)
  pad      <- diff(range(all_vals)) * 0.12
  c(min(all_vals) - pad, max(all_vals) + pad)
}

# ══════════════════════════════════════════════════════════════════
# PANEL — TRAJECTORY PER GROUP (strip label on top)
# ══════════════════════════════════════════════════════════════════

make_trajectory_group <- function(score_col, group_filter,
                                  strip_label, ylim,
                                  show_y_axis  = TRUE,
                                  show_legend  = FALSE) {
  
  df <- pbmc_meta %>%
    mutate(
      group_plot = case_when(
        m_group == "MC"                      ~ "MC",
        m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
        TRUE                                  ~ NA_character_
      ),
      mco_label = case_when(
        m_group == "MC"       ~ "MC",
        m_group == "MC+MS"    ~ "MC+MS",
        m_group == "MC+MS+MG" ~ "MC+MS+MG",
        TRUE                   ~ NA_character_
      )
    ) %>%
    filter(!is.na(group_plot), group_plot == group_filter) %>%
    group_by(donor, timepoint_cat, group_plot, mco_label) %>%
    summarise(mean_score = mean(.data[[score_col]]), .groups = "drop") %>%
    mutate(timepoint_cat = factor(timepoint_cat, levels = required_tp))
  
  col_use  <- col_vals_mco[group_filter]
  leg_pos  <- if (show_legend) "right" else "none"
  
  # dummy facet column for strip label
  df$panel_label <- strip_label
  
  p <- ggplot(df, aes(x     = timepoint_cat,
                      y     = mean_score,
                      group = donor,
                      color = group_plot)) +
    geom_line(linewidth = 0.8, alpha = 0.6) +
    geom_point(aes(shape = mco_label), size = 3.5, alpha = 0.9) +
    facet_wrap(~ panel_label) +
    coord_cartesian(ylim = ylim) +
    scale_color_manual(values = col_vals_mco) +
    scale_shape_manual(values = shape_vals) +
    labs(x     = NULL,
         y     = "Mean Module Score",
         color = "Group",
         shape = "Subgroup") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = leg_pos,
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line        = element_blank(),
      strip.background = element_rect(fill = col_use),
      strip.text       = element_text(color = "white", face = "bold", size = 11)
    )
  
  if (!show_y_axis) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  p
}

# ══════════════════════════════════════════════════════════════════
# PANEL — LME FITTED LINES (strip label on top, shared y)
# ══════════════════════════════════════════════════════════════════

make_lme_plot <- function(lme_model, ylim, strip_label) {
  
  coef_table      <- as.data.frame(coef(summary(lme_model)))
  coef_table$term <- rownames(coef_table)
  int_row         <- coef_table[grepl("tp_num:group_plotMCO", coef_table$term), ]
  estimate        <- int_row[["Estimate"]]
  two_sided_p     <- int_row[["Pr(>|t|)"]]
  ci              <- confint(lme_model, parm = "tp_num:group_plotMCO", method = "Wald")
  
  pred_grid <- expand.grid(
    tp_num     = 1:4,
    group_plot = factor(c("MC", "MCO"), levels = c("MC", "MCO"))
  ) %>%
    mutate(
      timepoint_cat = factor(
        c("pre_ici", "on_ici", "pre_steroid", "post_steroid")[tp_num],
        levels = required_tp
      )
    )
  pred_grid$fit    <- predict(lme_model, newdata = pred_grid, re.form = NA)
  pred_grid$panel_label <- strip_label
  
  # annotation y position near top of shared ylim
  y_annot <- ylim[2] * 0.97
  
  ggplot() +
    geom_line(data  = pred_grid,
              aes(x = timepoint_cat, y = fit,
                  group = group_plot, color = group_plot),
              linewidth = 1.8) +
    geom_point(data = pred_grid,
               aes(x = timepoint_cat, y = fit, color = group_plot),
               size = 4, shape = 18) +
    facet_wrap(~ panel_label) +
    coord_cartesian(ylim = ylim) +
    annotate("text",
             x     = 1, y = y_annot,
             label = paste0("β = ",       round(estimate,    4),
                            "  95% CI [", round(ci[1],       4),
                            ", ",         round(ci[2],       4), "]",
                            "  p = ",     round(two_sided_p, 3)),
             size  = 3, fontface = "bold", color = "black", hjust = 0) +
    scale_color_manual(values = col_vals_mco, name = "Group") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.y     = element_blank(),
      legend.position  = "right",
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line        = element_blank(),
      strip.background = element_rect(fill = "grey30"),
      strip.text       = element_text(color = "white", face = "bold", size = 11)
    )
}

# ══════════════════════════════════════════════════════════════════
# PANEL — LOG2 FOLD CHANGE
# ══════════════════════════════════════════════════════════════════

make_log2fc <- function(score_col, title_label) {
  
  fc_df <- pbmc_meta %>%
    mutate(
      group_plot = case_when(
        m_group == "MC"                      ~ "MC",
        m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
        TRUE                                  ~ NA_character_
      )
    ) %>%
    filter(!is.na(group_plot),
           timepoint_cat %in% c("pre_ici", "post_steroid")) %>%
    group_by(donor, group_plot, timepoint_cat) %>%
    summarise(mean_score = mean(.data[[score_col]]), .groups = "drop") %>%
    pivot_wider(names_from  = timepoint_cat,
                values_from = mean_score) %>%
    filter(!is.na(pre_ici), !is.na(post_steroid)) %>%
    mutate(
      log2_fc    = log2(post_steroid / pre_ici),
      group_plot = factor(group_plot, levels = c("MC", "MCO"))
    )
  
  stat_fc <- fc_df %>%
    wilcox_test(log2_fc ~ group_plot) %>%
    add_significance() %>%
    add_xy_position(x = "group_plot")
  
  ggplot(fc_df, aes(x = group_plot, y = log2_fc, fill = group_plot)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = group_plot), width = 0.1, size = 2.5) +
    stat_pvalue_manual(stat_fc,
                       label      = "p = {round(p, 3)}",
                       tip.length = 0.01,
                       size       = 3.5) +
    scale_fill_manual(values  = col_vals_mco) +
    scale_color_manual(values = col_vals_mco) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(title    = title_label,
         subtitle = paste0("MC n=", sum(fc_df$group_plot == "MC"),
                           " | MCO n=", sum(fc_df$group_plot == "MCO")),
         x        = NULL,
         y        = "Log2 FC\n(post_steroid / pre_ici)") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 11, face = "bold"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line       = element_blank()
    )
}

# ══════════════════════════════════════════════════════════════════
# FIT LME MODELS
# ══════════════════════════════════════════════════════════════════

prep_lme <- function(score_col) {
  pbmc_meta %>%
    mutate(
      group_plot = case_when(
        m_group == "MC"                      ~ "MC",
        m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
        TRUE                                  ~ NA_character_
      )
    ) %>%
    filter(!is.na(group_plot)) %>%
    group_by(donor, group_plot, timepoint_cat) %>%
    summarise(score = mean(.data[[score_col]]), .groups = "drop") %>%
    mutate(
      tp_num        = case_when(
        timepoint_cat == "pre_ici"      ~ 1,
        timepoint_cat == "on_ici"       ~ 2,
        timepoint_cat == "pre_steroid"  ~ 3,
        timepoint_cat == "post_steroid" ~ 4
      ),
      timepoint_cat = factor(timepoint_cat, levels = required_tp),
      group_plot    = factor(group_plot, levels = c("MC", "MCO"))
    )
}

lme_data_il17 <- prep_lme("IL17_sig")
lme_data_n4   <- prep_lme("N4_combined_sig")

lme_il17 <- lmer(score ~ tp_num * group_plot + (1 | donor),
                 data = lme_data_il17, REML = FALSE)
lme_n4   <- lmer(score ~ tp_num * group_plot + (1 | donor),
                 data = lme_data_n4,   REML = FALSE)

# ══════════════════════════════════════════════════════════════════
# SHARED Y LIMITS
# ══════════════════════════════════════════════════════════════════

ylim_il17 <- get_ylim("IL17_sig",        lme_il17)
ylim_n4   <- get_ylim("N4_combined_sig", lme_n4)

# ══════════════════════════════════════════════════════════════════
# BUILD ALL PANELS
# ══════════════════════════════════════════════════════════════════

# ── IL-17 ──
p_il17_mc  <- make_trajectory_group("IL17_sig", "MC",
                                    "MC — IL-17 Trajectory",
                                    ylim_il17,
                                    show_y_axis = TRUE,
                                    show_legend = FALSE)
p_il17_mco <- make_trajectory_group("IL17_sig", "MCO",
                                    "MCO — IL-17 Trajectory",
                                    ylim_il17,
                                    show_y_axis = FALSE,
                                    show_legend = TRUE)
p_il17_lme <- make_lme_plot(lme_il17, ylim_il17,
                            "LME — IL-17 Fitted Trajectories")
p_il17_fc  <- make_log2fc("IL17_sig",
                          "IL-17\npre_ici → post_steroid")

# ── N4 ──
p_n4_mc    <- make_trajectory_group("N4_combined_sig", "MC",
                                    "MC — N4 Trajectory",
                                    ylim_n4,
                                    show_y_axis = TRUE,
                                    show_legend = FALSE)
p_n4_mco   <- make_trajectory_group("N4_combined_sig", "MCO",
                                    "MCO — N4 Trajectory",
                                    ylim_n4,
                                    show_y_axis = FALSE,
                                    show_legend = TRUE)
p_n4_lme   <- make_lme_plot(lme_n4, ylim_n4,
                            "LME — N4 Fitted Trajectories")
p_n4_fc    <- make_log2fc("N4_combined_sig",
                          "N4 Combined\npre_ici → post_steroid")

# ══════════════════════════════════════════════════════════════════
# ASSEMBLE
# ══════════════════════════════════════════════════════════════════

p_combined <-
  (p_il17_mc | p_il17_mco | p_il17_lme | p_il17_fc) /
  (p_n4_mc   | p_n4_mco   | p_n4_lme   | p_n4_fc) +
  plot_layout(widths = c(2, 2, 2.5, 1)) +
  plot_annotation(
    title    = "IL-17 & N4 Immune Signatures — PBMC Longitudinal Analysis",
    subtitle = "MC vs MCO  |  Myeloid + T/NK lymphocytes",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold",  hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave("IL17_N4_Combined_Figure.pdf",
       p_combined, width = 16, height = 10)

cat("✅ Saved: IL17_N4_Combined_Figure.pdf\n")


## end of Part 6
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 7 - Fig S4C: pathogenic vs homeostatic Th17 scores, N4 vs other niches
##  (Clean/fig4/th17_homeo_patho.R)
##############################################################################


# ─── 0. Setup ─────────────────────────────────────────────────────────────────
setwd(file.path(TMOS_DIR, "Clean/fig4"))

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(ggpubr)

tnk_myeloid_heart_st <- fix_fov(readRDS("tnk_myeloid_heart_st.rds"))

# ─── 1. Subset ────────────────────────────────────────────────────────────────
th17_obj <- subset(
  tnk_myeloid_heart_st,
  subset = SingleR_label %in% c("Th17 cells", "Th1/Th17 cells") &
    niches %in% c("N1", "N2", "N3", "N4")
)

th17_obj$niche_group <- factor(
  ifelse(th17_obj$niches == "N4", "N4 (IE)", "Other (N1/N2/N3)"),
  levels = c("Other (N1/N2/N3)", "N4 (IE)")
)

th17_obj$SingleR_label <- factor(th17_obj$SingleR_label,
                                 levels = c("Th17 cells", "Th1/Th17 cells"))

# ─── 2. Genes ─────────────────────────────────────────────────────────────────
scrna_genes       <- rownames(th17_obj)
pathogenic_genes  <- c("CSF2", "IFNG", "TBX21", "IL23R", "GZMB")
homeostatic_genes <- c("IL10", "AHR", "MAF")
path_genes        <- intersect(pathogenic_genes,  scrna_genes)
homeo_genes       <- intersect(homeostatic_genes, scrna_genes)

# ─── 3. Score ─────────────────────────────────────────────────────────────────
th17_obj <- AddModuleScore(th17_obj, features = list(path_genes),
                           name = "Pathogenic_Th17", ctrl = 50, seed = 42)
th17_obj <- AddModuleScore(th17_obj, features = list(homeo_genes),
                           name = "Homeostatic_Th17", ctrl = 50, seed = 42)
th17_obj$Path_Homeo_Balance <- th17_obj$Pathogenic_Th171 - th17_obj$Homeostatic_Th171

# ─── 4. Metadata ──────────────────────────────────────────────────────────────
meta <- th17_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, SingleR_label, niche_group,
                Pathogenic_Th171, Homeostatic_Th171,
                Path_Homeo_Balance, IL17_Full)

# ─── 5. Summary ───────────────────────────────────────────────────────────────
bar_df <- meta %>%
  group_by(SingleR_label, niche_group) %>%
  summarise(
    N                = n(),
    Mean_Pathogenic  = mean(Pathogenic_Th171,   na.rm = TRUE),
    SE_Pathogenic    = sd(Pathogenic_Th171,     na.rm = TRUE) / sqrt(n()),
    Mean_Homeostatic = mean(Homeostatic_Th171,  na.rm = TRUE),
    SE_Homeostatic   = sd(Homeostatic_Th171,    na.rm = TRUE) / sqrt(n()),
    Mean_Balance     = mean(Path_Homeo_Balance, na.rm = TRUE),
    SE_Balance       = sd(Path_Homeo_Balance,   na.rm = TRUE) / sqrt(n()),
    Mean_IL17        = mean(IL17_Full,           na.rm = TRUE),
    .groups          = "drop"
  )

# ─── 6. P-values from per-cell Wilcoxon ───────────────────────────────────────
add_pval <- function(df, score_col) {
  df %>%
    group_by(SingleR_label) %>%
    summarise(
      p = wilcox.test(
        .data[[score_col]][niche_group == "N4 (IE)"],
        .data[[score_col]][niche_group == "Other (N1/N2/N3)"]
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        p < 0.001 ~ "p<0.001",
        p < 0.01  ~ "p<0.01",
        p < 0.05  ~ paste0("p=", round(p, 3)),
        TRUE      ~ paste0("p=", round(p, 3))
      )
    )
}

pval_path  <- add_pval(meta, "Pathogenic_Th171")
pval_homeo <- add_pval(meta, "Homeostatic_Th171")
pval_bal   <- add_pval(meta, "Path_Homeo_Balance")

cat("P-values Pathogenic:\n");  print(pval_path)
cat("P-values Homeostatic:\n"); print(pval_homeo)
cat("P-values Balance:\n");     print(pval_bal)

# ─── 7. Y positions for brackets ─────────────────────────────────────────────
y_path  <- max(bar_df$Mean_Pathogenic  + bar_df$SE_Pathogenic,  na.rm = TRUE) * 1.25
y_homeo <- max(abs(bar_df$Mean_Homeostatic) + bar_df$SE_Homeostatic, na.rm = TRUE) * 1.25
y_bal   <- max(abs(bar_df$Mean_Balance) + bar_df$SE_Balance,    na.rm = TRUE) * 1.25

# ─── 8. Shared aesthetics ─────────────────────────────────────────────────────
il17_min   <- min(bar_df$Mean_IL17, na.rm = TRUE)
il17_max   <- max(bar_df$Mean_IL17, na.rm = TRUE)

il17_scale <- scale_fill_gradient(
  low    = "#f7f7f7",
  high   = "#008080",
  limits = c(il17_min, il17_max),
  name   = "Mean\nIL17 score"
)

bar_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey95"),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    legend.position    = "right"
  )

# ─── 9. Panel A: Pathogenic ───────────────────────────────────────────────────
p_path <- ggplot(bar_df,
                 aes(x = niche_group, y = Mean_Pathogenic,
                     fill = Mean_IL17)) +
  geom_col(width = 0.55, alpha = 0.9, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean_Pathogenic - SE_Pathogenic,
                    ymax = Mean_Pathogenic + SE_Pathogenic),
                width = 0.15, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.35) +
  geom_segment(data = pval_path,
               aes(x = 1, xend = 2,
                   y = y_path * 0.88, yend = y_path * 0.88),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = pval_path,
            aes(x = 1.5, y = y_path, label = p_label),
            inherit.aes = FALSE, size = 3.2, fontface = "italic") +
  facet_wrap(~ SingleR_label, nrow = 1) +
  il17_scale +
  labs(title = "Pathogenic Th17 score — N4 vs other niches",
       y = "Mean module score ± SE", x = NULL) +
  bar_theme

# ─── 10. Panel B: Homeostatic ─────────────────────────────────────────────────
p_homeo <- ggplot(bar_df,
                  aes(x = niche_group, y = Mean_Homeostatic,
                      fill = Mean_IL17)) +
  geom_col(width = 0.55, alpha = 0.9, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean_Homeostatic - SE_Homeostatic,
                    ymax = Mean_Homeostatic + SE_Homeostatic),
                width = 0.15, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.35) +
  geom_segment(data = pval_homeo,
               aes(x = 1, xend = 2,
                   y = y_homeo * 0.88, yend = y_homeo * 0.88),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = pval_homeo,
            aes(x = 1.5, y = y_homeo, label = p_label),
            inherit.aes = FALSE, size = 3.2, fontface = "italic") +
  facet_wrap(~ SingleR_label, nrow = 1) +
  il17_scale +
  labs(title = "Homeostatic Th17 score — N4 vs other niches",
       y = "Mean module score ± SE", x = NULL) +
  bar_theme

# ─── 11. Panel C: Balance ─────────────────────────────────────────────────────
p_balance <- ggplot(bar_df,
                    aes(x = niche_group, y = Mean_Balance,
                        fill = Mean_IL17)) +
  geom_col(width = 0.55, alpha = 0.9, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean_Balance - SE_Balance,
                    ymax = Mean_Balance + SE_Balance),
                width = 0.15, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.35) +
  geom_segment(data = pval_bal,
               aes(x = 1, xend = 2,
                   y = y_bal * 0.88, yend = y_bal * 0.88),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = pval_bal,
            aes(x = 1.5, y = y_bal, label = p_label),
            inherit.aes = FALSE, size = 3.2, fontface = "italic") +
  facet_wrap(~ SingleR_label, nrow = 1) +
  il17_scale +
  labs(title = "Pathogenic − Homeostatic balance — N4 vs other niches",
       y = "Mean score difference ± SE", x = NULL) +
  bar_theme

# ─── 12. Save ─────────────────────────────────────────────────────────────────
pdf("PathogenicTh17_N4vsOther_BarPlots.pdf", width = 5, height = 8)
print(p_path / p_homeo / p_balance)
dev.off()
cat("✅ Saved: PathogenicTh17_N4vsOther_BarPlots.pdf\n")

# ─── 13. Excel ────────────────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", bar_df)
addWorksheet(wb, "Pvalues")
writeData(wb, "Pvalues", bind_rows(
  pval_path  %>% mutate(Score = "Pathogenic"),
  pval_homeo %>% mutate(Score = "Homeostatic"),
  pval_bal   %>% mutate(Score = "Balance")
))
addWorksheet(wb, "PerCell_Scores")
writeData(wb, "PerCell_Scores", meta %>% dplyr::select(-cell_id))
saveWorkbook(wb, "PathogenicTh17_N4vsOther_Results.xlsx", overwrite = TRUE)
cat("✅ Saved: PathogenicTh17_N4vsOther_Results.xlsx\n")

cat("\n=== Summary ===\n")
print(bar_df)


## end of Part 7
rm(list = setdiff(ls(), c("TMOS_DIR", "fix_fov"))); invisible(gc())


##############################################################################
##  PART 8 - Fig S4F: IL-17 receptor expression in CMs by proximity tier
##  (Clean/fig4/il17_receptor_ligand.R)
##############################################################################


# ─── 0. Setup ─────────────────────────────────────────────────────────────────
setwd(file.path(TMOS_DIR, "Clean/fig4"))

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(openxlsx)

xo_heart_cleaned <- fix_fov(readRDS("xo_heart_processed_withspatialil17labels.rds"))

# ─── 1. Subset CMs ────────────────────────────────────────────────────────────
cm_obj <- subset(
  xo_heart_cleaned,
  subset = updated_cell_labels == "Cardiomyocyte" &
    cm_IL17_group %in% c("Close_to_IL17_hub",
                         "Intermediate",
                         "Far_from_IL17_hub")
)

cm_obj$cm_IL17_group <- factor(
  cm_obj$cm_IL17_group,
  levels = c("Close_to_IL17_hub", "Intermediate", "Far_from_IL17_hub")
)

cat("CM cells per proximity tier:\n")
print(table(cm_obj$cm_IL17_group))

# ─── 2. Receptor genes ────────────────────────────────────────────────────────
receptor_genes <- c("IL17RA", "TRAF3IP2")

# ─── 3. Extract expression ────────────────────────────────────────────────────
expr_mat <- GetAssayData(cm_obj, assay = "Xenium", layer = "data")
expr_sub  <- as.data.frame(
  t(as.matrix(expr_mat[receptor_genes, , drop = FALSE]))
)
expr_sub$cell_id       <- rownames(expr_sub)
expr_sub$cm_IL17_group <- cm_obj$cm_IL17_group

# ─── 4. Summary stats ─────────────────────────────────────────────────────────
summary_df <- expr_sub %>%
  tidyr::pivot_longer(cols      = all_of(receptor_genes),
                      names_to  = "gene",
                      values_to = "expression") %>%
  group_by(gene, cm_IL17_group) %>%
  summarise(
    N          = n(),
    Mean_exp   = round(mean(expression,       na.rm = TRUE), 4),
    Median_exp = round(median(expression,     na.rm = TRUE), 4),
    Pct_pos    = round(mean(expression > 0) * 100, 2),
    .groups    = "drop"
  )

cat("\nSummary:\n")
print(summary_df)

# ─── 5. Shared aesthetics ─────────────────────────────────────────────────────
tier_cols <- c(
  "Close_to_IL17_hub" = "#b30000",
  "Intermediate"      = "#d9d9d9",
  "Far_from_IL17_hub" = "#4682B4"
)

tier_labels <- c(
  "Close_to_IL17_hub" = "Close\n(≤Q25)",
  "Intermediate"      = "Intermediate",
  "Far_from_IL17_hub" = "Far\n(≥Q75)"
)

plot_df <- expr_sub %>%
  tidyr::pivot_longer(cols      = all_of(receptor_genes),
                      names_to  = "gene",
                      values_to = "expression") %>%
  mutate(
    cm_IL17_group = factor(cm_IL17_group,
                           levels = c("Close_to_IL17_hub",
                                      "Intermediate",
                                      "Far_from_IL17_hub")),
    gene = factor(gene, levels = c("IL17RA", "TRAF3IP2"))
  )

pct_df <- summary_df %>%
  mutate(
    cm_IL17_group = factor(cm_IL17_group,
                           levels = c("Close_to_IL17_hub",
                                      "Intermediate",
                                      "Far_from_IL17_hub")),
    gene = factor(gene, levels = c("IL17RA", "TRAF3IP2"))
  )

# ─── Combined: % expressed + mean expression colored dots ────────────────────
p_combined <- ggplot(pct_df,
                     aes(x = cm_IL17_group, y = Pct_pos,
                         color = Mean_exp, group = gene)) +
  geom_line(linewidth = 0.8, color = "grey70") +
  geom_point(size = 6) +
  facet_wrap(~ gene, nrow = 1) +
  scale_color_gradient(low = "lightgrey", high = "#b30000",
                       name = "Mean\nexpression") +
  scale_x_discrete(labels = tier_labels) +
  scale_y_continuous(limits = c(0, NA),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "IL-17 receptor competence by proximity tier",
    subtitle = "Dot size fixed; color = mean expression; y = % cells expressing",
    y        = "% cells expressing",
    x        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(face = "bold.italic", size = 11),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(color = "grey50", size = 9)
  )

ggsave("IL17R_CM_ProximityTier_Combined.pdf",
       plot = p_combined, width = 5, height = 4, device = cairo_pdf)
cat("✅ Saved: IL17R_CM_ProximityTier_Combined.pdf\n")

# ─── 10. Save Excel ───────────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Summary_by_Tier")
writeData(wb, "Summary_by_Tier", summary_df)
saveWorkbook(wb, "IL17R_CM_ProximityTier_Results.xlsx", overwrite = TRUE)
cat("✅ Saved: IL17R_CM_ProximityTier_Results.xlsx\n")

cat("\n=== Summary ===\n")
print(summary_df)



## end of Part 8
