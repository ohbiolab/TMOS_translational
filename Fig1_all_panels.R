# ============================================================================
# Fig1_all_panels.R
# Figure 1 and Supplementary Figure S1 - all panels in one script.
#
#   Panel(s)            Part   Main inputs
#   1A                  -      schematic (BioRender), no code
#   1B                  1      Clean/fig1/TMOS_trans_outcome.xlsx
#   1C, 1D, 1E          2      Clean/fig1/prev_objects/heart_minor.rds,
#   S1B-S1F             2        skeletal_322026.rds, pbmc_2192026.rds
#   1F, S1G, S1H        3      Clean/reviewer_comments/R2.1_TNK_states/
#                                (state tables + paired CellChat objects made by
#                                 fig2/Fig2_all_panels.R (Part 3))
#   S1A (Sankey)        -      cohort Sankey from sample metadata; not included
# ============================================================================
if (!exists("TMOS_DIR")) source(file.path("config", "paths.R"))



##############################################################################
##  PART 1 - Fig 1B: fatal myocarditis, MC alone vs MCO (Cohort 1)
##############################################################################

# ============================================================================
# Part 1
# Fatal myocarditis rate, MC alone vs MCO (Cohort 1, Fisher's exact test)
#
# Figure panels : Fig 1B
# ============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)

# ── Load data ─────────────────────────────────────────────────────────────────
mgh   <- read_excel(file.path(TMOS_DIR, "Clean/fig1/TMOS_trans_outcome.xlsx"),
                    sheet = "MGH")   # adjust sheet name if needed

mdacc <- read_excel(file.path(TMOS_DIR, "Clean/fig1/TMOS_trans_outcome.xlsx"),
                    sheet = "MDACC")

# ── Recode MGH ────────────────────────────────────────────────────────────────
mgh <- mgh %>%
  mutate(
    Group     = if_else(M_group == "MC", "MC alone", "MC+MS /\nMC+MS+MG"),
    Fatal_bin = if_else(`Fatal due to myocarditis` == "Yes", 1, 0)
  )

mgh_summary <- mgh %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    events = sum(Fatal_bin),
    rate   = events / n * 100,
    se     = sqrt((rate/100 * (1 - rate/100)) / n) * 100,
    lower  = pmax(rate - 1.96 * se, 0),
    upper  = pmin(rate + 1.96 * se, 100),
    .groups = "drop"
  )

# Fisher's exact test for MGH
ct_mgh  <- table(mgh$Group, mgh$Fatal_bin)
p_val   <- fisher.test(ct_mgh)$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001",
                  ifelse(p_val < 0.05,  sprintf("p = %.3f", p_val),
                         sprintf("p = %.3f", p_val)))

# ── Recode MDACC ──────────────────────────────────────────────────────────────
mdacc <- mdacc %>%
  mutate(
    Group      = if_else(m_group == "MC", "MC alone", "MC+MS /\nMC+MS+MG"),
    MACE_bin   = if_else(MACE == "Yes", 1, 0)
  )

mdacc_summary <- mdacc %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    events = sum(MACE_bin),
    rate   = events / n * 100,
    se     = sqrt((rate/100 * (1 - rate/100)) / n) * 100,
    lower  = pmax(rate - 1.96 * se, 0),
    upper  = pmin(rate + 1.96 * se, 100),
    .groups = "drop"
  )

# ── MGH plot ──────────────────────────────────────────────────────────────────
y_bracket <- max(mgh_summary$upper) + 8

p1 <- ggplot(mgh_summary, aes(x = Group, y = rate, fill = Group)) +
  geom_col(width = 0.5, color = NA) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.7) +
  geom_text(aes(label = paste0(events, "/", n)),
            vjust = -1.8, size = 4.2, fontface = "bold") +
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) +
  annotate("segment", x = 1, xend = 1, y = y_bracket - 1.5, yend = y_bracket) +
  annotate("segment", x = 2, xend = 2, y = y_bracket - 1.5, yend = y_bracket) +
  annotate("text", x = 1.5, y = y_bracket + 3,
           label = p_label, size = 4, fontface = "italic") +
  scale_fill_manual(values = c("MC alone" = "#4682B4", "MC+MS /\nMC+MS+MG" = "#b30000")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title   = "MGH Cohort",
    subtitle = "Outcome: Fatal Myocarditis",
    x       = NULL,
    y       = "Event Rate (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.x      = element_text(face = "bold", size = 11)
  )

# ── MDACC plot ────────────────────────────────────────────────────────────────
p2 <- ggplot(mdacc_summary, aes(x = Group, y = rate, fill = Group)) +
  geom_col(width = 0.5, color = NA) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.7,
                data = filter(mdacc_summary, Group != "MC alone")) +
  geom_text(aes(label = paste0(events, "/", n)),
            vjust = -1.5, size = 4.2, fontface = "bold") +
  annotate("text", x = 1.5, y = max(mdacc_summary$rate) + 15,
           label = "No statistical comparison\n(MC alone: n = 1)",
           size = 3.5, fontface = "italic", color = "gray40", hjust = 0.5) +
  scale_fill_manual(values = c("MC alone" = "#4682B4", "MC+MS /\nMC+MS+MG" = "#b30000")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title    = "MDACC Cohort",
    subtitle = "Outcome: MACE",
    x        = NULL,
    y        = "Event Rate (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.x      = element_text(face = "bold", size = 11)
  )

# ── Combine with patchwork ────────────────────────────────────────────────────
combined <- p1 + p2 +
  plot_annotation(
    title   = "Fatality Rate and MACE by Myocarditis Group",
    caption = "Error bars represent 95% CI",
    theme   = theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.caption = element_text(color = "gray40")
    )
  )

combined

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(file.path(TMOS_DIR, "Clean/fig1/fatality_MACE_comparison.pdf"),
       plot = combined, width = 9, height = 6, dpi = 300)



ggsave(
  file.path(TMOS_DIR, "Clean/fig1/fatality_MGH.pdf"),
  plot   = p1,
  width  = 6,
  height = 5,
  device = cairo_pdf
)
cat("✅ Saved: fatality_MGH.pdf\n")


##############################################################################
##  PART 2 - Fig 1C, 1D, 1E; Fig S1B-S1F: cross-tissue scRNA-seq, paired heart vs skeletal DE and GSEA
##############################################################################

# ============================================================================
# Part 2
# TMOS cross-tissue scRNA-seq: UMAPs, marker dot plots, paired heart vs skeletal DE, Hallmark/Reactome GSEA, BRAF/NGF
#
# Figure panels : Fig 1C, 1D, 1E; Fig S1B, S1C, S1D, S1E, S1F
# ============================================================================
library(tidyr)   # used below but not loaded in the original script

setwd(file.path(TMOS_DIR, "Clean/fig1"))

heart <- readRDS("prev_objects/heart_minor.rds")
skeletal <- readRDS("prev_objects/skeletal_322026.rds")
pbmc <- readRDS("prev_objects/pbmc_2192026.rds")


######################################## Main: UMAPS #######################################################

# ── Unified cell type mapping ──────────────────────────────────────────────────

# HEART  (source column: major_cell_type)
heart_map <- c(
  "T_NK_lymphocytes"     = "T/NK Lymphocytes",
  "Myeloid_cells"        = "Myeloid",
  "B_cells"              = "B Cells",
  "Fibroblast_stromal"   = "Fibroblast/Stromal",
  "Smooth_muscle_mural"  = "Smooth Muscle",
  "Endothelial_cells"    = "Endothelial",
  "Cardiomyocytes"       = "Cardiomyocytes",      # heart-specific
  "Erythroid_RBC"        = "RBC"
)

heart@meta.data$cell_type_unified <-
  heart_map[heart@meta.data$major_cell_type]


skeletal_map <- c(
  "T/NK Lymphocytes"    = "T/NK Lymphocytes",
  "Myeloid"             = "Myeloid",
  "B cells"             = "B Cells",
  "Fibroblast/Stromal"  = "Fibroblast/Stromal",
  "Smooth Muscle"       = "Smooth Muscle",
  "Endothelial"         = "Endothelial",
  "RBC"                 = "RBC",
  "Skeletal Muscle"     = "Skeletal Muscle",
  "Myogenic Progenitor" = "Myogenic Progenitor",
  "Cycling"             = "T/NK Lymphocytes"      # ← merged with T/NK
)

skeletal@meta.data$cell_type_unified <-
  skeletal_map[skeletal@meta.data$major_match]

# Verify
table(skeletal@meta.data$cell_type_unified)



# PBMC  (source column: lineage)
pbmc_map <- c(
  "B and plasma" = "B Cells",
  "CD4"          = "T/NK Lymphocytes",
  "CD8 and NK"   = "T/NK Lymphocytes",
  "MNP"          = "Myeloid",           # Mono-Nuclear Phagocytes → Myeloid
  "cDCs"         = "Myeloid",           # dendritic cells → Myeloid
  "pDCs"         = "Myeloid"            # plasmacytoid DCs → Myeloid
)

pbmc@meta.data$cell_type_unified <-
  pbmc_map[pbmc@meta.data$lineage]


# ── Quick sanity checks ────────────────────────────────────────────────────────
table(heart@meta.data$cell_type_unified)
table(skeletal@meta.data$cell_type_unified)
table(pbmc@meta.data$cell_type_unified)

# ── Remove RBC from all objects ────────────────────────────────────────────────

heart    <- subset(heart,    cell_type_unified != "RBC")
skeletal <- subset(skeletal, cell_type_unified != "RBC")

# PBMC has no RBC in lineage but just in case:
pbmc     <- subset(pbmc,     cell_type_unified != "RBC")

# ── Verify ─────────────────────────────────────────────────────────────────────
table(heart@meta.data$cell_type_unified)
table(skeletal@meta.data$cell_type_unified)
table(pbmc@meta.data$cell_type_unified)

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(heart,    "heart_3232026.rds")
saveRDS(skeletal, "skeletal_3232026.rds")
saveRDS(pbmc,     "pbmc_3232026.rds")




library(Seurat)
library(ggplot2)

unified_colors <- c(
  "Endothelial"         = "#1b9e77",
  "Fibroblast/Stromal"  = "#7570b3",
  "Smooth Muscle"       = "#66a61e",
  "Myogenic Progenitor" = "#b2df8a",
  "T/NK Lymphocytes"    = "#1f78b4",
  "B Cells"             = "#d95f02",
  "Myeloid"             = "#e6ab02",
  "Skeletal Muscle"     = "#A50021",
  "Cardiomyocytes"      = "#d73027"   # ← red, added
)




# DimPlots
DimPlot(heart,    group.by = "cell_type_unified", cols = unified_colors)
DimPlot(skeletal, group.by = "cell_type_unified", cols = unified_colors)
DimPlot(pbmc,     group.by = "cell_type_unified", cols = unified_colors)


# ── Subset MC+MS+MG from each object ──────────────────────────────────────────

heart_tmos    <- subset(heart,    m_group == "MC+MS+MG")
skeletal_tmos <- subset(skeletal, m_group == "MC+MS+MG")
pbmc_tmos     <- subset(pbmc,     m_group == "MC+MS+MG")

# ── Verify ─────────────────────────────────────────────────────────────────────
table(heart_tmos@meta.data$m_group)
table(skeletal_tmos@meta.data$m_group)
table(pbmc_tmos@meta.data$m_group)

# expected counts from your original tables:
# heart_tmos    → 10,742 cells
# skeletal_tmos → 10,997 cells
# pbmc_tmos     → 18,703 cells

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(heart_tmos,    "heart_tmos_3232026.rds")
saveRDS(skeletal_tmos, "skeletal_tmos_3232026.rds")
saveRDS(pbmc_tmos,     "pbmc_tmos_3232026.rds")



DimPlot(heart_tmos,    group.by = "cell_type_unified", cols = unified_colors)
DimPlot(skeletal_tmos, group.by = "cell_type_unified", cols = unified_colors)
DimPlot(pbmc_tmos,     group.by = "cell_type_unified", cols = unified_colors)






library(Seurat)
library(harmony)
library(dplyr)


# ── Step 1: Merge all 3 objects ───────────────────────────────────────────────
combined_tmos <- merge(
  heart_tmos,
  y = list(skeletal_tmos, pbmc_tmos),
  add.cell.ids = c("Heart", "Skeletal", "PBMC"),
  project = "TMOS"
)

# ── Step 2: Standard preprocessing ───────────────────────────────────────────
combined_tmos <- NormalizeData(combined_tmos)
combined_tmos <- FindVariableFeatures(combined_tmos, nfeatures = 3000)
combined_tmos <- ScaleData(combined_tmos)
combined_tmos <- RunPCA(combined_tmos, npcs = 30)


# ── Assign tissue from cell barcode prefixes ──────────────────────────────────
combined_tmos$tissue <- case_when(
  grepl("^Heart_",    rownames(combined_tmos@meta.data)) ~ "Heart",
  grepl("^Skeletal_", rownames(combined_tmos@meta.data)) ~ "Skeletal",
  grepl("^PBMC_",     rownames(combined_tmos@meta.data)) ~ "PBMC",
  TRUE ~ NA_character_
)

# Verify
table(combined_tmos$tissue, useNA = "always")

# Verify — should give Heart=10736, Skeletal=10997, PBMC=18703
table(combined_tmos$tissue, useNA = "always")



# ── Re-run Harmony with tissue ────────────────────────────────────────────────
combined_tmos <- RunHarmony(combined_tmos, group.by.vars = "tissue")
combined_tmos <- RunUMAP(combined_tmos, reduction = "harmony", dims = 1:30)

# ── Plot ──────────────────────────────────────────────────────────────────────
p1 <- DimPlot(combined_tmos,
              group.by = "cell_type_unified",
              cols     = unified_colors,
              pt.size  = 0.3) +
  ggtitle("MC+MS+MG — all tissues")

p2 <- DimPlot(combined_tmos,
              group.by = "cell_type_unified",
              split.by = "tissue",
              cols     = unified_colors,
              pt.size  = 0.3) +
  theme(legend.position = "none")

library(patchwork)
p1 / p2 + plot_layout(heights = c(1, 1))



# ── Extract axis limits from combined UMAP ────────────────────────────────────
umap_coords <- Embeddings(combined_tmos, reduction = "umap")
x_limits    <- range(umap_coords[, 1]) + c(-0.5, 0.5)  # small buffer
y_limits    <- range(umap_coords[, 2]) + c(-0.5, 0.5)

# ── Shared theme with fixed axis limits ───────────────────────────────────────
umap_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title       = element_text(hjust = 0.5),
    axis.text        = element_text(color = "black"),
    legend.position  = "none"
  )

# ── Plot 1: Combined ──────────────────────────────────────────────────────────
p_merged <- DimPlot(
  combined_tmos,
  reduction  = "umap",
  group.by   = "cell_type_unified",
  cols       = unified_colors,
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("Combined") +
  xlim(x_limits) + ylim(y_limits) +
  umap_theme

# ── Plot 2: Heart ─────────────────────────────────────────────────────────────
p_heart <- DimPlot(
  subset(combined_tmos, tissue == "Heart"),
  reduction  = "umap",
  group.by   = "cell_type_unified",
  cols       = unified_colors,
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("Heart") +
  xlim(x_limits) + ylim(y_limits) +
  umap_theme

# ── Plot 3: PBMC ──────────────────────────────────────────────────────────────
p_pbmc <- DimPlot(
  subset(combined_tmos, tissue == "PBMC"),
  reduction  = "umap",
  group.by   = "cell_type_unified",
  cols       = unified_colors,
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("PBMC") +
  xlim(x_limits) + ylim(y_limits) +
  umap_theme

# ── Plot 4: Skeletal ──────────────────────────────────────────────────────────
p_skeletal <- DimPlot(
  subset(combined_tmos, tissue == "Skeletal"),
  reduction  = "umap",
  group.by   = "cell_type_unified",
  cols       = unified_colors,
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("Skeletal") +
  xlim(x_limits) + ylim(y_limits) +
  umap_theme

# ── Extract single shared legend ──────────────────────────────────────────────
shared_legend <- cowplot::get_legend(
  DimPlot(combined_tmos,
          group.by = "cell_type_unified",
          cols     = unified_colors) +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 12),
      legend.title    = element_blank(),
      legend.key.size = unit(0.5, "cm")
    )
)

# ── Arrange & save ────────────────────────────────────────────────────────────
library(cowplot)
library(patchwork)

four_panels <- (p_merged | p_heart) / (p_pbmc | p_skeletal)

final_plot <- plot_grid(
  four_panels,
  shared_legend,
  ncol       = 2,
  rel_widths = c(4, 0.6)
)

print(final_plot)

ggsave(
  filename    = "TMOS_4panel_UMAP.pdf",
  plot        = final_plot,
  device      = "pdf",
  path        = file.path(TMOS_DIR, "Clean/fig1"),
  width       = 14,
  height      = 12,
  useDingbats = FALSE
)




######################################## Main: Heart Dot plot #######################################################


library(Seurat)
library(ggplot2)

# =========================================================
# 1. Canonical markers (Top 3 per cell type)
# =========================================================
marker_list <- list(
  T_NK_Lymphocytes    = c("CD3D", "CD3E", "NKG7"),
  Myeloid             = c("LST1", "TYROBP", "FCGR3A"),
  B_Cells             = c("MS4A1", "CD79A", "IGHM"),
  Fibroblast_Stromal  = c("COL1A1", "DCN", "LUM"),
  Smooth_Muscle       = c("ACTA2", "TAGLN", "MYH11"),
  Endothelial         = c("CDH5", "VWF", "CLDN5"),
  Cardiomyocytes      = c("ACTC1", "TNNT2", "TNNI3")
)

# =========================================================
# 2. Flatten markers and keep only present genes
# =========================================================
markers_ordered <- unlist(marker_list, use.names = FALSE)
markers_ordered <- markers_ordered[markers_ordered %in% rownames(heart)]

# =========================================================
# 3. Order cell types (Y-axis) — rev() puts first at top
# =========================================================
celltype_order <- c(
  "T/NK Lymphocytes",
  "Myeloid",
  "B Cells",
  "Fibroblast/Stromal",
  "Smooth Muscle",
  "Endothelial",
  "Cardiomyocytes"
)

heart$cell_type_unified <- factor(
  heart$cell_type_unified,
  levels = rev(celltype_order)
)

Idents(heart) <- "cell_type_unified"

# =========================================================
# 4. DotPlot
# =========================================================
byr_cols <- c("#4FA3D9", "#FFF2A6", "#F15A4A")

dotplot_heart <- DotPlot(
  heart,
  features = markers_ordered,
  assay    = "RNA",
  dot.scale = 6
) +
  scale_color_gradientn(colors = byr_cols) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    axis.line       = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Average expression",
    size  = "% cells expressing"
  )

dotplot_heart

# =========================================================
# 5. Save
# =========================================================
ggsave(
  filename    = "DotPlot_Heart_CellType_BYR.pdf",
  plot        = dotplot_heart,
  device      = cairo_pdf,
  path        = file.path(TMOS_DIR, "Clean/fig1"),
  width       = 7.5,
  height      = 3
)




######################################## Suppl: Dot plots #######################################################


# =========================================================
# SKELETAL DOT PLOT — full skeletal object
# =========================================================

marker_list_skeletal <- list(
  `T/NK Lymphocytes`    = c("CD3D", "CD3E", "NKG7"),
  Myeloid               = c("LST1", "TYROBP", "FCGR3A"),
  `B Cells`             = c("MS4A1", "CD79A", "IGHM"),
  `Fibroblast/Stromal`  = c("COL1A1", "DCN", "LUM"),
  `Smooth Muscle`       = c("ACTA2", "TAGLN", "MYH11"),
  Endothelial           = c("CDH5", "VWF", "CLDN5"),
  `Skeletal Muscle`     = c("TTN", "MB", "TNNT3"),
  `Myogenic Progenitor` = c("PAX7", "MYF5", "MYOD1"),
  Cycling               = c("MKI67", "TOP2A", "PCNA")
)

# note: RBC excluded since we removed them from the object already

markers_skeletal <- unlist(marker_list_skeletal, use.names = FALSE)
markers_skeletal <- markers_skeletal[markers_skeletal %in% rownames(skeletal)]

celltype_order_skeletal <- c(
  "T/NK Lymphocytes",
  "Myeloid",
  "B Cells",
  "Fibroblast/Stromal",
  "Smooth Muscle",
  "Endothelial",
  "Skeletal Muscle",
  "Myogenic Progenitor",
  "Cycling"
)

skeletal$cell_type_unified <- factor(
  skeletal$cell_type_unified,
  levels = rev(celltype_order_skeletal)
)

Idents(skeletal) <- "cell_type_unified"

byr_cols <- c("#4FA3D9", "#FFF2A6", "#F15A4A")

dotplot_skeletal <- DotPlot(
  skeletal,
  features  = markers_skeletal,
  assay     = "RNA",
  dot.scale = 6
) +
  scale_color_gradientn(colors = byr_cols) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    axis.line       = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Average expression",
    size  = "% cells expressing"
  )

dotplot_skeletal

ggsave(
  filename = "DotPlot_Skeletal_CellType_BYR.pdf",
  plot     = dotplot_skeletal,
  device   = cairo_pdf,
  path     = file.path(TMOS_DIR, "Clean/fig1"),
  width    = 10,
  height   = 3
)


# =========================================================
# PBMC DOT PLOT — full pbmc object
# =========================================================

marker_list_pbmc <- list(
  `T/NK Lymphocytes` = c("CD3D", "CD3E", "NKG7"),
  Myeloid            = c("LYZ", "S100A8", "CD14"),
  `B Cells`          = c("MS4A1", "CD79A", "IGHM")
)

markers_pbmc <- unlist(marker_list_pbmc, use.names = FALSE)
markers_pbmc <- markers_pbmc[markers_pbmc %in% rownames(pbmc)]

celltype_order_pbmc <- c(
  "T/NK Lymphocytes",
  "Myeloid",
  "B Cells"
)

pbmc$cell_type_unified <- factor(
  pbmc$cell_type_unified,
  levels = rev(celltype_order_pbmc)
)

Idents(pbmc) <- "cell_type_unified"

dotplot_pbmc <- DotPlot(
  pbmc,
  features  = markers_pbmc,
  assay     = "RNA",
  dot.scale = 6
) +
  scale_color_gradientn(colors = byr_cols) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    axis.line       = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Average expression",
    size  = "% cells expressing"
  )

dotplot_pbmc

ggsave(
  filename = "DotPlot_PBMC_CellType_BYR.pdf",
  plot     = dotplot_pbmc,
  device   = cairo_pdf,
  path     = file.path(TMOS_DIR, "Clean/fig1"),
  width    = 6,
  height   = 2
)




########################################### Main: Paired Analysis ########################################################

# Heart vs Skeletal — only truly paired MDACC samples (4352 and 4420)
heart_paired_hs <- subset(heart_tmos,
                          sample_id %in% c("4352_myocarditis_S12_L001",
                                           "4420_myocarditis_S3_L001"))

skeletal_paired_hs <- subset(skeletal_tmos,
                             sample_id %in% c("4352_myositis_S11_L001",
                                              "4420_myositis_S4_L001"))

# Heart vs PBMC — MGH paired samples (SIC_197 and SIC_237)
heart_paired_hp <- subset(heart_tmos,
                          sample_id %in% c("SIC_197", "SIC_237"))

pbmc_paired_hp <- subset(pbmc_tmos,
                         donor %in% c("SIC_197", "SIC_237"))

# ── Verify ────────────────────────────────────────────────────────────────────
table(heart_paired_hs@meta.data$sample_id)
table(skeletal_paired_hs@meta.data$sample_id)
table(heart_paired_hp@meta.data$sample_id)
table(pbmc_paired_hp@meta.data$donor)





library(DESeq2)
library(dplyr)

# ── Step 1: Merge paired heart and skeletal into one object ───────────────────
heart_paired_hs$tissue    <- "Heart"
skeletal_paired_hs$tissue <- "Skeletal"

# add a clean patient ID column to both
heart_paired_hs$patient <- ifelse(
  grepl("4352", heart_paired_hs@meta.data$sample_id), "4352", "4420"
)
skeletal_paired_hs$patient <- ifelse(
  grepl("4352", skeletal_paired_hs@meta.data$sample_id), "4352", "4420"
)

hs_paired <- merge(
  heart_paired_hs,
  skeletal_paired_hs,
  add.cell.ids = c("Heart", "Skeletal")
)

# verify
table(hs_paired$tissue, hs_paired$patient)
table(hs_paired$tissue, hs_paired$cell_type_unified)



# ══════════════════════════════════════════════════════════════════════════════
# PSEUDOBULK DE — Heart vs Skeletal (paired MDACC samples)
# Convention: positive log2FC = higher in Heart
# ══════════════════════════════════════════════════════════════════════════════

library(DESeq2)
library(dplyr)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(openxlsx)

celltypes_to_test <- c(
  "T/NK Lymphocytes",
  "Myeloid",
  "Endothelial",
  "Fibroblast/Stromal",
  "Smooth Muscle",
  "B Cells"
)

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — PSEUDOBULK DE
# ══════════════════════════════════════════════════════════════════════════════

run_pseudobulk_de <- function(seurat_obj, celltype) {
  
  cat("\nRunning:", celltype, "\n")
  
  sub <- subset(seurat_obj, cell_type_unified == celltype)
  
  counts <- AggregateExpression(
    sub,
    group.by      = c("patient", "tissue"),
    assays        = "RNA",
    return.seurat = FALSE
  )$RNA
  
  col_meta <- data.frame(
    sample    = colnames(counts),
    patient   = sub("_Heart|_Skeletal", "", colnames(counts)),
    tissue    = ifelse(grepl("Heart", colnames(counts)), "Heart", "Skeletal"),
    row.names = colnames(counts)
  )
  
  cat("Pseudobulk samples:\n")
  print(col_meta)
  
  if (sum(col_meta$tissue == "Heart") < 2 |
      sum(col_meta$tissue == "Skeletal") < 2) {
    cat("Skipping", celltype, "— insufficient replicates\n")
    return(NULL)
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = col_meta,
    design    = ~ patient + tissue
  )
  
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  dds <- DESeq(dds)
  
  # ── contrast: Heart vs Skeletal (Heart positive) ──────────────────────────
  res <- results(
    dds,
    contrast = c("tissue", "Heart", "Skeletal"),  # Heart is numerator → positive = Heart
    alpha    = 0.05
  )
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(
      cell_type = celltype,
      direction = ifelse(log2FoldChange > 0, "Heart", "Skeletal")  # positive = Heart
    )
  
  return(res_df)
}

# ── Run ───────────────────────────────────────────────────────────────────────
de_results <- lapply(celltypes_to_test, function(ct) {
  tryCatch(
    run_pseudobulk_de(hs_paired, ct),
    error = function(e) { cat("Error in", ct, ":", conditionMessage(e), "\n"); return(NULL) }
  )
})
names(de_results) <- celltypes_to_test

# ── Summary ───────────────────────────────────────────────────────────────────
lapply(celltypes_to_test, function(ct) {
  res <- de_results[[ct]]
  if (is.null(res)) return(NULL)
  sig <- res %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
  cat(ct, "— significant DEGs:", nrow(sig),
      "(", sum(sig$direction == "Heart"), "up in Heart,",
      sum(sig$direction == "Skeletal"), "up in Skeletal )\n")
})

# ── Save DE results ───────────────────────────────────────────────────────────
de_path <- file.path(TMOS_DIR, "Clean/fig1/paired_de/")

for (ct in celltypes_to_test) {
  if (!is.null(de_results[[ct]])) {
    write.csv(
      de_results[[ct]],
      file      = paste0(de_path, "Pseudobulk_DE_", gsub("/| ", "_", ct), "_Heart_vs_Skeletal.csv"),
      row.names = FALSE
    )
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — GSEA (positive NES = Heart)
# ══════════════════════════════════════════════════════════════════════════════

hallmark_sets  <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
cat("Total Hallmark gene sets:", length(gene_sets_list), "\n")

run_gsea_fgsea_full <- function(res_df, celltype) {
  
  cat("\nRunning GSEA for:", celltype, "\n")
  
  ranked <- res_df %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue)) %>%  # positive = Heart
    arrange(desc(rank_metric)) %>%
    distinct(gene, .keep_all = TRUE)
  
  gene_ranks <- setNames(ranked$rank_metric, ranked$gene)
  
  set.seed(42)
  suppressWarnings(
    capture.output(
      fgsea_res <- fgsea(
        pathways    = gene_sets_list,
        stats       = gene_ranks,
        minSize     = 10,
        maxSize     = 500,
        nPermSimple = 10000,
        nproc       = 1
      )
    )
  )
  
  fgsea_res <- as.data.frame(fgsea_res) %>%
    mutate(
      cell_type     = celltype,
      direction     = ifelse(NES > 0, "Heart", "Skeletal"),  # positive NES = Heart
      sig_label     = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      ),
      pathway_clean = gsub("HALLMARK_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean)
    )
  
  return(fgsea_res)
}

# ── Run ───────────────────────────────────────────────────────────────────────
gsea_results_full <- lapply(celltypes_to_test, function(ct) {
  tryCatch(
    run_gsea_fgsea_full(de_results[[ct]], ct),
    error = function(e) { cat("Error in", ct, ":", conditionMessage(e), "\n"); return(NULL) }
  )
})
names(gsea_results_full) <- celltypes_to_test

# ── Save Excel — all Hallmark per cell type ───────────────────────────────────
wb <- createWorkbook()

for (ct in celltypes_to_test) {
  if (!is.null(gsea_results_full[[ct]])) {
    sheet_name <- gsub("/", "_", ct)
    addWorksheet(wb, sheet_name)
    gsea_results_full[[ct]] %>%
      select(-leadingEdge) %>%
      arrange(padj) %>%
      writeData(wb, sheet_name, .)
  }
}

saveWorkbook(
  wb,
  file      = paste0(de_path, "GSEA_Hallmark_AllCellTypes_Heart_vs_Skeletal.xlsx"),
  overwrite = TRUE
)
cat("Excel saved\n")

# ── Combine results ───────────────────────────────────────────────────────────
gsea_combined_full <- bind_rows(gsea_results_full) %>%
  filter(!is.na(NES))

# ── Filter to TMOS-relevant pathways significant in 2+ cell types ─────────────
tmos_relevant <- c(
  "INTERFERON GAMMA RESPONSE",
  "INTERFERON ALPHA RESPONSE",
  "INFLAMMATORY RESPONSE",
  "TNFA SIGNALING VIA NFKB",
  "IL6 JAK STAT3 SIGNALING",
  "IL2 STAT5 SIGNALING",
  "ALLOGRAFT REJECTION",
  "COMPLEMENT",
  "APOPTOSIS",
  "P53 PATHWAY",
  "HYPOXIA",
  "OXIDATIVE PHOSPHORYLATION",
  "REACTIVE OXYGEN SPECIES PATHWAY",
  "EPITHELIAL MESENCHYMAL TRANSITION",
  "COAGULATION",
  "MTORC1 SIGNALING",
  "UNFOLDED PROTEIN RESPONSE",
  "ANGIOGENESIS",
  "E2F TARGETS",
  "G2M CHECKPOINT"
)

sig_pathways_2plus <- gsea_combined_full %>%
  filter(padj < 0.05) %>%
  group_by(pathway_clean) %>%
  summarise(n_sig = n()) %>%
  filter(n_sig >= 2) %>%
  pull(pathway_clean)

final_pathways <- intersect(sig_pathways_2plus, tmos_relevant)
cat("\nFinal pathways for heatmap:", length(final_pathways), "\n")
print(final_pathways)

# ── Build plot dataframe ──────────────────────────────────────────────────────
gsea_plot_df <- gsea_combined_full %>%
  filter(pathway_clean %in% final_pathways) %>%
  mutate(
    direction = ifelse(NES > 0, "Heart", "Skeletal"),
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

gsea_plot_df$sig_label[is.na(gsea_plot_df$sig_label)] <- ""

# ── Shared objects ────────────────────────────────────────────────────────────
immune_celltypes  <- c("T/NK Lymphocytes", "Myeloid", "B Cells")
stromal_celltypes <- c("Endothelial", "Fibroblast/Stromal", "Smooth Muscle")

pathway_order <- gsea_plot_df %>%
  group_by(pathway_clean) %>%
  summarise(mean_NES = mean(NES, na.rm = TRUE)) %>%
  arrange(desc(mean_NES)) %>%
  pull(pathway_clean)

gsea_theme <- theme_bw(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y     = element_text(size = 11),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, size = 13),
    plot.caption    = element_text(size = 9, hjust = 0),
    plot.margin     = margin(8, 8, 8, 8)
  )

gsea_fill <- scale_fill_gradient2(
  low      = "#4FA3D9",   # enriched in Skeletal
  mid      = "white",
  high     = "#F15A4A",   # enriched in Heart
  midpoint = 0,
  limits   = c(-3, 3),
  oob      = scales::squish,
  name     = "NES"
)

# ── Immune heatmap ────────────────────────────────────────────────────────────
p_immune <- gsea_plot_df %>%
  filter(cell_type %in% immune_celltypes) %>%
  mutate(
    cell_type     = factor(cell_type, levels = immune_celltypes),
    pathway_clean = factor(pathway_clean, levels = pathway_order)
  ) %>%
  ggplot(aes(x = cell_type, y = pathway_clean, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4.5, vjust = 0.75, color = "black") +
  gsea_fill +
  gsea_theme +
  labs(
    x       = NULL,
    y       = NULL,
    title   = "Immune compartment",
    caption = "* FDR < 0.05   ** FDR < 0.01   *** FDR < 0.001\nRed = enriched in Heart; Blue = enriched in Skeletal"
  )

# ── Stromal heatmap ───────────────────────────────────────────────────────────
p_stromal <- gsea_plot_df %>%
  filter(cell_type %in% stromal_celltypes) %>%
  mutate(
    cell_type     = factor(cell_type, levels = stromal_celltypes),
    pathway_clean = factor(pathway_clean, levels = pathway_order)
  ) %>%
  ggplot(aes(x = cell_type, y = pathway_clean, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4.5, vjust = 0.75, color = "black") +
  gsea_fill +
  gsea_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x       = NULL,
    y       = NULL,
    title   = "Stromal compartment",
    caption = ""
  )

# ── Combine and save heatmap ──────────────────────────────────────────────────
p_combined <- p_immune + p_stromal +
  plot_layout(widths = c(3, 3), guides = "collect") &
  theme(legend.position = "right")

print(p_combined)

ggsave(
  filename = "GSEA_Hallmark_NES_Immune_vs_Stromal_Split.pdf",
  plot     = p_combined,
  device   = cairo_pdf,
  path     = de_path,
  width    = 7,
  height   = 6.5
)

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — Individual enrichment plot: TNFA NFKB in T/NK
# ══════════════════════════════════════════════════════════════════════════════

ranked_tnk <- de_results[["T/NK Lymphocytes"]] %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
  arrange(desc(rank_metric)) %>%
  distinct(gene, .keep_all = TRUE)

gene_ranks_tnk <- setNames(ranked_tnk$rank_metric, ranked_tnk$gene)

pathway_genes <- gene_sets_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]

set.seed(42)
suppressWarnings(
  capture.output(
    gsea_tnk <- fgsea(
      pathways    = gene_sets_list["HALLMARK_TNFA_SIGNALING_VIA_NFKB"],
      stats       = gene_ranks_tnk,
      minSize     = 10,
      maxSize     = 500,
      nPermSimple = 10000,
      nproc       = 1
    )
  )
)

cat("TNFA NF-kB in T/NK — NES:", round(gsea_tnk$NES, 3),
    "| padj:", round(gsea_tnk$padj, 4), "\n")

p_enrich <- plotEnrichment(
  pathway = pathway_genes,
  stats   = gene_ranks_tnk
) +
  labs(
    title    = "TNFA Signaling via NF-kB — T/NK Lymphocytes",
    subtitle = paste0(
      "NES = ", round(gsea_tnk$NES, 3),
      " | FDR = ", formatC(gsea_tnk$padj, format = "e", digits = 2),
      " | Heart vs Skeletal (paired)"
    ),
    x = "Rank in gene list\n(left = enriched in Heart | right = enriched in Skeletal)",
    y = "Enrichment score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid    = element_blank(),
    axis.text     = element_text(color = "black"),
    plot.margin   = margin(8, 8, 8, 8)
  )

print(p_enrich)

ggsave(
  filename = "GSEA_Enrichment_TNFA_NFkB_TNK_Heart_vs_Skeletal.pdf",
  plot     = p_enrich,
  device   = cairo_pdf,
  path     = de_path,
  width    = 7,
  height   = 5
)



















# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — GSEA (Reactome Pathways)
# ══════════════════════════════════════════════════════════════════════════════

# Pull Reactome gene sets instead of Hallmark
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
gene_sets_list <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)
cat("Total Reactome gene sets:", length(gene_sets_list), "\n")

run_gsea_fgsea_full <- function(res_df, celltype) {
  
  cat("\nRunning GSEA (Reactome) for:", celltype, "\n")
  
  ranked <- res_df %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
    arrange(desc(rank_metric)) %>%
    distinct(gene, .keep_all = TRUE)
  
  gene_ranks <- setNames(ranked$rank_metric, ranked$gene)
  
  set.seed(42)
  suppressWarnings(
    capture.output(
      fgsea_res <- fgsea(
        pathways    = gene_sets_list,
        stats       = gene_ranks,
        minSize     = 15,    # Slightly increased for Reactome specificity
        maxSize     = 500,
        nPermSimple = 10000,
        nproc       = 1
      )
    )
  )
  
  fgsea_res <- as.data.frame(fgsea_res) %>%
    mutate(
      cell_type     = celltype,
      direction     = ifelse(NES > 0, "Heart", "Skeletal"),
      sig_label     = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      ),
      # Clean up Reactome naming (e.g., REACTOME_TCF_DEPENDENT_WNT_SIGNALING -> TCF Dependent Wnt Signaling)
      pathway_clean = gsub("REACTOME_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = stringr::str_to_title(pathway_clean)
    )
  
  return(fgsea_res)
}

# ── Run (Assuming de_results exists from Part 1) ──────────────────────────────
gsea_results_full <- lapply(celltypes_to_test, function(ct) {
  tryCatch(
    run_gsea_fgsea_full(de_results[[ct]], ct),
    error = function(e) { cat("Error in", ct, ":", conditionMessage(e), "\n"); return(NULL) }
  )
})
names(gsea_results_full) <- celltypes_to_test

# ── Save Excel ───────────────────────────────────────────────────────────────
# ── Save Excel ───────────────────────────────────────────────────────────────
wb <- createWorkbook()

for (ct in celltypes_to_test) {
  if (!is.null(gsea_results_full[[ct]])) {
    # 1. Clean the sheet name (Excel doesn't like / or long names)
    sheet_name <- substr(gsub("/", "_", ct), 1, 31) 
    
    addWorksheet(wb, sheet_name)
    
    # 2. Use dplyr::select explicitly to avoid the "inherited method" error
    data_to_write <- gsea_results_full[[ct]] %>%
      as.data.frame() %>% 
      dplyr::select(-leadingEdge) %>%
      arrange(padj)
    
    writeData(wb, sheet_name, data_to_write)
  }
}

saveWorkbook(wb, file = paste0(de_path, "GSEA_Reactome_AllCellTypes.xlsx"), overwrite = TRUE)
# ── Combine results ───────────────────────────────────────────────────────────
gsea_combined_full <- bind_rows(gsea_results_full) %>% filter(!is.na(NES))




# ══════════════════════════════════════════════════════════════════════════════
# PART 4 — Individual GSEA Plots (Top 3 Heart vs Skeletal in T/NK)
# ══════════════════════════════════════════════════════════════════════════════

# 1. Prepare Ranks for T/NK specifically
res_tnk <- de_results[["T/NK Lymphocytes"]]

ranked_tnk <- res_tnk %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
  arrange(desc(rank_metric)) %>%
  distinct(gene, .keep_all = TRUE)

gene_ranks_tnk <- setNames(ranked_tnk$rank_metric, ranked_tnk$gene)

# 2. Define Plotting Function for consistency
make_gsea_plot <- function(p_id, p_clean, nes, fdr) {
  plotEnrichment(gene_sets_list[[p_id]], gene_ranks_tnk) +
    labs(
      title    = p_clean, 
      subtitle = paste0("NES: ", round(nes, 3), " | FDR: ", formatC(fdr, format = "e", digits = 2)),
      x = "Rank in Gene List", y = "Enrichment Score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face="bold", size=11, hjust=0.5),
      plot.subtitle = element_text(size=9, hjust=0.5, color="grey30"),
      panel.grid    = element_blank()
    )
}

# ── A. HEART-HIGH PLOTS (Top 3 Significant) ───────────────────────────────────
top3_heart <- gsea_results_full[["T/NK Lymphocytes"]] %>%
  filter(padj < 0.05, NES > 0) %>%
  slice_max(order_by = NES, n = 3)

heart_plots <- list()
for (i in 1:nrow(top3_heart)) {
  heart_plots[[i]] <- make_gsea_plot(top3_heart$pathway[i], top3_heart$pathway_clean[i], 
                                     top3_heart$NES[i], top3_heart$padj[i])
}

graphics.off() # Clear any stuck connections
if(length(heart_plots) > 0) {
  p_heart_final <- wrap_plots(heart_plots, ncol = 1) + 
    plot_annotation(title = "Top 3 Significant Heart-Enriched Pathways (T/NK)",
                    theme = theme(plot.title = element_text(face="bold", size=14, hjust=0.5)))
  
  pdf(file.path(de_path, "GSEA_TNK_Top3_Heart_High.pdf"), width = 4, height = 11)
  print(p_heart_final)
  dev.off()
  cat("Heart-high plots saved.\n")
}

# ── B. SKELETAL-HIGH PLOTS (Top 3 Significant) ────────────────────────────────
top3_skeletal <- gsea_results_full[["T/NK Lymphocytes"]] %>%
  filter(padj < 0.05, NES < 0) %>%
  slice_min(order_by = NES, n = 3)

sk_plots <- list()
for (i in 1:nrow(top3_skeletal)) {
  sk_plots[[i]] <- make_gsea_plot(top3_skeletal$pathway[i], top3_skeletal$pathway_clean[i], 
                                  top3_skeletal$NES[i], top3_skeletal$padj[i])
}

if(length(sk_plots) > 0) {
  p_sk_final <- wrap_plots(sk_plots, ncol = 1) + 
    plot_annotation(title = "Top 3 Significant Skeletal-Enriched Pathways (T/NK)",
                    theme = theme(plot.title = element_text(face="bold", size=14, hjust=0.5)))
  
  pdf(file.path(de_path, "GSEA_TNK_Top3_Skeletal_High.pdf"), width = 4, height = 11)
  print(p_sk_final)
  dev.off()
  cat("Skeletal-high plots saved.\n")
}



















library(dplyr)
library(ggplot2)
library(ggrepel)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(writexl)

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — DEG SUMMARY MIRROR BAR PLOT
# ══════════════════════════════════════════════════════════════════════════════

padj_cut <- 0.05
lfc_cut  <- 0.58

col_heart    <- "#D95F02"
col_skeletal <- "#1B9E77"

# ── Build summary counts per cell type ───────────────────────────────────────
summary_counts <- lapply(celltypes_to_test, function(ct) {
  res <- de_results[[ct]]
  if (is.null(res)) return(NULL)
  data.frame(
    CellType    = ct,
    Heart_Up    = sum(res$padj < padj_cut & res$log2FoldChange >  lfc_cut, na.rm = TRUE),
    Skeletal_Up = sum(res$padj < padj_cut & res$log2FoldChange < -lfc_cut, na.rm = TRUE)
  )
})

plot_df <- bind_rows(summary_counts) %>%
  pivot_longer(
    cols      = c(Heart_Up, Skeletal_Up),
    names_to  = "Tissue",
    values_to = "Count"
  ) %>%
  mutate(
    DisplayCount = if_else(Tissue == "Skeletal_Up", -Count, Count),
    Tissue       = factor(Tissue, levels = c("Skeletal_Up", "Heart_Up"))
  )

# ── Mirror bar plot ───────────────────────────────────────────────────────────
p_mirror <- ggplot(
  plot_df,
  aes(x = reorder(CellType, abs(DisplayCount)), y = DisplayCount, fill = Tissue)
) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(
    aes(label = abs(DisplayCount)),
    hjust    = ifelse(plot_df$DisplayCount >= 0, -0.3, 1.3),
    size     = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.25, 0.25))
  ) +
  scale_fill_manual(
    values = c("Heart_Up" = col_heart, "Skeletal_Up" = col_skeletal),
    labels = c("Heart_Up" = "Heart",   "Skeletal_Up" = "Skeletal")
  ) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line            = element_blank(),
    axis.text.y          = element_text(face = "bold", color = "black"),
    legend.position      = "top",
    legend.title         = element_text(face = "bold"),
    panel.grid.major.x   = element_line(color = "grey92", linetype = "dashed")
  ) +
  labs(
    title    = "DE Summary: Heart vs Skeletal",
    subtitle = paste0("padj < ", padj_cut, " | |log2FC| > ", lfc_cut),
    x        = NULL,
    y        = "Number of enriched genes",
    fill     = "Enriched in:"
  )

print(p_mirror)

ggsave(
  filename = "DEG_Summary_Mirror_Plot_Heart_vs_Skeletal.pdf",
  plot     = p_mirror,
  device   = cairo_pdf,
  path     = de_path,
  width    = 5,
  height   = 5
)

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — VOLCANO PLOT: T/NK LYMPHOCYTES
# ══════════════════════════════════════════════════════════════════════════════

de_tnk <- de_results[["T/NK Lymphocytes"]] %>%
  mutate(
    neglog10_padj = -log10(padj),
    direction     = case_when(
      padj < padj_cut & log2FoldChange >  lfc_cut ~ "Heart",
      padj < padj_cut & log2FoldChange < -lfc_cut ~ "Skeletal",
      TRUE ~ "NS"
    ),
    sig = padj < padj_cut & abs(log2FoldChange) >= lfc_cut
  )

de_tnk$direction <- factor(de_tnk$direction, levels = c("Skeletal", "Heart", "NS"))

vol_colors <- c(
  "Heart"    = col_heart,
  "Skeletal" = col_skeletal,
  "NS"       = "grey82"
)

# top labels — top 10 by padj + top 10 by |log2FC| per direction
top_labels_p <- de_tnk %>%
  filter(direction %in% c("Heart", "Skeletal")) %>%
  group_by(direction) %>%
  slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  ungroup()

top_labels_fc <- de_tnk %>%
  filter(direction %in% c("Heart", "Skeletal"), sig) %>%
  group_by(direction) %>%
  slice_max(order_by = abs(log2FoldChange), n = 10, with_ties = FALSE) %>%
  ungroup()

top_labels <- bind_rows(top_labels_p, top_labels_fc) %>%
  distinct(direction, gene, .keep_all = TRUE)

p_volcano <- ggplot(de_tnk, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(aes(color = direction), size = 2.2, alpha = 1) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut),
             linetype = "dotted", color = "black") +
  geom_hline(yintercept = -log10(padj_cut),
             linetype = "dotted", color = "black") +
  geom_text_repel(
    data          = top_labels,
    aes(label     = gene),
    size          = 3.5,
    fontface      = "bold",
    max.overlaps  = Inf,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right"
  ) +
  labs(
    title    = "T/NK Lymphocytes: Heart vs Skeletal",
    subtitle = "MC+MS+MG | paired pseudobulk DESeq2 | positive log2FC = Heart",
    x        = "log2FoldChange (Heart vs Skeletal)",
    y        = expression(-log[10]("adjusted p-value")),
    color    = "Enriched in:"
  )

print(p_volcano)

ggsave(
  filename = "Volcano_TNK_Heart_vs_Skeletal.pdf",
  plot     = p_volcano,
  device   = cairo_pdf,
  path     = de_path,
  width    = 8,
  height   = 7
)








# ── 1. Get the Stats from your DE Results (Explicit dplyr) ────────────────────
tnk_stats <- de_results[["T/NK Lymphocytes"]] %>%
  filter(gene %in% c("NGF", "BRAF")) %>%
  # Use dplyr:: to bypass the mask error
  dplyr::select(gene, log2FoldChange, padj)

# ── 2. Get the VST Expression ────────────────────────────────────────────────
sub_tnk <- subset(hs_paired, cell_type_unified == "T/NK Lymphocytes")

# Aggregate
counts_tnk <- AggregateExpression(
  sub_tnk, 
  group.by = c("patient", "tissue"), 
  assays = "RNA", 
  return.seurat = FALSE
)$RNA

col_meta_tnk <- data.frame(
  sample  = colnames(counts_tnk),
  patient = sub("_Heart|_Skeletal", "", colnames(counts_tnk)),
  tissue  = ifelse(grepl("Heart", colnames(counts_tnk)), "Heart", "Skeletal"),
  row.names = colnames(counts_tnk)
)

# Create dds and run VST
dds_temp <- DESeqDataSetFromMatrix(
  countData = round(counts_tnk), 
  colData = col_meta_tnk, 
  design = ~patient + tissue
)
vsd_tnk <- vst(dds_temp, blind = FALSE)
vst_counts <- assay(vsd_tnk)

# ── 3. Join Stats with Expression ───────────────────────────────────────────
plot_df <- as.data.frame(t(vst_counts[c("NGF", "BRAF"), ])) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(col_meta_tnk, by = "sample") %>%
  tidyr::pivot_longer(cols = c("NGF", "BRAF"), names_to = "gene", values_to = "Expression") %>%
  left_join(tnk_stats, by = "gene") %>%
  mutate(
    tissue = factor(tissue, levels = c("Skeletal", "Heart")),
    # Format labels for the plot
    stat_label = paste0("log2FC: ", round(log2FoldChange, 2), 
                        "\nFDR: ", formatC(padj, format = "e", digits = 1))
  )

# ── 4. Final Plot ───────────────────────────────────────────────────────────
p_tnk_final <- ggplot(plot_df, aes(x = tissue, y = Expression, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5, lwd = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.4, size = 1.5) +
  facet_wrap(~gene, scales = "free_y") +
  # Add Stats Annotation
  geom_text(data = distinct(plot_df, gene, stat_label), 
            aes(label = stat_label), 
            x = 1.5, y = Inf, vjust = 1.5, size = 3.5, 
            fontface = "italic", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Heart" = "#D95F02", "Skeletal" = "#1B9E77")) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text        = element_text(color = "black"),
    axis.title.x     = element_blank(),
    legend.position  = "none"
  ) +
  labs(y = "VST Normalized Expression", title = "T/NK Lymphocytes: NGF & BRAF")

print(p_tnk_final)

# ── 5. Save ─────────────────────────────────────────────────────────────────
ggsave(
  filename = "Boxplot_TNK_NGF_BRAF_Final.pdf",
  plot     = p_tnk_final,
  device   = cairo_pdf,
  path     = de_path,
  width    = 4.5,
  height   = 3.5
)


##############################################################################
##  PART 3 - Fig 1F; Fig S1G, S1H: paired T/NK states and CellChat
##############################################################################

# ============================================================================
# Part 3
# Paired heart vs skeletal muscle (same TMOS patients 4352 and 4420):
#   Fig 1F  - T/NK composition, cytotoxicity and exhaustion, heart - skeletal
#   Fig S1G - CellChat change in signalling role per cell group
#   Fig S1H - CellChat differential interaction strength (heart - muscle)
#
# Figure panels : Fig 1F; Fig S1G, S1H
# Uses the T/NK state tables and paired CellChat objects made by
# fig2/Fig2_all_panels.R (Part 3).
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)
library(patchwork)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")

col_heart <- "#D95F02"
col_skel  <- "#1B9E77"

short_names <- c(
  "Effector memory CD8 T cells"   = "EM CD8",
  "Terminal effector CD8 T cells" = "TE CD8",
  "Central memory CD8 T cells"    = "CM CD8",
  "Naive CD8 T cells"             = "Naive CD8",
  "Terminal effector CD4 T cells" = "TE CD4",
  "Naive CD4 T cells"             = "Naive CD4",
  "T regulatory cells"            = "Treg",
  "Th1 cells"                     = "Th1",
  "Th1/Th17 cells"                = "Th1/Th17",
  "Th17 cells"                    = "Th17",
  "Th2 cells"                     = "Th2",
  "Follicular helper T cells"     = "Tfh",
  "MAIT cells"                    = "MAIT",
  "Natural killer cells"          = "NK",
  "Vd2 gd T cells"                = "Vδ2 γδ T",
  "Non-Vd2 gd T cells"            = "Non-Vδ2 γδ T"
)

relab <- function(x) {
  out <- unname(short_names[x])
  ifelse(is.na(out), x, out)
}

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(),
        strip.text.y = element_text(face = "bold", size = 9, angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 5))

## ---- State scores: cytotoxicity, exhaustion ----
pdelta <- read_excel(file.path(out_dir, "R2.1_TNK_states_summary.xlsx"),
                     sheet = "Paired_heart_minus_skel") %>%
  dplyr::filter(signature %in% c("Cytotoxicity", "Exhaustion")) %>%
  dplyr::mutate(monaco_fine = relab(monaco_fine)) %>%
  dplyr::select(patient, monaco_fine, signature, heart, skeletal, delta)

ct_order <- pdelta %>% dplyr::filter(signature == "Cytotoxicity") %>%
  dplyr::group_by(monaco_fine) %>% dplyr::summarise(m = mean(delta), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(m)) %>% dplyr::pull(monaco_fine)
lim <- max(abs(pdelta$delta))

p_hm <- pdelta %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = ct_order),
                signature = factor(signature, levels = c("Cytotoxicity", "Exhaustion"))) %>%
  ggplot(aes(monaco_fine, patient, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(signature ~ .) +
  scale_fill_gradient2(low = col_skel, mid = "white", high = col_heart, midpoint = 0,
                       limits = c(-lim, lim), name = "Δ score\n(heart −\nskeletal)") +
  labs(x = NULL, y = NULL) +
  bt

## ---- Composition ----
comp <- read_excel(file.path(out_dir, "TNK_state_scores_heart_paired.xlsx"),
                   sheet = "Paired_composition") %>%
  dplyr::mutate(patient = paste("Patient", patient_id),
                monaco_fine = relab(monaco_fine)) %>%
  dplyr::select(patient, tissue, monaco_fine, prop) %>%
  tidyr::complete(patient, tissue, monaco_fine, fill = list(prop = 0)) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = prop) %>%
  dplyr::mutate(delta = heart - skeletal) %>%
  dplyr::filter(monaco_fine %in% ct_order)

lim_c <- max(abs(comp$delta))

p_comp <- comp %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = ct_order)) %>%
  ggplot(aes(monaco_fine, patient, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid("Composition" ~ .) +
  scale_fill_gradient2(low = col_skel, mid = "white", high = col_heart, midpoint = 0,
                       limits = c(-lim_c, lim_c), name = "Δ % of T/NK\n(heart −\nskeletal)") +
  labs(x = NULL, y = NULL) +
  bt +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## ---- Combine ----
fig_paired <- (p_comp / p_hm) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(title = "Heart minus skeletal muscle T/NK cells (paired TMOS, n = 2)",
                  theme = theme(plot.title = element_text(face = "bold", size = 11)))
fig_paired

ggsave(file.path(out_dir, "Paired_TNK_composition_states_heart_minus_skeletal.pdf"), fig_paired,
       width = 6.5, height = 3, device = cairo_pdf)

write_xlsx(list(Composition = comp, State_scores = pdelta),
           file.path(out_dir, "Paired_TNK_composition_states_heart_minus_skeletal.xlsx"))



library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")

col_heart <- "#D95F02"
col_skel  <- "#1B9E77"

short_names <- c(
  "Effector memory CD8 T cells"   = "EM CD8",
  "Terminal effector CD8 T cells" = "TE CD8",
  "Central memory CD8 T cells"    = "CM CD8",
  "Naive CD8 T cells"             = "Naive CD8",
  "Terminal effector CD4 T cells" = "TE CD4",
  "Naive CD4 T cells"             = "Naive CD4",
  "T regulatory cells"            = "Treg",
  "Th1 cells"                     = "Th1",
  "Th1/Th17 cells"                = "Th1/Th17",
  "Th17 cells"                    = "Th17",
  "Th2 cells"                     = "Th2",
  "Follicular helper T cells"     = "Tfh",
  "MAIT cells"                    = "MAIT",
  "Natural killer cells"          = "NK",
  "Vd2 gd T cells"                = "Vδ2 γδ T",
  "Non-Vd2 gd T cells"            = "Non-Vδ2 γδ T"
)
relab <- function(x) {
  out <- unname(short_names[x])
  ifelse(is.na(out), x, out)
}

non_immune <- c("Endothelial", "Fibroblast/Stromal", "Smooth muscle/Mural", "Parenchymal")
compartment <- function(x) dplyr::case_when(
  x %in% non_immune ~ "Stromal/vascular",
  x == "Myeloid"    ~ "Myeloid",
  x == "B cells"    ~ "B cells",
  TRUE              ~ "T/NK")
comp_cols <- c(`Stromal/vascular` = "#7f7f7f", Myeloid = "#E69F00",
               `B cells` = "#0072B2", `T/NK` = "#b30000")

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold", size = 11),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

cc_heart <- readRDS(file.path(out_dir, "CellChat_heart_paired.rds"))
cc_skel  <- readRDS(file.path(out_dir, "CellChat_skeletal_paired.rds"))

common <- setdiff(intersect(levels(cc_heart@idents), levels(cc_skel@idents)), "Parenchymal")

comm <- dplyr::bind_rows(
  subsetCommunication(cc_heart) %>% dplyr::mutate(tissue = "Heart"),
  subsetCommunication(cc_skel)  %>% dplyr::mutate(tissue = "Skeletal muscle")) %>%
  dplyr::mutate(source = as.character(source), target = as.character(target)) %>%
  dplyr::filter(source %in% common, target %in% common)

grp_order <- relab(c(intersect(c(non_immune, "Myeloid", "B cells"), common),
                     sort(setdiff(common, c(non_immune, "Myeloid", "B cells")))))

## ---- F: differential interaction strength ----
pair_w <- comm %>% dplyr::group_by(tissue, source, target) %>%
  dplyr::summarise(w = sum(prob), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue, values_from = w, values_fill = 0) %>%
  tidyr::complete(source = common, target = common, fill = list(Heart = 0, `Skeletal muscle` = 0)) %>%
  dplyr::mutate(delta  = Heart - `Skeletal muscle`,
                source = factor(relab(source), levels = grp_order),
                target = factor(relab(target), levels = grp_order))
lim <- max(abs(pair_w$delta))

pF <- ggplot(pair_w, aes(target, source, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  coord_fixed() +
  scale_fill_gradient2(low = col_skel, mid = "white", high = col_heart, midpoint = 0,
                       limits = c(-lim, lim), name = "Heart −\nmuscle") +
  labs(x = "Receiver", y = "Sender",
       title = "Differential interaction strength (heart − skeletal muscle)") +
  bt + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
             axis.text.y = element_text(size = 9))
pF

## ---- G: change in signalling role ----
roles <- dplyr::full_join(
  comm %>% dplyr::group_by(tissue, group = source) %>% dplyr::summarise(out = sum(prob), .groups = "drop"),
  comm %>% dplyr::group_by(tissue, group = target) %>% dplyr::summarise(inc = sum(prob), .groups = "drop"),
  by = c("tissue", "group")) %>%
  dplyr::mutate(dplyr::across(c(out, inc), ~ tidyr::replace_na(.x, 0))) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = c(out, inc), values_fill = 0) %>%
  dplyr::mutate(d_out = out_Heart - `out_Skeletal muscle`,
                d_in  = inc_Heart - `inc_Skeletal muscle`,
                comp  = compartment(group),
                group = relab(group))

pG <- ggplot(roles, aes(d_in, d_out, colour = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  geom_text_repel(aes(label = group), size = 3, seed = 1, show.legend = FALSE,
                  max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                  min.segment.length = 0, segment.colour = "grey60", force = 2) +
  scale_colour_manual(values = comp_cols, name = NULL) +
  labs(x = "Δ incoming strength (heart − muscle)", y = "Δ outgoing strength (heart − muscle)",
       title = "Change in signalling role per cell group") +
  bt + theme(aspect.ratio = 1, legend.position = "bottom")
pG

ggsave(file.path(out_dir, "CellChat_F_diff_heatmap.pdf"), pF, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "CellChat_G_roles.pdf"),        pG, width = 5, height = 5, device = cairo_pdf)
