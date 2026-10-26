# ============================================================================
# Fig2_all_panels.R
# Figure 2 and Supplementary Figure S2 - all panels in one script.
#
#   Panel(s)                 Part
#   2A                       -    schematic, no code
#   2B, 2C, 2D               1
#   2E; S2A, S2B, S2C        2
#   2F, 2G, 2H               3    (also makes inputs for Fig 1F, S1G, S1H)
#   2I, 2J                   4
#   2K, 2L, 2M; S2I          5
#   2N                       6
#   S2D, S2E                 7
#   S2F                      8a (CellChat objects) + 8b (plots)
#   S2G                      -    table, no code
#   S2H                      -    serum cytokines, analysed by collaborators
#   S2J                      9
# ============================================================================
if (!exists("TMOS_DIR")) source(file.path("config", "paths.R"))


##############################################################################
##  PART 1 - Fig 2B, 2C, 2D: heart T/NK TMOS vs MC, cell-level MAST DE and Reactome ORA
##############################################################################

# ============================================================================
# Heart T/NK TMOS vs MC: cell-level MAST DE, volcano, Reactome ORA (ReactomePA)
#
# Figure panels : Fig 2B, 2C, 2D
# ============================================================================

############################################################
## TNK CLEAN – DE: MC vs MC+MS+MG (MAST, institution-adjusted)
## Adds: enriched_in (which group higher)
############################################################

library(Seurat)
library(dplyr)
library(tibble)
library(writexl)
library(MAST)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
obj_path <- file.path(TMOS_DIR, "Clean/fig2/tnk_clean.rds")

outdir <- file.path(TMOS_DIR, "Clean/fig2/Heart_TNK cells")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, "TNK_DE_MAST_MCMSMG_vs_MC_institutionAdjusted.xlsx")

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
tnk_clean <- readRDS(obj_path)
DefaultAssay(tnk_clean) <- "RNA"

# Sanity checks
stopifnot(all(c("m_group", "institution") %in% colnames(tnk_clean@meta.data)))

# ------------------------------------------------------------
# Subset to only MC and MC+MS+MG
# ------------------------------------------------------------
keep_groups <- c("MC", "MC+MS+MG")
tnk_de <- subset(tnk_clean, subset = m_group %in% keep_groups)

tnk_de$m_group <- factor(tnk_de$m_group, levels = keep_groups)
Idents(tnk_de) <- "m_group"

# Optional: confirm cells per group
print(table(tnk_de$m_group))
print(table(tnk_de$institution, tnk_de$m_group))

# ------------------------------------------------------------
# Differential expression (MAST, institution adjusted)
# NOTE: avg_log2FC is (ident.1 - ident.2) on log-normalized scale
# Here: ident.1 = "MC+MS+MG" vs ident.2 = "MC"
# Positive avg_log2FC => higher in MC+MS+MG
# ------------------------------------------------------------
de <- FindMarkers(
  tnk_de,
  ident.1 = "MC+MS+MG",
  ident.2 = "MC",
  test.use = "MAST",
  latent.vars = c("institution"),
  min.pct = 0,
  logfc.threshold = 0
)

# ------------------------------------------------------------
# Clean + annotate results
# ------------------------------------------------------------
de_out <- de %>%
  rownames_to_column("gene") %>%
  mutate(
    comparison = "MC+MS+MG_vs_MC",
    group1 = "MC+MS+MG",
    group2 = "MC",
    enriched_in = case_when(
      avg_log2FC > 0  ~ "MC+MS+MG",
      avg_log2FC < 0  ~ "MC",
      TRUE            ~ "Equal"
    )
  ) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

# Optional: quick counts of enriched direction (all genes)
print(table(de_out$enriched_in))

# Optional: "significant" subset
de_sig <- de_out %>%
  filter(p_val_adj < 0.05)

# Optional: significant + |log2FC| threshold (pick your favorite)
de_sig_fc <- de_out %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.585)

# ------------------------------------------------------------
# Save Excel (multiple sheets)
# ------------------------------------------------------------
write_xlsx(
  list(
    All_Genes = de_out,
    Significant_padj_lt_0_05 = de_sig,
    Significant_padj_lt_0_05_abslog2FC_ge_0_585 = de_sig_fc
  ),
  path = outfile
)

message("✅ Saved TNK DE Excel to:/n", outfile)










############################################################
## TNK MAST DE (institution-adjusted) — downstream analysis
## Input: TNK_DE_MAST_MCMSMG_vs_MC_institutionAdjusted.xlsx
## Output: volcano + gene lists + Reactome enrichment + plots
############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(writexl)

# Reactome enrichment (ORA)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
indir <- file.path(TMOS_DIR, "Clean/fig2/Heart_TNK cells")

infile <- file.path(indir, "TNK_DE_MAST_MCMSMG_vs_MC_institutionAdjusted.xlsx")

outdir <- indir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# thresholds (match your prior choices)
padj_cut <- 0.05
lfc_cut  <- 0.585

# ------------------------------------------------------------
# 1) Load DE
# ------------------------------------------------------------
de_all <- read_excel(infile, sheet = "All_Genes")

# If gene symbols sometimes have trailing "...", clean it
de_all <- de_all %>%
  mutate(gene = str_replace(gene, "//.//.//..*$", ""))

# Safety: ensure columns exist
stopifnot(all(c("gene","avg_log2FC","p_val_adj","enriched_in") %in% colnames(de_all)))

# ------------------------------------------------------------
# 2) Define DE sets (Up in MC+MS+MG vs Up in MC)
# ------------------------------------------------------------
de_all <- de_all %>%
  mutate(
    neglog10_padj = -log10(p_val_adj),
    sig = p_val_adj < padj_cut & abs(avg_log2FC) >= lfc_cut,
    direction = case_when(
      sig & avg_log2FC > 0 ~ "Up in MC+MS+MG",
      sig & avg_log2FC < 0 ~ "Up in MC",
      TRUE ~ "Not significant"
    )
  )

genes_up_mcmsmg <- de_all %>%
  filter(direction == "Up in MC+MS+MG") %>%
  pull(gene) %>% unique()

genes_up_mc <- de_all %>%
  filter(direction == "Up in MC") %>%
  pull(gene) %>% unique()

cat("Up in MC+MS+MG:", length(genes_up_mcmsmg), "/n")
cat("Up in MC:", length(genes_up_mc), "/n")

# ------------------------------------------------------------
# 3) Volcano plot
# ------------------------------------------------------------
# ------------------------------------------------------------
# Volcano plot with EXACT group colors
# ------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

# Force filtering
de_plot <- de_all %>%
  mutate(
    gene = str_replace(gene, "//.//.//..*$", ""),
    neglog10_padj = -log10(p_val_adj)
  ) %>%
  filter(pct.1 > 0.25 | pct.2 > 0.25) %>%
  mutate(
    direction = case_when(
      p_val_adj < padj_cut & avg_log2FC >  lfc_cut  ~ "MC+MS+MG",
      p_val_adj < padj_cut & avg_log2FC < -lfc_cut  ~ "MC",
      TRUE ~ "NS"
    )
  )

# Make direction a factor with controlled order
de_plot$direction <- factor(
  de_plot$direction,
  levels = c("MC", "MC+MS+MG", "NS")
)

# Exact colors requested
vol_colors <- c(
  "MC"        = "#4682B4",
  "MC+MS+MG"  = "#b30000",
  "NS"        = "grey75"
)

# Top labels
top_labels <- de_plot %>%
  filter(direction %in% c("MC","MC+MS+MG")) %>%
  group_by(direction) %>%
  slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
  ungroup()

p_volcano <- ggplot(de_plot, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = direction), size = 1.8, alpha = 0.8) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    title = "TNK MAST DE: MC+MS+MG vs MC",
    subtitle = "Filtered: pct > 0.10 in ≥1 group",
    x = "avg_log2FC (MC+MS+MG − MC)",
    y = expression(-log[10]("adjusted p-value")),
    color = "Enriched in"
  )

ggsave(
  file.path(outdir, "TNK_MAST_Volcano_MCMSMG_vs_MC_pct0.10.pdf"),
  plot = p_volcano,
  width = 7.5,
  height = 6,
  device = cairo_pdf
)


library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

# Force filtering
de_plot <- de_all %>%
  mutate(
    gene = str_replace(gene, "//.//.//..*$", ""),
    neglog10_padj = -log10(p_val_adj)
  ) %>%
  filter(pct.1 > 0.25 | pct.2 > 0.25) %>%   # you set 0.25 here
  mutate(
    direction = case_when(
      p_val_adj < padj_cut & avg_log2FC >  lfc_cut  ~ "MC+MS+MG",
      p_val_adj < padj_cut & avg_log2FC < -lfc_cut  ~ "MC",
      TRUE ~ "NS"
    ),
    sig = p_val_adj < padj_cut & abs(avg_log2FC) >= lfc_cut
  )

# Make direction a factor with controlled order
de_plot$direction <- factor(
  de_plot$direction,
  levels = c("MC", "MC+MS+MG", "NS")
)

# Exact colors requested
vol_colors <- c(
  "MC"        = "#4682B4",
  "MC+MS+MG"  = "#b30000",
  "NS"        = "grey75"
)

# ------------------------------------------------------------
# Labels set A: Top 10 by adjusted p-value (per direction)
# ------------------------------------------------------------
top_labels_p <- de_plot %>%
  filter(direction %in% c("MC","MC+MS+MG")) %>%
  group_by(direction) %>%
  slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label_type = "Top10_padj")

# ------------------------------------------------------------
# Labels set B: Top 10 by |log2FC| among significant genes (per direction)
# ------------------------------------------------------------
top_labels_fc <- de_plot %>%
  filter(direction %in% c("MC","MC+MS+MG"), sig) %>%
  group_by(direction) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label_type = "Top10_abslog2FC")

# Combine + deduplicate genes (avoid double labels)
top_labels <- bind_rows(top_labels_p, top_labels_fc) %>%
  distinct(direction, gene, .keep_all = TRUE)

p_volcano <- ggplot(de_plot, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = direction), size = 2.6, alpha = 0.82) +  # bigger dots
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25
  ) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    title = "TNK MAST DE: MC+MS+MG vs MC",
    subtitle = "Filtered: pct > 0.25 in ≥1 group | Labels: top10 padj + top10 |log2FC| (sig only)",
    x = "avg_log2FC (MC+MS+MG − MC)",
    y = expression(-log[10]("adjusted p-value")),
    color = "Enriched in"
  )

ggsave(
  file.path(outdir, "TNK_MAST_Volcano_MCMSMG_vs_MC_pct0.25_Top10padj_Top10abslog2FC.pdf"),
  plot = p_volcano,
  width = 7.5,
  height = 6,
  device = cairo_pdf
)


# ------------------------------------------------------------
# 4) Save gene lists to Excel
# ------------------------------------------------------------
gene_lists <- list(
  Up_in_MCMSMG = data.frame(gene = genes_up_mcmsmg),
  Up_in_MC     = data.frame(gene = genes_up_mc)
)

write_xlsx(
  gene_lists,
  path = file.path(outdir, "TNK_MAST_GeneLists_padj0.05_logFC0.585.xlsx")
)

# ------------------------------------------------------------
# 5) Reactome enrichment (ORA) for each list
# ------------------------------------------------------------
convert_to_entrez <- function(genes) {
  bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  ) %>%
    distinct(ENTREZID)
}

run_reactome_ora <- function(entrez_ids) {
  enrichPathway(
    gene          = entrez_ids,
    organism      = "human",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

react_mcmsmg <- NULL
react_mc     <- NULL

if (length(genes_up_mcmsmg) >= 10) {
  entrez_mcmsmg <- convert_to_entrez(genes_up_mcmsmg)
  if (nrow(entrez_mcmsmg) > 0) react_mcmsmg <- run_reactome_ora(entrez_mcmsmg$ENTREZID)
}

if (length(genes_up_mc) >= 10) {
  entrez_mc <- convert_to_entrez(genes_up_mc)
  if (nrow(entrez_mc) > 0) react_mc <- run_reactome_ora(entrez_mc$ENTREZID)
}

react_df_mcmsmg <- if (!is.null(react_mcmsmg)) as.data.frame(react_mcmsmg) else data.frame()
react_df_mc     <- if (!is.null(react_mc))     as.data.frame(react_mc)     else data.frame()

write_xlsx(
  list(
    Reactome_Up_in_MCMSMG = react_df_mcmsmg,
    Reactome_Up_in_MC     = react_df_mc
  ),
  path = file.path(outdir, "TNK_MAST_Reactome_ORA.xlsx")
)

# ------------------------------------------------------------
# 6) Simple Reactome dotplots (top 15)
# ------------------------------------------------------------
library(enrichplot)

if (!is.null(react_mcmsmg)) {
  p1 <- dotplot(react_mcmsmg, showCategory = 15) +
    ggtitle("Reactome ORA — Up in MC+MS+MG")
  ggsave(file.path(outdir, "Reactome_ORA_Up_in_MCMSMG.pdf"),
         plot = p1, width = 8, height = 6, device = cairo_pdf)
}

if (!is.null(react_mc)) {
  p2 <- dotplot(react_mc, showCategory = 15) +
    ggtitle("Reactome ORA — Up in MC")
  ggsave(file.path(outdir, "Reactome_ORA_Up_in_MC.pdf"),
         plot = p2, width = 8, height = 6, device = cairo_pdf)
}

message("✅ Done. Outputs saved to:/n", outdir)






library(dplyr)
library(stringr)
library(ggplot2)
library(writexl)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1) Gene lists per group (SIGNIFICANT + direction)
# Uses same thresholds as volcano: padj_cut + lfc_cut
# ------------------------------------------------------------
genes_up_mcmsmg <- de_all %>%
  mutate(gene = str_replace(gene, "//.//.//..*$", "")) %>%
  filter(p_val_adj < padj_cut, avg_log2FC >  lfc_cut) %>%
  pull(gene) %>% unique()

genes_up_mc <- de_all %>%
  mutate(gene = str_replace(gene, "//.//.//..*$", "")) %>%
  filter(p_val_adj < padj_cut, avg_log2FC < -lfc_cut) %>%
  pull(gene) %>% unique()

cat("Genes up in MC+MS+MG:", length(genes_up_mcmsmg), "/n")
cat("Genes up in MC:", length(genes_up_mc), "/n")

# ------------------------------------------------------------
# 2) Helper: SYMBOL -> ENTREZ
# ------------------------------------------------------------
to_entrez <- function(symbols) {
  bitr(symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
}

entrez_mcmsmg <- to_entrez(genes_up_mcmsmg)
entrez_mc     <- to_entrez(genes_up_mc)

# ------------------------------------------------------------
# 3) Reactome ORA
# ------------------------------------------------------------
run_reactome <- function(entrez_ids) {
  enrichPathway(
    gene          = entrez_ids,
    organism      = "human",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

react_mcmsmg <- if (length(entrez_mcmsmg) >= 10) run_reactome(entrez_mcmsmg) else NULL
react_mc     <- if (length(entrez_mc)     >= 10) run_reactome(entrez_mc)     else NULL

df_mcmsmg <- if (!is.null(react_mcmsmg)) as.data.frame(react_mcmsmg) else data.frame()
df_mc     <- if (!is.null(react_mc))     as.data.frame(react_mc)     else data.frame()

# ------------------------------------------------------------
# 4) Save Excel
# ------------------------------------------------------------
out_xlsx <- file.path(outdir, "TNK_MAST_Reactome_ORA_ByGroup.xlsx")

write_xlsx(
  list(
    Genes_Up_in_MCMSMG = data.frame(gene = genes_up_mcmsmg),
    Genes_Up_in_MC     = data.frame(gene = genes_up_mc),
    Reactome_Up_in_MCMSMG = df_mcmsmg,
    Reactome_Up_in_MC     = df_mc
  ),
  path = out_xlsx
)

message("✅ Saved Reactome ORA Excel to:/n", out_xlsx)

# ------------------------------------------------------------
# 5) Barplot function: top 10, pathway name printed on bar
# ------------------------------------------------------------
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

make_bar <- function(df, title_text, fill_col) {
  if (nrow(df) == 0) return(NULL)
  
  df2 <- prep_bar_df(df)
  
  ggplot(df2, aes(x = neglog10_padj, y = Pathway_wrap)) +
    geom_col(fill = fill_col, width = 0.86) +
    geom_text(
      aes(label = Pathway_wrap),
      x = 0.02, hjust = 0,
      size = 4.2
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
      x = expression(-log[10]("adjusted p-value")),
      y = NULL
    )
}

# Use your group colors
col_mc     <- "#4682B4"
col_mcmsmg <- "#b30000"

p_mcmsmg <- make_bar(df_mcmsmg, "Reactome ORA — Up in MC+MS+MG", col_mcmsmg)
p_mc     <- make_bar(df_mc,     "Reactome ORA — Up in MC",       col_mc)

# ------------------------------------------------------------
# 6) Save plots
# ------------------------------------------------------------
if (!is.null(p_mcmsmg)) {
  ggsave(
    filename = file.path(outdir, "Reactome_ORA_Top10_Up_in_MCMSMG.pdf"),
    plot = p_mcmsmg,
    width = 9,
    height = 6,
    device = cairo_pdf
  )
}

if (!is.null(p_mc)) {
  ggsave(
    filename = file.path(outdir, "Reactome_ORA_Top10_Up_in_MC.pdf"),
    plot = p_mc,
    width = 9,
    height = 6,
    device = cairo_pdf
  )
}

message("✅ Saved Reactome barplots to:/n", outdir)


library(scales)

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

make_bar_gradient <- function(df, title_text, base_color) {
  
  if (nrow(df) == 0) return(NULL)
  
  df2 <- prep_bar_df(df)
  
  # create lighter-to-strong gradient
  low_col  <- alpha(base_color, 0.3)
  high_col <- base_color
  
  ggplot(df2, aes(x = neglog10_padj, y = Pathway_wrap, fill = FoldEnrichment)) +
    geom_col(width = 0.86) +
    geom_text(
      aes(label = Pathway_wrap),
      x = 0.02, hjust = 0,
      size = 4.2
    ) +
    scale_fill_gradient(
      low  = low_col,
      high = high_col,
      name = "Fold Enrichment"
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
      x = expression(-log[10]("adjusted p-value")),
      y = NULL
    )
}



col_mc     <- "#4682B4"
col_mcmsmg <- "#b30000"

p_mcmsmg <- make_bar_gradient(
  df_mcmsmg,
  "Reactome ORA — Up in MC+MS+MG",
  col_mcmsmg
)

p_mc <- make_bar_gradient(
  df_mc,
  "Reactome ORA — Up in MC",
  col_mc
)




if (!is.null(p_mcmsmg)) {
  ggsave(
    filename = file.path(outdir, "Reactome_ORA_Top10_Up_in_MCMSMG_gradient.pdf"),
    plot = p_mcmsmg,
    width = 6,
    height = 5,
    device = cairo_pdf
  )
}

if (!is.null(p_mc)) {
  ggsave(
    filename = file.path(outdir, "Reactome_ORA_Top10_Up_in_MC_gradient.pdf"),
    plot = p_mc,
    width = 6,
    height = 5,
    device = cairo_pdf
  )
}


######################## GSEA ##########################################


library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)



# ------------------------------------------------------------
# 1) Prepare the Ranked Gene List
# ------------------------------------------------------------
# We use: -log10(p-value) * sign(log2FC)
# This ensures highly significant UP genes are at the top (+) 
# and highly significant DOWN genes are at the bottom (-)
de_all_ranked <- de_all %>%
  mutate(gene = str_replace(gene, "//.//.//..*$", "")) %>%
  # Handle p-values of 0 to avoid Infinity
  mutate(p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj)) %>%
  mutate(rank_metric = -log10(p_val_adj) * sign(avg_log2FC)) %>%
  arrange(desc(rank_metric)) %>%
  distinct(gene, .keep_all = TRUE)

# Create a named vector for fgsea
ranks <- de_all_ranked$rank_metric
names(ranks) <- de_all_ranked$gene

# ------------------------------------------------------------
# 2) Fetch Hallmark Gene Sets from MSigDB
# ------------------------------------------------------------
h_df <- msigdbr(species = "human", category = "H")
h_list <- split(x = h_df$gene_symbol, f = h_df$gs_name)

# ------------------------------------------------------------
# 3) Run GSEA
# ------------------------------------------------------------
fgsea_res <- fgsea(pathways = h_list, 
                   stats = ranks,
                   minSize = 15,
                   maxSize = 500)

# Clean results for saving
fgsea_out <- fgsea_res %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(enriched_in = ifelse(NES > 0, "MC+MS+MG", "MC"))

# Save to Excel
write_xlsx(list(GSEA_Hallmark = fgsea_out), 
           path = file.path(outdir, "TNK_GSEA_Hallmark_MCMSMG_vs_MC.xlsx"))

# ------------------------------------------------------------
# 4) Plotting: GSEA Normalized Enrichment Score (NES)
# ------------------------------------------------------------
# Filter for significant pathways for the plot
sig_pathways <- fgsea_out %>% 
  filter(padj < 0.05) %>% 
  slice_max(order_by = abs(NES), n = 20)

p_gsea <- ggplot(sig_pathways, aes(x = NES, y = reorder(pathway, NES), fill = enriched_in)) +
  geom_col() +
  scale_fill_manual(values = c("MC" = "#fdae61", "MC+MS+MG" = "#b30000")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Hallmark GSEA: MC+MS+MG vs MC (T/NK)",
    subtitle = "Ranked by -log10(padj) * sign(logFC)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "Enriched in"
  ) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))

ggsave(file.path(outdir, "TNK_GSEA_Hallmark_NES_Plot.pdf"), 
       plot = p_gsea, width = 10, height = 7, device = cairo_pdf)





setwd(file.path(TMOS_DIR, "Single_cell/Heart/Seurat_by_CellType/TNK cells"))
tnk_clean <- readRDS("tnk_clean.rds")





library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(rstatix)

# --------------------------------------------------
# Settings
# --------------------------------------------------

outdir <- file.path(TMOS_DIR, "Spatial/Heart/MAPK_ReactomePA/MAP_kinase_activation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(tnk_clean) <- "RNA"

genes <- c(
  "ATF1","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","IKBKB",
  "IRAK1","IRAK2","MEF2A","MEF2C","NFKB1","PPP2CA","PPP2R1A",
  "MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
  "MAP2K1","MAP2K4","SKP1","MAP3K7","TRAF6","CUL1","IKBKG",
  "RIPK2","BTRC","RPS6KA5","MAPKAPK2","NOD1","TAB1","TAB2",
  "FBXW11","NOD2"
)

genes <- intersect(genes, rownames(tnk_clean))

m_group_colors <- c(
  "control"   = "#bdbdbd",
  "MC"        = "#fdae61",
  "MC+MS"     = "#e34a33",
  "MC+MS+MG"  = "#b30000"
)

group_levels <- c("control","MC","MC+MS","MC+MS+MG")

# --------------------------------------------------
# Extract expression
# --------------------------------------------------

expr_mat <- GetAssayData(
  tnk_clean,
  slot = "data"
)[genes, ]

meta <- tnk_clean@meta.data

df_long <- as.data.frame(t(expr_mat)) %>%
  mutate(m_group = meta$m_group) %>%
  filter(m_group != "MC+MG") %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(m_group = factor(m_group, levels = group_levels))

# --------------------------------------------------
# Summary stats
# --------------------------------------------------

sum_df <- df_long %>%
  group_by(Gene, m_group) %>%
  summarise(
    n_cells = n(),
    mean = mean(Expression, na.rm = TRUE),
    sem  = sd(Expression, na.rm = TRUE) / sqrt(n_cells),
    .groups = "drop"
  )

# Save summary
write.xlsx(
  sum_df,
  file.path(outdir, "MAPK_genes_TNK_mgroup_meanSEM_no_MCMG.xlsx"),
  overwrite = TRUE
)

# --------------------------------------------------
# Pairwise statistics per gene
# --------------------------------------------------

wilcox_all <- df_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ m_group, p.adjust.method = "BH") %>%
  ungroup()

ttest_all <- df_long %>%
  group_by(Gene) %>%
  pairwise_t_test(
    Expression ~ m_group,
    p.adjust.method = "BH",
    pool.sd = FALSE
  ) %>%
  ungroup()

# Save stats
wb <- createWorkbook()

addWorksheet(wb, "Wilcox_pairwise_BH")
writeData(wb, "Wilcox_pairwise_BH", wilcox_all)

addWorksheet(wb, "Ttest_pairwise_BH")
writeData(wb, "Ttest_pairwise_BH", ttest_all)

saveWorkbook(
  wb,
  file.path(outdir, "MAPK_genes_TNK_mgroup_pairwise_stats.xlsx"),
  overwrite = TRUE
)

cat("✅ Saved pairwise stats/n")

# --------------------------------------------------
# Order genes by overall mean
# --------------------------------------------------

gene_order <- sum_df %>%
  group_by(Gene) %>%
  summarise(overall_mean = mean(mean)) %>%
  arrange(desc(overall_mean)) %>%
  pull(Gene)

sum_df$Gene <- factor(sum_df$Gene, levels = gene_order)

# --------------------------------------------------
# 1️⃣ BARPLOT (mean ± SEM)
# --------------------------------------------------

p_bar <- ggplot(sum_df, aes(x = m_group, y = mean, fill = m_group)) +
  geom_col(width = 0.75) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(~Gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = m_group_colors) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "MAPK Activation Genes in T/NK cells",
    subtitle = "Mean ± SEM",
    x = "m_group",
    y = "Mean expression"
  )

ggsave(
  file.path(outdir, "MAPK_genes_TNK_mgroup_barplot_no_MCMG.pdf"),
  p_bar,
  width = 16,
  height = 18
)

# --------------------------------------------------
# 2️⃣ LINE PLOT (mean trend)
# --------------------------------------------------

p_line <- ggplot(sum_df, aes(x = m_group, y = mean, group = 1)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(aes(color = m_group), size = 2.5) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.15
  ) +
  facet_wrap(~Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = m_group_colors) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "MAPK Activation Genes in T/NK cells",
    subtitle = "Mean expression trend across groups",
    x = "m_group",
    y = "Mean expression"
  )

ggsave(
  file.path(outdir, "MAPK_genes_TNK_mgroup_lineplot_no_MCMG.pdf"),
  p_line,
  width = 16,
  height = 18
)

cat("✅ All outputs saved successfully/n")










library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------
# Settings
# --------------------------------------------------

outdir <- file.path(TMOS_DIR, "Spatial/Heart/MAPK_ReactomePA/MAP_kinase_activation")

DefaultAssay(tnk_clean) <- "RNA"

genes <- c(
  "ATF1","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","IKBKB",
  "IRAK1","IRAK2","MEF2A","MEF2C","NFKB1","PPP2CA","PPP2R1A",
  "MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
  "MAP2K1","MAP2K4","SKP1","MAP3K7","TRAF6","CUL1","IKBKG",
  "RIPK2","BTRC","RPS6KA5","MAPKAPK2","NOD1","TAB1","TAB2",
  "FBXW11","NOD2"
)

genes <- intersect(genes, rownames(tnk_clean))

m_group_colors <- c(
  "control"   = "#bdbdbd",
  "MC"        = "#fdae61",
  "MC+MS"     = "#e34a33",
  "MC+MS+MG"  = "#b30000"
)

# --------------------------------------------------
# Extract expression (remove MC+MG globally)
# --------------------------------------------------

expr_mat <- GetAssayData(tnk_clean, layer = "data")[genes, ]
meta <- tnk_clean@meta.data

df_long <- as.data.frame(t(expr_mat)) %>%
  mutate(m_group = meta$m_group) %>%
  filter(m_group != "MC+MG") %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Expression"
  )

# --------------------------------------------------
# Function to summarize + plot
# --------------------------------------------------

make_plots <- function(df, groups_keep, label_suffix) {
  
  df_sub <- df %>%
    filter(m_group %in% groups_keep) %>%
    mutate(m_group = factor(m_group, levels = groups_keep))
  
  sum_df <- df_sub %>%
    group_by(Gene, m_group) %>%
    summarise(
      mean = mean(Expression, na.rm = TRUE),
      sem  = sd(Expression, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    )
  
  # Order genes by overall mean
  gene_order <- sum_df %>%
    group_by(Gene) %>%
    summarise(overall_mean = mean(mean)) %>%
    arrange(desc(overall_mean)) %>%
    pull(Gene)
  
  sum_df$Gene <- factor(sum_df$Gene, levels = gene_order)
  
  # -------- BARPLOT --------
  p_bar <- ggplot(sum_df, aes(x = m_group, y = mean, fill = m_group)) +
    geom_col(width = 0.75) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = 0.2,
      linewidth = 0.6
    ) +
    facet_wrap(~Gene, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = m_group_colors) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste("MAPK Activation Genes in T/NK cells", label_suffix),
      x = "m_group",
      y = "Mean expression (± SEM)"
    )
  
  ggsave(
    file.path(outdir, paste0("MAPK_genes_TNK_barplot_", label_suffix, ".pdf")),
    p_bar,
    width = 16,
    height = 18
  )
  
  # -------- LINE PLOT --------
  p_line <- ggplot(sum_df, aes(x = m_group, y = mean, group = 1)) +
    geom_line(color = "black", linewidth = 0.7) +
    geom_point(aes(color = m_group), size = 2.5) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = 0.15
    ) +
    facet_wrap(~Gene, scales = "free_y", ncol = 6) +
    scale_color_manual(values = m_group_colors) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste("MAPK Activation Genes in T/NK cells", label_suffix),
      x = "m_group",
      y = "Mean expression"
    )
  
  ggsave(
    file.path(outdir, paste0("MAPK_genes_TNK_lineplot_", label_suffix, ".pdf")),
    p_line,
    width = 16,
    height = 18
  )
}

# --------------------------------------------------
# 1️⃣ Without control
# --------------------------------------------------

make_plots(
  df_long,
  groups_keep = c("MC","MC+MS","MC+MS+MG"),
  label_suffix = "no_control"
)

# --------------------------------------------------
# 2️⃣ Without control AND without MC+MS
# --------------------------------------------------

make_plots(
  df_long,
  groups_keep = c("MC","MC+MS+MG"),
  label_suffix = "no_control_no_MCMS"
)

cat("✅ Both additional figure sets saved./n")













# ==============================================================================
# Pseudobulk MAPK Analysis: MC vs MC+MS+MG (Patient-level Comparison)
# ==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

# 1. Setup and Paths
outdir <- file.path(TMOS_DIR, "Spatial/Heart/MAPK_ReactomePA/MAP_kinase_activation")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 2. Define Genes and Colors
genes <- c(
  "ATF1","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","IKBKB",
  "IRAK1","IRAK2","MEF2A","MEF2C","NFKB1","PPP2CA","PPP2R1A",
  "MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
  "MAP2K1","MAP2K4","SKP1","MAP3K7","TRAF6","CUL1","IKBKG",
  "RIPK2","BTRC","RPS6KA5","MAPKAPK2","NOD1","TAB1","TAB2",
  "FBXW11","NOD2"
)

# Ensure genes exist in the object
genes <- intersect(genes, rownames(tnk_clean))

# Define colors for the specific groups
m_group_colors <- c(
  "MC"        = "#4682B4",
  "MC+MS+MG"  = "#b30000"
)

# 3. Pseudobulk Aggregation
# Extract RNA data and metadata
expr_mat <- GetAssayData(tnk_clean, assay = "RNA", layer = "data")[genes, ]
meta <- tnk_clean@meta.data

# Create long-format pseudobulk dataframe
# We average expression of all cells within a sample_id
pb_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  mutate(sample_id = meta$sample_id,
         m_group = meta$m_group) %>%
  # Filter for only the two groups of interest
  filter(m_group %in% c("MC", "MC+MS+MG")) %>%
  group_by(sample_id, m_group) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(m_group = factor(m_group, levels = c("MC", "MC+MS+MG")))

# 4. Statistical Testing (Wilcoxon Rank-Sum at Sample Level)
stats_df <- pb_df %>%
  group_by(Gene) %>%
  do(w_test = wilcox.test(Expression ~ m_group, data = .)) %>%
  summarise(Gene = Gene, 
            p_value = w_test$p.value) %>%
  arrange(p_value)

# Save stats to Excel
write.xlsx(stats_df, file.path(outdir, "MAPK_Pseudobulk_Wilcoxon_Stats.xlsx"))

# 5. Create "Dots Only" Plot
# We use geom_jitter for patient points and a thin crossbar for the mean
p_dots <- ggplot(pb_df, aes(x = m_group, y = Expression, color = m_group)) +
  # Jittered points (Patient level)
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
  # Horizontal line representing the mean
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.4) +
  facet_wrap(~Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = m_group_colors) +
  # Statistical comparison with Wilcoxon
  stat_compare_means(
    comparisons = list(c("MC", "MC+MS+MG")),
    method = "wilcox.test",
    label = "p.format",
    size = 3,
    label.y.npc = 0.85 # Moves p-value down slightly from the very top
  ) +
  # Styling
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Pseudobulk Analysis: MAPK Activation Pathway",
    subtitle = "Points represent individual patient means | Wilcoxon Rank-Sum Test",
    x = NULL,
    y = "Average Normalized Expression (Per Sample)",
    color = "Group"
  ) +
  # Add breathing room at the top so p-values aren't cut off
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

# 6. Save Figure
pdf(file.path(outdir, "MAPK_Pseudobulk_DotPlots_MC_vs_Triple.pdf"), width = 16, height = 18)
print(p_dots)
dev.off()

cat("✅ Analysis complete. Stats and PDF saved to outdir./n")









# ==============================================================================
# Pseudobulk MAPK Analysis: Significant Genes Only (MC vs MC+MS+MG)
# ==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

# 1. Setup and Paths
outdir <- file.path(TMOS_DIR, "Spatial/Heart/MAPK_ReactomePA/MAP_kinase_activation")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 2. Define Genes and Colors
genes <- c(
  "ATF1","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","IKBKB",
  "IRAK1","IRAK2","MEF2A","MEF2C","NFKB1","PPP2CA","PPP2R1A",
  "MAPK1","MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10",
  "MAP2K1","MAP2K4","SKP1","MAP3K7","TRAF6","CUL1","IKBKG",
  "RIPK2","BTRC","RPS6KA5","MAPKAPK2","NOD1","TAB1","TAB2",
  "FBXW11","NOD2"
)

genes <- intersect(genes, rownames(tnk_clean))
m_group_colors <- c("MC" = "#4682B4", "MC+MS+MG" = "#b30000")

# 3. Pseudobulk Aggregation
expr_mat <- GetAssayData(tnk_clean, assay = "RNA", layer = "data")[genes, ]
meta <- tnk_clean@meta.data

pb_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  mutate(sample_id = meta$sample_id, m_group = meta$m_group) %>%
  filter(m_group %in% c("MC", "MC+MS+MG")) %>%
  group_by(sample_id, m_group) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(m_group = factor(m_group, levels = c("MC", "MC+MS+MG")))

# 4. Statistical Testing (Wilcoxon Rank-Sum)
stats_df <- pb_df %>%
  group_by(Gene) %>%
  do(w_test = wilcox.test(Expression ~ m_group, data = .)) %>%
  summarise(Gene = Gene, p_value = w_test$p.value) %>%
  arrange(p_value)

# Identify Significant Genes (p < 0.05)
sig_genes <- stats_df %>% filter(p_value < 0.05) %>% pull(Gene)

# Save stats to Excel
write.xlsx(stats_df, file.path(outdir, "MAPK_Pseudobulk_Wilcoxon_Stats.xlsx"))

# 5. Plotting Function (To maintain consistency)
make_dot_plot <- function(data_to_plot, title_suffix, num_cols) {
  ggplot(data_to_plot, aes(x = m_group, y = Expression, color = m_group)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.5) +
    facet_wrap(~Gene, scales = "free_y", ncol = num_cols) +
    scale_color_manual(values = m_group_colors) +
    stat_compare_means(
      comparisons = list(c("MC", "MC+MS+MG")),
      method = "wilcox.test",
      label = "p.format",
      size = 3.5,
      label.y.npc = 0.85
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = paste("MAPK Pathway:", title_suffix),
         subtitle = "Pseudobulk (Sample-level) | Wilcoxon Test",
         x = NULL, y = "Avg Expression")
}

# 6. Save BOTH Figures

# Figure A: All Genes
pdf(file.path(outdir, "MAPK_Pseudobulk_ALL_GENES.pdf"), width = 13, height = 18)
print(make_dot_plot(pb_df, "All Target Genes", 6))
dev.off()

# Figure B: ONLY SIGNIFICANT GENES
if(length(sig_genes) > 0) {
  pb_sig <- pb_df %>% filter(Gene %in% sig_genes)
  
  # Adjust dimensions based on number of significant genes
  plot_width <- ifelse(length(sig_genes) < 4, 8, 10)
  plot_height <- ceiling(length(sig_genes) / 4) * 4
  
  pdf(file.path(outdir, "MAPK_Pseudobulk_SIGNIFICANT_ONLY.pdf"), width = plot_width, height = plot_height)
  print(make_dot_plot(pb_sig, "Statistically Significant Genes (p < 0.05)", 4))
  dev.off()
  
  cat("✅ Significant genes found:", paste(sig_genes, collapse=", "), "/n")
} else {
  cat("❌ No genes reached significance at p < 0.05/n")
}

cat("✅ All outputs saved successfully to outdir./n")

## end of Part 1: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 2 - Fig 2E; Fig S2A-S2C: donor-level pseudobulk limma-voom and NF-kB/MAPK core program
##############################################################################

# ============================================================================
# Donor-level pseudobulk limma-voom (TMOS vs MC), Reactome ORA, 14-gene NF-kB/MAPK core program (site-control referenced)
#
# Figure panels : Fig 2E; Fig S2A, S2B, S2C
# ============================================================================
library(patchwork)   # used below but not loaded in the original script

library(Seurat)
library(dplyr)
library(tibble)
library(edgeR)
library(limma)
library(writexl)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(stringr)
library(scales)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cut <- 0.05
lfc_cut  <- 0.585
col_mc   <- "#4682B4"
col_tmos <- "#b30000"

## ---- Data ----
heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"

tnk <- subset(heart, subset = major_cell_type == "T_NK_lymphocytes" &
                m_group %in% c("MC", "MC+MS+MG"))

## ---- Pseudobulk per sample ----
pb <- AggregateExpression(tnk, assays = "RNA", group.by = "sample_id",
                          slot = "counts", return.seurat = FALSE)$RNA

meta <- tnk@meta.data %>%
  dplyr::count(sample_id, m_group, institution, name = "n_cells") %>%
  dplyr::mutate(sample_id = as.character(sample_id),
                pb_id     = gsub("_", "-", sample_id),
                pb_id     = ifelse(grepl("^[0-9]", pb_id), paste0("g", pb_id), pb_id),
                grp       = factor(ifelse(m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS"))) %>%
  dplyr::filter(n_cells >= 10)
rownames(meta) <- meta$pb_id

stopifnot(all(meta$pb_id %in% colnames(pb)))
table(meta$grp, meta$institution)

counts <- round(as.matrix(pb[, meta$pb_id]))

## ---- limma-voom ----
dge <- DGEList(counts, samples = meta)
dge <- dge[filterByExpr(dge, group = meta$grp), , keep.lib.sizes = FALSE]
dge <- normLibSizes(dge)

design <- model.matrix(~ grp, data = meta)

pdf(file.path(out_dir, "TNK_voom_mean_variance_TMOS_vs_MC.pdf"), width = 6, height = 5)
v <- voom(dge, design, plot = TRUE)
dev.off()

fit <- eBayes(lmFit(v, design))

res <- topTable(fit, coef = "grpTMOS", number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

sig <- res %>%
  dplyr::filter(padj < padj_cut, abs(log2FoldChange) >= lfc_cut) %>%
  dplyr::mutate(enriched_in = ifelse(log2FoldChange > 0, "TMOS", "MC"))

table(sig$enriched_in)

write_xlsx(list(Sig     = sig,
                Up_TMOS = dplyr::filter(sig, enriched_in == "TMOS"),
                Up_MC   = dplyr::filter(sig, enriched_in == "MC"),
                All     = res,
                Samples = meta),
           file.path(out_dir, "TNK_pseudobulk_limma_voom_TMOS_vs_MC.xlsx"))

## ---- Volcano ----
de_plot <- res %>%
  dplyr::mutate(neglog10_padj = -log10(padj),
                direction = dplyr::case_when(
                  padj < padj_cut & log2FoldChange >  lfc_cut ~ "TMOS",
                  padj < padj_cut & log2FoldChange < -lfc_cut ~ "MC",
                  TRUE ~ "NS"))
de_plot$direction <- factor(de_plot$direction, levels = c("MC", "TMOS", "NS"))

vol_colors <- c("MC" = col_mc, "TMOS" = col_tmos, "NS" = "grey75")

top_labels <- de_plot %>%
  dplyr::filter(direction %in% c("MC", "TMOS")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

p_volcano <- ggplot(de_plot, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(aes(color = direction), size = 1.8, alpha = 0.8) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dotted") +
  geom_text_repel(data = top_labels, aes(label = gene), size = 3.2, max.overlaps = Inf) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
  labs(title = "TNK pseudobulk DE (limma-voom): TMOS vs MC",
       subtitle = "Sample-level; padj < 0.05, |log2FC| ≥ 0.585",
       x = "log2FC (TMOS − MC)",
       y = expression(-log[10]("adjusted p-value")),
       color = "Enriched in")

ggsave(file.path(out_dir, "TNK_limma_Volcano_TMOS_vs_MC.pdf"),
       plot = p_volcano, width = 7.5, height = 6, device = cairo_pdf)

## ---- Reactome ORA (ReactomePA, as in original Methods) ----
to_entrez <- function(g) {
  unique(bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
}

run_ora <- function(genes) {
  if (length(genes) < 10) return(NULL)
  enrichPathway(gene = to_entrez(genes), organism = "human",
                pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
}

up_tmos <- sig$gene[sig$enriched_in == "TMOS"]
up_mc   <- sig$gene[sig$enriched_in == "MC"]

ora_tmos <- run_ora(up_tmos)
ora_mc   <- run_ora(up_mc)

add_fe <- function(df) {
  if (nrow(df) == 0 || "FoldEnrichment" %in% colnames(df)) return(df)
  gr <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  br <- sapply(strsplit(df$BgRatio,   "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df$FoldEnrichment <- gr / br
  df
}

df_tmos <- add_fe(if (!is.null(ora_tmos)) as.data.frame(ora_tmos) else data.frame())
df_mc   <- add_fe(if (!is.null(ora_mc))   as.data.frame(ora_mc)   else data.frame())

nrow(df_tmos); nrow(df_mc)

write_xlsx(list(Reactome_Up_TMOS = df_tmos,
                Reactome_Up_MC   = df_mc),
           file.path(out_dir, "TNK_Reactome_ORA_limma_TMOS_vs_MC.xlsx"))

## ---- Reactome barplots (gradient style, top 15) ----
prep_bar_df <- function(df, wrap_width = 55, top_n = 15) {
  df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(neglog10_padj = -log10(p.adjust),
                  Pathway_wrap  = str_wrap(Description, width = wrap_width),
                  Pathway_wrap  = factor(Pathway_wrap, levels = rev(Pathway_wrap)))
}

make_bar_gradient <- function(df, title_text, base_color) {
  if (nrow(df) == 0) return(NULL)
  df2 <- prep_bar_df(df)
  ggplot(df2, aes(x = neglog10_padj, y = Pathway_wrap, fill = FoldEnrichment)) +
    geom_col(width = 0.86) +
    geom_text(aes(label = Pathway_wrap), x = 0.02, hjust = 0, size = 4.2) +
    scale_fill_gradient(low = alpha(base_color, 0.3), high = base_color,
                        name = "Fold Enrichment") +
    theme_classic(base_size = 14) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line    = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(title = title_text,
         x = expression(-log[10]("adjusted p-value")),
         y = NULL)
}

p_tmos <- make_bar_gradient(df_tmos, "Reactome ORA — Up in TMOS", col_tmos)
p_mc   <- make_bar_gradient(df_mc,   "Reactome ORA — Up in MC",   col_mc)

if (!is.null(p_tmos)) ggsave(file.path(out_dir, "Reactome_ORA_Top15_Up_in_TMOS_gradient.pdf"),
                             plot = p_tmos, width = 6, height = 6.5, device = cairo_pdf)
if (!is.null(p_mc))   ggsave(file.path(out_dir, "Reactome_ORA_Top15_Up_in_MC_gradient.pdf"),
                             plot = p_mc, width = 6, height = 6.5, device = cairo_pdf)

p_volcano; p_tmos; p_mc










terms <- c("MAP kinase activation", "Interleukin-17 signaling")

## ---- Get both terms from TMOS ORA (unfiltered if needed) ----
ora_tmos_all <- enrichPathway(gene = to_entrez(up_tmos), organism = "human",
                              pvalueCutoff = 1, qvalueCutoff = 1,
                              pAdjustMethod = "BH", readable = TRUE)
tm_all <- as.data.frame(ora_tmos_all)

term_tab <- tm_all %>%
  dplyr::filter(Description %in% terms) %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
term_tab

de_genes <- setNames(lapply(terms, function(t)
  strsplit(term_tab$geneID[term_tab$Description == t], "/")[[1]]), terms)

## ---- Overlap of DE genes driving each term ----
shared   <- intersect(de_genes[[1]], de_genes[[2]])
only_mapk <- setdiff(de_genes[["MAP kinase activation"]], de_genes[["Interleukin-17 signaling"]])
only_il17 <- setdiff(de_genes[["Interleukin-17 signaling"]], de_genes[["MAP kinase activation"]])

list(shared = shared, only_MAPK_activation = only_mapk, only_IL17 = only_il17)

## ---- Overlap of the full Reactome gene sets ----
pw_genes <- lapply(term_tab$ID, function(id) {
  eg <- ReactomePA:::get_Reactome_DATA("human")$PATHID2EXTID[[id]]
  unique(bitr(eg, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL)
})
names(pw_genes) <- term_tab$Description

full_overlap <- data.frame(
  set = c("MAP kinase activation (full)", "Interleukin-17 signaling (full)", "Shared (full sets)",
          "MAPK activation genes contained in IL-17 set (%)"),
  n = c(length(pw_genes[["MAP kinase activation"]]),
        length(pw_genes[["Interleukin-17 signaling"]]),
        length(intersect(pw_genes[[1]], pw_genes[[2]])),
        round(100 * mean(pw_genes[["MAP kinase activation"]] %in% pw_genes[["Interleukin-17 signaling"]]), 1)))
full_overlap


## ---- Per-gene table with DE stats ----
gene_tab <- res %>%
  dplyr::filter(gene %in% union(de_genes[[1]], de_genes[[2]])) %>%
  dplyr::mutate(in_MAPK_activation = gene %in% de_genes[["MAP kinase activation"]],
                in_IL17_signaling  = gene %in% de_genes[["Interleukin-17 signaling"]]) %>%
  dplyr::select(gene, in_MAPK_activation, in_IL17_signaling, log2FoldChange, padj) %>%
  dplyr::arrange(dplyr::desc(in_MAPK_activation & in_IL17_signaling), padj)
gene_tab

write_xlsx(list(Terms        = term_tab,
                DE_gene_overlap = gene_tab,
                Full_set_overlap = full_overlap),
           file.path(out_dir, "TNK_IL17_vs_MAPKactivation_overlap_TMOS.xlsx"))


setdiff(pw_genes[["Interleukin-17 signaling"]], pw_genes[["MAP kinase activation"]])









nm_core <- c("TAB2","MAPK8","UBE2V1","NFKB1","MAP2K4","BTRC","CREB1",
             "PPP2CB","RPS6KA5","MEF2A","RPS6KA3","FBXW11","MAP2K1","DUSP4")

## ---- Per-sample expression (voom log2-CPM) ----
expr_df <- as.data.frame(t(v$E[nm_core, meta$pb_id])) %>%
  tibble::rownames_to_column("pb_id") %>%
  dplyr::left_join(meta[, c("pb_id", "sample_id", "grp")], by = "pb_id") %>%
  tidyr::pivot_longer(dplyr::all_of(nm_core), names_to = "gene", values_to = "log2CPM")

## ---- limma stats, ordered by log2FC (strongest first) ----
stat_df <- res %>%
  dplyr::filter(gene %in% nm_core) %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::mutate(label = sprintf("log2FC = %.2f\npadj = %.3f", log2FoldChange, padj))

gene_order   <- stat_df$gene
expr_df$gene <- factor(expr_df$gene, levels = gene_order)
stat_df$gene <- factor(stat_df$gene, levels = gene_order)

lab_pos <- expr_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(log2CPM) + 0.12 * diff(range(log2CPM)), .groups = "drop")
stat_df <- dplyr::left_join(stat_df, lab_pos, by = "gene")

## ---- Plot ----
p_core <- ggplot(expr_df, aes(x = grp, y = log2CPM, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.6) +
  geom_point(aes(colour = grp), size = 2.2,
             position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = stat_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 2.9, lineheight = 0.9) +
  facet_wrap(~ gene, scales = "free_y", ncol = 7) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Expression (voom log2-CPM)",
       title = "NF-κB/MAPK core genes in heart T/NK cells (sample level)",
       subtitle = "Pseudobulk limma-voom, TMOS vs MC") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank())

p_core
ggsave(file.path(out_dir, "TNK_NFkB_MAPK_core_genes_boxplots_TMOS_vs_MC.pdf"),
       p_core, width = 14, height = 4.7, device = cairo_pdf)









library(Seurat)
library(dplyr)
library(ggplot2)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_mc   <- "#4682B4"
col_tmos <- "#b30000"

nm_core <- c("TAB2","MAPK8","UBE2V1","NFKB1","MAP2K4","BTRC","CREB1",
             "PPP2CB","RPS6KA5","MEF2A","RPS6KA3","FBXW11","MAP2K1","DUSP4")

## ---- Data ----
heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"

tnk  <- subset(heart, subset = major_cell_type == "T_NK_lymphocytes" &
                 m_group %in% c("MC", "MC+MS+MG"))
tnk3 <- subset(heart, subset = major_cell_type == "T_NK_lymphocytes" &
                 m_group %in% c("control", "MC", "MC+MS+MG"))
rm(heart); gc()

## ---- Helper: two-group plot with bracket ----
plot_two <- function(df, yvar, ylab, pval, title_text) {
  rng   <- diff(range(df[[yvar]]))
  y_top <- max(df[[yvar]]) + 0.12 * rng
  tick  <- 0.03 * rng
  ggplot(df, aes(x = grp, y = .data[[yvar]], fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
    geom_point(aes(colour = grp), size = 3,
               position = position_jitter(width = 0.12, seed = 1)) +
    annotate("segment", x = 1, xend = 2, y = y_top, yend = y_top) +
    annotate("segment", x = 1, xend = 1, y = y_top, yend = y_top - tick) +
    annotate("segment", x = 2, xend = 2, y = y_top, yend = y_top - tick) +
    annotate("text", x = 1.5, y = y_top + 1.5 * tick,
             label = sprintf("p = %.3f", pval), size = 4) +
    scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
    scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
    labs(x = NULL, y = ylab, title = title_text) +
    theme_classic(base_size = 13) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
}

## ---- Non-adjusted (MC + TMOS) ----
tnk <- NormalizeData(tnk)
tnk <- AddModuleScore(tnk, features = list(nm_core), name = "NMcore_", seed = 1)

core_nonadj <- tnk@meta.data %>%
  dplyr::mutate(grp = factor(ifelse(m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS"))) %>%
  dplyr::group_by(sample_id, grp) %>%
  dplyr::summarise(n_cells = dplyr::n(), score = mean(NMcore_1), .groups = "drop") %>%
  dplyr::filter(n_cells >= 10)

w_nonadj <- wilcox.test(score ~ grp, data = core_nonadj, exact = FALSE)
w_nonadj

p_nonadj <- plot_two(core_nonadj, "score", "NF-κB/MAPK core score (mean per sample)",
                     w_nonadj$p.value, "NF-κB/MAPK core (non-adjusted)")
p_nonadj
ggsave(file.path(out_dir, "TNK_NFkB_MAPK_core_score_nonAdjusted_TMOS_vs_MC.pdf"),
       p_nonadj, width = 3.5, height = 3.5, device = cairo_pdf)

## ---- Inst-adjusted (control-referenced) ----
tnk3 <- NormalizeData(tnk3)
tnk3 <- AddModuleScore(tnk3, features = list(nm_core), name = "NMcore_", seed = 1)

core_df <- tnk3@meta.data %>%
  dplyr::mutate(grp = factor(dplyr::recode(as.character(m_group), "MC+MS+MG" = "TMOS"),
                             levels = c("control", "MC", "TMOS"))) %>%
  dplyr::group_by(sample_id, grp, institution) %>%
  dplyr::summarise(n_cells = dplyr::n(), score = mean(NMcore_1), .groups = "drop") %>%
  dplyr::filter(n_cells >= 10)

ctrl_med <- core_df %>%
  dplyr::filter(grp == "control") %>%
  dplyr::group_by(institution) %>%
  dplyr::summarise(n_ctrl = dplyr::n(), ctrl = median(score), .groups = "drop")
ctrl_med

core_df <- core_df %>%
  dplyr::left_join(ctrl_med, by = "institution") %>%
  dplyr::mutate(score_adj = score - ctrl)

core_instadj <- dplyr::filter(core_df, grp %in% c("MC", "TMOS")) %>% droplevels()

w_instadj <- wilcox.test(score_adj ~ grp, data = core_instadj, exact = FALSE)
delta     <- with(core_instadj, mean(score_adj[grp == "TMOS"]) - mean(score_adj[grp == "MC"]))

set.seed(1)
ctrl_df <- dplyr::filter(core_df, grp == "control")
boot_delta <- replicate(5000, {
  off <- ctrl_df %>% dplyr::group_by(institution) %>%
    dplyr::summarise(o = median(sample(score, replace = TRUE)), .groups = "drop")
  d <- core_instadj %>% dplyr::select(institution, grp, score) %>%
    dplyr::left_join(off, by = "institution") %>% dplyr::mutate(a = score - o)
  mean(d$a[d$grp == "TMOS"]) - mean(d$a[d$grp == "MC"])
})
ci <- quantile(boot_delta, c(0.025, 0.975))

w_instadj
round(c(delta = delta, ci), 3)
mean(boot_delta > 0)

p_instadj <- plot_two(core_instadj, "score_adj",
                      "NF-κB/MAPK core score\n(relative to site controls)",
                      w_instadj$p.value, "NF-κB/MAPK core (inst-adjusted)")
p_instadj
ggsave(file.path(out_dir, "TNK_NFkB_MAPK_core_score_instAdjusted_TMOS_vs_MC.pdf"),
       p_instadj, width = 3.8, height = 3.5, device = cairo_pdf)

## ---- Save ----
write_xlsx(list(Scores_nonAdjusted  = core_nonadj,
                Scores_instAdjusted = core_df,
                Ctrl_reference      = ctrl_med,
                Stats = data.frame(
                  analysis = c("nonAdjusted", "instAdjusted"),
                  wilcox_p = c(w_nonadj$p.value, w_instadj$p.value),
                  delta    = c(with(core_nonadj, mean(score[grp == "TMOS"]) - mean(score[grp == "MC"])), delta),
                  ci_low   = c(NA, ci[1]),
                  ci_high  = c(NA, ci[2]),
                  pct_boot_TMOS_gt_MC = c(NA, mean(boot_delta > 0)))),
           file.path(out_dir, "TNK_NFkB_MAPK_core_score_TMOS_vs_MC.xlsx"))










## Response-letter figure (R2.4), not a manuscript panel. It combines plot
## objects (fA, fB, fE, bt) that were created interactively, so it is skipped
## unless they exist in the session.
if (all(vapply(c("fA", "fB", "fE", "bt"), exists, logical(1)))) {

## ---- C: core genes per donor, all genes in one panel ----
g_stat <- res %>% dplyr::filter(gene %in% nm_core) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  dplyr::mutate(star = dplyr::case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"))

z <- t(scale(t(v$E[nm_core, meta$pb_id])))
gx_df <- as.data.frame(z) %>% tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "pb_id", values_to = "z") %>%
  dplyr::left_join(meta[, c("pb_id", "grp")], by = "pb_id") %>%
  dplyr::mutate(gene = factor(gene, levels = g_stat$gene))

star_df <- gx_df %>% dplyr::group_by(gene) %>% dplyr::summarise(y = max(z) + 0.35, .groups = "drop") %>%
  dplyr::left_join(g_stat[, c("gene", "star")], by = "gene") %>%
  dplyr::mutate(gene = factor(gene, levels = g_stat$gene))

pd <- position_dodge(width = 0.75)
fC <- ggplot(gx_df, aes(gene, z, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.65, position = pd) +
  geom_point(aes(colour = grp), size = 1.4,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75, seed = 1)) +
  geom_text(data = star_df, aes(gene, y, label = star), inherit.aes = FALSE, size = 3.5) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos),
                    labels = paste0(c("MC", "TMOS"), " (n = ", table(meta$grp)[c("MC", "TMOS")], ")")) +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(title = paste0("Core genes of the ", prog, " (per donor)"),
       x = NULL, y = "Expression\n(z-scored voom log2-CPM)", fill = NULL,
       caption = "limma padj: * < 0.05, ** < 0.01, *** < 0.001") +
  bt + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
             legend.position = "top")

## ---- Assemble ----
fig <- (fA | fB) / (fC | fE) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A",
                  title = paste0("Donor-level analyses of the ", prog, " in cardiac T/NK cells"),
                  caption = "Each donor contributes one pseudobulk profile/score (≥10 T/NK cells). Two-sided Wilcoxon rank-sum test (D). DE: padj < 0.05, |log2FC| ≥ 0.585.",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))) &
  theme(plot.tag = element_text(face = "bold", size = 15))
fig[[2]] <- fig[[2]] + plot_layout(widths = c(2.2, 1))

ggsave(file.path(out_dir, "R2.4_donor_level_combined_figure.pdf"), fig, width = 14, height = 11.5, device = cairo_pdf)
ggsave(file.path(out_dir, "R2.4_donor_level_combined_figure.png"), fig, width = 14, height = 11.5, dpi = 300, bg = "white")

}  # end R2.4 response figure

## end of Part 2: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 3 - Fig 2F, 2G, 2H: T/NK Monaco subtypes, gd T expansion, IL-17 program by subtype (also makes the state tables and paired CellChat objects used for Fig 1F, S1G, S1H)
##############################################################################

# ============================================================================
# Monaco T/NK subtypes, gd T expansion, program delta by subtype, paired heart-skeletal states, paired CellChat
#
# Figure panels : Fig 1F; Fig 2F, 2G, 2H; Fig S1G, S1H
# ============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(SingleR)
library(celldex)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

monaco <- celldex::MonacoImmuneData()

## ---- Helper: subset T/NK, Monaco fine labels (as is), save ----
label_tnk <- function(rds, tissue, ct_col, ct_value, groups = NULL) {
  obj <- readRDS(file.path(base_dir, rds))
  DefaultAssay(obj) <- "RNA"
  
  keep <- obj@meta.data[[ct_col]] == ct_value
  if (!is.null(groups)) keep <- keep & obj$m_group %in% groups
  tnk <- subset(obj, cells = colnames(obj)[keep])
  rm(obj); gc()
  
  tnk <- NormalizeData(tnk)
  sr  <- SingleR(test = GetAssayData(tnk, layer = "data"),
                 ref = monaco, labels = monaco$label.fine)
  tnk$monaco_fine <- sr$pruned.labels
  tnk$tissue      <- tissue
  if (!"patient_id" %in% colnames(tnk@meta.data))
    tnk$patient_id <- sub("_.*", "", tnk$sample_id)
  
  meta <- tnk@meta.data[, intersect(c("sample_id", "patient_id", "m_group", "institution",
                                      "tissue", "monaco_fine"), colnames(tnk@meta.data))]
  saveRDS(meta, file.path(out_dir, paste0(tissue, "_TNK_monaco_fine_meta.rds")))
  saveRDS(tnk,  file.path(out_dir, paste0(tissue, "_TNK_monaco.rds")))
  rm(tnk, sr); gc()
  meta
}

## ---- Run ----
meta_heart <- label_tnk("heart_3232026.rds", "heart",
                        ct_col = "major_cell_type", ct_value = "T_NK_lymphocytes",
                        groups = c("MC", "MC+MS+MG"))
meta_skel  <- label_tnk("sc_skeletal_322026.rds", "skeletal",
                        ct_col = "major_match", ct_value = "T/NK Lymphocytes")

## ---- Checks ----
table(meta_heart$sample_id, meta_heart$m_group)
table(meta_skel$sample_id,  meta_skel$m_group)


paired <- intersect(as.character(meta_heart$patient_id), as.character(meta_skel$patient_id))

## paired heart–skeletal patients
meta_heart$patient_id <- as.character(meta_heart$patient_id)
meta_skel$patient_id  <- as.character(meta_skel$patient_id)

dplyr::bind_rows(meta_heart, meta_skel) %>%
  dplyr::filter(patient_id %in% paired) %>%
  dplyr::count(patient_id, tissue, m_group)

## ---- Coverage: donors with >=10 cells per label ----
coverage <- function(meta) {
  meta %>%
    dplyr::filter(!is.na(monaco_fine)) %>%
    dplyr::count(sample_id, m_group, monaco_fine, name = "n_cells") %>%
    dplyr::filter(n_cells >= 10) %>%
    dplyr::count(monaco_fine, m_group, name = "n_donors") %>%
    tidyr::pivot_wider(names_from = m_group, values_from = n_donors, values_fill = 0)
}
cov_heart <- coverage(meta_heart); cov_heart
cov_skel  <- coverage(meta_skel);  cov_skel

write_xlsx(list(Heart_label_counts       = as.data.frame(table(meta_heart$monaco_fine, useNA = "ifany")),
                Skeletal_label_counts    = as.data.frame(table(meta_skel$monaco_fine,  useNA = "ifany")),
                Heart_coverage_ge10      = cov_heart,
                Skeletal_coverage_ge10   = cov_skel,
                Heart_cells_per_donor    = dplyr::count(meta_heart, sample_id, m_group, monaco_fine),
                Skeletal_cells_per_donor = dplyr::count(meta_skel,  sample_id, m_group, monaco_fine),
                Paired_patients          = data.frame(patient_id = paired)),
           file.path(out_dir, "TNK_monaco_fine_coverage_heart_skeletal.xlsx"))





library(ggplot2)

col_mc   <- "#4682B4"
col_tmos <- "#b30000"

tnk_labels <- c("Effector memory CD8 T cells", "Terminal effector CD8 T cells", "Central memory CD8 T cells",
                "Naive CD8 T cells", "Vd2 gd T cells", "Non-Vd2 gd T cells", "MAIT cells",
                "Natural killer cells", "T regulatory cells", "Th1 cells", "Th1/Th17 cells", "Th17 cells",
                "Th2 cells", "Follicular helper T cells", "Naive CD4 T cells", "Terminal effector CD4 T cells")

sigs <- list(
  IL17_program      = c("TAB2","MAPK8","UBE2V1","NFKB1","MAP2K4","BTRC","CREB1",
                        "PPP2CB","RPS6KA5","MEF2A","RPS6KA3","FBXW11","MAP2K1","DUSP4"),
  Cytotoxicity      = c("GZMA","GZMB","GZMH","PRF1","GNLY","NKG7"),
  Exhaustion        = c("PDCD1","HAVCR2","LAG3","TIGIT","TOX","CTLA4"),
  Tissue_residency  = c("ITGAE","CXCR6","ITGA1","CD69","ZNF683"),
  Effector_terminal = c("KLRG1","CX3CR1","TBX21","ZEB2","FGFBP2")
)

score_obj <- function(obj) {
  obj <- NormalizeData(obj)
  for (s in names(sigs)) {
    g <- intersect(sigs[[s]], rownames(obj))
    obj <- AddModuleScore(obj, features = list(g), name = paste0(s, "_"), seed = 1)
  }
  obj
}

sum_scores <- function(md, by) {
  md %>%
    tidyr::pivot_longer(dplyr::all_of(paste0(names(sigs), "_1")), names_to = "signature", values_to = "score") %>%
    dplyr::mutate(signature = sub("_1$", "", signature)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by)), monaco_fine, signature) %>%
    dplyr::summarise(n_cells = dplyr::n(), mean_score = mean(score), .groups = "drop")
}

## ---- 1. Heart: TMOS vs MC ----
heart <- readRDS(file.path(out_dir, "heart_TNK_monaco.rds"))
heart <- subset(heart, cells = colnames(heart)[!is.na(heart$monaco_fine) & heart$monaco_fine %in% tnk_labels])
heart$grp <- factor(ifelse(heart$m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS"))
heart <- score_obj(heart)
md_h  <- heart@meta.data
rm(heart); gc()

heart_group <- sum_scores(md_h, "grp")
heart_donor <- sum_scores(md_h, c("sample_id", "grp")) %>% dplyr::filter(n_cells >= 10)

## ---- 2. Paired heart vs skeletal: scored together ----
h <- readRDS(file.path(out_dir, "heart_TNK_monaco.rds"))
h <- subset(h, cells = colnames(h)[as.character(h$patient_id) %in% paired])
s <- readRDS(file.path(out_dir, "skeletal_TNK_monaco.rds"))
s <- subset(s, cells = colnames(s)[as.character(s$patient_id) %in% paired])
h$patient_id <- as.character(h$patient_id)
s$patient_id <- as.character(s$patient_id)

pr <- merge(h, s)
rm(h, s); gc()
pr <- JoinLayers(pr)
pr <- subset(pr, cells = colnames(pr)[!is.na(pr$monaco_fine) & pr$monaco_fine %in% tnk_labels])
pr <- score_obj(pr)
md_p <- pr@meta.data
rm(pr); gc()

paired_scores <- sum_scores(md_p, c("patient_id", "tissue")) %>% dplyr::filter(n_cells >= 10)
paired_comp   <- md_p %>% dplyr::count(patient_id, tissue, monaco_fine) %>%
  dplyr::group_by(patient_id, tissue) %>% dplyr::mutate(prop = 100 * n / sum(n)) %>% dplyr::ungroup()

## ---- Save ----
saveRDS(md_h, file.path(out_dir, "heart_TNK_state_scores_meta.rds"))
saveRDS(md_p, file.path(out_dir, "paired_TNK_state_scores_meta.rds"))
write_xlsx(list(Heart_by_group     = heart_group,
                Heart_by_donor     = heart_donor,
                Paired_scores      = paired_scores,
                Paired_composition = paired_comp,
                Signatures = data.frame(signature = rep(names(sigs), lengths(sigs)), gene = unlist(sigs))),
           file.path(out_dir, "TNK_state_scores_heart_paired.xlsx"))



## =====================================================================
## Composition (heart, donor level) — unadjusted Wilcoxon, plot significant only
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(patchwork)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
col_mc   <- "#4682B4"
col_tmos <- "#b30000"

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold", size = 11),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

meta_heart <- readRDS(file.path(out_dir, "heart_TNK_monaco_fine_meta.rds"))

tnk_labels <- c("Effector memory CD8 T cells", "Terminal effector CD8 T cells", "Central memory CD8 T cells",
                "Naive CD8 T cells", "Vd2 gd T cells", "Non-Vd2 gd T cells", "MAIT cells",
                "Natural killer cells", "T regulatory cells", "Th1 cells", "Th1/Th17 cells", "Th17 cells",
                "Th2 cells", "Follicular helper T cells", "Naive CD4 T cells", "Terminal effector CD4 T cells")

comp <- meta_heart %>%
  dplyr::filter(!is.na(monaco_fine), monaco_fine %in% tnk_labels) %>%
  dplyr::mutate(grp = factor(ifelse(m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS"))) %>%
  dplyr::group_by(sample_id) %>% dplyr::filter(dplyr::n() >= 10) %>% dplyr::ungroup() %>%
  dplyr::count(sample_id, grp, monaco_fine, name = "n") %>%
  tidyr::complete(tidyr::nesting(sample_id, grp), monaco_fine, fill = list(n = 0)) %>%
  dplyr::group_by(sample_id) %>% dplyr::mutate(prop = 100 * n / sum(n)) %>% dplyr::ungroup()

table(dplyr::distinct(comp, sample_id, grp)$grp)

comp_stats <- comp %>%
  dplyr::group_by(monaco_fine) %>%
  dplyr::summarise(median_MC   = median(prop[grp == "MC"]),
                   median_TMOS = median(prop[grp == "TMOS"]),
                   p = wilcox.test(prop ~ grp, exact = FALSE)$p.value, .groups = "drop") %>%
  dplyr::arrange(p)
comp_stats

sig_labels <- comp_stats$monaco_fine[comp_stats$p < 0.05]
sig_labels

comp_sig <- comp %>% dplyr::filter(monaco_fine %in% sig_labels) %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = sig_labels))
p_lab <- comp_stats %>% dplyr::filter(monaco_fine %in% sig_labels) %>%
  dplyr::left_join(comp_sig %>% dplyr::group_by(monaco_fine) %>%
                     dplyr::summarise(y = max(prop) * 1.12, .groups = "drop"), by = "monaco_fine") %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = sig_labels))

p_comp <- ggplot(comp_sig, aes(grp, prop, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2.2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_segment(data = p_lab, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = p_lab, aes(x = 1.5, y = y * 1.06, label = sprintf("p = %.3f", p)),
            inherit.aes = FALSE, size = 3.3) +
  facet_wrap(~ monaco_fine, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "% of T/NK cells per donor",
       title = "Heart T/NK composition (cell types with p < 0.05)") +
  bt
p_comp
ggsave(file.path(out_dir, "Heart_TNK_composition_sig_TMOS_vs_MC.pdf"), p_comp, width = 5.5, height = 4, device = cairo_pdf)

write_xlsx(list(Composition = comp, Stats = comp_stats),
           file.path(out_dir, "Heart_TNK_composition_TMOS_vs_MC.xlsx"))


## =====================================================================
## IL-17 program by cell type (heart): TMOS signal and delta vs MC
## =====================================================================
## load saved scores if the scoring block above was not re-run in this session
if (!exists("heart_group")) {
  sc_file <- file.path(out_dir, "TNK_state_scores_heart_paired.xlsx")
  heart_group   <- readxl::read_excel(sc_file, sheet = "Heart_by_group")
  heart_donor   <- readxl::read_excel(sc_file, sheet = "Heart_by_donor")
  paired_scores <- readxl::read_excel(sc_file, sheet = "Paired_scores")
}
heart_group$grp <- factor(heart_group$grp, levels = c("MC", "TMOS"))
heart_donor$grp <- factor(heart_donor$grp, levels = c("MC", "TMOS"))

## pooled (labels with >=10 cells in both groups)
prog_rank <- heart_group %>%
  dplyr::filter(signature == "IL17_program") %>%
  dplyr::select(grp, monaco_fine, n_cells, mean_score) %>%
  tidyr::pivot_wider(names_from = grp, values_from = c(mean_score, n_cells)) %>%
  dplyr::filter(n_cells_MC >= 10, n_cells_TMOS >= 10) %>%
  dplyr::mutate(delta = mean_score_TMOS - mean_score_MC)

## donor level (donors with >=10 cells of that label)
prog_donor <- heart_donor %>%
  dplyr::filter(signature == "IL17_program", monaco_fine %in% prog_rank$monaco_fine)
prog_donor_delta <- prog_donor %>%
  dplyr::group_by(monaco_fine, grp) %>%
  dplyr::summarise(mean_donor = mean(mean_score), n_donors = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = grp, values_from = c(mean_donor, n_donors)) %>%
  dplyr::mutate(delta_donor = mean_donor_TMOS - mean_donor_MC)

prog_rank <- prog_rank %>% dplyr::left_join(prog_donor_delta, by = "monaco_fine") %>%
  dplyr::arrange(dplyr::desc(delta))
as.data.frame(prog_rank)

ct_order <- rev(prog_rank$monaco_fine)          # largest delta on top
top3     <- prog_rank$monaco_fine[1:3]

## Panel: delta vs MC (lollipop), colour = TMOS score
p_delta <- prog_rank %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = ct_order),
                top = monaco_fine %in% top3) %>%
  ggplot(aes(delta, monaco_fine)) +
  geom_segment(aes(x = 0, xend = delta, yend = monaco_fine), colour = "grey60") +
  geom_point(aes(fill = mean_score_TMOS, stroke = ifelse(top, 1.2, 0.4)), shape = 21, size = 4.5, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", delta)), hjust = -0.6, size = 3) +
  scale_fill_gradient(low = "white", high = col_tmos, name = "Score in\nTMOS") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Δ program score (TMOS − MC)", y = NULL,
       title = "IL-17 signaling–related program by T/NK cell type") +
  bt

## Panel: per-donor scores per cell type
p_donor <- prog_donor %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = rev(ct_order))) %>%
  ggplot(aes(monaco_fine, mean_score, colour = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(size = n_cells), alpha = 0.85,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), name = NULL) +
  scale_size_continuous(range = c(1.5, 5), name = "Cells (n)") +
  labs(x = NULL, y = "Program score\n(mean per donor)",
       title = "Program score per donor and cell type (donors with ≥10 cells of that type)") +
  bt + theme(axis.text.x = element_text(angle = 40, hjust = 1))

p_delta; p_donor
ggsave(file.path(out_dir, "Heart_IL17program_delta_by_celltype.pdf"), p_delta, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Heart_IL17program_per_donor_by_celltype.pdf"), p_donor, width = 9, height = 4.5, device = cairo_pdf)


## =====================================================================
## Paired heart vs skeletal muscle (same patient)
## =====================================================================
sig_show <- c(IL17_program = "IL-17 program", Cytotoxicity = "Cytotoxicity", Exhaustion = "Exhaustion")

pdelta <- paired_scores %>%
  dplyr::mutate(patient = paste("Patient", patient_id)) %>%
  dplyr::select(patient, tissue, monaco_fine, signature, mean_score) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = mean_score) %>%
  dplyr::filter(!is.na(heart), !is.na(skeletal)) %>%
  dplyr::mutate(delta = heart - skeletal)

p_ct_order <- pdelta %>% dplyr::filter(signature == "IL17_program") %>%
  dplyr::group_by(monaco_fine) %>% dplyr::summarise(m = mean(delta), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(monaco_fine)

## Panel: program heart vs skeletal (dumbbell), per patient
p_dumb <- pdelta %>%
  dplyr::filter(signature == "IL17_program") %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = p_ct_order)) %>%
  ggplot(aes(y = monaco_fine)) +
  geom_segment(aes(x = skeletal, xend = heart, yend = monaco_fine), colour = "grey60") +
  geom_point(aes(x = skeletal, colour = "Skeletal muscle"), size = 3) +
  geom_point(aes(x = heart, colour = "Heart"), size = 3) +
  facet_wrap(~ patient) +
  scale_colour_manual(values = c(Heart = col_tmos, `Skeletal muscle` = "grey45"), name = NULL) +
  labs(x = "IL-17 program score (mean)", y = NULL,
       title = "Same patient: heart vs skeletal muscle T/NK cells") +
  bt + theme(legend.position = "top")

## Panel: heart − skeletal for program, cytotoxicity, exhaustion
p_hm <- pdelta %>%
  dplyr::filter(signature %in% names(sig_show)) %>%
  dplyr::mutate(signature = factor(sig_show[signature], levels = sig_show),
                monaco_fine = factor(monaco_fine, levels = p_ct_order)) %>%
  ggplot(aes(patient, monaco_fine, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", delta)), size = 2.8) +
  facet_wrap(~ signature, nrow = 1) +
  scale_fill_gradient2(low = col_mc, mid = "white", high = col_tmos, midpoint = 0,
                       name = "Heart −\nskeletal") +
  labs(x = NULL, y = NULL, title = "State scores: heart minus skeletal muscle") +
  bt

p_dumb; p_hm
ggsave(file.path(out_dir, "Paired_IL17program_heart_vs_skeletal.pdf"), p_dumb, width = 8, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Paired_state_delta_heart_minus_skeletal.pdf"), p_hm, width = 9, height = 5, device = cairo_pdf)


## =====================================================================
## Combined figure + tables
## =====================================================================
fig <- (p_comp | p_delta) / p_donor / (p_dumb | p_hm) +
  plot_layout(heights = c(1, 0.9, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(file.path(out_dir, "R2.1_TNK_states_combined.pdf"), fig, width = 17, height = 17, device = cairo_pdf)
ggsave(file.path(out_dir, "R2.1_TNK_states_combined.png"), fig, width = 17, height = 17, dpi = 300, bg = "white")

write_xlsx(list(Composition_stats      = comp_stats,
                Program_rank_celltype  = prog_rank,
                Program_per_donor      = prog_donor,
                Paired_heart_minus_skel = pdelta),
           file.path(out_dir, "R2.1_TNK_states_summary.xlsx"))








library(Seurat)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
paired   <- c("4352", "4420")

map_heart <- c(Myeloid_cells = "Myeloid", Endothelial_cells = "Endothelial",
               Fibroblast_stromal = "Fibroblast/Stromal", Smooth_muscle_mural = "Smooth muscle/Mural",
               B_cells = "B cells", Cardiomyocytes = "Parenchymal")

heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"

cells <- colnames(heart)[sub("_.*", "", heart$sample_id) %in% paired]
ct    <- as.character(heart$major_cell_type[cells])
cnt   <- LayerData(heart, assay = "RNA", layer = "counts")[, cells]
rm(heart); gc()

tnk_meta <- readRDS(file.path(out_dir, "heart_TNK_monaco_fine_meta.rds"))
lab  <- unname(map_heart[ct])
is_t <- ct == "T_NK_lymphocytes"
lab[is_t] <- tnk_meta[cells[is_t], "monaco_fine"]

keep <- !is.na(lab)
obj  <- CreateSeuratObject(counts = cnt[, keep])
rm(cnt); gc()
obj <- NormalizeData(obj)

heart_x <- list(data = GetAssayData(obj, layer = "data"),
                meta = data.frame(labels = lab[keep], row.names = colnames(obj)))
rm(obj); gc()

print(table(heart_x$meta$labels))
saveRDS(heart_x, file.path(out_dir, "CellChat_input_heart.rds"))






library(Seurat)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
paired   <- c("4352", "4420")

map_skel <- c(Myeloid = "Myeloid", Endothelial = "Endothelial",
              `Fibroblast/Stromal` = "Fibroblast/Stromal", `Smooth Muscle` = "Smooth muscle/Mural",
              `B cells` = "B cells", `Skeletal Muscle` = "Parenchymal")

skel <- readRDS(file.path(base_dir, "sc_skeletal_322026.rds"))
DefaultAssay(skel) <- "RNA"

cells <- colnames(skel)[as.character(skel$patient_id) %in% paired]
ct    <- as.character(skel$major_match[cells])
cnt   <- LayerData(skel, assay = "RNA", layer = "counts")[, cells]
rm(skel); gc()

tnk_meta <- readRDS(file.path(out_dir, "skeletal_TNK_monaco_fine_meta.rds"))
lab  <- unname(map_skel[ct])
is_t <- ct == "T/NK Lymphocytes"
lab[is_t] <- tnk_meta[cells[is_t], "monaco_fine"]

keep <- !is.na(lab)
obj  <- CreateSeuratObject(counts = cnt[, keep])
rm(cnt); gc()
obj <- NormalizeData(obj)

skel_x <- list(data = GetAssayData(obj, layer = "data"),
               meta = data.frame(labels = lab[keep], row.names = colnames(obj)))
rm(obj); gc()

print(table(skel_x$meta$labels))
saveRDS(skel_x, file.path(out_dir, "CellChat_input_skeletal.rds"))




library(CellChat)
library(future)
plan("sequential")

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")

tissue <- "skeletal"     # "skeletal" for script 2b

x <- readRDS(file.path(out_dir, paste0("CellChat_input_", tissue, ".rds")))

## drop groups with <10 cells
keep_groups <- names(which(table(x$meta$labels) >= 10))
keep   <- x$meta$labels %in% keep_groups
x$data <- x$data[, keep]
x$meta <- x$meta[keep, , drop = FALSE]
print(table(x$meta$labels))

cc <- createCellChat(object = x$data, meta = x$meta, group.by = "labels")
rm(x); gc()

cc@DB <- CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type = "triMean", nboot = 20)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

saveRDS(cc, file.path(out_dir, paste0("CellChat_", tissue, "_paired.rds")))







library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
col_tmos <- "#b30000"

t_focus <- c("Th1 cells", "Th1/Th17 cells", "Effector memory CD8 T cells",
             "Non-Vd2 gd T cells", "Vd2 gd T cells")

cc_heart <- readRDS(file.path(out_dir, "CellChat_heart_paired.rds"))
cc_skel  <- readRDS(file.path(out_dir, "CellChat_skeletal_paired.rds"))

## ---- Interactions involving the focus T subsets ----
get_comm <- function(cc, tissue) {
  subsetCommunication(cc) %>%
    dplyr::filter(target %in% t_focus | source %in% t_focus) %>%
    dplyr::mutate(tissue = tissue,
                  direction = ifelse(target %in% t_focus, "to T subsets", "from T subsets"))
}
comm <- dplyr::bind_rows(get_comm(cc_heart, "heart"), get_comm(cc_skel, "skeletal"))

## ---- Pathway strength per tissue and direction ----
pw <- comm %>%
  dplyr::group_by(pathway_name, direction, tissue) %>%
  dplyr::summarise(strength = sum(prob), n_pairs = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue, values_from = c(strength, n_pairs), values_fill = 0) %>%
  dplyr::mutate(delta = strength_heart - strength_skeletal) %>%
  dplyr::arrange(dplyr::desc(delta))
as.data.frame(head(pw, 20))

## ---- Plot: top 10 heart-higher and top 5 skeletal-higher per direction ----
pw_plot <- pw %>%
  dplyr::group_by(direction) %>%
  dplyr::arrange(dplyr::desc(delta), .by_group = TRUE) %>%
  dplyr::filter(dplyr::row_number() <= 10 | dplyr::row_number() > dplyr::n() - 5) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(higher_in = ifelse(delta > 0, "Heart", "Skeletal muscle"),
                pathway_name = tidytext::reorder_within(pathway_name, delta, direction))

p_lr <- ggplot(pw_plot, aes(delta, pathway_name, fill = higher_in)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ direction, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = c(Heart = col_tmos, `Skeletal muscle` = "grey55"), name = NULL) +
  labs(x = "Communication strength (heart − skeletal muscle)", y = NULL,
       title = "Signalling to and from Th1, Th1/Th17, effector memory CD8 and γδ T cells",
       subtitle = "CellChat, paired TMOS patients (n = 2), heart vs skeletal muscle") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), legend.position = "top")
p_lr
ggsave(file.path(out_dir, "CellChat_heart_vs_skeletal_Tfocus_pathways.pdf"), p_lr,
       width = 11, height = 6, device = cairo_pdf)

## ---- Save ----
write_xlsx(list(Pathways = pw, Interactions = comm),
           file.path(out_dir, "CellChat_heart_vs_skeletal_Tfocus.xlsx"))







## =====================================================================
## R2.1 – CellChat: global heart vs skeletal muscle comparison
## Uses the two saved CellChat objects only (light; no Seurat objects)
## =====================================================================
library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")

col_heart <- "#b30000"
col_skel  <- "grey45"
tis_cols  <- c(Heart = col_heart, `Skeletal muscle` = col_skel)

non_immune <- c("Endothelial", "Fibroblast/Stromal", "Smooth muscle/Mural", "Parenchymal")
compartment <- function(x) dplyr::case_when(
  x %in% non_immune ~ "Stromal/vascular/parenchymal",
  x == "Myeloid"    ~ "Myeloid",
  x == "B cells"    ~ "B cells",
  TRUE              ~ "T/NK")
comp_cols <- c(`Stromal/vascular/parenchymal` = "#7f7f7f", Myeloid = "#E69F00",
               `B cells` = "#009E73", `T/NK` = col_heart)

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold", size = 11),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

cc_heart <- readRDS(file.path(out_dir, "CellChat_heart_paired.rds"))
cc_skel  <- readRDS(file.path(out_dir, "CellChat_skeletal_paired.rds"))

## ---- All significant interactions, restricted to groups present in both tissues ----
common <- intersect(levels(cc_heart@idents), levels(cc_skel@idents))
common <- setdiff(common, "Parenchymal")   # cardiomyocytes too few to assess
common
setdiff(union(levels(cc_heart@idents), levels(cc_skel@idents)), common)   # groups dropped

comm <- dplyr::bind_rows(
  subsetCommunication(cc_heart) %>% dplyr::mutate(tissue = "Heart"),
  subsetCommunication(cc_skel)  %>% dplyr::mutate(tissue = "Skeletal muscle")) %>%
  dplyr::mutate(source = as.character(source), target = as.character(target)) %>%
  dplyr::filter(source %in% common, target %in% common) %>%
  dplyr::mutate(tissue = factor(tissue, levels = names(tis_cols)),
                src_comp = compartment(source), tgt_comp = compartment(target))

## =====================================================================
## A. Overall number and strength of interactions
## =====================================================================
overall <- comm %>% dplyr::group_by(tissue) %>%
  dplyr::summarise(`Number of interactions` = dplyr::n(),
                   `Interaction strength`   = sum(prob), .groups = "drop") %>%
  tidyr::pivot_longer(-tissue, names_to = "metric", values_to = "value")

pA <- ggplot(overall, aes(tissue, value, fill = tissue)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(value > 100, round(value), sprintf("%.1f", value))), vjust = -0.4, size = 3.3) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = tis_cols, guide = "none") +
  labs(x = NULL, y = NULL, title = "Overall communication") +
  bt + theme(axis.text.x = element_text(angle = 20, hjust = 1))

## =====================================================================
## B. Information flow per pathway (all cell groups)
## =====================================================================
flow <- comm %>% dplyr::group_by(pathway_name, tissue) %>%
  dplyr::summarise(flow = sum(prob), .groups = "drop") %>%
  tidyr::complete(pathway_name, tissue, fill = list(flow = 0)) %>%
  dplyr::group_by(pathway_name) %>%
  dplyr::mutate(total = sum(flow), rel = flow / total) %>% dplyr::ungroup()

flow_wide <- flow %>% dplyr::select(pathway_name, tissue, flow) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = flow) %>%
  dplyr::mutate(rel_heart = Heart / (Heart + `Skeletal muscle`),
                log2_ratio = log2((Heart + 1e-3) / (`Skeletal muscle` + 1e-3))) %>%
  dplyr::arrange(dplyr::desc(rel_heart))

top_pw <- flow %>% dplyr::distinct(pathway_name, total) %>%
  dplyr::slice_max(total, n = 30) %>% dplyr::pull(pathway_name)
pw_order <- flow_wide %>% dplyr::filter(pathway_name %in% top_pw) %>%
  dplyr::arrange(rel_heart) %>% dplyr::pull(pathway_name)

pB <- flow %>% dplyr::filter(pathway_name %in% top_pw) %>%
  dplyr::mutate(pathway_name = factor(pathway_name, levels = pw_order)) %>%
  ggplot(aes(rel, pathway_name, fill = tissue)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "white") +
  scale_fill_manual(values = tis_cols, name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Relative information flow", y = NULL,
       title = "Signalling pathways: heart vs skeletal muscle (top 30 by total flow)") +
  bt + theme(legend.position = "top")

## =====================================================================
## C. Differential interaction strength between cell groups (heart − muscle)
## =====================================================================
grp_order <- c(intersect(c(non_immune, "Myeloid", "B cells"), common),
               sort(setdiff(common, c(non_immune, "Myeloid", "B cells"))))

pair_w <- comm %>% dplyr::group_by(tissue, source, target) %>%
  dplyr::summarise(w = sum(prob), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue, values_from = w, values_fill = 0) %>%
  tidyr::complete(source = common, target = common, fill = list(Heart = 0, `Skeletal muscle` = 0)) %>%
  dplyr::mutate(delta = Heart - `Skeletal muscle`,
                source = factor(source, levels = grp_order),
                target = factor(target, levels = grp_order))
lim <- max(abs(pair_w$delta))

pC <- ggplot(pair_w, aes(target, source, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = col_skel, mid = "white", high = col_heart, midpoint = 0,
                       limits = c(-lim, lim), name = "Heart −\nmuscle") +
  labs(x = "Receiver", y = "Sender",
       title = "Differential interaction strength between cell groups") +
  bt + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 8),
             axis.text.y = element_text(size = 8))

## =====================================================================
## D. Signalling roles: change in outgoing vs incoming strength per cell group
## =====================================================================
roles <- dplyr::full_join(
  comm %>% dplyr::group_by(tissue, group = source) %>% dplyr::summarise(out = sum(prob), .groups = "drop"),
  comm %>% dplyr::group_by(tissue, group = target) %>% dplyr::summarise(inc = sum(prob), .groups = "drop"),
  by = c("tissue", "group")) %>%
  dplyr::mutate(dplyr::across(c(out, inc), ~ tidyr::replace_na(.x, 0))) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = c(out, inc), values_fill = 0) %>%
  dplyr::mutate(d_out = out_Heart - `out_Skeletal muscle`,
                d_in  = inc_Heart - `inc_Skeletal muscle`,
                comp  = compartment(group))

pD <- ggplot(roles, aes(d_in, d_out, colour = comp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  geom_text_repel(aes(label = group), size = 2.8, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = comp_cols, name = NULL) +
  labs(x = "Δ incoming strength (heart − muscle)", y = "Δ outgoing strength (heart − muscle)",
       title = "Change in signalling role per cell group") +
  bt + theme(legend.position = "bottom")

## =====================================================================
## E. Key ligand–receptor pairs by sender compartment
## =====================================================================
lr <- comm %>% dplyr::group_by(interaction_name_2, pathway_name, tissue) %>%
  dplyr::summarise(s = sum(prob), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tissue, values_from = s, values_fill = 0) %>%
  dplyr::mutate(delta = Heart - `Skeletal muscle`) %>% dplyr::arrange(dplyr::desc(delta))

top_lr <- c(head(lr$interaction_name_2, 10), tail(lr$interaction_name_2, 6))

lr_bub <- comm %>%
  dplyr::filter(interaction_name_2 %in% top_lr, tgt_comp == "T/NK",
                src_comp %in% c("Stromal/vascular/parenchymal", "T/NK")) %>%
  dplyr::mutate(route = ifelse(src_comp == "T/NK", "T/NK → T/NK", "Stromal/vascular → T/NK")) %>%
  dplyr::group_by(interaction_name_2, tissue, route) %>%
  dplyr::summarise(s = sum(prob), .groups = "drop") %>%
  tidyr::complete(interaction_name_2 = top_lr, tissue, route, fill = list(s = 0)) %>%
  dplyr::mutate(interaction_name_2 = factor(interaction_name_2, levels = rev(top_lr)),
                tissue = factor(tissue, levels = names(tis_cols)))

pE <- ggplot(dplyr::filter(lr_bub, s > 0), aes(tissue, interaction_name_2, size = s, colour = tissue)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = 6.5, linetype = "dashed", colour = "grey60") +
  facet_wrap(~ route) +
  scale_colour_manual(values = tis_cols, guide = "none") +
  scale_size_continuous(range = c(1, 8), name = "Summed\nprobability") +
  labs(x = NULL, y = NULL,
       title = "Key ligand–receptor pairs to T/NK cells",
       subtitle = "Above dashed line: heart-higher; below: muscle-higher") +
  bt + theme(axis.text.x = element_text(angle = 20, hjust = 1),
             plot.subtitle = element_text(size = 8.5, colour = "grey30"))
## =====================================================================
## Save
## =====================================================================
ggsave(file.path(out_dir, "CellChat_A_overall.pdf"),        pA, width = 5,  height = 4,   device = cairo_pdf)
ggsave(file.path(out_dir, "CellChat_B_pathway_flow.pdf"),   pB, width = 6,  height = 7,   device = cairo_pdf)
ggsave(file.path(out_dir, "CellChat_C_diff_heatmap.pdf"),   pC, width = 8,  height = 7,   device = cairo_pdf)
ggsave(file.path(out_dir, "CellChat_D_roles.pdf"),          pD, width = 7,  height = 6,   device = cairo_pdf)
ggsave(file.path(out_dir, "CellChat_E_LR_bubble.pdf"),      pE, width = 12, height = 6,   device = cairo_pdf)

fig <- (pB | pE) / (pC | pD) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = list(c("F", "I", "G", "H")),
                  title = "Cell–cell communication: heart vs skeletal muscle (paired TMOS patients, n = 2)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))) &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(file.path(out_dir, "R2.1_CellChat_FtoI.pdf"), fig, width = 17, height = 15, device = cairo_pdf)
ggsave(file.path(out_dir, "R2.1_CellChat_FtoI.png"), fig, width = 17, height = 15, dpi = 300, bg = "white")

write_xlsx(list(Overall       = overall,
                Pathway_flow  = flow_wide,
                Pair_strength = pair_w,
                Roles         = roles,
                LR_pairs      = lr),
           file.path(out_dir, "R2.1_CellChat_global_summary.xlsx"))




library(magick)
library(patchwork)
library(ggplot2)
library(grid)

out_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments/R2.1_TNK_states")

## read a saved PDF panel as an image panel
img_panel <- function(file) {
  img <- magick::image_read_pdf(file.path(out_dir, file), density = 300)
  wrap_elements(full = rasterGrob(as.raster(img), interpolate = TRUE))
}

p1 <- img_panel("Heart_TNK_composition_sig_TMOS_vs_MC.pdf")
p2 <- img_panel("Heart_IL17program_delta_by_celltype.pdf")
p3 <- img_panel("Paired_IL17program_heart_vs_skeletal.pdf")
p4 <- img_panel("Paired_state_delta_heart_minus_skeletal.pdf")

## CellChat panels already in the session
p5 <- pB + labs(title = "Signalling pathways: heart vs skeletal muscle")

p7 <- pD + labs(title = "Change in signalling role per cell group")


p6 <- img_panel("CellChat_C_diff_heatmap.pdf")    # or: pC + labs(title = "Differential interaction strength (heart − muscle)")

row1 <- (p1 | p2)      + plot_layout(widths = c(1.375, 1.3))   # A 5.5x4, B 6.5x5
row2 <- (p3 | p4)      + plot_layout(widths = c(1.6, 1.8))     # C 8x5,   D 9x5
row3 <- (p5 | p6 | p7) + plot_layout(widths = c(5, 8, 6))      # E, F (8x7), G

fig_r21 <- row1 / row2 / row3 +
  plot_layout(heights = c(7.1, 5.6, 7)) +
  plot_annotation(tag_levels = "A",
                  title = "T/NK subsets, states and local interactions in TMOS",
                  theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))) &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(file.path(out_dir, "Response_Figure_R2.1.pdf"), fig_r21, width = 19, height = 21, device = cairo_pdf)
ggsave(file.path(out_dir, "Response_Figure_R2.1.png"), fig_r21, width = 19, height = 21, dpi = 300, bg = "white")






library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")
col_mc   <- "#4682B4"
col_tmos <- "#b30000"

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold", size = 11))

tnk_labels <- c("Effector memory CD8 T cells", "Terminal effector CD8 T cells", "Central memory CD8 T cells",
                "Naive CD8 T cells", "Vd2 gd T cells", "Non-Vd2 gd T cells", "MAIT cells",
                "Natural killer cells", "T regulatory cells", "Th1 cells", "Th1/Th17 cells", "Th17 cells",
                "Th2 cells", "Follicular helper T cells", "Naive CD4 T cells", "Terminal effector CD4 T cells")
gd_labels <- c("Vd2 gd T cells", "Non-Vd2 gd T cells")

meta_heart <- readRDS(file.path(out_dir, "heart_TNK_monaco_fine_meta.rds"))

## ---- % γδ T (Vd2 + Non-Vd2) of T/NK per donor ----
gd <- meta_heart %>%
  dplyr::filter(!is.na(monaco_fine), monaco_fine %in% tnk_labels) %>%
  dplyr::mutate(grp = factor(ifelse(m_group == "MC", "MC", "TMOS"), levels = c("MC", "TMOS"))) %>%
  dplyr::group_by(sample_id, grp) %>%
  dplyr::summarise(n_total = dplyr::n(),
                   n_gd    = sum(monaco_fine %in% gd_labels),
                   prop    = 100 * n_gd / n_total, .groups = "drop") %>%
  dplyr::filter(n_total >= 10)

table(gd$grp)

## ---- Stats ----
wt <- wilcox.test(prop ~ grp, data = gd, exact = FALSE)
gd_stats <- gd %>% dplyr::group_by(grp) %>%
  dplyr::summarise(n_donors = dplyr::n(), median = median(prop),
                   IQR_low = quantile(prop, 0.25), IQR_high = quantile(prop, 0.75), .groups = "drop") %>%
  dplyr::mutate(wilcox_p = wt$p.value)
gd_stats

## ---- Plot ----
y_bar <- max(gd$prop) * 1.12

p_gd <- ggplot(gd, aes(grp, prop, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2.2, position = position_jitter(width = 0.12, seed = 1)) +
  annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar) +
  annotate("text", x = 1.5, y = y_bar * 1.06, label = sprintf("p = %.3f", wt$p.value), size = 3.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "% of T/NK cells per donor",
       title = "Heart γδ T cells (Vd2 + non-Vd2)") +
  bt
p_gd

ggsave(file.path(out_dir, "Heart_gdT_combined_TMOS_vs_MC.pdf"), p_gd, width = 3, height = 3, device = cairo_pdf)
write_xlsx(list(Per_donor = gd, Stats = gd_stats),
           file.path(out_dir, "Heart_gdT_combined_TMOS_vs_MC.xlsx"))






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







library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(patchwork)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.1_TNK_states")

col_tmos  <- "#b30000"
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
        plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

## ---- Lollipop: TMOS − MC by cell type ----
prog_rank <- read_excel(file.path(out_dir, "R2.1_TNK_states_summary.xlsx"),
                        sheet = "Program_rank_celltype") %>%
  dplyr::mutate(monaco_fine = relab(monaco_fine)) %>%
  dplyr::arrange(dplyr::desc(delta))

ct_order <- rev(prog_rank$monaco_fine)          # largest delta on top
top3     <- prog_rank$monaco_fine[1:3]

p_delta <- prog_rank %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = ct_order),
                top = monaco_fine %in% top3) %>%
  ggplot(aes(delta, monaco_fine)) +
  geom_segment(aes(x = 0, xend = delta, yend = monaco_fine), colour = "grey60") +
  geom_point(aes(fill = mean_score_TMOS, stroke = ifelse(top, 1.2, 0.4)),
             shape = 21, size = 3.5, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", delta)), hjust = -0.6, size = 2.8) +
  scale_fill_gradient(low = "white", high = col_tmos, name = "Score in\nTMOS") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Δ program score (TMOS − MC)", y = NULL, title = "Heart: TMOS vs MC") +
  bt

## ---- Heatmap: heart − skeletal, same rows ----
pd_il17 <- read_excel(file.path(out_dir, "R2.1_TNK_states_summary.xlsx"),
                      sheet = "Paired_heart_minus_skel") %>%
  dplyr::filter(signature == "IL17_program") %>%
  dplyr::mutate(monaco_fine = relab(monaco_fine)) %>%
  dplyr::filter(monaco_fine %in% ct_order) %>%
  dplyr::mutate(monaco_fine = factor(monaco_fine, levels = ct_order))

lim_i <- max(abs(pd_il17$delta))

p_il17_hm <- ggplot(pd_il17, aes(patient, monaco_fine, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient2(low = col_skel, mid = "white", high = col_heart, midpoint = 0,
                       limits = c(-lim_i, lim_i), name = "Heart −\nskeletal") +
  labs(x = NULL, y = NULL, title = "Paired TMOS:\nheart vs skeletal") +
  bt +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ---- Combine ----
fig_il17 <- (p_delta | p_il17_hm) +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(title = "IL-17 signaling–related program in T/NK cells",
                  theme = theme(plot.title = element_text(face = "bold", size = 11)))
fig_il17

ggsave(file.path(out_dir, "Heart_IL17program_TMOSvsMC_and_paired.pdf"), fig_il17,
       width = 6.5, height = 4, device = cairo_pdf)

## end of Part 3: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 4 - Fig 2I, 2J: Reactome pathway module scores across tissues and groups
##############################################################################

# ============================================================================
# Reactome pathway module scores (IL-17, MAPK, NF-kB, TNF, IFNg-MAPK, RUNX1-FOXP3) across tissues/groups
#
# Figure panels : Fig 2I, 2J
# ============================================================================
library(tidyr)   # used below but not loaded in the original script

setwd(file.path(TMOS_DIR, "Clean/fig2/tissue biology analysis"))

library(Seurat)
library(reactome.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(tibble)

out_dir <- file.path(TMOS_DIR, "Clean/fig2/tissue biology analysis")

# ─── 1. Load Data ─────────────────────────────────────────────────────────────
heart_sc    <- readRDS("sc_heart_minor.rds")
skeletal_sc <- readRDS("sc_skeletal_322026.rds")

# ─── 2. Define Pathways ───────────────────────────────────────────────────────
reac_to_entrez <- as.list(reactomePATHID2EXTID)

pull_pathway <- function(pathway_id, pathway_name) {
  entrez <- reac_to_entrez[[pathway_id]]
  if (is.null(entrez)) { cat("WARNING: Not found:", pathway_id, "\n"); return(NULL) }
  syms <- AnnotationDbi::select(org.Hs.eg.db,
                                keys    = entrez,
                                columns = "SYMBOL",
                                keytype = "ENTREZID")$SYMBOL
  syms <- unique(na.omit(syms))
  cat(pathway_name, ":", length(syms), "genes\n")
  return(syms)
}

pathway_list <- list(
  MAPK        = pull_pathway("R-HSA-450294",  "MAPK Activation"),
  IL17        = pull_pathway("R-HSA-448424",  "IL-17 Signaling"),
  RUNX1_FOXP3 = pull_pathway("R-HSA-8877330", "RUNX1_FOXP3_Treg"),
  TNF         = pull_pathway("R-HSA-75893",   "TNF Signaling"),
  NFkB        = pull_pathway("R-HSA-9758274", "NFkB Regulation"),
  IFNg_MAPK   = pull_pathway("R-HSA-9732724", "IFNg_MAPK")
)

# ─── 3. Score Objects ─────────────────────────────────────────────────────────
add_scores <- function(obj, assay_use) {
  valid <- rownames(obj)
  for (path_name in names(pathway_list)) {
    genes_filtered <- pathway_list[[path_name]][pathway_list[[path_name]] %in% valid]
    cat(" ", path_name, "genes in panel:", length(genes_filtered), "\n")
    tmp_name <- paste0(path_name, "_tmp_")
    obj <- AddModuleScore(obj,
                          features = list(genes_filtered),
                          name     = tmp_name,
                          assay    = assay_use)
    obj@meta.data[[path_name]] <- obj@meta.data[[paste0(tmp_name, "1")]]
  }
  obj@meta.data <- obj@meta.data[, !grepl("_tmp_", colnames(obj@meta.data))]
  return(obj)
}

cat("\n--- Scoring Heart SC ---\n");    heart_sc    <- add_scores(heart_sc,    "RNA")
cat("\n--- Scoring Skeletal SC ---\n"); skeletal_sc <- add_scores(skeletal_sc, "RNA")

# ─── 4. Matrix Builder ────────────────────────────────────────────────────────
build_zscore_matrix_named <- function(meta, group_col) {
  path_cols <- names(pathway_list)
  meta      <- meta %>% filter(!is.na(.data[[group_col]]))
  for (p in path_cols) {
    meta[[paste0(p, "_Z")]] <- (meta[[p]] - mean(meta[[p]], na.rm = TRUE)) /
      sd(meta[[p]],   na.rm = TRUE)
  }
  mat <- meta %>%
    group_by(.data[[group_col]]) %>%
    summarise(across(all_of(paste0(path_cols, "_Z")), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    column_to_rownames(group_col) %>%
    as.matrix() %>%
    t()
  rownames(mat) <- path_cols
  return(mat)
}

# ─── 5. Filter & Build Matrix ─────────────────────────────────────────────────
h_filt <- subset(heart_sc,
                 subset = m_group         %in% c("MC", "MC+MS+MG") &
                   major_cell_type %in% c("T_NK_lymphocytes", "Myeloid_cells", "B_cells"))

s_filt <- subset(skeletal_sc,
                 subset = m_group     %in% c("MC", "MC+MS+MG") &
                   major_match %in% c("T/NK Lymphocytes", "Myeloid", "B cells"))

# Use full label internally for unique column IDs
h_filt$celltype_group <- paste0("H | ", h_filt$major_cell_type, " | ", h_filt$m_group)
s_filt$celltype_group <- paste0("S | ", s_filt$major_match,     " | ", s_filt$m_group)

mat_h        <- build_zscore_matrix_named(h_filt@meta.data, group_col = "celltype_group")
mat_s        <- build_zscore_matrix_named(s_filt@meta.data, group_col = "celltype_group")
shared_rows  <- intersect(rownames(mat_h), rownames(mat_s))
mat_combined <- cbind(mat_h[shared_rows, ], mat_s[shared_rows, ])

# ─── 6. Column Ordering ───────────────────────────────────────────────────────
col_hclust     <- hclust(dist(t(mat_combined)))
clustered_cols <- colnames(mat_combined)[col_hclust$order]
cols_mcmsmg    <- clustered_cols[grepl("MC\\+MS\\+MG$", clustered_cols)]
cols_mc        <- clustered_cols[grepl("MC$",            clustered_cols)]
mat_ordered    <- mat_combined[, c(cols_mcmsmg, cols_mc)]

# ─── 7. Annotations & clean labels ───────────────────────────────────────────
col_names  <- colnames(mat_ordered)
tissue_vec <- ifelse(grepl("^H \\|", col_names), "Heart", "Skeletal")
treat_vec  <- ifelse(grepl("MC\\+MS\\+MG$", col_names), "MC+MS+MG", "MC")

# Strip everything except the cell type (middle part between first and last " | ")
display_labels <- gsub("^[HS] \\| (.+) \\| MC(\\+MS\\+MG)?$", "\\1", col_names)

cat("Display labels:\n"); print(display_labels)  # verify before plotting

top_anno <- HeatmapAnnotation(
  Treatment = treat_vec,
  Tissue    = tissue_vec,
  col = list(
    Treatment = c("MC" = "#4682B4", "MC+MS+MG" = "#b30000"),
    Tissue    = c("Heart" = "#D95F02", "Skeletal" = "#1B9E77")
  ),
  annotation_name_side = "right",
  simple_anno_size     = unit(3, "mm")
)

# ─── 8. Build & Save Heatmap ──────────────────────────────────────────────────
n_mcmsmg <- length(cols_mcmsmg)
limit    <- max(abs(mat_ordered), na.rm = TRUE)
col_fun  <- colorRamp2(c(-limit, 0, limit), c("#4682B4", "white", "#b30000"))

ht <- Heatmap(
  mat_ordered,
  name                = "Z-score",
  col                 = col_fun,
  top_annotation      = top_anno,
  cluster_rows        = TRUE,
  cluster_columns     = TRUE,
  column_order        = seq_len(ncol(mat_ordered)),
  column_dend_reorder = FALSE,
  column_split        = c(rep(1, n_mcmsmg), rep(2, ncol(mat_ordered) - n_mcmsmg)),
  column_title_gp     = gpar(fontsize = 0),
  column_gap          = unit(3, "mm"),
  column_labels       = display_labels,
  row_names_gp        = gpar(fontsize = 10),
  column_names_gp     = gpar(fontsize = 7),
  column_names_rot    = 45,
  border              = FALSE,
  rect_gp             = gpar(col = NA),
  width               = unit(ncol(mat_ordered) * 10, "mm"),
  height              = unit(nrow(mat_ordered) * 9,  "mm")
)

# ─── Save with tight dimensions ───────────────────────────────────────────────
pdf(file.path(out_dir, "Combined_SC_6Pathways_ImmuneOnly_Final.pdf"))
ht_drawn <- draw(ht, merge_legend = TRUE, padding = unit(c(5, 15, 5, 5), "mm"))
w <- convertX(ComplexHeatmap:::width(ht_drawn),  "inches", valueOnly = TRUE)
h <- convertY(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
dev.off()

pdf(file.path(out_dir, "Combined_SC_6Pathways_ImmuneOnly_Final.pdf"),
    width  = w + 1,
    height = h + 0.5)
draw(ht, merge_legend = TRUE, padding = unit(c(5, 15, 5, 5), "mm"))
dev.off()

cat("Done — dimensions:", round(w + 1, 1), "x", round(h + 0.5, 1), "inches\n")





library(rstatix)
library(openxlsx)

h_tnk <- subset(h_filt, subset = major_cell_type == "T_NK_lymphocytes")
s_tnk <- subset(s_filt, subset = major_match == "T/NK Lymphocytes")
path_cols <- names(pathway_list)

# ─── Sample-level scores ──────────────────────────────────────────────────────
h_samp <- h_tnk@meta.data %>%
  group_by(sample_id, m_group) %>%
  summarise(across(all_of(path_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Group = paste0("Heart\n", m_group))

s_samp <- s_tnk@meta.data %>%
  group_by(sample_id) %>%
  summarise(across(all_of(path_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(m_group = "MC+MS+MG",
         Group   = "Skeletal\nMC+MS+MG")

samp_df <- bind_rows(h_samp, s_samp) %>%
  pivot_longer(cols      = all_of(path_cols),
               names_to  = "Pathway",
               values_to = "Score") %>%
  mutate(Group = factor(Group, levels = c("Heart\nMC",
                                          "Heart\nMC+MS+MG",
                                          "Skeletal\nMC+MS+MG")))

# ─── Sample counts per group (for x-axis labels) ──────────────────────────────
n_per_group <- samp_df %>%
  filter(Pathway == path_cols[1]) %>%   # one pathway is enough to count
  group_by(Group) %>%
  summarise(n = n_distinct(sample_id), .groups = "drop")

cat("Sample counts per group:\n"); print(n_per_group)

# Build label: "Heart\nMC\n(n=12)"
group_labels <- setNames(
  paste0(n_per_group$Group, "\n(n=", n_per_group$n, ")"),
  n_per_group$Group
)

samp_df <- samp_df %>%
  mutate(Group_label = factor(group_labels[as.character(Group)],
                              levels = group_labels))

# ─── Pairwise Wilcoxon per pathway ────────────────────────────────────────────
stats_list <- lapply(path_cols, function(pw) {
  df_pw <- samp_df %>% filter(Pathway == pw)
  
  pw_test <- df_pw %>%
    wilcox_test(Score ~ Group, comparisons = list(
      c("Heart\nMC",       "Heart\nMC+MS+MG"),
      c("Heart\nMC",       "Skeletal\nMC+MS+MG"),
      c("Heart\nMC+MS+MG", "Skeletal\nMC+MS+MG")
    )) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    mutate(Pathway = pw)
  
  pw_test
})

stats_df <- bind_rows(stats_list) %>%
  dplyr::select(Pathway, group1, group2, n1, n2,
                statistic, p, p.adj, p.adj.signif)

cat("\nPairwise stats:\n"); print(stats_df)

# ─── Save Excel ───────────────────────────────────────────────────────────────
wb <- createWorkbook()

addWorksheet(wb, "Pairwise_Stats")
writeData(wb, "Pairwise_Stats", stats_df)
setColWidths(wb, "Pairwise_Stats", cols = 1:ncol(stats_df), widths = "auto")

# Also add descriptive summary
desc_df <- samp_df %>%
  group_by(Pathway, Group) %>%
  summarise(
    N        = n_distinct(sample_id),
    Mean     = round(mean(Score,   na.rm = TRUE), 4),
    Median   = round(median(Score, na.rm = TRUE), 4),
    SD       = round(sd(Score,     na.rm = TRUE), 4),
    SEM      = round(SD / sqrt(N), 4),
    .groups  = "drop"
  )

addWorksheet(wb, "Descriptive_Summary")
writeData(wb, "Descriptive_Summary", desc_df)
setColWidths(wb, "Descriptive_Summary",
             cols = 1:ncol(desc_df), widths = "auto")

saveWorkbook(wb,
             file.path(out_dir, "TNK_PathwayScores_SampleLevel_Stats.xlsx"),
             overwrite = TRUE)

cat("Excel stats saved.\n")

library(ggpubr)

# ─── Prepare p-value annotations per pathway ──────────────────────────────────
# Use raw p values, compute y position per pathway dynamically
pval_anno <- stats_df %>%
  mutate(
    group1 = group_labels[group1],
    group2 = group_labels[group2]
  ) %>%
  left_join(
    samp_df %>%
      group_by(Pathway) %>%
      summarise(y_max = max(Score, na.rm = TRUE),
                y_range = max(Score, na.rm = TRUE) - min(Score, na.rm = TRUE),
                .groups = "drop"),
    by = "Pathway"
  ) %>%
  group_by(Pathway) %>%
  arrange(Pathway, group1, group2) %>%
  mutate(
    step      = row_number(),
    y.position = y_max + (y_range * 0.1 * step)  # stack p-values above max
  ) %>%
  ungroup() %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "p<0.001",
      p < 0.01  ~ paste0("p=", round(p, 3)),
      p < 0.05  ~ paste0("p=", round(p, 3)),
      TRUE      ~ paste0("p=", round(p, 2))
    )
  )


group_levels <- c("Heart\nMC", "Heart\nMC+MS+MG", "Skeletal\nMC+MS+MG")

fill_colors <- c(
  "Heart\nMC\n(n=9)"          = "#4682B4",
  "Heart\nMC+MS+MG\n(n=5)"    = "#b30000",
  "Skeletal\nMC+MS+MG\n(n=2)" = "#1B9E77"
)
color_colors <- fill_colors

# ─── Plot with p-values ───────────────────────────────────────────────────────
p_box <- ggplot(samp_df, aes(x = Group_label, y = Score,
                             fill = Group_label, color = Group_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4,
               width = 0.6, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
  stat_pvalue_manual(
    pval_anno,
    label         = "p_label",
    xmin          = "group1",
    xmax          = "group2",
    y.position    = "y.position",
    tip.length    = 0.01,
    bracket.size  = 0.4,
    size          = 2.8
  ) +
  facet_wrap(~ Pathway, nrow = 2, scales = "free_y") +
  scale_fill_manual(values  = fill_colors) +
  scale_color_manual(values = color_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw(base_size = 12) +
  labs(title    = "T/NK Cell Pathway Activation — Sample Level",
       subtitle = "Each dot = one sample; box = median ± IQR; p-values from Wilcoxon test",
       x = NULL, y = "Mean Module Score") +
  theme(panel.grid         = element_blank(),
        axis.text.x        = element_text(angle = 30, hjust = 1),
        strip.text         = element_text(face = "bold", size = 11),
        legend.position    = "none",
        panel.spacing      = unit(0.5, "lines"))

ggsave(file.path(out_dir, "TNK_PathwayModuleScores_SampleLevel_Boxplot.pdf"),
       plot = p_box, width = 10, height = 7)

cat("Box plot with p-values saved.\n")

## end of Part 4: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 5 - Fig 2K, 2L, 2M; Fig S2I: pre-ICI PBMC pseudobulk DE and MCO signature
##############################################################################

# ============================================================================
# Pre-ICI PBMC pseudobulk DE (MCO vs MC), 4-gene MCO signature score, longitudinal trajectories
#
# Figure panels : Fig 2K, 2L, 2M; Fig S2I
# ============================================================================

setwd(file.path(TMOS_DIR, "Clean/fig2"))
pbmc <- readRDS("pbmc_642025.rds")

# Memory: this script only uses the RNA counts and the metadata, so the other
# layers, assays, reductions and graphs are dropped right after loading. The
# full ~330,000-cell object otherwise runs out of memory at subset() on a
# desktop. Results are unchanged.
pbmc <- Seurat::DietSeurat(pbmc, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
invisible(gc())

# =============================================================================
# 0. LIBRARIES & SETUP
# =============================================================================

library(Seurat)
library(DESeq2)
library(apeglm)
library(tidyverse)
library(glmnet)
library(pROC)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(scales)
library(ComplexHeatmap)
library(openxlsx)

set.seed(77)

out_dir <- "TMOS_MCO_Signature_2026"
dir.create(out_dir,                      showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"),  showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

save_table <- function(df, name) {
  write.csv(df,
            file.path(out_dir, "tables", paste0(name, ".csv")),
            row.names = FALSE)
  invisible(df)
}

grp_colors <- c("MC" = "#4682B4", "MCO" = "#b30000")
grp_fill   <- c("MC" = "#4682B4", "MCO" = "#b30000")

border_theme <- theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  axis.line    = element_blank()
)

# =============================================================================
# 1. METADATA ENGINEERING
# =============================================================================

pbmc@meta.data <- pbmc@meta.data %>%
  mutate(
    deg_group = case_when(
      m_group == "MC"                      ~ "MC",
      m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
      TRUE                                 ~ NA_character_
    )
  )

message("=== Group counts ===")
print(table(pbmc@meta.data$deg_group, useNA = "ifany"))

pbmc_deg <- subset(pbmc, subset = deg_group %in% c("MC", "MCO"))
pbmc_deg <- subset(pbmc_deg, subset = timepoint_cat %in% c("pre_ici", "on_ici"))

message("=== Timepoint distribution ===")
print(table(pbmc_deg@meta.data$timepoint_cat, pbmc_deg@meta.data$deg_group))

message("=== Samples ===")
print(table(pbmc_deg@meta.data$sample_id, pbmc_deg@meta.data$deg_group))

# =============================================================================
# 2. PSEUDOBULK AGGREGATION — SAMPLE LEVEL
# =============================================================================

aggregate_pseudobulk <- function(seurat_obj,
                                 group_by   = "sample_id",
                                 count_slot = "counts") {
  meta       <- seurat_obj@meta.data
  samples    <- unique(meta[[group_by]])
  counts_mat <- GetAssayData(seurat_obj, layer = count_slot, assay = "RNA")
  
  pb_list <- lapply(samples, function(s) {
    cells <- rownames(meta[meta[[group_by]] == s, ])
    if (length(cells) == 1) {
      Matrix::Matrix(counts_mat[, cells, drop = FALSE], sparse = TRUE)
    } else {
      Matrix::rowSums(counts_mat[, cells])
    }
  })
  
  pb_mat           <- do.call(cbind, pb_list)
  colnames(pb_mat) <- samples
  
  sample_meta <- meta %>%
    mutate(across(where(is.factor), as.character)) %>%
    group_by(.data[[group_by]]) %>%
    summarise(
      donor     = dplyr::first(donor),
      deg_group = dplyr::first(deg_group),
      m_group   = dplyr::first(m_group),
      n_cells   = n(),
      .groups   = "drop"
    ) %>%
    column_to_rownames(group_by)
  
  sample_meta <- sample_meta[colnames(pb_mat), ]
  list(counts = pb_mat, meta = sample_meta)
}

pb <- aggregate_pseudobulk(pbmc_deg, group_by = "sample_id")

message(sprintf("Pseudobulk matrix: %d genes x %d samples",
                nrow(pb$counts), ncol(pb$counts)))
print(table(pb$meta$deg_group))

keep_genes <- rowSums(pb$counts >= 1) >= 2
pb$counts  <- pb$counts[keep_genes, ]
message(sprintf("Genes after filter: %d", sum(keep_genes)))

# =============================================================================
# 3. DESeq2 DIFFERENTIAL EXPRESSION
# Note: donor is fully confounded with deg_group (each donor = one group only)
# so design = ~ deg_group only. Samples from same donor are pseudo-replicates;
# this limitation is acknowledged in methods.
# =============================================================================

counts_int <- as.matrix(pb$counts)
storage.mode(counts_int) <- "integer"

pb$meta$deg_group <- factor(pb$meta$deg_group, levels = c("MC", "MCO"))
pb$meta$donor     <- factor(pb$meta$donor)

dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = pb$meta,
  design    = ~ deg_group
)

message("Running DESeq2 ...")
dds <- DESeq(dds, test = "Wald", fitType = "parametric", parallel = FALSE)

res_shrunk <- lfcShrink(dds,
                        coef  = "deg_group_MCO_vs_MC",
                        type  = "apeglm",
                        quiet = TRUE)

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(
    direction = case_when(
      log2FoldChange > 0.58 & padj < 0.05  ~ "Up in MCO",
      log2FoldChange < -0.58 & padj < 0.05 ~ "Down in MCO",
      TRUE                                  ~ "NS"
    ),
    neg_log10_padj = -log10(padj)
  )

message("=== DESeq2 complete ===")
print(table(res_df$direction))

sig_degs <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.58) %>%
  arrange(padj)

save_table(sig_degs, "DESeq2_significant_DEGs")
message(sprintf("Significant DEGs: %d", nrow(sig_degs)))
print(sig_degs)

# Save to Excel
wb <- createWorkbook()
addWorksheet(wb, "Full_DESeq2_results")
writeDataTable(wb, "Full_DESeq2_results", res_df,
               tableStyle = "TableStyleMedium9")
setColWidths(wb, "Full_DESeq2_results",
             cols = seq_len(ncol(res_df)), widths = "auto")
addWorksheet(wb, "Significant_DEGs")
writeDataTable(wb, "Significant_DEGs", sig_degs,
               tableStyle = "TableStyleMedium9")
setColWidths(wb, "Significant_DEGs",
             cols = seq_len(ncol(sig_degs)), widths = "auto")
saveWorkbook(wb, file.path(out_dir, "tables", "DESeq2_results.xlsx"),
             overwrite = TRUE)
message("DESeq2 results saved")

# =============================================================================
# 4. NORMALIZATION
# =============================================================================

norm_counts  <- counts(dds, normalized = TRUE)
log_norm_mat <- log2(norm_counts + 1)

message(sprintf("Normalized matrix: %d genes x %d samples",
                nrow(log_norm_mat), ncol(log_norm_mat)))

write.csv(as.data.frame(log_norm_mat),
          file.path(out_dir, "tables", "log2_normalized_counts_all_genes.csv"))

# =============================================================================
# 5. SIGNATURE GENES
# sig_genes     = upregulated only → used for scoring + trajectory
# sig_genes_all = up + down       → used for heatmap
# =============================================================================

sig_genes <- res_df %>%
  filter(padj < 0.05, log2FoldChange > 0.58) %>%
  arrange(padj) %>%
  pull(gene)

sig_genes_all <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.58) %>%
  arrange(desc(log2FoldChange)) %>%   # up first, then down
  pull(gene)

cat("Upregulated sig genes:", length(sig_genes), "\n")
print(sig_genes)
cat("Total sig genes (up + down):", length(sig_genes_all), "\n")

# Expression matrix for scoring — upregulated only
sig_log_mat <- t(log_norm_mat[sig_genes, ])
y           <- as.numeric(pb$meta$deg_group == "MCO")

# Save full sig gene table
sig_gene_df <- res_df %>%
  filter(gene %in% sig_genes_all) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange, lfcSE, pvalue, padj, direction) %>%
  mutate(fold_change = round(2^abs(log2FoldChange), 2))

save_table(sig_gene_df, "MCO_signature_genes_all")

# =============================================================================
# 6. MCO SIGNATURE SCORE — mean z-score of upregulated genes
# =============================================================================

sig_scores <- rowMeans(scale(sig_log_mat))

score_df <- data.frame(
  sample_id = rownames(sig_log_mat),
  deg_group = pb$meta$deg_group,
  score     = sig_scores
) %>%
  mutate(deg_group = factor(deg_group, levels = c("MC", "MCO")))

print(score_df)

wt <- wilcox.test(score ~ deg_group, data = score_df, exact = FALSE)
message(sprintf("Wilcoxon p-value: %.4f", wt$p.value))

save_table(score_df, "MCO_signature_scores")

# =============================================================================
# 7. VOLCANO PLOT — label all significant genes
# =============================================================================

# Label all significant genes
label_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.58) %>%
  pull(gene)

volcano_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(label = ifelse(gene %in% label_genes, gene, NA))

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj,
                        color = direction)) +
  
  geom_point(size = 1.8, alpha = 0.6) +
  
  geom_text_repel(aes(label = label),
                  size          = 3.5,
                  fontface      = "bold",
                  box.padding   = 0.5,
                  point.padding = 0.3,
                  max.overlaps  = 50,
                  na.rm         = TRUE) +
  
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  
  scale_color_manual(
    values = c("Up in MCO"   = "#b30000",
               "Down in MCO" = "#4682B4",
               "NS"          = "grey70"),
    name = NULL
  ) +
  
  labs(
    title    = "Differential Expression: MCO vs MC",
    subtitle = "Pseudobulk DESeq2 | FDR < 0.05 | |LFC| > 0.58",
    x        = "Log2 Fold Change (MCO vs MC)",
    y        = "-Log10 Adjusted P-value"
  ) +
  
  theme_classic(base_size = 13) +
  border_theme +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "top"
  )

print(p_volcano)

ggsave(file.path(out_dir, "plots", "volcano_MCO_vs_MC.pdf"),
       plot = p_volcano, width = 4, height = 5, device = "pdf")

# =============================================================================
# 8. SCORE BOXPLOT
# =============================================================================

p_score <- ggplot(score_df,
                  aes(x = deg_group, y = score, color = deg_group)) +
  geom_jitter(width = 0.15, size = 3) +
  geom_boxplot(aes(fill = deg_group), alpha = 0.3,
               outlier.shape = NA, width = 0.4) +
  scale_color_manual(values = grp_colors) +
  scale_fill_manual(values  = grp_fill) +
  annotate("text", x = 1.5,
           y     = max(score_df$score) + 0.3,
           label = sprintf("p = %.3f", wt$p.value),
           size  = 4) +
  labs(
    title = "MCO Signature Score",
    x     = "Group",
    y     = "Signature Score (mean z-score)"
  ) +
  theme_classic(base_size = 13) +
  border_theme +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

print(p_score)

ggsave(file.path(out_dir, "plots", "MCO_signature_score.pdf"),
       plot = p_score, width = 4, height = 4, device = "pdf")

# =============================================================================
# 9. HEATMAP — all significant DEGs (up + down)
# =============================================================================

heat_mat   <- log_norm_mat[sig_genes_all, ]
heat_mat_z <- t(scale(t(heat_mat)))  # z-score per gene

# Column annotation: group
col_anno <- HeatmapAnnotation(
  Group = pb$meta$deg_group,
  col   = list(Group = grp_colors),
  annotation_name_side = "left"
)

# Row annotation: direction — matched to sig_genes_all order
row_dir <- res_df %>%
  filter(gene %in% sig_genes_all) %>%
  arrange(match(gene, sig_genes_all)) %>%
  pull(direction)

row_anno <- rowAnnotation(
  Direction = row_dir,
  col = list(Direction = c("Up in MCO"   = "#b30000",
                           "Down in MCO" = "#4682B4")),
  annotation_name_side = "top"
)

# Order columns: MC first, then MCO
col_order <- order(pb$meta$deg_group)

ht <- Heatmap(
  heat_mat_z[, col_order],
  name                 = "Z-score",
  col                  = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  top_annotation       = col_anno[col_order],
  right_annotation     = row_anno,
  cluster_columns      = FALSE,
  cluster_rows         = TRUE,
  show_column_names    = TRUE,
  show_row_names       = TRUE,
  row_names_gp         = gpar(fontsize = 10),
  column_names_gp      = gpar(fontsize = 9),
  column_names_rot     = 45,
  heatmap_legend_param = list(title = "Z-score")
)

pdf(file.path(out_dir, "plots", "heatmap_sig_DEGs.pdf"),
    width  = 7,
    height = max(4, nrow(heat_mat_z) * 0.4))
draw(ht)
dev.off()



# =============================================================================
# 10. TRAJECTORY PLOT — sample level, upregulated sig genes
# =============================================================================

tp_order <- c("pre_ici", "on_ici", "pre_steroid", "post_steroid")

pbmc_traj <- subset(pbmc,
                    subset = m_group %in% c("MC", "MC+MS", "MC+MS+MG") &
                      timepoint_cat %in% tp_order)

pbmc_traj@meta.data <- pbmc_traj@meta.data %>%
  mutate(
    traj_group = case_when(
      m_group == "MC"                      ~ "MC",
      m_group %in% c("MC+MS", "MC+MS+MG") ~ "MCO",
      TRUE                                 ~ NA_character_
    ),
    donor = sub("_[^_]+$", "", sample_id)
  )

# Pseudobulk full counts per sample then CPM
counts_full  <- GetAssayData(pbmc_traj, layer = "counts", assay = "RNA")
meta_traj    <- pbmc_traj@meta.data
samples_traj <- unique(meta_traj$sample_id)

pb_full_list <- lapply(samples_traj, function(s) {
  cells <- rownames(meta_traj[meta_traj$sample_id == s, ])
  if (length(cells) == 1) {
    Matrix::Matrix(counts_full[, cells, drop = FALSE], sparse = TRUE)[, 1]
  } else {
    Matrix::rowSums(counts_full[, cells])
  }
})

pb_full_mat           <- do.call(cbind, pb_full_list)
colnames(pb_full_mat) <- samples_traj

lib_sizes    <- colSums(pb_full_mat)
cpm_full     <- sweep(pb_full_mat, 2, lib_sizes, "/") * 1e6
log_cpm_full <- log2(cpm_full + 1)

samp_meta_traj <- meta_traj %>%
  group_by(sample_id) %>%
  summarise(
    donor         = dplyr::first(donor),
    traj_group    = dplyr::first(traj_group),
    timepoint_cat = dplyr::first(timepoint_cat),
    n_cells       = n(),
    .groups       = "drop"
  )

# Upregulated genes only for trajectory
sig_genes_traj <- sig_genes[sig_genes %in% rownames(log_cpm_full)]

traj_long <- as.data.frame(t(log_cpm_full[sig_genes_traj, ])) %>%
  rownames_to_column("sample_id") %>%
  left_join(samp_meta_traj, by = "sample_id") %>%
  pivot_longer(cols      = all_of(sig_genes_traj),
               names_to  = "gene",
               values_to = "log2_cpm") %>%
  mutate(
    timepoint_cat = factor(timepoint_cat, levels = tp_order),
    traj_group    = factor(traj_group,    levels = c("MC", "MCO"))
  ) %>%
  filter(!is.na(timepoint_cat), !is.na(traj_group),
         donor != "SIC_3")

traj_summary <- traj_long %>%
  group_by(gene, timepoint_cat, traj_group) %>%
  summarise(
    mean_exp = mean(log2_cpm, na.rm = TRUE),
    se_exp   = sd(log2_cpm,  na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

p_traj <- ggplot(traj_summary,
                 aes(x = timepoint_cat, y = mean_exp,
                     color = traj_group, group = traj_group)) +
  
  geom_jitter(data  = traj_long,
              aes(x = timepoint_cat, y = log2_cpm,
                  color = traj_group, group = traj_group),
              width = 0.1, size = 1.5, alpha = 0.4) +
  
  geom_line(linewidth = 1.2) +
  geom_point(aes(y = mean_exp), size = 3) +
  
  geom_errorbar(aes(ymin = mean_exp - se_exp,
                    ymax = mean_exp + se_exp),
                width = 0.15, linewidth = 0.8) +
  
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = grp_colors, name = "Group") +
  
  scale_x_discrete(labels = c("pre_ici"      = "Pre-ICI",
                              "on_ici"       = "On-ICI",
                              "pre_steroid"  = "Pre-steroid",
                              "post_steroid" = "Post-steroid")) +
  
  labs(
    title   = "MCO Signature Gene Expression Across Timepoints",
    x       = "Timepoint",
    y       = "Expression (log2 CPM)",
    caption = "Points = individual samples; Lines = group mean ± SE"
  ) +
  
  theme_classic(base_size = 13) +
  border_theme +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    strip.text      = element_text(face = "bold", size = 11),
    legend.position = "top",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

print(p_traj)

n_genes <- length(sig_genes_traj)
plot_w  <- max(7, ceiling(n_genes / 2) * 3.5)
plot_h  <- ceiling(n_genes / 2) * 4

ggsave(file.path(out_dir, "plots", "trajectory_signature_genes.pdf"),
       plot = p_traj, width = plot_w, height = plot_h, device = "pdf")






# =============================================================================
# SIGNATURE GENE ANNOTATION TABLE — Excel
# =============================================================================

table_df <- data.frame(
  Gene = c("GPRC5D-AS1", "NBPF10", "TMEM51", "ZNF578"),
  Gene_Type = c("lncRNA", "Protein-coding", "Protein-coding", "Protein-coding"),
  Log2FoldChange = c(2.93, 0.74, 1.61, 2.16),
  FDR = c(0.000278, 0.012036, 0.023114, 0.023385),
  Known_Biology = c(
    "Antisense lncRNA to plasma cell marker GPRC5D; regulates skeletal muscle factors MyoD/MyoG via miR-520d-5p sponging",
    "NBPF family member; contains DUF1220 protein domains; 1q21.1 CNVs linked to inflammatory/immune disease; NBPF10 variants associated with ME/CFS",
    "Transmembrane protein; TMEM51 polymorphisms associated with Raynaud's syndrome via vasculature and immunity loci",
    "Zinc finger transcription factor; related ZNF members regulate NF-κB/IKK signaling and inflammatory cytokine production"
  ),
  Reference = c(
    "Yu et al., Aging 2023. DOI: 10.18632/aging.205279",
    "Nguyen et al., 2024; GeneCards NBPF10",
    "Hartmann et al., Cell Genomics 2024. DOI: 10.1016/j.xgen.2024.100630",
    "Liu et al., Cells 2019. DOI: 10.3390/cells8121604"
  ),
  stringsAsFactors = FALSE
)

# Create workbook
wb <- createWorkbook()
addWorksheet(wb, "Signature_Gene_Annotations")

# Write data
writeData(wb, "Signature_Gene_Annotations", table_df, startRow = 2, startCol = 1)

# Style header row
header_style <- createStyle(
  fontName    = "Arial",
  fontSize    = 11,
  fontColour  = "white",
  fgFill      = "#2C3E50",
  halign      = "center",
  valign      = "center",
  textDecoration = "bold",
  border      = "TopBottomLeftRight",
  borderColour = "white",
  wrapText    = TRUE
)

addStyle(wb, "Signature_Gene_Annotations",
         style = header_style,
         rows  = 2,
         cols  = 1:ncol(table_df),
         gridExpand = TRUE)

# Style odd rows
odd_style <- createStyle(
  fontName  = "Arial",
  fontSize  = 10,
  fgFill    = "#F7F7F7",
  valign    = "center",
  border    = "TopBottomLeftRight",
  borderColour = "grey80",
  wrapText  = TRUE
)

# Style even rows
even_style <- createStyle(
  fontName  = "Arial",
  fontSize  = 10,
  fgFill    = "white",
  valign    = "center",
  border    = "TopBottomLeftRight",
  borderColour = "grey80",
  wrapText  = TRUE
)

for (i in 1:nrow(table_df)) {
  style <- if (i %% 2 == 1) odd_style else even_style
  addStyle(wb, "Signature_Gene_Annotations",
           style = style,
           rows  = i + 2,
           cols  = 1:ncol(table_df),
           gridExpand = TRUE)
}

# Gene column bold
gene_style <- createStyle(
  fontName       = "Arial",
  fontSize       = 10,
  textDecoration = "bold",
  valign         = "center",
  wrapText       = TRUE
)

addStyle(wb, "Signature_Gene_Annotations",
         style = gene_style,
         rows  = 3:(nrow(table_df) + 2),
         cols  = 1,
         gridExpand = TRUE)

# Add title above table
title_style <- createStyle(
  fontName       = "Arial",
  fontSize       = 13,
  textDecoration = "bold",
  fontColour     = "#2C3E50"
)

writeData(wb, "Signature_Gene_Annotations",
          "MCO Blood Transcriptomic Signature — Upregulated Genes (Pre/On-ICI PBMC)",
          startRow = 1, startCol = 1)

addStyle(wb, "Signature_Gene_Annotations",
         style = title_style,
         rows  = 1, cols  = 1)

# Set column widths
setColWidths(wb, "Signature_Gene_Annotations",
             cols   = 1:ncol(table_df),
             widths = c(14, 14, 14, 10, 55, 45))

# Set row heights
setRowHeights(wb, "Signature_Gene_Annotations",
              rows    = 3:(nrow(table_df) + 2),
              heights = 60)

# Save
saveWorkbook(
  wb,
  file.path(out_dir, "tables", "signature_gene_annotations.xlsx"),
  overwrite = TRUE
)

message("Annotation table saved as Excel")

## end of Part 5: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 6 - Fig 2N: leave-one-donor-out validation of the blood signature
##############################################################################

# ============================================================================
# Leave-one-donor-out Firth logistic validation of the pre-MC blood signature
#
# Figure panels : Fig 2N
# ============================================================================

## =====================================================================
## R3.4 – Internal validation of the pre-ICI / on-ICI MCO blood signature
##   Logistic regression on the 4 signature genes (TMEM51, NBPF10, GPRC5D-AS1, ZNF578)
##   (A) Leave-one-donor-out (LODO), Firth logistic    <- primary
##   (B) LODO, standard glm logistic                   <- original approach
##   (C) Leave-one-sample-out (LOSO), Firth logistic   <- donor's other sample stays in training
##   (D) LODO, genes re-selected in each fold (Firth)  <- nested, from lodo_checkpoint.rds
## Normalization: DESeq2 size factors from the training samples; held-out samples
##   normalized to the training geometric means; log2(norm + 1); genes standardized
##   on training mean/SD. No housekeeping normalization.
## Needs: pseudobulk_inputs.rds and lodo_checkpoint.rds from the R3.4 main script.
## =====================================================================
library(DESeq2)
library(tidyverse)
library(pROC)
library(openxlsx)
library(cowplot)
if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf")
library(logistf)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R3.4_preMC_signature")
grp_colors <- c(MC = "#4682B4", MCO = "#b30000")

inp <- readRDS(file.path(out_dir, "pseudobulk_inputs.rds"))
pb_counts <- inp$pb_counts; pb_meta <- inp$pb_meta; rm(inp); gc()
lodo_ckpt <- readRDS(file.path(out_dir, "lodo_checkpoint.rds"))

y         <- as.integer(pb_meta$deg_group == "MCO")
obs_genes <- c("TMEM51", "NBPF10", "GPRC5D-AS1", "ZNF578")
stopifnot(all(obs_genes %in% rownames(pb_counts)))

## ── 1. HELPERS ───────────────────────────────────────────────────────────
norm_split <- function(tr, te, genes) {
  sf_tr <- estimateSizeFactorsForMatrix(pb_counts[, tr, drop = FALSE])
  geo   <- exp(rowMeans(log(pb_counts[, tr, drop = FALSE])))
  sf_te <- estimateSizeFactorsForMatrix(pb_counts[, te, drop = FALSE], geoMeans = geo)
  list(train = t(log2(sweep(pb_counts[genes, tr, drop = FALSE], 2, sf_tr, "/") + 1)),
       test  = t(log2(sweep(pb_counts[genes, te, drop = FALSE], 2, sf_te, "/") + 1)))
}

prep <- function(Xtr, ytr, Xte) {
  mu <- colMeans(Xtr); s <- apply(Xtr, 2, sd); s[s == 0] <- 1
  dtr <- as.data.frame(scale(Xtr, mu, s)); dte <- as.data.frame(scale(Xte, mu, s))
  names(dtr) <- names(dte) <- make.names(colnames(Xtr))
  dtr$y <- ytr
  list(dtr = dtr, dte = dte)
}

fit_firth <- function(Xtr, ytr, Xte) {
  d <- prep(Xtr, ytr, Xte)
  fit <- tryCatch(logistf(y ~ ., data = d$dtr, plcontrol = logistpl.control(maxit = 500)),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(d$dte)))
  as.numeric(predict(fit, newdata = d$dte, type = "response"))
}

fit_glm <- function(Xtr, ytr, Xte) {
  d <- prep(Xtr, ytr, Xte)
  fit <- tryCatch(suppressWarnings(glm(y ~ ., data = d$dtr, family = binomial(link = "logit"),
                                       control = glm.control(maxit = 200))),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(d$dte)))
  as.numeric(predict(fit, newdata = d$dte, type = "response"))
}

## unit = "donor" (LODO) or "sample" (LOSO); genes = "fixed" or "nested"
run_cv <- function(unit, fitter, genes = "fixed") {
  ids <- if (unit == "donor") unique(pb_meta$donor) else colnames(pb_counts)
  map_dfr(ids, function(id) {
    te <- if (unit == "donor") which(pb_meta$donor == id) else which(colnames(pb_counts) == id)
    tr <- setdiff(seq_along(y), te)
    if (length(unique(y[tr])) < 2) { message("  skip ", id, ": one class in training"); return(NULL) }
    g <- obs_genes
    if (genes == "nested") {
      g <- strsplit(unique(lodo_ckpt$genes[lodo_ckpt$donor == id]), ", ")[[1]]
      g <- intersect(g[g != ""], rownames(pb_counts))
    }
    pred <- if (length(g) == 0) rep(NA_real_, length(te)) else {
      x <- norm_split(tr, te, g); fitter(x$train, y[tr], x$test) }
    tibble(left_out = id, donor = pb_meta$donor[te], sample_id = colnames(pb_counts)[te],
           group = pb_meta$deg_group[te], true_y = y[te], n_genes = length(g), pred_prob = pred)
  })
}

## AUC, 95% CI (DeLong), Youden operating point
## AUC, 95% CI (DeLong), Youden operating point; NA if a class is missing
evaluate <- function(d) {
  ok <- !is.na(d$pred_prob)
  base <- tibble(n_predicted = sum(ok), n_not_predicted = sum(!ok),
                 n_MC = sum(d$true_y[ok] == 0), n_MCO = sum(d$true_y[ok] == 1))
  if (base$n_MC == 0 || base$n_MCO == 0)
    return(list(roc = NULL, stats = mutate(base, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                                           threshold = NA_real_, sensitivity = NA_real_,
                                           specificity = NA_real_, youden_J = NA_real_)))
  r <- roc(d$true_y[ok], d$pred_prob[ok], levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  b <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"),
              transpose = FALSE)[1, ]
  list(roc = r,
       stats = mutate(base, AUC = as.numeric(auc(r)), CI_low = ci(r)[1], CI_high = ci(r)[3],
                      threshold = b$threshold, sensitivity = b$sensitivity, specificity = b$specificity,
                      youden_J = b$sensitivity + b$specificity - 1))
}

## ── 2. RUN ───────────────────────────────────────────────────────────────
cv_list <- list(
  "LODO - Firth (fixed genes)"         = run_cv("donor",  fit_firth),
  "LODO - glm (fixed genes)"           = run_cv("donor",  fit_glm),
  "LOSO - Firth (fixed genes)"         = run_cv("sample", fit_firth),
  "LODO - Firth (genes re-selected)"   = run_cv("donor",  fit_firth, genes = "nested"))

ev <- map(cv_list, evaluate)
perf <- imap_dfr(ev, ~ mutate(.x$stats, method = .y, .before = 1)) %>%
  mutate(across(c(AUC, CI_low, CI_high, sensitivity, specificity, youden_J), ~ round(.x, 3)),
         threshold = signif(threshold, 3))
print(perf, width = Inf)

## predicted class at each method's Youden threshold
cv_list <- imap(cv_list, function(d, nm) {
  th <- ev[[nm]]$stats$threshold
  d %>% mutate(pred_class = ifelse(is.na(pred_prob), NA, ifelse(pred_prob >= th, "MCO", "MC")),
               correct = pred_class == group)
})

primary <- cv_list[["LODO - Firth (fixed genes)"]]
conf_mat <- table(Predicted = factor(primary$pred_class, levels = c("MC", "MCO")),
                  Actual    = factor(primary$group,      levels = c("MC", "MCO")))
cat("\nConfusion matrix, LODO Firth (Youden threshold):\n"); print(conf_mat)

donor_summary <- primary %>%
  group_by(donor, group) %>%
  summarise(n_samples = n(), mean_pred_prob = round(mean(pred_prob), 3),
            pct_correct = round(mean(correct) * 100, 1), .groups = "drop") %>%
  arrange(group, donor)
print(donor_summary)

## ── 3. SAVE TABLES ───────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Performance");  writeData(wb, "Performance", perf)
addWorksheet(wb, "Confusion_LODO_Firth")
writeData(wb, "Confusion_LODO_Firth", as.data.frame.matrix(conf_mat) %>% rownames_to_column("Predicted"))
addWorksheet(wb, "Per_donor_LODO_Firth"); writeData(wb, "Per_donor_LODO_Firth", donor_summary)
sheet_names <- c("Pred_LODO_Firth", "Pred_LODO_glm", "Pred_LOSO_Firth", "Pred_LODO_nested")
walk2(cv_list, sheet_names, function(d, sh) { addWorksheet(wb, sh); writeData(wb, sh, d) })
saveWorkbook(wb, file.path(out_dir, "R3.4_LODO_LOSO_logistic.xlsx"), overwrite = TRUE)

## ── 4. FIGURE ────────────────────────────────────────────────────────────
bt <- theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold"))

show <- c("LODO - Firth (fixed genes)", "LODO - glm (fixed genes)", "LOSO - Firth (fixed genes)")
roc_df <- map_dfr(show, function(nm) {
  r <- ev[[nm]]$roc
  tibble(fpr = 1 - r$specificities, tpr = r$sensitivities,
         model = sprintf("%s: AUC %.2f (%.2f–%.2f)", nm,
                         ev[[nm]]$stats$AUC, ev[[nm]]$stats$CI_low, ev[[nm]]$stats$CI_high))
}) %>% arrange(model, fpr, tpr)
lvl <- unique(roc_df$model)[match(show, sub(":.*", "", unique(roc_df$model)))]
roc_df$model <- factor(roc_df$model, levels = lvl)

p_roc <- ggplot(roc_df, aes(fpr, tpr, colour = model, linetype = model)) +
  geom_abline(linetype = "dashed", colour = "grey70") +
  geom_step(direction = "vh", linewidth = 1) +
  scale_colour_manual(values = c("#b30000", "grey30", "grey60"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate", title = "Held-out ROC") +
  bt + theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 9))

th_primary <- ev[["LODO - Firth (fixed genes)"]]$stats$threshold
p_prob <- ggplot(primary, aes(group, pred_prob, colour = group)) +
  geom_hline(yintercept = th_primary, linetype = "dashed", colour = "grey50") +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 3) +
  scale_colour_manual(values = grp_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Held-out predicted P(MCO)", title = "LODO – Firth",
       subtitle = sprintf("Sens %.2f, Spec %.2f at Youden threshold",
                          ev[["LODO - Firth (fixed genes)"]]$stats$sensitivity,
                          ev[["LODO - Firth (fixed genes)"]]$stats$specificity)) + bt

fig <- plot_grid(p_roc, p_prob, nrow = 1, rel_widths = c(1.25, 0.75), labels = c("A", "B"), align = "h", axis = "b")
ggsave(file.path(out_dir, "R3.4_LODO_LOSO_logistic.pdf"), fig, width = 9.5, height = 5.5)
ggsave(file.path(out_dir, "R3.4_LODO_LOSO_logistic.png"), fig, width = 9.5, height = 5.5, dpi = 300, bg = "white")
cat("Saved to", out_dir, "\n")







p_prob <- ggplot(primary, aes(group, pred_prob, colour = group)) +
  geom_hline(yintercept = th_primary, linetype = "dashed", colour = "grey50") +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 3) +
  scale_colour_manual(values = grp_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Held-out predicted P(MCO)", title = "LODO – Firth",
       subtitle = sprintf("Sens %.2f, Spec %.2f\n(Youden threshold)",
                          ev[["LODO - Firth (fixed genes)"]]$stats$sensitivity,
                          ev[["LODO - Firth (fixed genes)"]]$stats$specificity)) + bt

fig <- plot_grid(p_roc, p_prob, nrow = 1, rel_widths = c(1.35, 0.65), labels = c("A", "B"), align = "h", axis = "b")
ggsave(file.path(out_dir, "R3.4_LODO_LOSO_logistic.pdf"), fig, width = 8, height = 5.5)
ggsave(file.path(out_dir, "R3.4_LODO_LOSO_logistic.png"), fig, width = 8, height = 5.5, dpi = 300, bg = "white")

## end of Part 6: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 7 - Fig S2D, S2E: heart broad cell-type composition, MC vs MC+MS+MG
##############################################################################

# ============================================================================
# Heart broad cell-type composition, MC vs MC+MS+MG (TMOS)
#
# Figure panels : Fig S2D (row "Total Tissue" of
#                   MCvsTMOS_Full_Atlas_Boxplots_PerfectSpacing.pdf; also saved
#                   on its own as MCvsTMOS_Major_CellType_Boxplots_TotalTissue.pdf)
#                 Fig S2E (MCvsTMOS_CellType_Heatmap.pdf)
# ============================================================================

library(Seurat)
library(dplyr)

setwd(file.path(TMOS_DIR, "Single_cell/Heart/Seurat_by_CellType/Combined"))
merged <- readRDS("heart_minor.rds")


# ---------------- Fig S2E: broad cell-type proportion heatmap ----------------
################ Heat map #########################################


# 1. Filter and set specific factor levels
groups_to_keep <- c("MC", "MC+MS+MG")

meta_sub <- merged@meta.data %>%
  filter(m_group %in% groups_to_keep) %>%
  mutate(m_group = factor(m_group, levels = groups_to_keep))

# 2. Calculate cell type percentages
count_table <- table(meta_sub$major_cell_type, meta_sub$sample_id)
prop_matrix <- prop.table(count_table, margin = 2) * 100

# 3. APPLY ROW-WISE Z-SCORE SCALING
# t() is used because scale() works on columns by default
scaled_matrix <- t(scale(t(prop_matrix)))

# 4. Prepare the Column Annotation and Order
# Ensure dplyr is loaded
library(dplyr)

# 4. Prepare the Column Annotation and Order
anno_df <- meta_sub %>%
  dplyr::select(sample_id, m_group) %>%  # This avoids the "inherited method" conflict
  distinct() %>%
  arrange(m_group, sample_id)


rownames(anno_df) <- anno_df$sample_id

# Reorder the scaled matrix columns to match the sorted annotation
scaled_matrix <- scaled_matrix[, rownames(anno_df)]


m_group_colors <- c(
  "control"   = "#bdbdbd",
  "MC"        = "#4682B4", # Using Teal for high differentiation
  "MC+MS+MG"  = "#b30000"
)


library(ComplexHeatmap)
library(circlize)

col_anno <- HeatmapAnnotation(
  m_group = anno_df$m_group,
  col = list(m_group = m_group_colors),
  show_legend = TRUE
)

# 6. Create the Heatmap
Heatmap(
  scaled_matrix,
  name = "Z-score",
  top_annotation = col_anno,
  cluster_rows = TRUE,          # Cluster cell types based on scaling
  cluster_columns = FALSE,      # Keep your manual group order
  column_split = anno_df$m_group, # Visual gaps between groups
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 10),
  # Color scale centered at 0 (White). Blue = Low, Red = High.
  col = colorRamp2(c(-2, 0, 2), c("#0571b0", "#f7f7f7", "#ca0020")),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)




# 1. Define the specific directory and file name
output_dir <- paste0(TMOS_DIR, "/Single_cell/Heart/Seurat_by_CellType/Combined/")
file_name <- "MCvsTMOS_CellType_Heatmap.pdf"
full_path <- paste0(output_dir, file_name)

# 2. Open PDF device 
# Width 12 and Height 7 is usually a good starting point for these sample sizes
pdf(full_path, width = 10, height = 5.5)

# 3. Generate the Heatmap
# Note: rect_gp is removed to ensure no internal grid lines/borders
ht <- Heatmap(
  scaled_matrix,
  name = "Z-score",
  top_annotation = col_anno,
  cluster_rows = TRUE,          
  cluster_columns = FALSE,      
  column_split = anno_df$m_group, 
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 10),
  
  # Your requested color palette: Blue (Low), Yellow (Mid), Red (High)
  col = colorRamp2(c(-2, 0, 2), c("#0571b0", "#f7f7f7", "#ca0020")),
  
  border = TRUE,
  use_raster = FALSE # Set to TRUE if the file size becomes too large
)

# 4. Draw the heatmap and close the file
draw(ht, merge_legend = TRUE)
dev.off()

message("Heatmap successfully saved to: ", full_path)

# ---------------- T/NK subset (used by the T/NK row of the boxplot figure) ----------------
library(ComplexHeatmap)
library(circlize)
tnk_subset <- subset(merged, subset = major_cell_type == "T_NK_lymphocytes" & 
                       m_group %in% c("MC", "MC+MS+MG"))

# ---------------- Fig S2D: proportion boxplots (Total Tissue row = S2D) ----------------
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(writexl)

# --- SETTINGS ---
output_dir <- paste0(TMOS_DIR, "/Single_cell/Heart/Seurat_by_CellType/Combined/")
m_group_colors <- c("MC" = "#4682B4", "MC+MS+MG" = "#b30000")
my_comparison <- list(c("MC", "MC+MS+MG"))

all_stats <- list()

# --- REFINED HELPER FUNCTION ---
process_compartment <- function(df, major_type = NULL, plot_title, y_label, group_col = "minor") {
  
  if (!is.null(major_type)) {
    sub_data <- df %>% filter(major_cell_type == major_type)
  } else {
    sub_data <- df 
    group_col <- "major_cell_type"
  }
  
  plot_data <- sub_data %>%
    filter(m_group %in% c("MC", "MC+MS+MG")) %>%
    group_by(sample_id, m_group, !!sym(group_col)) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(proportion = (count / sum(count)) * 100) %>%
    ungroup() %>%
    mutate(m_group = factor(m_group, levels = c("MC", "MC+MS+MG")))
  
  # Stats Calculation
  stats <- plot_data %>%
    group_by(!!sym(group_col)) %>%
    do(w_test = wilcox.test(proportion ~ m_group, data = .)) %>%
    summarise(CellType = !!sym(group_col), 
              p_value = w_test$p.value,
              Compartment = plot_title)
  all_stats[[plot_title]] <<- stats 
  
  # Create Plot
  p <- ggplot(plot_data, aes(x = m_group, y = proportion, fill = m_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) + 
    facet_wrap(as.formula(paste("~", group_col)), scales = "free_y", ncol = 4) +
    scale_fill_manual(values = m_group_colors) +
    # Adjusting stat_compare_means: 
    # vjust = 0 resets the text position to the default (above the bracket)
    # bracket.nudge.y pushes the whole bracket up slightly to clear the dots
    stat_compare_means(comparisons = my_comparison, 
                       method = "wilcox.test", 
                       label = "p.format", 
                       size = 3.2,
                       vjust = 0,
                       bracket.nudge.y = 1) + 
    # expand adds 20% space to the top of the axis so the bracket doesn't hit the ceiling
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) + 
    theme_bw() +
    labs(title = plot_title, y = y_label, x = NULL) +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.y = element_text(size = 9),
          plot.title = element_text(size = 12, face = "bold"),
          panel.spacing = unit(0.8, "lines"),
          legend.position = "none")
  
  return(p)
}

# --- GENERATE INDIVIDUAL PLOTS ---
p1 <- process_compartment(merged@meta.data, NULL, "Total Tissue", "% of Total Cells")
p2 <- process_compartment(tnk_subset@meta.data, "T_NK_lymphocytes", "T/NK Pool", "% of T/NK")
p3 <- process_compartment(merged@meta.data, "Myeloid_cells", "Myeloid Pool", "% of Myeloid")
p4 <- process_compartment(merged@meta.data, "Fibroblast_stromal", "Fibroblast Pool", "% of Fibroblast")

# --- COMBINE AND SAVE ---
combined_all <- (p1 / p2 / p3 / p4) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(paste0(output_dir, "MCvsTMOS_Full_Atlas_Boxplots_PerfectSpacing.pdf"), 
       combined_all, width = 11, height = 22, limitsize = FALSE)

# Save Excel Stats
final_stats_df <- bind_rows(all_stats)
write_xlsx(final_stats_df, paste0(output_dir, "MCvsTMOS_Statistics_Summary.xlsx"))

# Fig S2D on its own (same plot object as the "Total Tissue" row above)
ggsave(paste0(output_dir, "MCvsTMOS_Major_CellType_Boxplots_TotalTissue.pdf"), p1, width = 11, height = 6)

## end of Part 7: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 8a - Fig S2F (part a): CellChat objects, merged and by institution
##############################################################################

# ============================================================================
# CellChat MC+MS+MG vs MC: merged, MDACC-only and MGH-only
#
# Figure panels : Fig S2F
# ============================================================================

setwd(file.path(TMOS_DIR, "Single_cell/Heart/Seurat_by_CellType/Combined/updated_cellchat_institute"))
heart <- readRDS("heart_minor.rds")



# Subset heart to only MC and MC+MS+MG
heart_MC_subset <- subset(
  heart,
  subset = m_group %in% c("MC", "MC+MS+MG")
)

# Check
table(heart_MC_subset$m_group)
table(heart_MC_subset$minor)

heart <- heart_MC_subset
heart$minor <- NULL

saveRDS(heart, "heart_cellchat.rds")



###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

library(CellChat)

# Split by condition
heart_MC <- subset(heart, m_group == "MC")
heart_MCMSMG <- subset(heart, m_group == "MC+MS+MG")

# Create CellChat objects
cc_MC <- createCellChat(object = heart_MC, group.by = "major_cell_type")
cc_MCMSMG <- createCellChat(object = heart_MCMSMG, group.by = "major_cell_type")

# Set DB
CellChatDB <- CellChatDB.human
cc_MC@DB <- CellChatDB
cc_MCMSMG@DB <- CellChatDB

# Preprocessing
cc_MC <- subsetData(cc_MC)
cc_MC <- identifyOverExpressedGenes(cc_MC)
cc_MC <- identifyOverExpressedInteractions(cc_MC)

cc_MCMSMG <- subsetData(cc_MCMSMG)
cc_MCMSMG <- identifyOverExpressedGenes(cc_MCMSMG)
cc_MCMSMG <- identifyOverExpressedInteractions(cc_MCMSMG)


detach("package:CellChat", unload = TRUE)
library(CellChat)
library(future)
future::plan("sequential") 

# Infer communication
cc_MC <- computeCommunProb(cc_MC, type = "triMean", nboot = 1)
cc_MC <- filterCommunication(cc_MC, min.cells = 10)
cc_MC <- computeCommunProbPathway(cc_MC)
cc_MC <- aggregateNet(cc_MC)

detach("package:CellChat", unload = TRUE)
library(CellChat)
library(future)
future::plan("sequential") 

cc_MCMSMG <- computeCommunProb(cc_MCMSMG, type = "triMean", nboot = 1)
cc_MCMSMG <- filterCommunication(cc_MCMSMG, min.cells = 10)
cc_MCMSMG <- computeCommunProbPathway(cc_MCMSMG)
cc_MCMSMG <- aggregateNet(cc_MCMSMG)

# Merge
cc_list <- list(MC = cc_MC, MC_MS_MG = cc_MCMSMG)
cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))




# Number and strength of interactions - heatmap
compareInteractions(cc_merged, show.legend = FALSE, group = c(1,2), measure = "count")
compareInteractions(cc_merged, show.legend = FALSE, group = c(1,2), measure = "weight")

# Circle plots
par(mfrow = c(1,2))
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count", title.name = "Differential number of interactions")
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight", title.name = "Differential interaction strength")

gg1 <- netVisual_heatmap(cc_merged, measure = "count")
gg2 <- netVisual_heatmap(cc_merged, measure = "weight")
gg1 + gg2




# Save circle plots
object.list <- list(MC = cc_MC, MC_MS_MG = cc_MCMSMG)
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))


table(heart_MC@meta.data$institution)
table(heart_MCMSMG@meta.data$institution)






# Split by institution and condition
heart_MC_mgh <- subset(heart_MC, institution == "mgh")
heart_MC_mdacc <- subset(heart_MC, institution == "mdacc")
heart_MCMSMG_mgh <- subset(heart_MCMSMG, institution == "mgh")
heart_MCMSMG_mdacc <- subset(heart_MCMSMG, institution == "mdacc")

# (interactive checks, moved off: these objects are created further down)
# levels(cc_MC_mgh@idents)
# levels(cc_MC_mdacc@idents)
# levels(cc_MCMSMG_mgh@idents)
# levels(cc_MCMSMG_mdacc@idents)


keep_celltypes <- c("T_NK_lymphocytes", "Myeloid_cells", "Fibroblast_stromal", "Smooth_muscle_mural", "Endothelial_cells", "Cardiomyocytes")

heart_MC_mgh <- subset(heart_MC_mgh, major_cell_type %in% keep_celltypes)
heart_MC_mdacc <- subset(heart_MC_mdacc, major_cell_type %in% keep_celltypes)
heart_MCMSMG_mgh <- subset(heart_MCMSMG_mgh, major_cell_type %in% keep_celltypes)
heart_MCMSMG_mdacc <- subset(heart_MCMSMG_mdacc, major_cell_type %in% keep_celltypes)

# Drop the now-empty cell-type levels (newer CellChat versions stop on unused
# factor levels in the identities; the cells used are unchanged)
heart_MC_mgh$major_cell_type       <- droplevels(factor(heart_MC_mgh$major_cell_type))
heart_MC_mdacc$major_cell_type     <- droplevels(factor(heart_MC_mdacc$major_cell_type))
heart_MCMSMG_mgh$major_cell_type   <- droplevels(factor(heart_MCMSMG_mgh$major_cell_type))
heart_MCMSMG_mdacc$major_cell_type <- droplevels(factor(heart_MCMSMG_mdacc$major_cell_type))




# Recreate CellChat objects with consistent cell types
cc_MC_mgh <- createCellChat(object = heart_MC_mgh, group.by = "major_cell_type")
cc_MC_mdacc <- createCellChat(object = heart_MC_mdacc, group.by = "major_cell_type")
cc_MCMSMG_mgh <- createCellChat(object = heart_MCMSMG_mgh, group.by = "major_cell_type")
cc_MCMSMG_mdacc <- createCellChat(object = heart_MCMSMG_mdacc, group.by = "major_cell_type")

# Set DB
cc_MC_mgh@DB <- CellChatDB.human
cc_MC_mdacc@DB <- CellChatDB.human
cc_MCMSMG_mgh@DB <- CellChatDB.human
cc_MCMSMG_mdacc@DB <- CellChatDB.human

# Preprocess
cc_MC_mgh <- subsetData(cc_MC_mgh)
cc_MC_mgh <- identifyOverExpressedGenes(cc_MC_mgh)
cc_MC_mgh <- identifyOverExpressedInteractions(cc_MC_mgh)

cc_MC_mdacc <- subsetData(cc_MC_mdacc)
cc_MC_mdacc <- identifyOverExpressedGenes(cc_MC_mdacc)
cc_MC_mdacc <- identifyOverExpressedInteractions(cc_MC_mdacc)

cc_MCMSMG_mgh <- subsetData(cc_MCMSMG_mgh)
cc_MCMSMG_mgh <- identifyOverExpressedGenes(cc_MCMSMG_mgh)
cc_MCMSMG_mgh <- identifyOverExpressedInteractions(cc_MCMSMG_mgh)

cc_MCMSMG_mdacc <- subsetData(cc_MCMSMG_mdacc)
cc_MCMSMG_mdacc <- identifyOverExpressedGenes(cc_MCMSMG_mdacc)
cc_MCMSMG_mdacc <- identifyOverExpressedInteractions(cc_MCMSMG_mdacc)

# Infer communication
cc_MC_mgh <- computeCommunProb(cc_MC_mgh, type = "triMean", nboot = 1)
cc_MC_mgh <- filterCommunication(cc_MC_mgh, min.cells = 10)
cc_MC_mgh <- computeCommunProbPathway(cc_MC_mgh)
cc_MC_mgh <- aggregateNet(cc_MC_mgh)

cc_MC_mdacc <- computeCommunProb(cc_MC_mdacc, type = "triMean", nboot = 1)
cc_MC_mdacc <- filterCommunication(cc_MC_mdacc, min.cells = 10)
cc_MC_mdacc <- computeCommunProbPathway(cc_MC_mdacc)
cc_MC_mdacc <- aggregateNet(cc_MC_mdacc)

cc_MCMSMG_mgh <- computeCommunProb(cc_MCMSMG_mgh, type = "triMean", nboot = 1)
cc_MCMSMG_mgh <- filterCommunication(cc_MCMSMG_mgh, min.cells = 10)
cc_MCMSMG_mgh <- computeCommunProbPathway(cc_MCMSMG_mgh)
cc_MCMSMG_mgh <- aggregateNet(cc_MCMSMG_mgh)

cc_MCMSMG_mdacc <- computeCommunProb(cc_MCMSMG_mdacc, type = "triMean", nboot = 1)
cc_MCMSMG_mdacc <- filterCommunication(cc_MCMSMG_mdacc, min.cells = 10)
cc_MCMSMG_mdacc <- computeCommunProbPathway(cc_MCMSMG_mdacc)
cc_MCMSMG_mdacc <- aggregateNet(cc_MCMSMG_mdacc)

# Merge by institution
cc_merged_mgh <- mergeCellChat(list(MC = cc_MC_mgh, MC_MS_MG = cc_MCMSMG_mgh), add.names = c("MC", "MC_MS_MG"))
cc_merged_mdacc <- mergeCellChat(list(MC = cc_MC_mdacc, MC_MS_MG = cc_MCMSMG_mdacc), add.names = c("MC", "MC_MS_MG"))









# ── MERGED (all samples) ──────────────────────────────────────────────────────

object.list <- list(MC = cc_MC, MC_MS_MG = cc_MCMSMG)
weight.max.all <- getMaxWeight(object.list, attribute = c("idents", "count"))
weight.max.all.w <- getMaxWeight(object.list, attribute = c("idents", "weight"))

pdf("circle_count_all.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.all[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("circle_strength_all.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.all.w[2], edge.width.max = 12,
                   title.name = paste0("Interaction strength - ", names(object.list)[i]))
}
dev.off()

pdf("circle_diff_all.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength (MC+MS+MG vs MC)")
dev.off()

pdf("heatmaps_all.pdf", width = 12, height = 6)
gg1 <- netVisual_heatmap(cc_merged, measure = "count", title.name = "All - Count")
gg2 <- netVisual_heatmap(cc_merged, measure = "weight", title.name = "All - Strength")
print(gg1 + gg2)
dev.off()

# ── MGH ──────────────────────────────────────────────────────────────────────

object.list.mgh <- list(MC = cc_MC_mgh, MC_MS_MG = cc_MCMSMG_mgh)
weight.max.mgh <- getMaxWeight(object.list.mgh, attribute = c("idents", "count"))
weight.max.mgh.w <- getMaxWeight(object.list.mgh, attribute = c("idents", "weight"))

pdf("circle_count_mgh.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list.mgh)) {
  netVisual_circle(object.list.mgh[[i]]@net$count, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.mgh[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list.mgh)[i], " (MGH)"))
}
dev.off()

pdf("circle_strength_mgh.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list.mgh)) {
  netVisual_circle(object.list.mgh[[i]]@net$weight, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.mgh.w[2], edge.width.max = 12,
                   title.name = paste0("Interaction strength - ", names(object.list.mgh)[i], " (MGH)"))
}
dev.off()

pdf("circle_diff_mgh.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged_mgh, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions - MGH (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged_mgh, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength - MGH (MC+MS+MG vs MC)")
dev.off()

pdf("heatmaps_mgh.pdf", width = 12, height = 6)
gg1 <- netVisual_heatmap(cc_merged_mgh, measure = "count", title.name = "MGH - Count")
gg2 <- netVisual_heatmap(cc_merged_mgh, measure = "weight", title.name = "MGH - Strength")
print(gg1 + gg2)
dev.off()

# ── MDACC ─────────────────────────────────────────────────────────────────────

object.list.mdacc <- list(MC = cc_MC_mdacc, MC_MS_MG = cc_MCMSMG_mdacc)
weight.max.mdacc <- getMaxWeight(object.list.mdacc, attribute = c("idents", "count"))
weight.max.mdacc.w <- getMaxWeight(object.list.mdacc, attribute = c("idents", "weight"))

pdf("circle_count_mdacc.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list.mdacc)) {
  netVisual_circle(object.list.mdacc[[i]]@net$count, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.mdacc[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list.mdacc)[i], " (MDACC)"))
}
dev.off()

pdf("circle_strength_mdacc.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list.mdacc)) {
  netVisual_circle(object.list.mdacc[[i]]@net$weight, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.mdacc.w[2], edge.width.max = 12,
                   title.name = paste0("Interaction strength - ", names(object.list.mdacc)[i], " (MDACC)"))
}
dev.off()

pdf("circle_diff_mdacc.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged_mdacc, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions - MDACC (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged_mdacc, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength - MDACC (MC+MS+MG vs MC)")
dev.off()

pdf("heatmaps_mdacc.pdf", width = 12, height = 6)
gg1 <- netVisual_heatmap(cc_merged_mdacc, measure = "count", title.name = "MDACC - Count")
gg2 <- netVisual_heatmap(cc_merged_mdacc, measure = "weight", title.name = "MDACC - Strength")
print(gg1 + gg2)
dev.off()




# Get difference matrices directly
diff_count_all <- cc_merged@net[[2]]$count - cc_merged@net[[1]]$count
diff_count_mgh <- cc_merged_mgh@net[[2]]$count - cc_merged_mgh@net[[1]]$count
diff_count_mdacc <- cc_merged_mdacc@net[[2]]$count - cc_merged_mdacc@net[[1]]$count

diff_weight_all <- cc_merged@net[[2]]$weight - cc_merged@net[[1]]$weight
diff_weight_mgh <- cc_merged_mgh@net[[2]]$weight - cc_merged_mgh@net[[1]]$weight
diff_weight_mdacc <- cc_merged_mdacc@net[[2]]$weight - cc_merged_mdacc@net[[1]]$weight

# Get symmetric limits
count_limit <- max(abs(c(diff_count_all, diff_count_mgh, diff_count_mdacc)), na.rm = TRUE)
weight_limit <- max(abs(c(diff_weight_all, diff_weight_mgh, diff_weight_mdacc)), na.rm = TRUE)

count_limit
weight_limit


col_fun_count <- colorRamp2(c(-53, 0, 53), c("#2166AC", "#FFFFBF", "#D6604D"))
col_fun_weight <- colorRamp2(c(-1.682688, 0, 1.682688), c("#2166AC", "#FFFFBF", "#D6604D"))



make_heatmaps <- function(diff_count, diff_weight, title_prefix) {
  h1 <- Heatmap(diff_count, name = "count", col = col_fun_count,
                column_title = paste0(title_prefix, " - Count (MC+MS+MG vs MC)"),
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_side = "left", column_names_side = "bottom",
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 9),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(round(diff_count[i, j], 1), x, y, gp = gpar(fontsize = 7))
                })
  h2 <- Heatmap(diff_weight, name = "weight", col = col_fun_weight,
                column_title = paste0(title_prefix, " - Strength (MC+MS+MG vs MC)"),
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_side = "left", column_names_side = "bottom",
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 9),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(round(diff_weight[i, j], 2), x, y, gp = gpar(fontsize = 7))
                })
  draw(h1 + h2)
}

pdf("uni_heatmaps_all.pdf", width = 14, height = 7)
make_heatmaps(diff_count_all, diff_weight_all, "All")
dev.off()

pdf("uni_heatmaps_mgh.pdf", width = 14, height = 7)
make_heatmaps(diff_count_mgh, diff_weight_mgh, "MGH")
dev.off()

pdf("uni_heatmaps_mdacc.pdf", width = 14, height = 7)
make_heatmaps(diff_count_mdacc, diff_weight_mdacc, "MDACC")
dev.off()









# Get sign matrices for each comparison
sign_all <- sign(diff_count_all)
sign_mgh <- sign(diff_count_mgh)
sign_mdacc <- sign(diff_count_mdacc)

# Find cells where all three agree in direction (and are not zero)
# First align row/col names since all and institute have different cell types
common_cells <- intersect(rownames(diff_count_all), rownames(diff_count_mgh))

sign_all_sub <- sign_all[common_cells, common_cells]
sign_mgh_sub <- sign_mgh[common_cells, common_cells]
sign_mdacc_sub <- sign_mdacc[common_cells, common_cells]

# Consistent = all three same sign and none are zero
consistent_count <- (sign_all_sub == sign_mgh_sub) & 
  (sign_mgh_sub == sign_mdacc_sub) & 
  (sign_all_sub != 0)

# Same for weight
sign_all_w <- sign(diff_weight_all[common_cells, common_cells])
sign_mgh_w <- sign(diff_weight_mgh[common_cells, common_cells])
sign_mdacc_w <- sign(diff_weight_mdacc[common_cells, common_cells])

consistent_weight <- (sign_all_w == sign_mgh_w) & 
  (sign_mgh_w == sign_mdacc_w) & 
  (sign_all_w != 0)

# Print consistent interactions
cat("=== COUNT: Consistent direction across all 3 ===\n")
which(consistent_count, arr.ind = TRUE) |> 
  as.data.frame() |>
  dplyr::mutate(sender = rownames(diff_count_all[common_cells, common_cells])[row],
                receiver = colnames(diff_count_all[common_cells, common_cells])[col],
                direction = ifelse(sign_all_sub[consistent_count] > 0, "UP in MC+MS+MG", "DOWN in MC+MS+MG"),
                all = diff_count_all[common_cells, common_cells][consistent_count],
                mgh = diff_count_mgh[consistent_count],
                mdacc = diff_count_mdacc[consistent_count]) |>
  print()

cat("\n=== WEIGHT: Consistent direction across all 3 ===\n")
which(consistent_weight, arr.ind = TRUE) |>
  as.data.frame() |>
  dplyr::mutate(sender = rownames(diff_weight_all[common_cells, common_cells])[row],
                receiver = colnames(diff_weight_all[common_cells, common_cells])[col],
                direction = ifelse(sign_all_w[consistent_weight] > 0, "UP in MC+MS+MG", "DOWN in MC+MS+MG"),
                all = diff_weight_all[common_cells, common_cells][consistent_weight],
                mgh = diff_weight_mgh[consistent_weight],
                mdacc = diff_weight_mdacc[consistent_weight]) |>
  print()










library(openxlsx)

# Build count dataframe
count_df <- data.frame(
  sender = rep(rownames(diff_count_all[common_cells, common_cells]), times = length(common_cells)),
  receiver = rep(colnames(diff_count_all[common_cells, common_cells]), each = length(common_cells)),
  diff_count_all = as.vector(diff_count_all[common_cells, common_cells]),
  diff_count_mgh = as.vector(diff_count_mgh),
  diff_count_mdacc = as.vector(diff_count_mdacc),
  consistent = as.vector(consistent_count),
  direction = ifelse(as.vector(sign_all_sub) > 0, "UP in MC+MS+MG", 
                     ifelse(as.vector(sign_all_sub) < 0, "DOWN in MC+MS+MG", "no change"))
)

# Build weight dataframe
weight_df <- data.frame(
  sender = rep(rownames(diff_weight_all[common_cells, common_cells]), times = length(common_cells)),
  receiver = rep(colnames(diff_weight_all[common_cells, common_cells]), each = length(common_cells)),
  diff_weight_all = as.vector(diff_weight_all[common_cells, common_cells]),
  diff_weight_mgh = as.vector(diff_weight_mgh),
  diff_weight_mdacc = as.vector(diff_weight_mdacc),
  consistent = as.vector(consistent_weight),
  direction = ifelse(as.vector(sign_all_w) > 0, "UP in MC+MS+MG",
                     ifelse(as.vector(sign_all_w) < 0, "DOWN in MC+MS+MG", "no change"))
)

# Write to excel
wb <- createWorkbook()
addWorksheet(wb, "Count")
addWorksheet(wb, "Weight")
writeData(wb, "Count", count_df)
writeData(wb, "Weight", weight_df)

# Highlight consistent rows in yellow
consistent_count_rows <- which(count_df$consistent) + 1  # +1 for header
consistent_weight_rows <- which(weight_df$consistent) + 1

addStyle(wb, "Count", style = createStyle(fgFill = "#FFFF00"), 
         rows = consistent_count_rows, cols = 1:ncol(count_df), gridExpand = TRUE)
addStyle(wb, "Weight", style = createStyle(fgFill = "#FFFF00"), 
         rows = consistent_weight_rows, cols = 1:ncol(weight_df), gridExpand = TRUE)

saveWorkbook(wb, "cellchat_interactions.xlsx", overwrite = TRUE)









# ── BUBBLE PLOTS: Fibroblast → T_NK ──────────────────────────────────────────

pdf("bubble_fibro_tnk_all.pdf", width = 5, height = 8)
netVisual_bubble(cc_merged, sources.use = "Fibroblast_stromal", targets.use = "T_NK_lymphocytes",
                 comparison = c(1,2), angle.x = 45, remove.isolate = TRUE,
                 title.name = "Fibroblast → T_NK (All)")
dev.off()

pdf("bubble_fibro_tnk_mgh.pdf", width = 5, height = 8)
netVisual_bubble(cc_merged_mgh, sources.use = "Fibroblast_stromal", targets.use = "T_NK_lymphocytes",
                 comparison = c(1,2), angle.x = 45, remove.isolate = TRUE,
                 title.name = "Fibroblast → T_NK (MGH)")
dev.off()

pdf("bubble_fibro_tnk_mdacc.pdf", width = 5, height = 8)
netVisual_bubble(cc_merged_mdacc, sources.use = "Fibroblast_stromal", targets.use = "T_NK_lymphocytes",
                 comparison = c(1,2), angle.x = 45, remove.isolate = TRUE,
                 title.name = "Fibroblast → T_NK (MDACC)")
dev.off()





library(openxlsx)

# ── EXTRACT LR PAIRS ─────────────────────────────────────────────────────────

# Get LR communication dataframes
lr_MC <- subsetCommunication(cc_MC)
lr_MCMSMG <- subsetCommunication(cc_MCMSMG)
lr_MC_mgh <- subsetCommunication(cc_MC_mgh)
lr_MCMSMG_mgh <- subsetCommunication(cc_MCMSMG_mgh)
lr_MC_mdacc <- subsetCommunication(cc_MC_mdacc)
lr_MCMSMG_mdacc <- subsetCommunication(cc_MCMSMG_mdacc)

# Add condition label
lr_MC$condition <- "MC"
lr_MCMSMG$condition <- "MC+MS+MG"
lr_MC_mgh$condition <- "MC"
lr_MCMSMG_mgh$condition <- "MC+MS+MG"
lr_MC_mdacc$condition <- "MC"
lr_MCMSMG_mdacc$condition <- "MC+MS+MG"

# Combine per comparison
lr_all <- rbind(lr_MC, lr_MCMSMG)
lr_mgh <- rbind(lr_MC_mgh, lr_MCMSMG_mgh)
lr_mdacc <- rbind(lr_MC_mdacc, lr_MCMSMG_mdacc)

# Create LR pair identifier
lr_all$lr_pair <- paste(lr_all$source, lr_all$target, lr_all$ligand, lr_all$receptor, sep = "_")
lr_mgh$lr_pair <- paste(lr_mgh$source, lr_mgh$target, lr_mgh$ligand, lr_mgh$receptor, sep = "_")
lr_mdacc$lr_pair <- paste(lr_mdacc$source, lr_mdacc$target, lr_mdacc$ligand, lr_mdacc$receptor, sep = "_")

# Find LR pairs present in MC+MS+MG across all three comparisons
mcmsmg_pairs_all <- lr_all$lr_pair[lr_all$condition == "MC+MS+MG"]
mcmsmg_pairs_mgh <- lr_mgh$lr_pair[lr_mgh$condition == "MC+MS+MG"]
mcmsmg_pairs_mdacc <- lr_mdacc$lr_pair[lr_mdacc$condition == "MC+MS+MG"]

consistent_lr <- Reduce(intersect, list(mcmsmg_pairs_all, mcmsmg_pairs_mgh, mcmsmg_pairs_mdacc))
cat("Consistent LR pairs across all 3:", length(consistent_lr), "\n")

lr_all$consistent_all3 <- lr_all$lr_pair %in% consistent_lr
lr_mgh$consistent_all3 <- lr_mgh$lr_pair %in% consistent_lr
lr_mdacc$consistent_all3 <- lr_mdacc$lr_pair %in% consistent_lr



# ── WRITE EXCEL ──────────────────────────────────────────────────────────────

wb <- createWorkbook()
highlight <- createStyle(fgFill = "#FFFF00")

# LR sheets
addWorksheet(wb, "LR_All")
addWorksheet(wb, "LR_MGH")
addWorksheet(wb, "LR_MDACC")

writeData(wb, "LR_All", lr_all)
writeData(wb, "LR_MGH", lr_mgh)
writeData(wb, "LR_MDACC", lr_mdacc)

# Highlight consistent rows
consistent_rows_all <- which(lr_all$consistent_all3) + 1
consistent_rows_mgh <- which(lr_mgh$consistent_all3) + 1
consistent_rows_mdacc <- which(lr_mdacc$consistent_all3) + 1

if (length(consistent_rows_all) > 0)
  addStyle(wb, "LR_All", style = highlight, rows = consistent_rows_all, cols = 1:ncol(lr_all), gridExpand = TRUE)
if (length(consistent_rows_mgh) > 0)
  addStyle(wb, "LR_MGH", style = highlight, rows = consistent_rows_mgh, cols = 1:ncol(lr_mgh), gridExpand = TRUE)
if (length(consistent_rows_mdacc) > 0)
  addStyle(wb, "LR_MDACC", style = highlight, rows = consistent_rows_mdacc, cols = 1:ncol(lr_mdacc), gridExpand = TRUE)

saveWorkbook(wb, "cellchat_LR.xlsx", overwrite = TRUE)
cat("Saved cellchat_LR.xlsx\n")



library(ggplot2)
## CXCL12/CXCR4 bar plot (exploratory, not a manuscript panel). It needs the
## fibro and tnk subsets, which were created interactively, so it is skipped
## unless they exist in the session.
if (exists("fibro") && exists("tnk")) {

library(dplyr)
library(patchwork)

# Function to compute mean and SE
summarize_expr <- function(seurat_obj, gene, group_col = "m_group", inst_col = "institution") {
  expr <- GetAssayData(seurat_obj, layer = "data")[gene, ]
  df <- data.frame(
    expression = expr,
    m_group = seurat_obj@meta.data[[group_col]],
    institution = seurat_obj@meta.data[[inst_col]]
  )
  
  # Per institution
  per_inst <- df %>%
    group_by(institution, m_group) %>%
    summarize(mean = mean(expression),
              se = sd(expression) / sqrt(n()),
              .groups = "drop") %>%
    mutate(group_label = paste(institution, m_group, sep = "\n"))
  
  # Merged
  merged <- df %>%
    group_by(m_group) %>%
    summarize(mean = mean(expression),
              se = sd(expression) / sqrt(n()),
              .groups = "drop") %>%
    mutate(institution = "All", 
           group_label = paste("All", m_group, sep = "\n"))
  
  bind_rows(per_inst, merged) %>%
    mutate(institution = factor(institution, levels = c("mgh", "mdacc", "All")))
}

cxcl12_sum <- summarize_expr(fibro, "CXCL12")
cxcr4_sum <- summarize_expr(tnk, "CXCR4")

# Set x-axis order: mgh MC, mgh MC+MS+MG, mdacc MC, mdacc MC+MS+MG, All MC, All MC+MS+MG
x_order <- c("mgh\nMC", "mgh\nMC+MS+MG", "mdacc\nMC", "mdacc\nMC+MS+MG", "All\nMC", "All\nMC+MS+MG")
cxcl12_sum$group_label <- factor(cxcl12_sum$group_label, levels = x_order)
cxcr4_sum$group_label <- factor(cxcr4_sum$group_label, levels = x_order)

colors <- c("MC" = "#4393C3", "MC+MS+MG" = "#D6604D")

p1 <- ggplot(cxcl12_sum, aes(x = group_label, y = mean, fill = m_group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = colors) +
  ggtitle("CXCL12 in Fibroblast_stromal") +
  ylab("Mean normalized expression") + xlab("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "dashed", color = "grey50")

p2 <- ggplot(cxcr4_sum, aes(x = group_label, y = mean, fill = m_group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = colors) +
  ggtitle("CXCR4 in T_NK_lymphocytes") +
  ylab("Mean normalized expression") + xlab("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "dashed", color = "grey50")

pdf("barplot_CXCL12_CXCR4.pdf", width = 12, height = 6)
p1 + p2
dev.off()

}  # end CXCL12/CXCR4 bar plot







# Save all CellChat objects
saveRDS(cc_MC, "cc_MC.rds")
saveRDS(cc_MCMSMG, "cc_MCMSMG.rds")
saveRDS(cc_merged, "cc_merged.rds")

saveRDS(cc_MC_mgh, "cc_MC_mgh.rds")
saveRDS(cc_MCMSMG_mgh, "cc_MCMSMG_mgh.rds")
saveRDS(cc_merged_mgh, "cc_merged_mgh.rds")

saveRDS(cc_MC_mdacc, "cc_MC_mdacc.rds")
saveRDS(cc_MCMSMG_mdacc, "cc_MCMSMG_mdacc.rds")
saveRDS(cc_merged_mdacc, "cc_merged_mdacc.rds")






# (the S2F panel plots follow in Part 8b)

## end of Part 8a: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 8b - Fig S2F (part b): circle plots and differential-strength heatmaps
##############################################################################

# ============================================================================
# Fig S2F - CellChat, MC+MS+MG vs MC: differential interaction circle plots
# (merged, MDACC, MGH; top row) and differential interaction-strength
# heatmaps (bottom row).
#
# Figure panels : Fig S2F
# Uses the CellChat objects made in Part 8a.
# ============================================================================

setwd(file.path(TMOS_DIR, "Single_cell/Heart/Seurat_by_CellType/Combined/updated_cellchat_institute"))

library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ── Load merged CellChat objects ──────────────────────────────────────────────
cc_merged       <- readRDS("cc_merged.rds")
cc_merged_mgh   <- readRDS("cc_merged_mgh.rds")
cc_merged_mdacc <- readRDS("cc_merged_mdacc.rds")

# ---------------- top row: differential interaction circle plots ----------------

## Merged (MDACC + MGH)
pdf("circle_diff_all.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength (MC+MS+MG vs MC)")
dev.off()

## MDACC
pdf("circle_diff_mdacc.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged_mdacc, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions - MDACC (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged_mdacc, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength - MDACC (MC+MS+MG vs MC)")
dev.off()

## MGH
pdf("circle_diff_mgh.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cc_merged_mgh, weight.scale = TRUE, measure = "count", title.name = "Differential Number of Interactions - MGH (MC+MS+MG vs MC)")
netVisual_diffInteraction(cc_merged_mgh, weight.scale = TRUE, measure = "weight", title.name = "Differential Interaction Strength - MGH (MC+MS+MG vs MC)")
dev.off()

# ---------------- bottom row: differential interaction-strength heatmaps ----------------

# ── Compute difference matrices (MC+MS+MG minus MC) ──────────────────────────
diff_weight_all   <- cc_merged@net[[2]]$weight       - cc_merged@net[[1]]$weight
diff_weight_mgh   <- cc_merged_mgh@net[[2]]$weight   - cc_merged_mgh@net[[1]]$weight
diff_weight_mdacc <- cc_merged_mdacc@net[[2]]$weight - cc_merged_mdacc@net[[1]]$weight

# ── Restrict to cell types present in ALL THREE comparisons ──────────────────
common_cells <- Reduce(intersect, list(
  rownames(diff_weight_all),
  rownames(diff_weight_mgh),
  rownames(diff_weight_mdacc)
))
cat("Common cell types across all 3:", paste(common_cells, collapse = ", "), "\n")

diff_weight_all   <- diff_weight_all[common_cells, common_cells]
diff_weight_mgh   <- diff_weight_mgh[common_cells, common_cells]
diff_weight_mdacc <- diff_weight_mdacc[common_cells, common_cells]

# ── Shared symmetric color scale (blue - white - red) ────────────────────────
weight_limit <- max(abs(c(diff_weight_all, diff_weight_mgh, diff_weight_mdacc)), na.rm = TRUE)
cat("Color scale limit (+/-):", weight_limit, "\n")

col_fun <- colorRamp2(
  c(-weight_limit, 0, weight_limit),
  c("#2166AC", "white", "#D6604D")
)

# ── Helper to build each heatmap with a unique name ──────────────────────────
make_ht <- function(mat, title, ht_name) {
  Heatmap(
    matrix            = mat,
    name              = ht_name,          # unique per heatmap -> no duplicate warning
    col               = col_fun,
    column_title      = title,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_names_side    = "left",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = 9),
    column_names_gp   = gpar(fontsize = 9),
    show_heatmap_legend = FALSE,
    border              = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(round(mat[i, j], 2), x, y, gp = gpar(fontsize = 7))
    }
  )
}

ht_all   <- make_ht(diff_weight_all,   "All\n(MC+MS+MG vs MC)",   "All")
ht_mgh   <- make_ht(diff_weight_mgh,   "MGH\n(MC+MS+MG vs MC)",   "MGH")
ht_mdacc <- make_ht(diff_weight_mdacc, "MDACC\n(MC+MS+MG vs MC)", "MDACC")

# ── Shared legend (ASCII-safe title) ─────────────────────────────────────────
lgd <- Legend(
  col_fun   = col_fun,
  title     = "D Strength\n(MC+MS+MG - MC)",
  direction = "vertical"
)

# ── Save PDF ──────────────────────────────────────────────────────────────────
pdf("diff_strength_heatmaps_3way.pdf", width = 18, height = 7)
draw(
  ht_all + ht_mgh + ht_mdacc,
  annotation_legend_list = list(lgd),
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  padding                = unit(c(5, 5, 5, 5), "mm")
)
dev.off()

cat("Saved: diff_strength_heatmaps_3way.pdf\n")

## end of Part 8b: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())


##############################################################################
##  PART 9 - Fig S2J: benchmark of the IL-17 programs vs other inflammatory programs (AUC)
##############################################################################

# ============================================================================
# Benchmark of IL-17 signaling-related programs vs other inflammatory programs (donor-level AUC); IL-17 source/receiver
#
# Figure panels : Fig S2J
# ============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(msigdbr)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.2_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
col_mc   <- "#4682B4"
col_tmos <- "#b30000"

nm_core <- c("TAB2","MAPK8","UBE2V1","NFKB1","MAP2K4","BTRC","CREB1",
             "PPP2CB","RPS6KA5","MEF2A","RPS6KA3","FBXW11","MAP2K1","DUSP4")

## ---- Gene sets (Hallmark sets as in R1.3, minus the 14 core genes) ----
h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
hall <- function(nm) setdiff(unique(h$gene_symbol[h$gs_name == nm]), nm_core)

progs <- list(
  `IL-17 signaling–related` = nm_core,
  `Cytotoxic CD8`           = c("GZMA","GZMB","GZMH","PRF1","GNLY","NKG7"),
  `IFN-γ response`          = hall("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  `TNF-α/NF-κB`             = hall("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  `IL-6/JAK-STAT3`          = hall("HALLMARK_IL6_JAK_STAT3_SIGNALING")
)
sapply(progs, length)

## ---- Data: heart T/NK, control + MC + TMOS ----
heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"
tnk <- subset(heart, subset = major_cell_type == "T_NK_lymphocytes" &
                m_group %in% c("control", "MC", "MC+MS+MG"))
rm(heart); gc()

tnk <- NormalizeData(tnk)
for (i in seq_along(progs)) {
  g <- intersect(progs[[i]], rownames(tnk))
  tnk <- AddModuleScore(tnk, features = list(g), name = paste0("P", i, "_"), seed = 1)
}
md <- tnk@meta.data
rm(tnk); gc()

## ---- Donor-level scores (>=10 T/NK cells) ----
sc_cols <- paste0("P", seq_along(progs), "_1")
donor <- md %>%
  dplyr::mutate(grp = factor(dplyr::recode(as.character(m_group), "MC+MS+MG" = "TMOS"),
                             levels = c("control", "MC", "TMOS"))) %>%
  dplyr::group_by(sample_id, grp, institution) %>%
  dplyr::summarise(n_cells = dplyr::n(), dplyr::across(dplyr::all_of(sc_cols), mean), .groups = "drop") %>%
  dplyr::filter(n_cells >= 10) %>%
  tidyr::pivot_longer(dplyr::all_of(sc_cols), names_to = "prog", values_to = "score") %>%
  dplyr::mutate(prog = names(progs)[as.integer(sub("P(\\d+)_1", "\\1", prog))])

## site-referenced: score minus median of same-institution controls (as in R2.4)
ctrl <- donor %>% dplyr::filter(grp == "control") %>%
  dplyr::group_by(prog, institution) %>% dplyr::summarise(ctrl = median(score), .groups = "drop")
donor <- donor %>% dplyr::left_join(ctrl, by = c("prog", "institution")) %>%
  dplyr::mutate(score_adj = score - ctrl)

dd <- donor %>% dplyr::filter(grp %in% c("MC", "TMOS")) %>% droplevels()
table(dplyr::distinct(dd, sample_id, grp)$grp)

## ---- Stats: Wilcoxon p, AUC (TMOS > MC), for raw and site-referenced ----
auc_fun <- function(x, g) {
  w <- wilcox.test(x[g == "TMOS"], x[g == "MC"], exact = FALSE)
  c(p = w$p.value, auc = unname(w$statistic) / (sum(g == "TMOS") * sum(g == "MC")))
}
bench <- dd %>% dplyr::group_by(prog) %>%
  dplyr::summarise(p_raw = auc_fun(score, grp)["p"],         auc_raw = auc_fun(score, grp)["auc"],
                   p_adj = auc_fun(score_adj, grp)["p"],     auc_adj = auc_fun(score_adj, grp)["auc"],
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(auc_adj))
as.data.frame(bench)

## ---- Plot: per-donor scores per program (site-referenced) ----
dd$prog <- factor(dd$prog, levels = bench$prog)
lab <- bench %>% dplyr::mutate(prog = factor(prog, levels = bench$prog),
                               lab = sprintf("AUC = %.2f\np = %.3f", auc_adj, p_adj))
p_bench <- ggplot(dd, aes(grp, score_adj, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.3, size = 3) +
  facet_wrap(~ prog, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Program score per donor\n(relative to site-matched controls)",
       title = "Heart T/NK cells: IL-17 signaling–related program vs other inflammatory programs") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
p_bench
ggsave(file.path(out_dir, "TNK_program_benchmark_TMOS_vs_MC.pdf"), p_bench, width = 13, height = 4.5, device = cairo_pdf)

write_xlsx(list(Stats = bench, Donor_scores = donor,
                Gene_sets = data.frame(prog = rep(names(progs), lengths(progs)), gene = unlist(progs))),
           file.path(out_dir, "TNK_program_benchmark_TMOS_vs_MC.xlsx"))






## ---- 43-gene IL-17 signaling–related module (Reactome IL-17 ∩ Xenium panel) ----
mapk38 <- c("NOD1","TAB1","CHUK","CREB1","ATF2","MAPK14","DUSP6","ELK1","TAB2","FBXW11",
            "IKBKB","IRAK1","IRAK2","MEF2A","MEF2C","ATF1","NFKB1","PPP2CA","PPP2R1A","MAPK1",
            "MAPK3","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10","MAP2K1","NOD2","MAP2K4","SKP1",
            "MAP3K7","TRAF6","CUL1","IKBKG","RIPK2","BTRC","RPS6KA5","MAPKAPK2")
il17_43 <- c(mapk38, "IL17A","IL17F","IL25","IL17RA","IL17RB")
length(il17_43)                       # 43
held_out <- setdiff(il17_43, nm_core) # 43 minus the 14 DE-derived core genes
length(held_out)                      # 33

heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"
tnk <- subset(heart, subset = major_cell_type == "T_NK_lymphocytes" &
                m_group %in% c("control", "MC", "MC+MS+MG"))
rm(heart); gc()

tnk <- NormalizeData(tnk)
g43 <- intersect(il17_43,  rownames(tnk)); cat(length(g43), "of 43 found\n")
gho <- intersect(held_out, rownames(tnk)); cat(length(gho), "of 33 found\n")

tnk <- AddModuleScore(tnk, features = list(g43), name = "M43_", seed = 1)
tnk <- AddModuleScore(tnk, features = list(gho), name = "MHO_", seed = 1)
md43 <- tnk@meta.data
rm(tnk); gc()

## ---- Donor-level for both new sets ----
d43 <- md43 %>%
  dplyr::mutate(grp = factor(dplyr::recode(as.character(m_group), "MC+MS+MG" = "TMOS"),
                             levels = c("control", "MC", "TMOS"))) %>%
  dplyr::group_by(sample_id, grp, institution) %>%
  dplyr::summarise(n_cells = dplyr::n(),
                   `IL-17 signaling–related (43-gene, spatial)`     = mean(M43_1),
                   `IL-17 signaling–related (33-gene, held-out)`    = mean(MHO_1), .groups = "drop") %>%
  dplyr::filter(n_cells >= 10) %>%
  tidyr::pivot_longer(dplyr::starts_with("IL-17"), names_to = "prog", values_to = "score")

ctrl43 <- d43 %>% dplyr::filter(grp == "control") %>%
  dplyr::group_by(prog, institution) %>% dplyr::summarise(ctrl = median(score), .groups = "drop")
d43 <- d43 %>% dplyr::left_join(ctrl43, by = c("prog", "institution")) %>%
  dplyr::mutate(score_adj = score - ctrl)





## ---- Combine with existing benchmark ----
donor_all <- dplyr::bind_rows(donor, d43)
dd_all <- donor_all %>% dplyr::filter(grp %in% c("MC", "TMOS")) %>% droplevels()

n_genes_all <- c(sapply(progs, length),
                 `IL-17 signaling–related (43-gene, spatial)`  = length(g43),
                 `IL-17 signaling–related (33-gene, held-out)` = length(gho))

bench_all <- dd_all %>% dplyr::group_by(prog) %>%
  dplyr::summarise(p_raw = auc_fun(score, grp)["p"],     auc_raw = auc_fun(score, grp)["auc"],
                   p_adj = auc_fun(score_adj, grp)["p"], auc_adj = auc_fun(score_adj, grp)["auc"],
                   .groups = "drop") %>%
  dplyr::mutate(n_genes = n_genes_all[prog]) %>%
  dplyr::arrange(dplyr::desc(auc_adj))
as.data.frame(bench_all)

## ---- Plot ----
prog_order <- c("IL-17 signaling–related (43-gene, spatial)",
                "IL-17 signaling–related (33-gene, held-out)",
                "Cytotoxic CD8", "IFN-γ response", "TNF-α/NF-κB", "IL-6/JAK-STAT3",
                "IL-17 signaling–related")
dd_all$prog <- factor(dd_all$prog, levels = prog_order)
lab_all <- bench_all %>% dplyr::mutate(prog = factor(prog, levels = prog_order),
                                       lab = sprintf("AUC = %.2f\np = %.3f", auc_adj, p_adj))

p_bench2 <- ggplot(dd_all, aes(grp, score_adj, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab_all, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.3, size = 3) +
  facet_wrap(~ prog, nrow = 1, scales = "free_y",
             labeller = labeller(prog = function(x) gsub("IL-17 signaling–related", "IL-17 sig.–related", x))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Program score per donor\n(relative to site-matched controls)",
       title = "Heart T/NK cells: IL-17 signaling–related program vs other inflammatory programs") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold"))
p_bench2
ggsave(file.path(out_dir, "TNK_program_benchmark_with43_TMOS_vs_MC.pdf"), p_bench2,
       width = 18, height = 4.5, device = cairo_pdf)

## ---- Save ----
write_xlsx(list(Stats = bench_all, Donor_scores = donor_all,
                Gene_sets = data.frame(
                  prog = c(rep(names(progs), lengths(progs)),
                           rep("IL-17 signaling–related (43-gene, spatial)",  length(g43)),
                           rep("IL-17 signaling–related (33-gene, held-out)", length(gho))),
                  gene = c(unlist(progs), g43, gho))),
           file.path(out_dir, "TNK_program_benchmark_with43_TMOS_vs_MC.xlsx"))








library(patchwork)

drop_ho <- "IL-17 signaling–related (33-gene, held-out)"
dd_all    <- dd_all    %>% dplyr::filter(prog != drop_ho) %>% droplevels()
bench_all <- bench_all %>% dplyr::filter(prog != drop_ho)

nice <- c(`IL-17 signaling–related`                    = "IL-17 sig.–related (14-gene core)",
          `IL-17 signaling–related (43-gene, spatial)` = "IL-17 sig.–related (43-gene module)",
          `Cytotoxic CD8` = "Cytotoxic CD8", `IFN-γ response` = "IFN-γ response",
          `TNF-α/NF-κB` = "TNF-α/NF-κB", `IL-6/JAK-STAT3` = "IL-6/JAK-STAT3")

## ---- A: AUC per program (raw and site-referenced) ----
auc_df <- bench_all %>%
  dplyr::select(prog, auc_raw, auc_adj) %>%
  tidyr::pivot_longer(c(auc_raw, auc_adj), names_to = "type", values_to = "auc") %>%
  dplyr::mutate(type = ifelse(type == "auc_raw", "Raw", "Site-referenced"),
                prog = factor(nice[prog], levels = rev(nice[bench_all$prog[order(bench_all$auc_raw)]])),
                il17 = grepl("IL-17", prog))

pA_auc <- ggplot(auc_df, aes(auc, prog)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_line(aes(group = prog), colour = "grey70") +
  geom_point(aes(shape = type, colour = il17), size = 3.5) +
  scale_colour_manual(values = c(`TRUE` = col_tmos, `FALSE` = "grey40"), guide = "none") +
  scale_shape_manual(values = c(Raw = 16, `Site-referenced` = 17), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "AUC (TMOS vs MC, donor level)", y = NULL,
       title = "Separation of TMOS from MC by inflammatory program") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold"), legend.position = "top")

## ---- B: per-donor scores (raw), AUC labels only ----
prog_order <- names(nice)[order(match(names(nice), bench_all$prog[order(-bench_all$auc_raw)]))]
dd_all$prog <- factor(dd_all$prog, levels = prog_order)
lab_all <- bench_all %>% dplyr::mutate(prog = factor(prog, levels = prog_order),
                                       lab = sprintf("AUC = %.2f", auc_raw))

pB_box <- ggplot(dd_all, aes(grp, score, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab_all, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.5, size = 3.2) +
  facet_wrap(~ prog, nrow = 1, scales = "free_y", labeller = labeller(prog = nice)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Program score\n(mean per donor)", title = "Program scores per donor") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold"))

fig_bench <- (pA_auc | plot_spacer()) / pB_box +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
fig_bench[[1]] <- fig_bench[[1]] + plot_layout(widths = c(1.2, 1))

ggsave(file.path(out_dir, "TNK_program_benchmark_AUC.pdf"), fig_bench, width = 16, height = 8.5, device = cairo_pdf)
ggsave(file.path(out_dir, "TNK_program_benchmark_AUC.png"), fig_bench, width = 16, height = 8.5, dpi = 300, bg = "white")

write_xlsx(list(Stats = bench_all, Donor_scores = dplyr::filter(donor_all, prog != drop_ho)),
           file.path(out_dir, "TNK_program_benchmark_AUC.xlsx"))








library(patchwork)

drop_ho <- "IL-17 signaling–related (33-gene, held-out)"
dd_all    <- dd_all    %>% dplyr::filter(prog != drop_ho) %>% droplevels()
bench_all <- bench_all %>% dplyr::filter(prog != drop_ho)

nice <- c(`IL-17 signaling–related`                    = "IL-17 sig.–related (14-gene core)",
          `IL-17 signaling–related (43-gene, spatial)` = "IL-17 sig.–related (43-gene module)",
          `IL-6/JAK-STAT3` = "IL-6/JAK-STAT3", `TNF-α/NF-κB` = "TNF-α/NF-κB",
          `IFN-γ response` = "IFN-γ response", `Cytotoxic CD8` = "Cytotoxic CD8")

## order programs by site-referenced AUC
ord <- bench_all %>% dplyr::arrange(dplyr::desc(auc_adj)) %>% dplyr::pull(prog)
leg <- setNames(sprintf("%s (AUC %.2f)", nice[ord], bench_all$auc_adj[match(ord, bench_all$prog)]), ord)

## ---- A: ROC curves (site-referenced scores) ----
roc_df <- dd_all %>%
  dplyr::group_by(prog) %>%
  dplyr::group_modify(~ {
    s <- .x$score_adj; y <- .x$grp == "TMOS"
    th <- c(Inf, sort(unique(s), decreasing = TRUE))
    data.frame(fpr = sapply(th, function(t) mean(s[!y] >= t)),
               tpr = sapply(th, function(t) mean(s[y]  >= t)))
  }) %>% dplyr::ungroup() %>%
  dplyr::mutate(prog = factor(prog, levels = ord))

roc_cols <- setNames(c("#b30000", "#e8743b", "grey15", "grey40", "grey60", "grey75"), ord)
roc_lty  <- setNames(c("solid", "solid", "dashed", "dotdash", "dotted", "longdash"), ord)
roc_cols[ord[!grepl("IL-17", ord)]] <- c("grey15", "grey40", "grey60", "grey75")
roc_cols[grepl("IL-17", ord)] <- c("#b30000", "#e8743b")[seq_len(sum(grepl("IL-17", ord)))]
roc_lty[grepl("IL-17", ord)] <- "solid"

pA_roc <- ggplot(roc_df, aes(fpr, tpr, colour = prog, linetype = prog)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80") +
  geom_step(direction = "vh", linewidth = 1) +
  scale_colour_manual(values = roc_cols, labels = leg, name = NULL) +
  scale_linetype_manual(values = roc_lty, labels = leg, name = NULL) +
  coord_equal() +
  labs(x = "False positive rate (MC donors)", y = "True positive rate (TMOS donors)",
       title = "TMOS vs MC donor classification by program",
       subtitle = "Program scores relative to site-matched controls") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "right", legend.key.width = unit(1.4, "cm"))

## ---- B: per-donor site-referenced scores ----
dd_all$prog <- factor(dd_all$prog, levels = ord)
lab_all <- bench_all %>% dplyr::mutate(prog = factor(prog, levels = ord),
                                       lab = sprintf("AUC = %.2f", auc_adj))

pB_box <- ggplot(dd_all, aes(grp, score_adj, fill = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab_all, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.5, size = 3.2) +
  facet_wrap(~ prog, nrow = 1, scales = "free_y", labeller = labeller(prog = nice)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Program score per donor\n(relative to site-matched controls)",
       title = "Program scores per donor") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold"))

fig_bench <- pA_roc / pB_box +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(file.path(out_dir, "TNK_program_benchmark_ROC.pdf"), fig_bench, width = 15, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "TNK_program_benchmark_ROC.png"), fig_bench, width = 15, height = 10, dpi = 300, bg = "white")

write_xlsx(list(Stats = bench_all, ROC = roc_df,
                Donor_scores = dplyr::filter(donor_all, prog != drop_ho)),
           file.path(out_dir, "TNK_program_benchmark_ROC.xlsx"))













library(patchwork)

drop_ho <- "IL-17 signaling–related (33-gene, held-out)"
dd_all    <- dd_all    %>% dplyr::filter(prog != drop_ho) %>% droplevels()
bench_all <- bench_all %>% dplyr::filter(prog != drop_ho)

nice <- c(`IL-17 signaling–related`                    = "IL-17 sig.–related\n(14-gene core)",
          `IL-17 signaling–related (43-gene, spatial)` = "IL-17 sig.–related\n(43-gene module)",
          `IL-6/JAK-STAT3` = "IL-6/JAK-STAT3", `TNF-α/NF-κB` = "TNF-α/NF-κB",
          `IFN-γ response` = "IFN-γ response", `Cytotoxic CD8` = "Cytotoxic CD8")

ord <- bench_all %>% dplyr::arrange(dplyr::desc(auc_adj)) %>% dplyr::pull(prog)

## ---- A: AUC bar plot (control-corrected scores) ----
auc_bar <- bench_all %>%
  dplyr::mutate(label = factor(nice[prog], levels = rev(nice[ord])),
                il17  = grepl("IL-17", prog))

pA_auc <- ggplot(auc_bar, aes(auc_adj, label, fill = il17)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.2f", auc_adj)), hjust = -0.2, size = 3.5) +
  scale_fill_manual(values = c(`TRUE` = col_tmos, `FALSE` = "grey60"), guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "AUC, TMOS vs MC (donor level; 0.5 = no separation)", y = NULL,
       title = "Separation of TMOS from MC by program") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold"))

## ---- B: per-donor control-corrected scores ----
dd_all$prog <- factor(dd_all$prog, levels = ord)
lab_all <- bench_all %>% dplyr::mutate(prog = factor(prog, levels = ord),
                                       lab = sprintf("AUC = %.2f", auc_adj))

pB_box <- ggplot(dd_all, aes(grp, score_adj, fill = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab_all, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.5, size = 3.2) +
  facet_wrap(~ prog, nrow = 2, scales = "free_y", labeller = labeller(prog = nice)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "Program score per donor\n(relative to site-matched controls)",
       title = "Program scores per donor") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold"))

fig_bench <- (pA_auc | pB_box) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(file.path(out_dir, "TNK_program_benchmark.pdf"), fig_bench, width = 15, height = 6.5, device = cairo_pdf)
ggsave(file.path(out_dir, "TNK_program_benchmark.png"), fig_bench, width = 15, height = 6.5, dpi = 300, bg = "white")

write_xlsx(list(Stats = bench_all, Donor_scores = dplyr::filter(donor_all, prog != drop_ho)),
           file.path(out_dir, "TNK_program_benchmark.xlsx"))
















library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(writexl)

base_dir <- file.path(TMOS_DIR, "Clean/reviewer_comments")
out_dir  <- file.path(base_dir, "R2.2_source_receiver")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
col_mc   <- "#4682B4"
col_tmos <- "#b30000"

src_genes <- c("IL17A", "IL17F", "RORC", "IL23R", "CCR6", "KLRB1")
rec_genes <- c("IL17RA", "IL17RC", "TRAF3IP2")
resp_set  <- c("CXCL1", "CXCL2", "CXCL8", "CCL20", "IL6", "CSF3", "LCN2", "S100A8", "S100A9", "MMP3")
receivers <- c("Myeloid_cells", "Fibroblast_stromal", "Endothelial_cells", "Smooth_muscle_mural", "Cardiomyocytes")

bt <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        plot.title = element_text(face = "bold"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

heart <- readRDS(file.path(base_dir, "heart_3232026.rds"))
DefaultAssay(heart) <- "RNA"
heart <- subset(heart, subset = m_group %in% c("control", "MC", "MC+MS+MG"))
heart$grp <- factor(dplyr::recode(as.character(heart$m_group), "MC+MS+MG" = "TMOS"),
                    levels = c("control", "MC", "TMOS"))
heart <- NormalizeData(heart)

## =====================================================================
## A. Source: type-17 machinery by Monaco T/NK subset (MC + TMOS)
## =====================================================================
tnk_meta <- readRDS(file.path(base_dir, "R2.1_TNK_states", "heart_TNK_monaco_fine_meta.rds"))
tnk <- subset(heart, cells = intersect(rownames(tnk_meta), colnames(heart)))
tnk$monaco_fine <- tnk_meta[colnames(tnk), "monaco_fine"]
keep_lab <- names(which(table(tnk$monaco_fine) >= 20))
tnk <- subset(tnk, cells = colnames(tnk)[tnk$monaco_fine %in% keep_lab])
tnk$lab_grp <- paste(tnk$monaco_fine, tnk$grp, sep = " | ")

src_present <- intersect(src_genes, rownames(tnk))
src_missing <- setdiff(src_genes, rownames(tnk)); src_missing

src_tab <- FetchData(tnk, vars = c(src_present, "monaco_fine", "grp")) %>%
  tidyr::pivot_longer(dplyr::all_of(src_present), names_to = "gene", values_to = "expr") %>%
  dplyr::group_by(monaco_fine, grp, gene) %>%
  dplyr::summarise(pct = 100 * mean(expr > 0), avg = mean(expm1(expr)), n = dplyr::n(), .groups = "drop")

## add genes absent from the dataset as 0% / 0 expression
if (length(src_missing) > 0) {
  src_tab <- dplyr::bind_rows(
    src_tab,
    tidyr::expand_grid(dplyr::distinct(src_tab, monaco_fine, grp), gene = src_missing) %>%
      dplyr::mutate(pct = 0, avg = 0, n = NA_integer_))
}

src_tab <- src_tab %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(avg_z = if (sd(avg) > 0) as.numeric(scale(avg)) else 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = factor(gene, levels = src_genes))
rm(tnk); gc()

pA_src <- ggplot(src_tab, aes(gene, monaco_fine, size = pct, colour = avg_z)) +
  geom_point() +
  facet_wrap(~ grp) +
  scale_colour_gradient2(low = col_mc, mid = "grey90", high = col_tmos, name = "Scaled\nexpression") +
  scale_size_continuous(range = c(0, 6), name = "% cells") +
  labs(x = NULL, y = NULL, title = "Type-17 machinery by T/NK subset (IL-17 source capacity)") +
  bt + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

## =====================================================================
## B. Receiver: IL-17 receptor genes across heart cell types (MC + TMOS)
## =====================================================================
rec_tab <- FetchData(heart, vars = c(rec_genes, "major_cell_type", "grp")) %>%
  dplyr::filter(grp != "control") %>%
  tidyr::pivot_longer(dplyr::all_of(rec_genes), names_to = "gene", values_to = "expr") %>%
  dplyr::group_by(major_cell_type, grp, gene) %>%
  dplyr::summarise(pct = 100 * mean(expr > 0), avg = mean(expm1(expr)), .groups = "drop") %>%
  dplyr::group_by(gene) %>% dplyr::mutate(avg_z = as.numeric(scale(avg))) %>% dplyr::ungroup()

pB_rec <- ggplot(rec_tab, aes(gene, major_cell_type, size = pct, colour = avg_z)) +
  geom_point() +
  facet_wrap(~ grp) +
  scale_colour_gradient2(low = col_mc, mid = "grey90", high = col_tmos, name = "Scaled\nexpression") +
  scale_size_continuous(range = c(0, 6), name = "% cells") +
  labs(x = NULL, y = NULL, title = "IL-17 receptor machinery by heart cell type (receiver capacity)") +
  bt + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

## =====================================================================
## C. Receiver response: IL-17 target-gene score per donor, per receiver type
## =====================================================================
g_resp <- intersect(resp_set, rownames(heart)); g_resp
rx <- subset(heart, subset = major_cell_type %in% receivers)
rm(heart); gc()
rx <- AddModuleScore(rx, features = list(g_resp), name = "IL17resp_", seed = 1)
md_rx <- rx@meta.data
rm(rx); gc()

resp <- md_rx %>%
  dplyr::group_by(major_cell_type, sample_id, grp, institution) %>%
  dplyr::summarise(n_cells = dplyr::n(), score = mean(IL17resp_1), .groups = "drop") %>%
  dplyr::filter(n_cells >= 10)
ctrl <- resp %>% dplyr::filter(grp == "control") %>%
  dplyr::group_by(major_cell_type, institution) %>% dplyr::summarise(ctrl = median(score), .groups = "drop")
resp <- resp %>% dplyr::left_join(ctrl, by = c("major_cell_type", "institution")) %>%
  dplyr::mutate(score_adj = score - ctrl)

rr <- resp %>% dplyr::filter(grp %in% c("MC", "TMOS"), !is.na(score_adj)) %>% droplevels()
resp_stats <- rr %>% dplyr::group_by(major_cell_type) %>%
  dplyr::summarise(n_MC = sum(grp == "MC"), n_TMOS = sum(grp == "TMOS"),
                   auc = if (n_MC > 0 && n_TMOS > 0)
                     unname(wilcox.test(score_adj[grp == "TMOS"], score_adj[grp == "MC"], exact = FALSE)$statistic) / (n_MC * n_TMOS)
                   else NA_real_, .groups = "drop")
as.data.frame(resp_stats)

lab_r <- resp_stats %>% dplyr::mutate(lab = sprintf("AUC = %.2f\n(%d vs %d)", auc, n_TMOS, n_MC))
pC_resp <- ggplot(rr, aes(grp, score_adj, fill = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.55) +
  geom_point(aes(colour = grp), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = lab_r, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.3, size = 3) +
  facet_wrap(~ major_cell_type, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  scale_colour_manual(values = c(MC = col_mc, TMOS = col_tmos), guide = "none") +
  labs(x = NULL, y = "IL-17 response score per donor\n(relative to site-matched controls)",
       title = "IL-17 target-gene response in receiver cell types") +
  bt

## ---- Save ----
fig <- (pA_src | pB_rec) / pC_resp + plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
ggsave(file.path(out_dir, "R2.2_source_receiver.pdf"), fig, width = 17, height = 12, device = cairo_pdf)
ggsave(file.path(out_dir, "R2.2_source_receiver.png"), fig, width = 17, height = 12, dpi = 300, bg = "white")

write_xlsx(list(Source = src_tab, Receptor = rec_tab, Response_donor = resp,
                Response_stats = resp_stats,
                Gene_sets = data.frame(set = c(rep("source", length(src_genes)), rep("receptor", length(rec_genes)),
                                               rep("response", length(g_resp))),
                                       gene = c(src_genes, rec_genes, g_resp))),
           file.path(out_dir, "R2.2_source_receiver.xlsx"))

## end of Part 9: clear the workspace before the next part
rm(list = setdiff(ls(), "TMOS_DIR")); invisible(gc())
