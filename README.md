# IL-17/MAPK Spatial Axis in Triple M Overlap Syndrome

Analysis code for:

> **Single-Cell and Spatial Transcriptomics Identify an IL-17/MAPK Immune Program in ICI-Associated Triple M Overlap Syndrome**
> Abushukair H, *et al.* (Johnson DB, Naqash AR, Oh TG, co-senior). *Manuscript under review.*

Triple M Overlap Syndrome (TMOS) is immune checkpoint inhibitor (ICI) myocarditis together with myositis and myasthenia gravis. This repository contains the R code behind the main and supplementary figures. It covers reanalysis of public heart, skeletal muscle and PBMC scRNA-seq, Xenium 5K spatial transcriptomics of TMOS heart and skeletal muscle, and serum cytokine profiling.

## Group definitions

| Label | Definition |
|---|---|
| MC | ICI myocarditis alone |
| MCO | MC with myositis and/or myasthenia gravis (MC+MS, MC+MG, TMOS) |
| TMOS | MC + MS + MG (labelled `MC+MS+MG` in some object metadata) |

## Repository layout

```
config/paths.R          set the project folder once (TMOS_DIR)
install_packages.R      CRAN / Bioconductor / GitHub dependencies
fig1/ … fig4/           scripts that produce the figure panels (main + supplementary)
```

Each figure script starts with a header listing the panels it produces, the part that makes each one, the `.rds` objects it reads and writes, and the original analysis files.

## Running

1. Install dependencies: `Rscript install_packages.R`
2. Set `TMOS_DIR` in `config/paths.R` (or as an environment variable) to the folder holding the processed objects.
3. Run the figure scripts from the repository root. Each script is split into parts (`##  PART n`) that are self-contained; the objects are large, so run one part at a time in a fresh R session.

The figure scripts start from processed Seurat and Xenium objects (QC-filtered, integrated and annotated), listed below. Upstream processing (import, QC, integration, clustering and annotation) is not included.

## Figure-to-script map

| Panel(s) | Script |
|---|---|
| **Figure 1** | |
| 1A | Schematic (BioRender), no code |
| 1B | `fig1/Fig1_all_panels.R` (Part 1) |
| 1C, 1D, 1E | `fig1/Fig1_all_panels.R` (Part 2) |
| 1F | `fig1/Fig1_all_panels.R` (Part 3) |
| **Figure S1** | |
| S1A | Cohort Sankey diagram from sample metadata; code not included |
| S1B–S1F | `fig1/Fig1_all_panels.R` (Part 2) |
| S1G, S1H | `fig1/Fig1_all_panels.R` (Part 3; inputs made by `fig2/Fig2_all_panels.R` Part 3) |
| **Figure 2** | |
| 2A | Schematic, no code |
| 2B, 2C, 2D | `fig2/Fig2_all_panels.R` (Part 1) |
| 2E | `fig2/Fig2_all_panels.R` (Part 2) |
| 2F, 2G, 2H | `fig2/Fig2_all_panels.R` (Part 3) |
| 2I, 2J | `fig2/Fig2_all_panels.R` (Part 4) |
| 2K, 2L, 2M | `fig2/Fig2_all_panels.R` (Part 5) |
| 2N | `fig2/Fig2_all_panels.R` (Part 6) |
| **Figure S2** | |
| S2A, S2B, S2C | `fig2/Fig2_all_panels.R` (Part 2) |
| S2D, S2E | `fig2/Fig2_all_panels.R` (Part 7) |
| S2F | `fig2/Fig2_all_panels.R` (Part 8a CellChat objects, Part 8b plots) |
| S2G | Table, no code |
| S2H | Serum cytokine analysis performed by collaborators; data and code not included |
| S2I | `fig2/Fig2_all_panels.R` (Part 5) |
| S2J | `fig2/Fig2_all_panels.R` (Part 9) |
| **Figure 3** | |
| 3A, 3B, 3J | Xenium Explorer images, no code |
| 3C | `fig3/Fig3_all_panels.R` (Part 1) |
| 3H, 3I | `fig3/Fig3_all_panels.R` (Part 3) |
| 3K | `fig3/Fig3_all_panels.R` (Part 2) |
| 3D, 3F | `fig3/Fig3_all_panels.R` (Part 5) |
| 3E, 3G | `fig3/Fig3_all_panels.R` (Part 4) |
| **Figure S3** | |
| S3A | Niche UMAP made interactively; no script |
| S3B | `fig3/Fig3_all_panels.R` (Part 3) |
| S3C–S3F | `fig3/Fig3_all_panels.R` (Part 6) |
| **Figure 4** | |
| 4A, 4G | Xenium Explorer images, no code |
| 4B, 4C | `fig4/Fig4_all_panels.R` (Part 1) |
| 4D, 4E, 4H | `fig4/Fig4_all_panels.R` (Part 2a) |
| 4F | `fig4/Fig4_all_panels.R` (Part 3) |
| 4I, 4J | `fig4/Fig4_all_panels.R` (Part 4b) |
| 4K | `fig4/Fig4_all_panels.R` (Part 5) |
| 4L | `fig4/Fig4_all_panels.R` (Part 6) |
| **Figure S4** | |
| S4A | `fig4/Fig4_all_panels.R` (Part 1) |
| S4B, S4D | `fig4/Fig4_all_panels.R` (Part 2a) |
| S4I | `fig4/Fig4_all_panels.R` (Part 2b) |
| S4C | `fig4/Fig4_all_panels.R` (Part 7) |
| S4E, S4G | `fig4/Fig4_all_panels.R` (Part 6) |
| S4F | `fig4/Fig4_all_panels.R` (Part 8) |
| S4H | `fig4/Fig4_all_panels.R` (Part 5) |
| **Figure S5** | |
| S5A, S5B | `fig4/Fig4_all_panels.R` (Part 4a) |
| S5C–S5G | `fig4/Fig4_all_panels.R` (Part 4b) |

## Input objects

| Object | Content | Used by |
|---|---|---|
| `heart_minor.rds` (copy: `sc_heart_minor.rds`) | Combined heart scRNA-seq (both cohorts) with minor cell-type labels | Fig 1, Fig 2 |
| `skeletal_322026.rds` (copy: `sc_skeletal_322026.rds`) | Annotated skeletal muscle scRNA-seq | Fig 1, Fig 2 |
| `tnk_clean.rds` | Heart T/NK sub-clusters | Fig 2 |
| `pbmc_2192026.rds`, `pbmc_642025.rds` | Annotated PBMC scRNA-seq | Fig 1, Fig 2, Fig 4 |
| `xo_heart_processed.rds` | Xenium heart with cell types and niches | Fig 3, Fig 4 |
| `skeletal_combined_352026.rds` | Xenium skeletal muscle with cell types | Fig 3 |

Objects made by the figure scripts and used later:

| Object | Created by | Content |
|---|---|---|
| `heart_3232026.rds` | `fig1/Fig1_all_panels.R` (Part 2) | Final heart object used for Fig 2 and Fig 4 |
| `tnk_myeloid_heart_st.rds` | `fig3/Fig3_all_panels.R` (Part 6a) | Xenium T/NK and myeloid cells with SingleR subtypes and the IL-17 program score |
| `xo_heart_processed_withspatialil17labels.rds` | `fig4/Fig4_all_panels.R` (Part 2a) | Xenium heart with IL-17 hotspot and proximity labels |

## Data availability

| Dataset | Source |
|---|---|
| Cohort 1 scRNA-seq (MGH; Blum *et al.*, *Nature* 2024) | GEO [GSE228597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228597) |
| Cohort 2 scRNA-seq (MDACC; Siddiqui *et al.*) | EGA EGAS00001007511 (controlled access) |
| Xenium 5K spatial transcriptomics (VUMC) | *Accession to be added* |
| Serum cytokines (Luminex 71-plex) | Generated and analysed by collaborators; available from the corresponding authors on reasonable request |

Processed objects are not included in this repository. Patient-level clinical data are not shared.

## Contact

Hassan Abushukair, University of Oklahoma Health Sciences Center.
