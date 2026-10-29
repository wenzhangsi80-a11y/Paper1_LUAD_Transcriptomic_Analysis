## 00_setup.R

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

work_dir <- input_dir

setwd(work_dir)

getwd()


load("input_dir/Paper1_PRE_STEP19_WORKSPACE_20260914.RData")

ls()

# ============================================================
# Paper 1
# Supplementary Figure 4
# Tumor-only WGCNA sensitivity/refinement analysis
#
# Plotting only
# NO WGCNA recalculation
# NO GS recalculation
# NO MM recalculation
# NO correlation recalculation
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE 4 — PLOTTING WORKFLOW\n")
cat("============================================================\n")

##Supplementary Figure 4A
#50 modules × 3 TME traits heatmap

##SF4A-1：Prepare data
# ============================================================
# SF4A
# Tumor-only WGCNA module–TME association heatmap
# ============================================================

cat("\n============================================================\n")
cat("SF4A: Module–TME association heatmap\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check locked object
# ------------------------------------------------------------

stopifnot(
  exists("module_tme_results")
)

cat(
  "module_tme_results dimensions:",
  paste(dim(module_tme_results), collapse = " × "),
  "\n"
)

print(
  names(module_tme_results)
)


# ------------------------------------------------------------
# 2. Preserve locked data
# ------------------------------------------------------------

sf4a_df <- module_tme_results


# ------------------------------------------------------------
# 3. Check expected traits
# ------------------------------------------------------------

expected_traits <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)

missing_traits <- setdiff(
  expected_traits,
  unique(sf4a_df$Trait)
)

if (length(missing_traits) > 0) {
  
  stop(
    paste(
      "Missing expected traits:",
      paste(missing_traits, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Module ordering
# ------------------------------------------------------------

module_order <- unique(
  sf4a_df$Module
)

cat(
  "Number of modules:",
  length(module_order),
  "\n"
)


# ------------------------------------------------------------
# 5. Trait ordering
# ------------------------------------------------------------

sf4a_df$Trait <- factor(
  sf4a_df$Trait,
  levels = expected_traits
)

sf4a_df$Module <- factor(
  sf4a_df$Module,
  levels = rev(module_order)
)


# ------------------------------------------------------------
# 6. Significance label
# ------------------------------------------------------------

sf4a_df <- sf4a_df %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )


##SF4A-2：Mark ME2 / ME9
#ME2 = blue
#ME9 = magenta
# ------------------------------------------------------------
# 7. Create display labels
# ------------------------------------------------------------

sf4a_df <- sf4a_df %>%
  mutate(
    Module_label = case_when(
      as.character(Module) == "ME2" ~ "ME2 (blue)",
      as.character(Module) == "ME9" ~ "ME9 (magenta)",
      TRUE ~ as.character(Module)
    )
  )


# Preserve ordering after relabeling
module_display_order <- rev(
  c(
    "ME2 (blue)",
    "ME9 (magenta)"
  )
)


# ------------------------------------------------------------
# 8. Build heatmap
# ------------------------------------------------------------

p_sf4a <- ggplot(
  sf4a_df,
  aes(
    x = Trait,
    y = Module,
    fill = Spearman_rho
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = Sig),
    size = 3
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    name = expression(rho)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 7
    ),
    axis.title = element_blank(),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.margin = margin(
      5, 5, 5, 5
    )
  )

p_sf4a


# ============================================================
# Paper 1
# Supplementary Figure 4A
#
# FINAL VERSION
# Tumor-only WGCNA module–TME association heatmap
#
# Design:
#   1. Heatmap fill = locked Spearman's rho
#   2. Significance symbols = locked BH-adjusted FDR
#   3. Left-side strip = highlight ME2 and ME9
#   4. Y-axis = WGCNA module names
#   5. No additional WGCNA color squares beside Y-axis labels
#   6. Right-side legend:
#        FDR explanation
#        rho colorbar
#
# IMPORTANT:
#   Plotting only.
#   No WGCNA statistics are recalculated.
#   No module–trait correlation is recalculated.
#   No P value is recalculated.
#   No FDR is recalculated.
# ============================================================


cat("\n============================================================\n")
cat("SF4A FINAL — WGCNA–TME HEATMAP\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(showtext)
library(sysfonts)
library(ragg)


# ------------------------------------------------------------
# 1. Register Arial font
# ------------------------------------------------------------

font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/arial.ttf",
  bold = "C:/Windows/Fonts/arialbd.ttf"
)

showtext_auto()


# ------------------------------------------------------------
# 2. Validate sf4a_df
# ------------------------------------------------------------

required_cols <- c(
  "Module",
  "Trait",
  "Spearman_rho",
  "FDR"
)

missing_cols <- setdiff(
  required_cols,
  names(sf4a_df)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns in sf4a_df:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 3. Create plotting copy
#
# IMPORTANT:
#   No statistical values are recalculated.
# ------------------------------------------------------------

sf4a_plot_df <- sf4a_df %>%
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  )


# ------------------------------------------------------------
# 4. Significance annotation
#
# Uses the already locked FDR values.
#
# IMPORTANT:
#   This is only a visual categorization of FDR.
#   No FDR is recalculated.
# ------------------------------------------------------------

sf4a_plot_df <- sf4a_plot_df %>%
  dplyr::mutate(
    Sig = dplyr::case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )


# ------------------------------------------------------------
# 5. Validate module / trait structure
# ------------------------------------------------------------

cat("Number of modules: ",
    length(unique(sf4a_plot_df$Module)),
    "\n",
    sep = "")

cat("Number of traits: ",
    length(unique(sf4a_plot_df$Trait)),
    "\n",
    sep = "")

cat("\nModules:\n")
print(unique(sf4a_plot_df$Module))

cat("\nTraits:\n")
print(unique(sf4a_plot_df$Trait))


# ------------------------------------------------------------
# 6. Lock trait order
# ------------------------------------------------------------

sf4a_plot_df <- sf4a_plot_df %>%
  dplyr::mutate(
    Trait = factor(
      Trait,
      levels = c(
        "ImmuneScore",
        "StromalScore",
        "ESTIMATEScore"
      )
    )
  )


# ------------------------------------------------------------
# 7. Lock module order
#
# IMPORTANT:
#   This does NOT reorder according to correlation.
#   It preserves the existing module order in sf4a_df.
# ------------------------------------------------------------

module_order <- unique(
  sf4a_plot_df$Module
)

sf4a_plot_df <- sf4a_plot_df %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(module_order)
    )
  )


# ------------------------------------------------------------
# 8. Create highlight annotation
#
# ME2  -> blue
# ME9  -> magenta
#
# Other modules -> blank
#
# This is a visual annotation only.
# ------------------------------------------------------------

sf4a_plot_df <- sf4a_plot_df %>%
  dplyr::mutate(
    Highlight = dplyr::case_when(
      as.character(Module) == "ME2" ~ "ME2",
      as.character(Module) == "ME9" ~ "ME9",
      TRUE ~ "Other"
    )
  )


# ------------------------------------------------------------
# 9. Create separate row annotation dataframe
#
# One row per WGCNA module.
# ------------------------------------------------------------

sf4a_annotation_df <- data.frame(
  Module = module_order,
  Highlight = dplyr::case_when(
    module_order == "ME2" ~ "ME2",
    module_order == "ME9" ~ "ME9",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

sf4a_annotation_df$Module <- factor(
  sf4a_annotation_df$Module,
  levels = rev(module_order)
)


# ------------------------------------------------------------
# 10. Main heatmap
# ------------------------------------------------------------

p_sf4a_heatmap <- ggplot(
  sf4a_plot_df,
  aes(
    x = Trait,
    y = Module,
    fill = Spearman_rho
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  
  geom_text(
    aes(
      label = Sig
    ),
    family = "Arial",
    size = 4.0,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_fill_gradient2(
    name = "ρ",
    low = "#82B5EC",
    mid = "white",
    high = "#F57C6E",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  scale_x_discrete(
    labels = c(
      "ImmuneScore" = "ImmuneScore",
      "StromalScore" = "StromalScore",
      "ESTIMATEScore" = "ESTIMATEScore"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +
  
  theme(
    
    axis.text.x = element_text(
      family = "Arial",
      size = 10,
      face = "bold",
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.text.y = element_text(
      family = "Arial",
      size = 8.5,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    
    axis.line = element_blank(),
    
    legend.position = "none",
    
    panel.background = element_blank(),
    
    plot.background = element_blank(),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )


# ------------------------------------------------------------
# 11. Left-side ME2 / ME9 highlight strip
#
# This is deliberately separate from the heatmap.
# ------------------------------------------------------------

p_sf4a_strip <- ggplot(
  sf4a_annotation_df,
  aes(
    x = 1,
    y = Module,
    fill = Highlight
  )
) +
  
  geom_tile(
    width = 1,
    height = 0.85
  ) +
  
  scale_fill_manual(
    values = c(
      "ME2" = "#4C78A8",
      "ME9" = "#B279A2",
      "Other" = "white"
    )
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  scale_x_continuous(
    limits = c(0.5, 1.5),
    expand = c(0, 0)
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "none",
    plot.margin = margin(
      t = 5,
      r = 2,
      b = 5,
      l = 0
    )
  )


# ------------------------------------------------------------
# 12. Combine strip + heatmap
# ------------------------------------------------------------

sf4a_main <- cowplot::plot_grid(
  p_sf4a_strip,
  p_sf4a_heatmap,
  nrow = 1,
  rel_widths = c(
    0.035,
    1
  ),
  align = "h",
  axis = "tb"
)


# ------------------------------------------------------------
# 13. Create FDR explanation legend
#
# IMPORTANT:
#   This explicitly refers to BH-adjusted FDR.
# ------------------------------------------------------------

p_sf4a_fdr <- ggplot() +
  
  annotate(
    "text",
    x = 0,
    y = 4,
    label = "FDR",
    hjust = 0,
    family = "Arial",
    size = 8 / .pt,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 3,
    label = "***  FDR < 0.001",
    hjust = 0,
    family = "Arial",
    size = 7 / .pt
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 2,
    label = "**   FDR < 0.01",
    hjust = 0,
    family = "Arial",
    size = 7 / .pt
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 1,
    label = "*    FDR < 0.05",
    hjust = 0,
    family = "Arial",
    size = 7 / .pt
  ) +
  
  coord_cartesian(
    xlim = c(
      -0.05,
      1.8
    ),
    ylim = c(
      0.7,
      4.3
    ),
    expand = FALSE,
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )


# ------------------------------------------------------------
# 14. Extract rho colorbar
# ------------------------------------------------------------

rho_legend_sf4a <- cowplot::get_legend(
  
  p_sf4a_heatmap +
    
    theme(
      legend.position = "right",
      
      legend.title = element_text(
        family = "Arial",
        size = 8,
        face = "bold",
        hjust = 0
      ),
      
      legend.text = element_text(
        family = "Arial",
        size = 7
      ),
      
      legend.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    ) +
    
    guides(
      fill = guide_colorbar(
        title = "ρ",
        title.position = "top",
        title.hjust = 0,
        
        barwidth = unit(
          0.35,
          "cm"
        ),
        
        barheight = unit(
          2.2,
          "cm"
        )
      )
    )
)


# ------------------------------------------------------------
# 15. Stack FDR explanation above rho colorbar
# ------------------------------------------------------------

right_legend_sf4a <- cowplot::ggdraw() +
  
  cowplot::draw_plot(
    p_sf4a_fdr,
    x = 0,
    y = 0.48,
    width = 1,
    height = 0.24
  ) +
  
  cowplot::draw_grob(
    rho_legend_sf4a,
    x = 0,
    y = 0.03,
    width = 1,
    height = 0.42
  )


# ------------------------------------------------------------
# 16. Final SF4A
#
# Main heatmap remains separate from right-side legend.
# ------------------------------------------------------------

SF4A_final <- cowplot::plot_grid(
  sf4a_main,
  right_legend_sf4a,
  nrow = 1,
  rel_widths = c(
    1,
    0.17
  ),
  align = "h",
  axis = "tb"
)


# ------------------------------------------------------------
# 17. Display final SF4A
# ------------------------------------------------------------

SF4A_final


# ------------------------------------------------------------
# 18. Create output directory
# ------------------------------------------------------------

sf4a_output_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4A"
)

if (
  !dir.exists(sf4a_output_dir)
) {
  dir.create(
    sf4a_output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 19. Export 600 dpi PNG
# ------------------------------------------------------------

sf4a_png <- file.path(
  sf4a_output_dir,
  "Supplementary_Figure_4A_WGCNA_TME_heatmap_600dpi.png"
)

ggsave(
  filename = sf4a_png,
  plot = SF4A_final,
  width = 6.8,
  height = 7.2,
  units = "in",
  dpi = 600,
  device = ragg::agg_png,
  bg = "white"
)


# ------------------------------------------------------------
# 20. Export 600 dpi TIFF
# ------------------------------------------------------------

sf4a_tiff <- file.path(
  sf4a_output_dir,
  "Supplementary_Figure_4A_WGCNA_TME_heatmap_600dpi.tiff"
)

ggsave(
  filename = sf4a_tiff,
  plot = SF4A_final,
  width = 6.8,
  height = 7.2,
  units = "in",
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 21. Export vector PDF
# ------------------------------------------------------------

sf4a_pdf <- file.path(
  sf4a_output_dir,
  "Supplementary_Figure_4A_WGCNA_TME_heatmap.pdf"
)

ggsave(
  filename = sf4a_pdf,
  plot = SF4A_final,
  width = 6.8,
  height = 7.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ------------------------------------------------------------
# 22. Export plotting dataframe for audit
# ------------------------------------------------------------

sf4a_table_file <- file.path(
  sf4a_output_dir,
  "Supplementary_Figure_4A_plotting_data_audit.csv"
)

write.csv(
  sf4a_plot_df,
  sf4a_table_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Output QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("SF4A OUTPUT QC\n")
cat("------------------------------------------------------------\n")

cat(
  "PNG exists: ",
  file.exists(sf4a_png),
  "\n",
  sep = ""
)

cat(
  "TIFF exists: ",
  file.exists(sf4a_tiff),
  "\n",
  sep = ""
)

cat(
  "PDF exists: ",
  file.exists(sf4a_pdf),
  "\n",
  sep = ""
)

cat(
  "Audit table exists: ",
  file.exists(sf4a_table_file),
  "\n",
  sep = ""
)


if (
  !all(
    file.exists(
      c(
        sf4a_png,
        sf4a_tiff,
        sf4a_pdf,
        sf4a_table_file
      )
    )
  )
) {
  stop(
    "SF4A output QC FAILED."
  )
}


# ------------------------------------------------------------
# 24. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SF4A FINAL STATUS: PASS\n")
cat("============================================================\n")

cat(
  "Heatmap uses locked Spearman's rho values.\n"
)

cat(
  "Significance symbols use locked BH-adjusted FDR values.\n"
)

cat(
  "ME2 and ME9 are highlighted by a left-side annotation strip.\n"
)

cat(
  "No additional WGCNA module-color square was added beside Y-axis labels.\n"
)

cat(
  "FDR explanation is positioned above the rho colorbar.\n"
)

cat(
  "SF4A exported separately for external figure assembly.\n"
)

cat("============================================================\n")


# ============================================================
# Arial font diagnostic
# Plotting-only diagnostic; does not modify SF4A results
# ============================================================

cat("\n============================================================\n")
cat("Arial font diagnostic\n")
cat("============================================================\n")

cat("\nR version:\n")
print(R.version.string)

cat("\nWindows Arial font files:\n")

arial_candidates <- c(
  "C:/Windows/Fonts/arial.ttf",
  "C:/Windows/Fonts/arialbd.ttf",
  "C:/Windows/Fonts/ariali.ttf",
  "C:/Windows/Fonts/arialbi.ttf"
)

print(
  data.frame(
    File = arial_candidates,
    Exists = file.exists(arial_candidates)
  )
)

cat("\nRegistered font families in sysfonts:\n")

if (requireNamespace("sysfonts", quietly = TRUE)) {
  print(sysfonts::font_families())
} else {
  cat("sysfonts is not installed.\n")
}

cat("\nshowtext status:\n")

if (requireNamespace("showtext", quietly = TRUE)) {
  cat("showtext installed: YES\n")
} else {
  cat("showtext installed: NO\n")
}

cat("\n============================================================\n")

# ============================================================
# Paper 1
# Supplementary Figure 4A
#SF4A FINAL v4 — CAIRO PDF MASTER ONLY

# ============================================================
#
# Paper 1
#
# Supplementary Figure 4A
#
# SF4A FINAL v4 — CAIRO PDF MASTER ONLY
#
# WGCNA–TME HEATMAP
#
# ============================================================
#
# PURPOSE
#
#   Final vector-PDF master for Supplementary Figure 4A.
#
#   This version is based on the validated SF4A FINAL v3
#   plotting parameters, with the following final changes:
#
#   1. Cairo PDF only
#   2. Arial font for final vector PDF
#   3. Main heatmap explicitly contains NO legend
#   4. A separate legend-only plot is used to extract
#      exactly ONE Spearman's rho colorbar
#   5. FDR explanation is placed above the single rho colorbar
#   6. ME2 / ME9 highlight colors use the locked visual colors:
#
#        ME2 = #4472C4
#        ME9 = #C000C0
#        Other = white
#
# IMPORTANT
#
#   - No statistical values are recalculated.
#   - No locked WGCNA results are modified.
#   - sf4a_df remains unchanged.
#   - Only plotting / rendering / export settings are used.
#   - No ragg.
#   - No showtext.
#   - No sysfonts.
#   - No Ghostscript.
#   - No PNG/TIFF export from R.
#   - PDF only.
#
# LOCKED CONTENT PRESERVED
#
#   - Spearman's rho
#   - BH-adjusted FDR
#   - ME2 / ME9 highlighting
#   - module order
#   - trait order
#   - heatmap layout
#
# EXPORT
#
#   PDF : 6.8 × 7.2 in, vector PDF
#
# OUTPUT
#
#   Paper1_STEP19_PREP_20260914/
#       audit/
#           SF4A/
#               FINAL_v4_CAIRO_PDF_MASTER/
#
# ============================================================


cat("\n============================================================\n")
cat("SF4A FINAL v4 — CAIRO PDF MASTER ONLY\n")
cat("WGCNA–TME HEATMAP\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "cowplot",
  "grid"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste(
      "Missing required packages:",
      paste(missing_packages, collapse = ", ")
    )
  )
  
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)


# ------------------------------------------------------------
# 1. Final PDF rendering strategy
# ------------------------------------------------------------

#
# Cairo PDF has already been verified to render Arial correctly
# in the current Windows/R environment.
#
# Therefore:
#
#   PDF master = Cairo PDF + Arial
#
# No showtext/sysfonts/ragg are used.
#

sf4a_font_family <- "Arial"

cat(
  "SF4A rendering font family: ",
  sf4a_font_family,
  "\n",
  sep = ""
)

cat(
  "Rendering device: Cairo PDF\n"
)

cat(
  "ragg/showtext/sysfonts intentionally disabled.\n\n"
)


# ------------------------------------------------------------
# 2. Validate required plotting object
# ------------------------------------------------------------

if (!exists("sf4a_df", inherits = TRUE)) {
  
  stop(
    "Object 'sf4a_df' was not found in the current R session."
  )
  
}


# ------------------------------------------------------------
# 3. Validate sf4a_df required columns
# ------------------------------------------------------------

required_cols <- c(
  "Module",
  "Trait",
  "Spearman_rho",
  "FDR"
)

missing_cols <- setdiff(
  required_cols,
  names(sf4a_df)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      "Missing required columns in sf4a_df:",
      paste(missing_cols, collapse = ", ")
    )
  )
  
}

cat(
  "Required sf4a_df columns: PASS\n"
)


# ------------------------------------------------------------
# 4. Preserve locked reference data
# ------------------------------------------------------------

sf4a_locked_reference <- sf4a_df %>%
  
  dplyr::select(
    all_of(required_cols)
  ) %>%
  
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  )


# ------------------------------------------------------------
# 5. Basic locked-data validation
# ------------------------------------------------------------

if (nrow(sf4a_locked_reference) != 150) {
  
  stop(
    paste(
      "Unexpected number of SF4A rows:",
      nrow(sf4a_locked_reference),
      "(expected 150)"
    )
  )
  
}


if (
  length(unique(sf4a_locked_reference$Module)) != 50
) {
  
  stop(
    paste(
      "Unexpected number of modules:",
      length(unique(sf4a_locked_reference$Module)),
      "(expected 50)"
    )
  )
  
}


expected_traits <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)


if (
  !all(
    expected_traits %in%
    unique(sf4a_locked_reference$Trait)
  )
) {
  
  stop(
    "SF4A trait structure does not contain all three expected TME traits."
  )
  
}


cat(
  "Locked SF4A structure: PASS\n"
)

cat(
  "  Rows    : ",
  nrow(sf4a_locked_reference),
  "\n",
  sep = ""
)

cat(
  "  Modules : ",
  length(unique(sf4a_locked_reference$Module)),
  "\n",
  sep = ""
)

cat(
  "  Traits  : ",
  length(unique(sf4a_locked_reference$Trait)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Create plotting copy
# ------------------------------------------------------------

#
# sf4a_df itself is NOT modified.
#

sf4a_plot_df <- sf4a_df %>%
  
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  )


# ------------------------------------------------------------
# 7. Verify plotting copy preserves locked values
# ------------------------------------------------------------

sf4a_plot_reference <- sf4a_plot_df %>%
  
  dplyr::select(
    all_of(required_cols)
  ) %>%
  
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  )


if (
  !isTRUE(
    all.equal(
      sf4a_plot_reference,
      sf4a_locked_reference,
      check.attributes = FALSE,
      tolerance = 1e-12
    )
  )
) {
  
  stop(
    "Plotting copy does not exactly preserve locked SF4A data."
  )
  
}


cat(
  "Locked plotting values preserved: PASS\n"
)


# ------------------------------------------------------------
# 8. Significance annotation
# ------------------------------------------------------------

#
# Uses locked BH-adjusted FDR values.
#
# No statistical test is performed here.
#

sf4a_plot_df <- sf4a_plot_df %>%
  
  dplyr::mutate(
    
    Sig = dplyr::case_when(
      
      FDR < 0.001 ~ "***",
      
      FDR < 0.01 ~ "**",
      
      FDR < 0.05 ~ "*",
      
      TRUE ~ ""
      
    )
    
  )


# ------------------------------------------------------------
# 9. Lock trait order
# ------------------------------------------------------------

sf4a_plot_df <- sf4a_plot_df %>%
  
  dplyr::mutate(
    
    Trait = factor(
      Trait,
      levels = expected_traits
    )
    
  )


if (
  !identical(
    levels(sf4a_plot_df$Trait),
    expected_traits
  )
) {
  
  stop(
    "Trait order lock failed."
  )
  
}


# ------------------------------------------------------------
# 10. Lock module order
# ------------------------------------------------------------

module_order <- unique(
  sf4a_plot_df$Module
)


if (length(module_order) != 50) {
  
  stop(
    paste(
      "Module order validation failed:",
      length(module_order),
      "modules detected."
    )
  )
  
}


sf4a_plot_df <- sf4a_plot_df %>%
  
  dplyr::mutate(
    
    Module = factor(
      Module,
      levels = rev(module_order)
    )
    
  )


# ------------------------------------------------------------
# 11. Create ME2 / ME9 highlight annotation
# ------------------------------------------------------------

#
# LOCKED VISUAL COLORS
#
# ME2 = blue
# ME9 = magenta
# Other = white
#

sf4a_plot_df <- sf4a_plot_df %>%
  
  dplyr::mutate(
    
    Highlight = dplyr::case_when(
      
      as.character(Module) == "ME2" ~ "ME2",
      
      as.character(Module) == "ME9" ~ "ME9",
      
      TRUE ~ "Other"
      
    )
    
  )


# ------------------------------------------------------------
# 12. Create row annotation dataframe
# ------------------------------------------------------------

sf4a_annotation_df <- data.frame(
  
  Module = module_order,
  
  Highlight = dplyr::case_when(
    
    module_order == "ME2" ~ "ME2",
    
    module_order == "ME9" ~ "ME9",
    
    TRUE ~ "Other"
    
  ),
  
  stringsAsFactors = FALSE
  
)


sf4a_annotation_df$Module <- factor(
  
  sf4a_annotation_df$Module,
  
  levels = rev(module_order)
  
)


# ------------------------------------------------------------
# 13. Rendering warning capture
# ------------------------------------------------------------

sf4a_render_warnings <- character(0)


# ============================================================
#
# 14. MAIN HEATMAP
#
# IMPORTANT:
#
#   THE MAIN HEATMAP HAS NO LEGEND.
#
#   Both:
#
#       legend.position = "none"
#
#   and:
#
#       guides(fill = "none")
#
#   are explicitly specified.
#
#   Therefore the main heatmap can NEVER contribute a rho
#   colorbar to the final cowplot.
#
# ============================================================

withCallingHandlers({
  
  p_sf4a_heatmap <- ggplot(
    
    sf4a_plot_df,
    
    aes(
      x = Trait,
      y = Module,
      fill = Spearman_rho
    )
    
  ) +
    
    geom_tile(
      
      color = "white",
      linewidth = 0.35
      
    ) +
    
    geom_text(
      
      aes(
        label = Sig
      ),
      
      family = sf4a_font_family,
      
      size = 3.2,
      
      fontface = "bold",
      
      color = "black"
      
    ) +
    
    scale_fill_gradient2(
      
      name = "ρ",
      
      low = "#82B5EC",
      
      mid = "white",
      
      high = "#F57C6E",
      
      midpoint = 0,
      
      limits = c(-1, 1)
      
    ) +
    
    scale_y_discrete(
      drop = FALSE
    ) +
    
    scale_x_discrete(
      
      labels = c(
        
        "ImmuneScore" = "ImmuneScore",
        
        "StromalScore" = "StromalScore",
        
        "ESTIMATEScore" = "ESTIMATEScore"
        
      )
      
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_classic(
      
      base_size = 11,
      
      base_family = sf4a_font_family
      
    ) +
    
    theme(
      
      axis.text.x = element_text(
        
        family = sf4a_font_family,
        
        size = 10,
        
        face = "bold",
        
        color = "black",
        
        angle = 0,
        
        hjust = 0.5
        
      ),
      
      axis.text.y = element_text(
        
        family = sf4a_font_family,
        
        size = 9,
        
        color = "black"
        
      ),
      
      axis.ticks = element_blank(),
      
      axis.line = element_blank(),
      
      #
      # CRITICAL:
      #
      legend.position = "none",
      
      panel.background = element_blank(),
      
      plot.background = element_blank(),
      
      plot.margin = margin(
        
        t = 5,
        
        r = 5,
        
        b = 5,
        
        l = 5
        
      )
      
    ) +
    
    #
    # CRITICAL:
    #
    # Explicitly remove ALL fill guides from the main plot.
    #
    guides(
      fill = "none"
    )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ============================================================
#
# 15. LEFT-SIDE ME2 / ME9 HIGHLIGHT STRIP
#
# ============================================================

withCallingHandlers({
  
  p_sf4a_strip <- ggplot(
    
    sf4a_annotation_df,
    
    aes(
      x = 1,
      y = Module,
      fill = Highlight
    )
    
  ) +
    
    geom_tile(
      
      width = 1,
      
      height = 0.85
      
    ) +
    
    scale_fill_manual(
      
      values = c(
        
        "ME2" = "#4472C4",
        
        "ME9" = "#C000C0",
        
        "Other" = "white"
        
      ),
      
      guide = "none"
      
    ) +
    
    scale_y_discrete(
      drop = FALSE
    ) +
    
    scale_x_continuous(
      
      limits = c(
        0.5,
        1.5
      ),
      
      expand = c(0, 0)
      
    ) +
    
    theme_void() +
    
    theme(
      
      legend.position = "none",
      
      plot.margin = margin(
        
        t = 5,
        
        r = 2,
        
        b = 5,
        
        l = 0
        
      )
      
    )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ------------------------------------------------------------
# 16. Combine strip + main heatmap
# ------------------------------------------------------------

withCallingHandlers({
  
  sf4a_main <- cowplot::plot_grid(
    
    p_sf4a_strip,
    
    p_sf4a_heatmap,
    
    nrow = 1,
    
    rel_widths = c(
      0.035,
      1
    ),
    
    align = "h",
    
    axis = "tb"
    
  )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ============================================================
#
# 17. FDR EXPLANATION
#
# ============================================================

withCallingHandlers({
  
  p_sf4a_fdr <- ggplot() +
    
    annotate(
      
      "text",
      
      x = 0,
      
      y = 4,
      
      label = "FDR",
      
      hjust = 0,
      
      family = sf4a_font_family,
      
      size = 9 / .pt,
      
      fontface = "bold"
      
    ) +
    
    annotate(
      
      "text",
      
      x = 0,
      
      y = 3,
      
      label = "***  FDR < 0.001",
      
      hjust = 0,
      
      family = sf4a_font_family,
      
      size = 8 / .pt
      
    ) +
    
    annotate(
      
      "text",
      
      x = 0,
      
      y = 2,
      
      label = "**   FDR < 0.01",
      
      hjust = 0,
      
      family = sf4a_font_family,
      
      size = 8 / .pt
      
    ) +
    
    annotate(
      
      "text",
      
      x = 0,
      
      y = 1,
      
      label = "*    FDR < 0.05",
      
      hjust = 0,
      
      family = sf4a_font_family,
      
      size = 8 / .pt
      
    ) +
    
    coord_cartesian(
      
      xlim = c(
        -0.05,
        1.8
      ),
      
      ylim = c(
        0.7,
        4.3
      ),
      
      expand = FALSE,
      
      clip = "off"
      
    ) +
    
    theme_void() +
    
    theme(
      
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
      
    )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ============================================================
#
# 18. CREATE A DEDICATED LEGEND-ONLY PLOT
#
# IMPORTANT:
#
#   This plot exists ONLY to extract the rho colorbar.
#
#   It is NEVER placed directly into the final figure.
#
#   This avoids accidental duplication of the colorbar.
#
# ============================================================

withCallingHandlers({
  
  p_sf4a_rho_legend_source <- ggplot(
    
    sf4a_plot_df,
    
    aes(
      x = Trait,
      y = Module,
      fill = Spearman_rho
    )
    
  ) +
    
    geom_tile() +
    
    scale_fill_gradient2(
      
      name = "ρ",
      
      low = "#82B5EC",
      
      mid = "white",
      
      high = "#F57C6E",
      
      midpoint = 0,
      
      limits = c(-1, 1)
      
    ) +
    
    guides(
      
      fill = guide_colorbar(
        
        title = "ρ",
        
        title.position = "top",
        
        title.hjust = 0,
        
        barwidth = unit(
          0.38,
          "cm"
        ),
        
        barheight = unit(
          2.4,
          "cm"
        )
        
      )
      
    ) +
    
    theme_void() +
    
    theme(
      
      legend.position = "right",
      
      legend.title = element_text(
        
        family = sf4a_font_family,
        
        size = 9,
        
        face = "bold",
        
        hjust = 0
        
      ),
      
      legend.text = element_text(
        
        family = sf4a_font_family,
        
        size = 8
        
      ),
      
      legend.margin = margin(
        0,
        0,
        0,
        0
      )
      
    )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ------------------------------------------------------------
# 19. Extract exactly ONE rho colorbar
# ------------------------------------------------------------

withCallingHandlers({
  
  rho_legend_sf4a <- cowplot::get_legend(
    
    p_sf4a_rho_legend_source
    
  )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ============================================================
#
# 20. STACK FDR ABOVE SINGLE RHO COLORBAR
#
# ============================================================

withCallingHandlers({
  
  right_legend_sf4a <- cowplot::ggdraw() +
    
    #
    # FDR explanation
    #
    cowplot::draw_plot(
      
      p_sf4a_fdr,
      
      x = 0,
      
      y = 0.48,
      
      width = 1,
      
      height = 0.24
      
    ) +
    
    #
    # SINGLE rho colorbar
    #
    cowplot::draw_grob(
      
      rho_legend_sf4a,
      
      x = 0,
      
      y = 0.10,
      
      width = 1,
      
      height = 0.42
      
    )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ============================================================
#
# 21. FINAL SF4A
#
# ============================================================

withCallingHandlers({
  
  SF4A_final <- cowplot::plot_grid(
    
    sf4a_main,
    
    right_legend_sf4a,
    
    nrow = 1,
    
    rel_widths = c(
      1,
      0.17
    ),
    
    align = "h",
    
    axis = "tb"
    
  )
  
},
warning = function(w) {
  
  sf4a_render_warnings <<-
    
    c(
      sf4a_render_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ------------------------------------------------------------
# 22. Display final SF4A
# ------------------------------------------------------------

SF4A_final


# ------------------------------------------------------------
# 23. Rendering warning checkpoint
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("SF4A RENDERING WARNING CHECKPOINT\n")
cat("------------------------------------------------------------\n")

if (length(sf4a_render_warnings) == 0) {
  
  cat(
    "Plot construction warnings: NONE\n"
  )
  
} else {
  
  cat(
    "Plot construction warnings detected: ",
    length(sf4a_render_warnings),
    "\n",
    sep = ""
  )
  
  cat("\nWarnings:\n")
  
  print(
    unique(sf4a_render_warnings)
  )
  
}


# ============================================================
#
# 24. Create v4 output directory
#
# ============================================================

sf4a_output_dir_v4 <- file.path(
  
  "Paper1_STEP19_PREP_20260914",
  
  "audit",
  
  "SF4A",
  
  "FINAL_v4_CAIRO_PDF_MASTER"
  
)


if (!dir.exists(sf4a_output_dir_v4)) {
  
  dir.create(
    
    sf4a_output_dir_v4,
    
    recursive = TRUE
    
  )
  
}


# ------------------------------------------------------------
# 25. Define output files
# ------------------------------------------------------------

sf4a_pdf_v4 <- file.path(
  
  sf4a_output_dir_v4,
  
  "Paper1_Supplementary_Figure4A_FINAL_Cairo.pdf"
  
)


sf4a_table_file_v4 <- file.path(
  
  sf4a_output_dir_v4,
  
  "Paper1_Supplementary_Figure4A_FINAL_Cairo_audit.csv"
  
)


sf4a_warning_file_v4 <- file.path(
  
  sf4a_output_dir_v4,
  
  "SF4A_FINAL_v4_Cairo_render_warning_audit.txt"
  
)


# ============================================================
#
# 26. Export vector PDF using Cairo
#
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Exporting FINAL Cairo PDF...\n")
cat("------------------------------------------------------------\n")


pdf_export_warnings <- character(0)


withCallingHandlers({
  
  grDevices::cairo_pdf(
    
    filename = sf4a_pdf_v4,
    
    width = 6.8,
    
    height = 7.2,
    
    family = sf4a_font_family,
    
    onefile = TRUE
    
  )
  
  print(SF4A_final)
  
  grDevices::dev.off()
  
},
warning = function(w) {
  
  pdf_export_warnings <<-
    
    c(
      pdf_export_warnings,
      conditionMessage(w)
    )
  
  invokeRestart("muffleWarning")
  
})


# ------------------------------------------------------------
# 27. Export plotting dataframe for audit
# ------------------------------------------------------------

write.csv(
  
  sf4a_plot_df,
  
  sf4a_table_file_v4,
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 28. Export warning audit
# ------------------------------------------------------------

all_export_warnings <- c(
  
  sf4a_render_warnings,
  
  pdf_export_warnings
  
)


writeLines(
  
  c(
    
    "SF4A FINAL v4 — CAIRO PDF RENDER WARNING AUDIT",
    
    "=================================================",
    
    "",
    
    paste(
      "Total warnings:",
      length(all_export_warnings)
    ),
    
    "",
    
    if (
      length(all_export_warnings) == 0
    ) {
      
      "NO WARNINGS DETECTED."
      
    } else {
      
      unique(all_export_warnings)
      
    }
    
  ),
  
  con = sf4a_warning_file_v4
  
)


# ============================================================
#
# 29. File existence QC
#
# ============================================================

cat("\n============================================================\n")
cat("SF4A FINAL v4 OUTPUT QC\n")
cat("============================================================\n")


cat(
  "PDF exists: ",
  file.exists(sf4a_pdf_v4),
  "\n",
  sep = ""
)


cat(
  "Audit table exists: ",
  file.exists(sf4a_table_file_v4),
  "\n",
  sep = ""
)


cat(
  "Warning audit exists: ",
  file.exists(sf4a_warning_file_v4),
  "\n",
  sep = ""
)


required_output_files <- c(
  
  sf4a_pdf_v4,
  
  sf4a_table_file_v4,
  
  sf4a_warning_file_v4
  
)


if (
  !all(
    file.exists(
      required_output_files
    )
  )
) {
  
  stop(
    "SF4A FINAL v4 output file QC FAILED."
  )
  
}


# ============================================================
#
# 30. File-size QC
#
# ============================================================

output_sizes <- file.info(
  
  required_output_files
  
)$size


cat("\nOutput file sizes:\n")

print(
  
  data.frame(
    
    File = basename(required_output_files),
    
    Size_bytes = output_sizes,
    
    stringsAsFactors = FALSE
    
  )
  
)


if (
  any(
    is.na(output_sizes) |
    output_sizes <= 0
  )
) {
  
  stop(
    "SF4A FINAL v4 file-size QC FAILED."
  )
  
}


# ============================================================
#
# 31. Locked-data audit QC
#
# ============================================================

sf4a_exported_audit <- read.csv(
  
  sf4a_table_file_v4,
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


export_reference <- sf4a_exported_audit %>%
  
  dplyr::select(
    all_of(required_cols)
  ) %>%
  
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  )


if (
  !isTRUE(
    all.equal(
      export_reference,
      sf4a_locked_reference,
      check.attributes = FALSE,
      tolerance = 1e-12
    )
  )
) {
  
  stop(
    "Exported SF4A audit data do not match locked reference data."
  )
  
}


cat(
  "\nLocked-data export audit: PASS\n"
)


# ============================================================
#
# 32. Explicit rho / FDR integrity checks
#
# ============================================================

rho_match <- isTRUE(
  
  all.equal(
    
    as.numeric(
      export_reference$Spearman_rho
    ),
    
    as.numeric(
      sf4a_locked_reference$Spearman_rho
    ),
    
    tolerance = 1e-12
    
  )
  
)


fdr_match <- isTRUE(
  
  all.equal(
    
    as.numeric(
      export_reference$FDR
    ),
    
    as.numeric(
      sf4a_locked_reference$FDR
    ),
    
    tolerance = 1e-12
    
  )
  
)


if (!rho_match) {
  
  stop(
    "Spearman's rho integrity QC FAILED."
  )
  
}


if (!fdr_match) {
  
  stop(
    "FDR integrity QC FAILED."
  )
  
}


cat(
  "Spearman's rho integrity: PASS\n"
)

cat(
  "BH-adjusted FDR integrity: PASS\n"
)


# ============================================================
#
# 33. Module / trait order QC
#
# ============================================================

export_module_order <- unique(
  
  export_reference$Module
  
)


export_trait_order <- unique(
  
  export_reference$Trait
  
)


if (
  !identical(
    export_module_order,
    unique(
      sf4a_locked_reference$Module
    )
  )
) {
  
  stop(
    "Module order integrity QC FAILED."
  )
  
}


if (
  !identical(
    export_trait_order,
    unique(
      sf4a_locked_reference$Trait
    )
  )
) {
  
  stop(
    "Trait order integrity QC FAILED."
  )
  
}


cat(
  "Module order integrity: PASS\n"
)

cat(
  "Trait order integrity: PASS\n"
)


# ============================================================
#
# 34. ME2 / ME9 annotation QC
#
# ============================================================

me2_present <- "ME2" %in%
  sf4a_plot_df$Module

me9_present <- "ME9" %in%
  sf4a_plot_df$Module


if (!me2_present) {
  
  stop(
    "ME2 highlight annotation missing."
  )
  
}


if (!me9_present) {
  
  stop(
    "ME9 highlight annotation missing."
  )
  
}


me2_highlight_ok <- all(
  
  sf4a_plot_df$Highlight[
    as.character(sf4a_plot_df$Module) == "ME2"
  ] == "ME2"
  
)


me9_highlight_ok <- all(
  
  sf4a_plot_df$Highlight[
    as.character(sf4a_plot_df$Module) == "ME9"
  ] == "ME9"
  
)


if (!me2_highlight_ok) {
  
  stop(
    "ME2 highlight integrity QC FAILED."
  )
  
}


if (!me9_highlight_ok) {
  
  stop(
    "ME9 highlight integrity QC FAILED."
  )
  
}


cat(
  "ME2 highlight integrity: PASS\n"
)

cat(
  "ME9 highlight integrity: PASS\n"
)


# ============================================================
#
# 35. Rendering warning QC
#
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Rendering warning QC\n")
cat("------------------------------------------------------------\n")


if (length(all_export_warnings) > 0) {
  
  cat(
    "WARNING: ",
    length(all_export_warnings),
    " rendering/export warning(s) detected.\n",
    sep = ""
  )
  
  cat("\nUnique warnings:\n")
  
  print(
    unique(all_export_warnings)
  )
  
  stop(
    "SF4A FINAL v4 rendering warning QC FAILED."
  )
  
} else {
  
  cat(
    "Rendering/export warnings: NONE\n"
  )
  
}


# ============================================================
#
# 36. Final QC summary
#
# ============================================================

cat("\n============================================================\n")
cat("SF4A FINAL v4 QC SUMMARY\n")
cat("============================================================\n")

cat(
  "Locked data preserved             : PASS\n"
)

cat(
  "Spearman's rho preserved          : PASS\n"
)

cat(
  "BH-adjusted FDR preserved         : PASS\n"
)

cat(
  "ME2 annotation preserved          : PASS\n"
)

cat(
  "ME9 annotation preserved          : PASS\n"
)

cat(
  "Module order preserved            : PASS\n"
)

cat(
  "Trait order preserved             : PASS\n"
)

cat(
  "Main heatmap rho legend removed   : PASS\n"
)

cat(
  "Single rho colorbar source        : PASS\n"
)

cat(
  "FDR explanation                   : PASS\n"
)

cat(
  "Cairo PDF generated               : PASS\n"
)

cat(
  "Audit table generated             : PASS\n"
)

cat(
  "Rendering warnings                : NONE\n"
)


# ============================================================
#
# 37. Final status
#
# ============================================================

cat("\n============================================================\n")
cat("SF4A FINAL v4 STATUS: PASS\n")
cat("============================================================\n")

cat(
  "Cairo vector PDF master completed.\n"
)

cat(
  "Main heatmap contains no colorbar.\n"
)

cat(
  "Exactly one rho colorbar is included in the right panel.\n"
)

cat(
  "FDR explanation is positioned above the rho colorbar.\n"
)

cat(
  "ME2 = #4472C4.\n"
)

cat(
  "ME9 = #C000C0.\n"
)

cat(
  "No statistical values were recalculated.\n"
)

cat(
  "No locked WGCNA results were modified.\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nFinal PDF:\n",
  normalizePath(
    sf4a_pdf_v4,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "\nAudit CSV:\n",
  normalizePath(
    sf4a_table_file_v4,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "\nWarning audit:\n",
  normalizePath(
    sf4a_warning_file_v4,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("============================================================\n")



##Supplementary Figure 4B
============================================================
# Supplementary Figure 4B
# SF4B-1 FINAL
# Locked core-gene module positioning table
# ============================================================

cat("\n============================================================\n")
cat("SF4B-1 FINAL: Core-gene module positioning\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Use locked positioning object
# ------------------------------------------------------------

if (exists("core_position_final")) {
  
  sf4b_df <- core_position_final
  
} else {
  
  stop(
    "Locked object 'core_position_final' is not available."
  )
  
}


# ------------------------------------------------------------
# 2. Validate required columns
# ------------------------------------------------------------

required_cols <- c(
  "Symbol",
  "Gene",
  "Module",
  "ME",
  "Module_Size",
  "MM",
  "MM_percentile",
  "MM_rank"
)

missing_cols <- setdiff(
  required_cols,
  names(sf4b_df)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      "Required columns are missing:",
      paste(missing_cols, collapse = ", ")
    )
  )
  
}


# ------------------------------------------------------------
# 3. Locked core genes
# ------------------------------------------------------------

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ------------------------------------------------------------
# 4. Extract core genes using verified Symbol column
# ------------------------------------------------------------

sf4b_core <- sf4b_df %>%
  dplyr::filter(
    Symbol %in% core_genes_locked
  )


cat(
  "Core genes found:",
  nrow(sf4b_core),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Validate exact three genes
# ------------------------------------------------------------

if (
  nrow(sf4b_core) != 3 ||
  !all(
    core_genes_locked %in% sf4b_core$Symbol
  )
) {
  
  stop(
    "SF4B validation failed: AOC3, MS4A7, and IL7R were not all found."
  )
  
}


# ------------------------------------------------------------
# 6. Check for duplicated core genes
# ------------------------------------------------------------

if (
  anyDuplicated(sf4b_core$Symbol) > 0
) {
  
  stop(
    "SF4B validation failed: duplicated core-gene entries detected."
  )
  
}


# ------------------------------------------------------------
# 7. Build locked plotting table
# ------------------------------------------------------------

sf4b_plot_df <- sf4b_core %>%
  dplyr::transmute(
    Gene = Symbol,
    Module = Module,
    ME = ME,
    Module_Size = Module_Size,
    MM = MM,
    MM_Rank = MM_rank,
    MM_Percentile = MM_percentile
  )


# ------------------------------------------------------------
# 8. Set biological display order
# ------------------------------------------------------------

sf4b_plot_df$Gene <- factor(
  sf4b_plot_df$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)


# ------------------------------------------------------------
# 9. Arrange in display order
# ------------------------------------------------------------

sf4b_plot_df <- sf4b_plot_df %>%
  dplyr::arrange(Gene)


# ------------------------------------------------------------
# 10. Print final locked table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SF4B FINAL LOCKED POSITIONING TABLE\n")
cat("============================================================\n\n")

print(
  sf4b_plot_df
)


# ------------------------------------------------------------
# 11. Compact validation summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("SF4B LOCKED VALUES\n")
cat("------------------------------------------------------------\n")

print(
  sf4b_plot_df %>%
    dplyr::select(
      Gene,
      Module,
      ME,
      Module_Size,
      MM_Rank,
      MM_Percentile
    )
)


cat("\n============================================================\n")
cat("SF4B-1 STATUS: PASS\n")
cat("Locked positioning values were used directly.\n")
cat("No MM, rank, percentile, or module assignment was recalculated.\n")
cat("============================================================\n")



# ============================================================
# Paper 1
# Supplementary Figure 4B
# SF4B-2 FINAL — Locked MM Positioning Plot
#
# Purpose:
#   Visualize the locked within-module positioning of
#   AOC3, MS4A7, and IL7R.
#
# IMPORTANT:
#   - Plotting only
#   - Uses ONLY sf4b_plot_df
#   - No WGCNA object is accessed
#   - No MM is recalculated
#   - No rank is recalculated
#   - No percentile is recalculated
#   - No module assignment is recalculated
#   - Locked values are preserved exactly
# ============================================================

# ============================================================
# Paper 1
# Supplementary Figure 4B
# SF4B-2 FINAL — Cairo PDF Master
# Locked MM Positioning Plot
# ============================================================

cat("\n============================================================\n")
cat("SF4B-2 FINAL — LOCKED MM POSITIONING PLOT\n")
cat("Cairo PDF Master / Arial\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(grid)


# ------------------------------------------------------------
# 2. Validate plotting input
# ------------------------------------------------------------

required_cols <- c(
  "Gene",
  "Module",
  "ME",
  "Module_Size",
  "MM",
  "MM_Rank",
  "MM_Percentile"
)

missing_cols <- setdiff(
  required_cols,
  names(sf4b_plot_df)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      "Missing required columns in sf4b_plot_df:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 3. Validate the locked three-gene dataset
# ------------------------------------------------------------

expected_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

if (
  !setequal(
    as.character(sf4b_plot_df$Gene),
    expected_genes
  )
) {
  
  stop(
    "sf4b_plot_df does not contain exactly AOC3, MS4A7, and IL7R."
  )
}

if (nrow(sf4b_plot_df) != 3) {
  
  stop(
    "sf4b_plot_df must contain exactly 3 rows."
  )
}

if (
  anyDuplicated(
    as.character(sf4b_plot_df$Gene)
  ) > 0
) {
  
  stop(
    "Duplicated gene detected in sf4b_plot_df."
  )
}


# ------------------------------------------------------------
# 4. Create plotting-only copy
# ------------------------------------------------------------

# IMPORTANT:
# No statistical quantity is modified or recalculated.

sf4b_final_plot_df <- sf4b_plot_df %>%
  
  dplyr::mutate(
    
    Gene = factor(
      as.character(Gene),
      levels = c(
        "AOC3",
        "MS4A7",
        "IL7R"
      )
    ),
    
    Module = as.character(Module),
    ME = as.character(ME)
    
  ) %>%
  
  dplyr::arrange(Gene)


# ------------------------------------------------------------
# 5. Locked-value validation
# ------------------------------------------------------------

expected_locked <- data.frame(
  
  Gene = c(
    "AOC3",
    "IL7R",
    "MS4A7"
  ),
  
  Module = c(
    "magenta",
    "blue",
    "blue"
  ),
  
  ME = c(
    "ME9",
    "ME2",
    "ME2"
  ),
  
  Module_Size = c(
    534,
    1798,
    1798
  ),
  
  MM = c(
    0.8655922,
    0.8163713,
    0.7498121
  ),
  
  MM_Rank = c(
    7,
    120,
    267
  ),
  
  MM_Percentile = c(
    98.87430,
    93.37785,
    85.19755
  ),
  
  stringsAsFactors = FALSE
)


check_df <- sf4b_final_plot_df %>%
  
  dplyr::mutate(
    Gene = as.character(Gene)
  ) %>%
  
  dplyr::select(
    Gene,
    Module,
    ME,
    Module_Size,
    MM,
    MM_Rank,
    MM_Percentile
  ) %>%
  
  dplyr::arrange(Gene)


# ------------------------------------------------------------
# 6. Exact locked-value assertions
# ------------------------------------------------------------

if (
  !identical(
    as.character(check_df$Gene),
    as.character(expected_locked$Gene)
  )
) {
  stop("Locked Gene values do not match.")
}


if (
  !identical(
    as.character(check_df$Module),
    as.character(expected_locked$Module)
  )
) {
  stop("Locked Module values do not match.")
}


if (
  !identical(
    as.character(check_df$ME),
    as.character(expected_locked$ME)
  )
) {
  stop("Locked ME values do not match.")
}


if (
  !isTRUE(
    all.equal(
      as.numeric(check_df$Module_Size),
      as.numeric(expected_locked$Module_Size),
      tolerance = 0
    )
  )
) {
  stop("Locked Module_Size values do not match.")
}


if (
  !isTRUE(
    all.equal(
      as.numeric(check_df$MM),
      as.numeric(expected_locked$MM),
      tolerance = 1e-7
    )
  )
) {
  stop("Locked MM values do not match.")
}


if (
  !isTRUE(
    all.equal(
      as.numeric(check_df$MM_Rank),
      as.numeric(expected_locked$MM_Rank),
      tolerance = 0
    )
  )
) {
  stop("Locked MM rank values do not match.")
}


if (
  !isTRUE(
    all.equal(
      as.numeric(check_df$MM_Percentile),
      as.numeric(expected_locked$MM_Percentile),
      tolerance = 1e-5
    )
  )
) {
  stop("Locked MM percentile values do not match.")
}


cat("------------------------------------------------------------\n")
cat("LOCKED-VALUE ASSERTION: PASS\n")
cat("------------------------------------------------------------\n")


# ------------------------------------------------------------
# 7. Display labels
# ------------------------------------------------------------

sf4b_final_plot_df <- sf4b_final_plot_df %>%
  
  dplyr::mutate(
    
    Position_label = paste0(
      "Rank ",
      MM_Rank,
      "/",
      Module_Size,
      "  |  ",
      sprintf("%.1f", MM_Percentile),
      "%"
    ),
    
    Module_label = paste0(
      ME,
      " / ",
      Module
    )
    
  )


# ------------------------------------------------------------
# 8. Module color mapping
# ------------------------------------------------------------

sf4b_module_colors <- c(
  "blue"    = "#4472C4",
  "magenta" = "#C000C0"
)

# ------------------------------------------------------------
# 9. Font setting
# ------------------------------------------------------------

# Same font family as the successful SF4A Cairo PDF master.

sf4b_font_family <- "Arial"


# ------------------------------------------------------------
# 10. Build final SF4B plot
# ------------------------------------------------------------

p_sf4b <- ggplot(
  
  sf4b_final_plot_df,
  
  aes(
    x = MM_Percentile,
    y = Gene
  )
  
) +
  
  
  # ----------------------------------------------------------
# Reference axis
# ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  linewidth = 0.45,
  linetype = "solid"
) +
  
  
  # ----------------------------------------------------------
# Light reference grid
# ----------------------------------------------------------

geom_segment(
  data = data.frame(
    x = c(20, 40, 60, 80, 100)
  ),
  aes(
    x = x,
    xend = x,
    y = 0.55,
    yend = 3.45
  ),
  inherit.aes = FALSE,
  linewidth = 0.25,
  linetype = "dotted"
) +
  
  
  # ----------------------------------------------------------
# Colored MM positioning bar
# ----------------------------------------------------------

geom_segment(
  aes(
    x = 0,
    xend = MM_Percentile,
    y = Gene,
    yend = Gene,
    color = Module
  ),
  linewidth = 4.2,
  lineend = "round"
) +
  
  
  # ----------------------------------------------------------
# MM percentile position point
# ----------------------------------------------------------

geom_point(
  aes(
    fill = Module
  ),
  size = 5.2,
  shape = 21,
  color = "white",
  stroke = 1.1
) +
  
  
  # ----------------------------------------------------------
# Gene labels
# ----------------------------------------------------------

geom_text(
  aes(
    x = -2,
    label = as.character(Gene)
  ),
  hjust = 1,
  fontface = "bold",
  family = sf4b_font_family,
  size = 4.3
) +
  
  
  # ------------------------------------------------------------
# Rank / percentile information
# ------------------------------------------------------------

geom_text(
  aes(
    x = 104,
    label = Position_label
  ),
  hjust = 0,
  family = sf4b_font_family,
  size = 3.6
) +
  
  
  # ------------------------------------------------------------
# Module label
# ------------------------------------------------------------

geom_label(
  aes(
    x = 148,
    label = Module_label,
    fill = Module
  ),
  hjust = 0,
  family = sf4b_font_family,
  fontface = "bold",
  size = 3.4,
  linewidth = 0,
  color = "white",
  label.padding = unit(
    c(
      0.20,
      0.38,
      0.20,
      0.38
    ),
    "lines"
  )
) +
  
  
  # ------------------------------------------------------------
# Module colors
# ------------------------------------------------------------

scale_color_manual(
  values = sf4b_module_colors,
  guide = "none"
) +
  
  scale_fill_manual(
    values = sf4b_module_colors,
    guide = "none"
  ) +
  
  
  # ------------------------------------------------------------
# X axis
# ------------------------------------------------------------

scale_x_continuous(
  limits = c(
    -18,
    178
  ),
  breaks = c(
    0,
    20,
    40,
    60,
    80,
    100
  ),
  expand = c(
    0,
    0
  )
) +
  
  
  # ----------------------------------------------------------
# Y axis
# ----------------------------------------------------------

scale_y_discrete(
  limits = rev(
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )
) +
  
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  x = "Module membership percentile (%)",
  y = NULL
) +
  
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 11,
  base_family = sf4b_font_family
) +
  
  
  theme(
    
    # X axis
    axis.title.x = element_text(
      family = sf4b_font_family,
      size = 11,
      margin = margin(
        t = 8
      )
    ),
    
    axis.text.x = element_text(
      family = sf4b_font_family,
      size = 10
    ),
    
    # Y axis
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    # Grid
    panel.grid = element_blank(),
    
    # Overall margins
    plot.margin = margin(
      t = 15,
      r = 8,
      b = 15,
      l = 12
    )
    
  )

# ------------------------------------------------------------
# 11. Display
# ------------------------------------------------------------

print(p_sf4b)


# ============================================================
# 12. Output directory
# ============================================================

sf4b_output_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4B",
  "FINAL_CAIRO_PDF_MASTER"
)

if (
  !dir.exists(sf4b_output_dir)
) {
  
  dir.create(
    sf4b_output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 13. Cairo PDF master
# ============================================================

sf4b_pdf <- file.path(
  sf4b_output_dir,
  "Paper1_Supplementary_Figure4B_FINAL_Cairo.pdf"
)


ggsave(
  
  filename = sf4b_pdf,
  
  plot = p_sf4b,
  
  width = 8.4,
  height = 3.9,
  units = "in",
  
  device = grDevices::cairo_pdf,
  
  bg = "white"
)


# ============================================================
# 14. Export locked plotting table
# ============================================================

sf4b_table_file <- file.path(
  sf4b_output_dir,
  "Paper1_Supplementary_Figure4B_locked_positioning_table.csv"
)


write.csv(
  check_df,
  sf4b_table_file,
  row.names = FALSE
)


# ============================================================
# 15. Plotting-data audit
# ============================================================

sf4b_audit_file <- file.path(
  sf4b_output_dir,
  "SF4B_plotting_data_audit.csv"
)


write.csv(
  check_df,
  sf4b_audit_file,
  row.names = FALSE
)


# ============================================================
# 16. Output QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SF4B FINAL OUTPUT QC\n")
cat("============================================================\n")


pdf_exists <- file.exists(sf4b_pdf)

table_exists <- file.exists(
  sf4b_table_file
)

audit_exists <- file.exists(
  sf4b_audit_file
)


if (pdf_exists) {
  
  pdf_size <- file.info(
    sf4b_pdf
  )$size
  
} else {
  
  pdf_size <- NA_real_
}


cat(
  "PDF exists: ",
  pdf_exists,
  "\n",
  sep = ""
)

cat(
  "PDF size (bytes): ",
  pdf_size,
  "\n",
  sep = ""
)

cat(
  "Locked table exists: ",
  table_exists,
  "\n",
  sep = ""
)

cat(
  "Audit table exists: ",
  audit_exists,
  "\n",
  sep = ""
)


if (
  !pdf_exists ||
  is.na(pdf_size) ||
  pdf_size < 10000 ||
  !table_exists ||
  !audit_exists
) {
  
  stop(
    "SF4B FINAL output QC FAILED."
  )
}


# ============================================================
# 17. Final status
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SF4B FINAL STATUS: PASS\n")
cat("============================================================\n")

cat(
  "Input: sf4b_plot_df only.\n"
)

cat(
  "Locked MM, rank, percentile, and module assignments were not recalculated.\n"
)

cat(
  "ME2 and ME9 use the locked display colors.\n"
)

cat(
  "Cairo PDF generated as the master figure.\n"
)

cat(
  "Arial font setting matches the successful SF4A Cairo PDF workflow.\n"
)

cat(
  "Locked positioning table and plotting audit were exported.\n"
)

cat("============================================================\n")

##Supplementary Figure 4C

##SF4C-1：check MM_all 和 GS_all
# ============================================================
# SF4C
# MM–GS integration
# ============================================================

cat("\n============================================================\n")
cat("SF4C: MM–GS relationships\n")
cat("============================================================\n")


stopifnot(
  exists("MM_all"),
  exists("GS_all")
)


cat("MM_all:\n")
cat(
  "Dimensions:",
  paste(dim(MM_all), collapse = " × "),
  "\n"
)

print(
  head(
    MM_all
  )
)


cat("\nGS_all:\n")
cat(
  "Dimensions:",
  paste(dim(GS_all), collapse = " × "),
  "\n"
)

print(
  head(
    GS_all
  )
)


cat("\nMM_all column names:\n")
print(names(MM_all))

cat("\nGS_all column names:\n")
print(names(GS_all))


# ============================================================
# SF4C-2
# Extract MM for ME2/blue and ME9/magenta
# ============================================================

cat("\n============================================================\n")
cat("SF4C-2: Extract MM and align with GS_all\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Confirm locked module correspondence
# ------------------------------------------------------------

stopifnot(
  "MMblue" %in% colnames(MM_all),
  "MMmagenta" %in% colnames(MM_all)
)

cat("ME2  -> MMblue    : confirmed\n")
cat("ME9  -> MMmagenta : confirmed\n")


# ------------------------------------------------------------
# 2. Extract MM using row names as Gene IDs
# ------------------------------------------------------------

mm_core <- data.frame(
  
  Gene = rownames(MM_all),
  
  ME2 = MM_all[, "MMblue"],
  
  ME9 = MM_all[, "MMmagenta"],
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


# ------------------------------------------------------------
# 3. Basic QC
# ------------------------------------------------------------

stopifnot(
  nrow(mm_core) == nrow(MM_all),
  all(mm_core$Gene == rownames(MM_all)),
  !anyDuplicated(mm_core$Gene)
)


cat(
  "MM core dimensions:",
  paste(dim(mm_core), collapse = " × "),
  "\n"
)

cat(
  "ME2/blue MM range:",
  paste(
    round(
      range(
        mm_core$ME2,
        na.rm = TRUE
      ),
      4
    ),
    collapse = " to "
  ),
  "\n"
)

cat(
  "ME9/magenta MM range:",
  paste(
    round(
      range(
        mm_core$ME9,
        na.rm = TRUE
      ),
      4
    ),
    collapse = " to "
  ),
  "\n"
)


# ------------------------------------------------------------
# 4. Extract GS information
# ------------------------------------------------------------

gs_core <- GS_all %>%
  
  dplyr::select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  )


# ------------------------------------------------------------
# 5. Confirm Gene IDs are compatible
# ------------------------------------------------------------

cat(
  "MM genes:",
  nrow(mm_core),
  "\n"
)

cat(
  "GS genes:",
  nrow(gs_core),
  "\n"
)

cat(
  "Common genes:",
  length(
    intersect(
      mm_core$Gene,
      gs_core$Gene
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. Merge MM and GS
# ------------------------------------------------------------

sf4c_base <- mm_core %>%
  
  inner_join(
    gs_core,
    by = "Gene"
  )


# ------------------------------------------------------------
# 7. Final merge QC
# ------------------------------------------------------------

stopifnot(
  nrow(sf4c_base) == length(
    intersect(
      mm_core$Gene,
      gs_core$Gene
    )
  )
)

stopifnot(
  all(
    sf4c_base$Gene %in% rownames(MM_all)
  )
)

stopifnot(
  all(
    sf4c_base$Gene %in% GS_all$Gene
  )
)


cat(
  "SF4C base dimensions:",
  paste(dim(sf4c_base), collapse = " × "),
  "\n"
)

cat("SF4C-2 STATUS: PASS\n")


print(
  head(
    sf4c_base
  )
)


# ============================================================
# STEP 19C-1
# SF4C: Add locked module assignment
# ============================================================

cat("\n============================================================\n")
cat("STEP 19C-1: Add locked module assignment\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Basic integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(sf4c_base) == nrow(MM_all),
  length(moduleColors) == nrow(MM_all),
  identical(
    rownames(MM_all),
    names(moduleColors)
  )
)


# ------------------------------------------------------------
# 2. Add module assignment
# ------------------------------------------------------------

sf4c_base$Module <- unname(moduleColors)


# ------------------------------------------------------------
# 3. Check module distribution
# ------------------------------------------------------------

cat("\nModule distribution:\n")
print(
  sort(
    table(sf4c_base$Module),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 4. Check blue and magenta module sizes
# ------------------------------------------------------------

cat("\nBlue module size:\n")
print(
  sum(sf4c_base$Module == "blue")
)

cat("\nMagenta module size:\n")
print(
  sum(sf4c_base$Module == "magenta")
)


# ------------------------------------------------------------
# 5. Create module-specific datasets
# ------------------------------------------------------------

sf4c_blue <- sf4c_base[
  sf4c_base$Module == "blue",
  ,
  drop = FALSE
]

sf4c_magenta <- sf4c_base[
  sf4c_base$Module == "magenta",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Basic checks
# ------------------------------------------------------------

cat("\nBlue dataset dimensions:\n")
print(dim(sf4c_blue))

cat("\nMagenta dataset dimensions:\n")
print(dim(sf4c_magenta))


# ------------------------------------------------------------
# 7. Confirm MM values are the correct module-specific MM
# ------------------------------------------------------------

cat("\nBlue module MM summary:\n")
print(
  summary(sf4c_blue$ME2)
)

cat("\nMagenta module MM summary:\n")
print(
  summary(sf4c_magenta$ME9)
)


# ------------------------------------------------------------
# 8. Check missing values
# ------------------------------------------------------------

cat("\nMissing values in blue dataset:\n")
print(
  colSums(
    is.na(
      sf4c_blue[
        ,
        c(
          "ME2",
          "GS_Immune",
          "GS_Stromal",
          "GS_ESTIMATE"
        )
      ]
    )
  )
)

cat("\nMissing values in magenta dataset:\n")
print(
  colSums(
    is.na(
      sf4c_magenta[
        ,
        c(
          "ME9",
          "GS_Immune",
          "GS_Stromal",
          "GS_ESTIMATE"
        )
      ]
    )
  )
)


cat("\nSTEP 19C-1 completed.\n")


# ============================================================
# STEP 19C-2
# SF4C: Within-module MM-GS concordance analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 19C-2: Within-module MM-GS concordance analysis\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define the six locked analyses
# ------------------------------------------------------------

sf4c_tests <- data.frame(
  Module = c(
    "blue",
    "blue",
    "blue",
    "magenta",
    "magenta",
    "magenta"
  ),
  Module_Label = c(
    "ME2 / blue",
    "ME2 / blue",
    "ME2 / blue",
    "ME9 / magenta",
    "ME9 / magenta",
    "ME9 / magenta"
  ),
  MM = c(
    "ME2",
    "ME2",
    "ME2",
    "ME9",
    "ME9",
    "ME9"
  ),
  GS = c(
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Function for one Spearman correlation
# ------------------------------------------------------------

run_sf4c_spearman <- function(
    df,
    mm_col,
    gs_col,
    module_name
) {
  
  x <- df[[mm_col]]
  y <- df[[gs_col]]
  
  keep <- complete.cases(x, y)
  
  x <- x[keep]
  y <- y[keep]
  
  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    Module = module_name,
    MM = mm_col,
    GS = gs_col,
    N = length(x),
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 3. Run all six tests
# ------------------------------------------------------------

result_list <- list()

result_list[[1]] <- run_sf4c_spearman(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_Immune",
  module_name = "blue"
)

result_list[[2]] <- run_sf4c_spearman(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_Stromal",
  module_name = "blue"
)

result_list[[3]] <- run_sf4c_spearman(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_ESTIMATE",
  module_name = "blue"
)

result_list[[4]] <- run_sf4c_spearman(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_Immune",
  module_name = "magenta"
)

result_list[[5]] <- run_sf4c_spearman(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_Stromal",
  module_name = "magenta"
)

result_list[[6]] <- run_sf4c_spearman(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_ESTIMATE",
  module_name = "magenta"
)


sf4c_results <- do.call(
  rbind,
  result_list
)


# ------------------------------------------------------------
# 4. BH correction across all six tests jointly
# ------------------------------------------------------------

sf4c_results$FDR_BH <- p.adjust(
  sf4c_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 5. Add module labels
# ------------------------------------------------------------

sf4c_results$Module_Label <- ifelse(
  sf4c_results$Module == "blue",
  "ME2 / blue",
  "ME9 / magenta"
)


# ------------------------------------------------------------
# 6. Add significance category
# ------------------------------------------------------------

sf4c_results$Significance <- ifelse(
  sf4c_results$FDR_BH < 0.001,
  "***",
  ifelse(
    sf4c_results$FDR_BH < 0.01,
    "**",
    ifelse(
      sf4c_results$FDR_BH < 0.05,
      "*",
      ""
    )
  )
)


# ------------------------------------------------------------
# 7. Reorder columns
# ------------------------------------------------------------

sf4c_results <- sf4c_results[
  ,
  c(
    "Module",
    "Module_Label",
    "MM",
    "GS",
    "N",
    "Spearman_rho",
    "P_value",
    "FDR_BH",
    "Significance"
  )
]


# ------------------------------------------------------------
# 8. Print complete results
# ------------------------------------------------------------

cat("\nSF4C MM-GS correlation results:\n\n")

print(
  sf4c_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Compact numerical output
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Compact SF4C results\n")
cat("------------------------------------------------------------\n")

print(
  sf4c_results[
    ,
    c(
      "Module_Label",
      "MM",
      "GS",
      "N",
      "Spearman_rho",
      "P_value",
      "FDR_BH"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Basic integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(sf4c_results) == 6,
  all(sf4c_results$N[c(1:3)] == 1798),
  all(sf4c_results$N[c(4:6)] == 534),
  all(is.finite(sf4c_results$Spearman_rho)),
  all(is.finite(sf4c_results$P_value)),
  all(is.finite(sf4c_results$FDR_BH))
)


cat("\nSTEP 19C-2 completed successfully.\n")

##STEP 19C-3：SF4C statistical audit + lock
# ============================================================
# STEP 19C-3
# SF4C: Statistical audit + lock
# ============================================================

cat("\n============================================================\n")
cat("STEP 19C-3: SF4C statistical audit + lock\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Create audit directory
# ------------------------------------------------------------

sf4c_audit_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4C"
)

dir.create(
  sf4c_audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Final integrity checks
# ------------------------------------------------------------

stopifnot(
  exists("sf4c_results"),
  exists("sf4c_blue"),
  exists("sf4c_magenta")
)

stopifnot(
  nrow(sf4c_results) == 6,
  nrow(sf4c_blue) == 1798,
  nrow(sf4c_magenta) == 534
)

stopifnot(
  all(sf4c_results$N[c(1:3)] == 1798),
  all(sf4c_results$N[c(4:6)] == 534)
)

stopifnot(
  all(sf4c_results$Module_Label == c(
    "ME2 / blue",
    "ME2 / blue",
    "ME2 / blue",
    "ME9 / magenta",
    "ME9 / magenta",
    "ME9 / magenta"
  ))
)

stopifnot(
  all(sf4c_results$GS == c(
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ))
)

stopifnot(
  all(is.finite(sf4c_results$Spearman_rho)),
  all(is.finite(sf4c_results$P_value)),
  all(is.finite(sf4c_results$FDR_BH))
)


# ------------------------------------------------------------
# 2. Verify locked source-object dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(MM_all) == 31609,
  ncol(MM_all) == 50,
  nrow(GS_all) == 31609,
  ncol(GS_all) == 13,
  length(moduleColors) == 31609
)


# ------------------------------------------------------------
# 3. Verify module sizes against locked WGCNA structure
# ------------------------------------------------------------

stopifnot(
  sum(moduleColors == "blue") == 1798,
  sum(moduleColors == "magenta") == 534
)


# ------------------------------------------------------------
# 4. Verify the six reported rho values
# ------------------------------------------------------------

expected_rho <- c(
  0.9433047,
  0.7068876,
  0.9081796,
  0.4023558,
  0.5221497,
  0.5088140
)

stopifnot(
  isTRUE(
    all.equal(
      sf4c_results$Spearman_rho,
      expected_rho,
      tolerance = 1e-7
    )
  )
)


# ------------------------------------------------------------
# 5. Verify BH-adjusted FDR values
# ------------------------------------------------------------

expected_fdr <- c(
  0,
  4.394097e-272,
  0,
  3.373681e-22,
  1.660438e-38,
  1.984393e-36
)

stopifnot(
  isTRUE(
    all.equal(
      sf4c_results$FDR_BH,
      expected_fdr,
      tolerance = 1e-7
    )
  )
)


# ------------------------------------------------------------
# 6. Create compact locked audit table
# ------------------------------------------------------------

sf4c_audit <- sf4c_results[
  ,
  c(
    "Module_Label",
    "MM",
    "GS",
    "N",
    "Spearman_rho",
    "P_value",
    "FDR_BH",
    "Significance"
  )
]

sf4c_audit$Analysis <- "Within-module Spearman MM-GS correlation"
sf4c_audit$Multiple_testing <- "BH across all six tests"

sf4c_audit <- sf4c_audit[
  ,
  c(
    "Analysis",
    "Module_Label",
    "MM",
    "GS",
    "N",
    "Spearman_rho",
    "P_value",
    "FDR_BH",
    "Significance",
    "Multiple_testing"
  )
]


# ------------------------------------------------------------
# 7. Save statistical audit table
# ------------------------------------------------------------

sf4c_audit_csv <- file.path(
  sf4c_audit_dir,
  "SF4C_statistical_audit_LOCKED_20260914.csv"
)

write.csv(
  sf4c_audit,
  sf4c_audit_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 8. Save locked SF4C objects
# ------------------------------------------------------------

sf4c_locked_rds <- file.path(
  sf4c_audit_dir,
  "SF4C_statistics_LOCKED_20260914.rds"
)

saveRDS(
  list(
    sf4c_results = sf4c_results,
    sf4c_blue = sf4c_blue,
    sf4c_magenta = sf4c_magenta
  ),
  sf4c_locked_rds
)


# ------------------------------------------------------------
# 9. Create concise audit summary
# ------------------------------------------------------------

sf4c_summary_file <- file.path(
  sf4c_audit_dir,
  "SF4C_statistical_audit_summary_20260914.txt"
)

summary_lines <- c(
  "Paper 1 — Supplementary Figure 4C",
  "SF4C statistical audit + lock",
  "Date: 2026-09-14",
  "",
  "Analysis:",
  "Within-module Spearman MM-GS concordance analysis.",
  "",
  "Source objects:",
  "MM_all: 31,609 genes x 50 modules",
  "GS_all: 31,609 genes x 13 columns",
  "moduleColors: 31,609 genes",
  "",
  "Module-specific analysis:",
  "ME2 / blue: 1,798 genes",
  "ME9 / magenta: 534 genes",
  "",
  "Multiple-testing correction:",
  "Benjamini-Hochberg correction applied jointly across all six tests.",
  "",
  "Locked results:",
  sprintf(
    "ME2 / blue - ImmuneScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[1],
    sf4c_results$FDR_BH[1]
  ),
  sprintf(
    "ME2 / blue - StromalScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[2],
    sf4c_results$FDR_BH[2]
  ),
  sprintf(
    "ME2 / blue - ESTIMATEScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[3],
    sf4c_results$FDR_BH[3]
  ),
  sprintf(
    "ME9 / magenta - ImmuneScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[4],
    sf4c_results$FDR_BH[4]
  ),
  sprintf(
    "ME9 / magenta - StromalScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[5],
    sf4c_results$FDR_BH[5]
  ),
  sprintf(
    "ME9 / magenta - ESTIMATEScore: rho = %.7f; FDR = %.6g",
    sf4c_results$Spearman_rho[6],
    sf4c_results$FDR_BH[6]
  ),
  "",
  "Interpretation:",
  "Strong positive MM-GS concordance was observed within the blue module.",
  "Moderate positive MM-GS concordance was observed within the magenta module.",
  "",
  "Lock status:",
  "PASS / LOCKED",
  "",
  "No WGCNA parameters, MM values, GS values, or module assignments were recalculated.",
  "SF4C statistics are derived from the locked WGCNA objects."
)

writeLines(
  summary_lines,
  sf4c_summary_file
)


# ------------------------------------------------------------
# 10. Final file check
# ------------------------------------------------------------

stopifnot(
  file.exists(sf4c_audit_csv),
  file.exists(sf4c_locked_rds),
  file.exists(sf4c_summary_file)
)


# ------------------------------------------------------------
# 11. Report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 19C-3 STATUS: PASS / LOCKED\n")
cat("============================================================\n")

cat("\nLocked statistical audit:\n")
print(sf4c_audit, row.names = FALSE)

cat("\nFiles created:\n")
cat("1.", sf4c_audit_csv, "\n")
cat("2.", sf4c_locked_rds, "\n")
cat("3.", sf4c_summary_file, "\n")

cat("\nSTEP 19C-3 completed successfully.\n")



#STEP 19C-4：SF4C Cairo PDF Final Plot
# ============================================================
# PAPER 1
# STEP 19C-4
# Supplementary Figure 4C
# Within-module MM-GS concordance
#
# FINAL CAIRO PDF MASTER
# ============================================================


cat("\n============================================================\n")
cat("STEP 19C-4: SF4C Final Cairo PDF plotting\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)
library(Cairo)


# ------------------------------------------------------------
# 1. Define locked input
# ------------------------------------------------------------

sf4c_locked_file <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4C",
  "SF4C_statistics_LOCKED_20260914.rds"
)


stopifnot(
  file.exists(sf4c_locked_file)
)


# ------------------------------------------------------------
# 2. Read ONLY locked SF4C RDS
# ------------------------------------------------------------

sf4c_locked <- readRDS(
  sf4c_locked_file
)


sf4c_results <- sf4c_locked$sf4c_results
sf4c_blue <- sf4c_locked$sf4c_blue
sf4c_magenta <- sf4c_locked$sf4c_magenta


# ------------------------------------------------------------
# 3. Integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(sf4c_results) == 6,
  nrow(sf4c_blue) == 1798,
  nrow(sf4c_magenta) == 534
)


stopifnot(
  identical(
    sf4c_results$Module_Label,
    c(
      "ME2 / blue",
      "ME2 / blue",
      "ME2 / blue",
      "ME9 / magenta",
      "ME9 / magenta",
      "ME9 / magenta"
    )
  )
)


stopifnot(
  identical(
    sf4c_results$GS,
    c(
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE",
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE"
    )
  )
)


# ------------------------------------------------------------
# 4. Fixed module colors
# ------------------------------------------------------------

sf4c_module_colors <- c(
  "blue" = "#4472C4",
  "magenta" = "#C000C0"
)


# ------------------------------------------------------------
# 5. Create plotting function
# ------------------------------------------------------------

plot_sf4c_panel <- function(
    df,
    mm_col,
    gs_col,
    module_name,
    panel_label,
    gs_label,
    result_row
) {
  
  rho_value <- sf4c_results$Spearman_rho[result_row]
  fdr_value <- sf4c_results$FDR_BH[result_row]
  n_value <- sf4c_results$N[result_row]
  
  module_color <- sf4c_module_colors[module_name]
  
  
  # ----------------------------------------------------------
  # FDR display
  # ----------------------------------------------------------
  
  fdr_text <- if (
    fdr_value == 0
  ) {
    "FDR < 1e-300"
  } else if (
    fdr_value < 0.001
  ) {
    "FDR < 0.001"
  } else {
    paste0(
      "FDR = ",
      formatC(
        fdr_value,
        format = "e",
        digits = 2
      )
    )
  }
  
  
  annotation_text <- paste0(
    "Spearman's \u03c1 = ",
    sprintf(
      "%.3f",
      rho_value
    ),
    "\n",
    fdr_text,
    "\n",
    "N = ",
    format(
      n_value,
      big.mark = ","
    )
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    df,
    aes(
      x = .data[[mm_col]],
      y = .data[[gs_col]]
    )
  ) +
    geom_point(
      color = module_color,
      alpha = 0.35,
      size = 0.55
    ) +
    
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.50,
      color = module_color
    )+
    
    annotate(
      "text",
      x = -0.53,
      y = 0.91,
      label = annotation_text,
      hjust = 0,
      vjust = 1,
      family = "Arial",
      fontface = "plain",
      size = 3.15,
      lineheight = 1.05
    ) +
    
    labs(
      title = paste0(
        panel_label,
        "  ",
        module_name,
        " module — ",
        gs_label
      ),
      x = "Module Membership (MM)",
      y = "Gene Significance (GS)"
    ) +
    
    coord_cartesian(
      xlim = c(-0.6, 1.0),
      ylim = c(-1.0, 1.0),
      expand = FALSE
    ) +
    
    theme_classic(
      base_family = "Arial"
    ) +
    
    theme(
      plot.title = element_text(
        family = "Arial",
        size = 9.8,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 6
        )
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 9.5,
        face = "plain"
      ),
      
      axis.text = element_text(
        family = "Arial",
        size = 8.5,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.45
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4
      ),
      
      panel.grid = element_blank(),
      
      plot.margin = margin(
        t = 5,
        r = 12,
        b = 5,
        l = 5
      )
    )
  
  return(p)
}


# ------------------------------------------------------------
# 6. Generate six panels
# ------------------------------------------------------------

p_sf4c_A <- plot_sf4c_panel(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_Immune",
  module_name = "blue",
  panel_label = "(a)",
  gs_label = "ImmuneScore",
  result_row = 1
)


p_sf4c_B <- plot_sf4c_panel(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_Stromal",
  module_name = "blue",
  panel_label = "(b)",
  gs_label = "StromalScore",
  result_row = 2
)


p_sf4c_C <- plot_sf4c_panel(
  df = sf4c_blue,
  mm_col = "ME2",
  gs_col = "GS_ESTIMATE",
  module_name = "blue",
  panel_label = "(c)",
  gs_label = "ESTIMATEScore",
  result_row = 3
)


p_sf4c_D <- plot_sf4c_panel(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_Immune",
  module_name = "magenta",
  panel_label = "(d)",
  gs_label = "ImmuneScore",
  result_row = 4
)


p_sf4c_E <- plot_sf4c_panel(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_Stromal",
  module_name = "magenta",
  panel_label = "(e)",
  gs_label = "StromalScore",
  result_row = 5
)


p_sf4c_F <- plot_sf4c_panel(
  df = sf4c_magenta,
  mm_col = "ME9",
  gs_col = "GS_ESTIMATE",
  module_name = "magenta",
  panel_label = "(f)",
  gs_label = "ESTIMATEScore",
  result_row = 6
)


# ------------------------------------------------------------
# 7. Assemble final six-panel figure
# ------------------------------------------------------------

sf4c_final <- (
  p_sf4c_A | p_sf4c_B | p_sf4c_C
) /
  (
    p_sf4c_D | p_sf4c_E | p_sf4c_F
  )


sf4c_final <- sf4c_final +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  )


# ------------------------------------------------------------
# 8. Output directory
# ------------------------------------------------------------

sf4c_plot_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4C"
)

dir.create(
  sf4c_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# STEP 19C-4
# SF4C — Final Cairo PDF export
# FIX: CairoPDF() uses unnamed file argument
# ============================================================


# ------------------------------------------------------------
# 1. Output file
# ------------------------------------------------------------

sf4c_pdf <- file.path(
  sf4c_plot_dir,
  "Paper1_Supplementary_Figure4C_FINAL_Cairo.pdf"
)


# ------------------------------------------------------------
# 2. Remove old incomplete file if present
# ------------------------------------------------------------

if (file.exists(sf4c_pdf)) {
  file.remove(sf4c_pdf)
}


# ------------------------------------------------------------
# 3. Open Cairo PDF device
# ------------------------------------------------------------

CairoPDF(
  sf4c_pdf,
  width = 8.6,
  height = 6.7,
  family = "Arial"
)


# ------------------------------------------------------------
# 4. Draw final six-panel figure
# ------------------------------------------------------------

print(
  sf4c_final
)


# ------------------------------------------------------------
# 5. Close device
# ------------------------------------------------------------

dev.off()


# ------------------------------------------------------------
# 6. File integrity check
# ------------------------------------------------------------

stopifnot(
  file.exists(sf4c_pdf)
)


sf4c_file_info <- file.info(
  sf4c_pdf
)


stopifnot(
  !is.na(sf4c_file_info$size),
  sf4c_file_info$size > 10000
)


# ------------------------------------------------------------
# 7. Final status
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "STEP 19C-4: SF4C Cairo PDF export\n",
  "============================================================\n",
  "\nOutput file:\n",
  sf4c_pdf,
  "\n\nFile size:",
  round(sf4c_file_info$size / 1024, 1),
  "KB\n",
  "\nSTATUS: PASS / LOCKED\n",
  "============================================================\n"
)


##Supplementary Figure 4D
##SF4D-1：Prepare three-gene GS data
# ============================================================
# ============================================================
# SF4D-1
# Gene-level TME associations
# Correct extraction using locked ENSEMBL IDs
# ============================================================

cat("\n============================================================\n")
cat("SF4D-1: Gene-level TME associations\n")
cat("Using locked GS_all values\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check GS_all
# ------------------------------------------------------------

stopifnot(
  exists("GS_all"),
  nrow(GS_all) == 31609
)


# ------------------------------------------------------------
# 2. Locked core-gene mapping
# ------------------------------------------------------------

sf4d_gene_map <- data.frame(
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  ENSEMBL = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  stringsAsFactors = FALSE
)


print(sf4d_gene_map)


# ------------------------------------------------------------
# 3. Extract locked GS values by ENSEMBL ID
# ------------------------------------------------------------

sf4d_df <- GS_all %>%
  filter(
    Gene %in% sf4d_gene_map$ENSEMBL
  ) %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  )


# ------------------------------------------------------------
# 4. Restore gene symbols
# ------------------------------------------------------------

sf4d_df <- sf4d_gene_map %>%
  left_join(
    sf4d_df,
    by = c(
      "ENSEMBL" = "Gene"
    )
  )


# ------------------------------------------------------------
# 5. Integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(sf4d_df) == 3
)

stopifnot(
  identical(
    sf4d_df$Gene,
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )
)

stopifnot(
  all(
    complete.cases(
      sf4d_df[
        ,
        c(
          "GS_Immune",
          "GS_Stromal",
          "GS_ESTIMATE"
        )
      ]
    )
  )
)


# ------------------------------------------------------------
# 6. Display final SF4D data
# ------------------------------------------------------------

print(
  sf4d_df
)


cat(
  "\nSF4D-1 STATUS: PASS / LOCKED\n"
)

### ============================================================
# SF4D-2
# Prepare locked plotting table
# Gene-level TME associations
# ============================================================

cat("\n============================================================\n")
cat("SF4D-2: Prepare locked SF4D plotting table\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check locked source
# ------------------------------------------------------------

stopifnot(
  nrow(sf4d_df) == 3
)

stopifnot(
  identical(
    sf4d_df$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )
)


# ------------------------------------------------------------
# 2. Convert to long format
# ------------------------------------------------------------

sf4d_plot_df <- sf4d_df %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      GS_Immune,
      GS_Stromal,
      GS_ESTIMATE
    ),
    names_to = "Trait",
    values_to = "GS"
  )


# ------------------------------------------------------------
# 3. Convert trait labels
# ------------------------------------------------------------

sf4d_plot_df$Trait <- factor(
  sf4d_plot_df$Trait,
  levels = c(
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ),
  labels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)


# ------------------------------------------------------------
# 4. Preserve gene order
# ------------------------------------------------------------

sf4d_plot_df$Gene <- factor(
  sf4d_plot_df$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)


# ------------------------------------------------------------
# 5. Absolute GS
# ------------------------------------------------------------

sf4d_plot_df$Abs_GS <- abs(
  sf4d_plot_df$GS
)


# ------------------------------------------------------------
# 6. Integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(sf4d_plot_df) == 9
)

stopifnot(
  all(
    complete.cases(
      sf4d_plot_df[, c("Gene", "Trait", "GS", "Abs_GS")]
    )
  )
)

stopifnot(
  all(sf4d_plot_df$GS > 0)
)

stopifnot(
  all(
    sf4d_plot_df$Abs_GS ==
      abs(sf4d_plot_df$GS)
  )
)


# ------------------------------------------------------------
# 7. Display plotting table
# ------------------------------------------------------------

print(
  sf4d_plot_df
)


cat(
  "\nSF4D-2 STATUS: PASS / LOCKED\n"
)

# ============================================================
# STEP 19D-3
# Supplementary Figure 4D
# Gene-level TME associations
# Final Cairo PDF
# ============================================================

cat("\n============================================================\n")
cat("STEP 19D-3: SF4D Cairo PDF plotting\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(ggplot2)
library(Cairo)


# ------------------------------------------------------------
# 2. Final plotting table integrity
# ------------------------------------------------------------

stopifnot(
  exists("sf4d_plot_df")
)

stopifnot(
  nrow(sf4d_plot_df) == 9
)

stopifnot(
  all(
    complete.cases(
      sf4d_plot_df[
        ,
        c(
          "Gene",
          "Trait",
          "GS",
          "Abs_GS"
        )
      ]
    )
  )
)


# ------------------------------------------------------------
# 3. Fixed factor order
# ------------------------------------------------------------

sf4d_plot_df$Gene <- factor(
  sf4d_plot_df$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

sf4d_plot_df$Trait <- factor(
  sf4d_plot_df$Trait,
  levels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)


# ------------------------------------------------------------
# 4. Create SF4D plot
# ------------------------------------------------------------

p_sf4d <- ggplot(
  sf4d_plot_df,
  aes(
    x = Trait,
    y = Gene
  )
) +
  
  geom_point(
    aes(
      size = Abs_GS,
      fill = GS
    ),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  
  scale_size_continuous(
    range = c(
      3.5,
      9.0
    ),
    limits = c(
      0.35,
      0.82
    ),
    breaks = c(
      0.4,
      0.6,
      0.8
    ),
    name = "|GS|"
  ) +
  
  scale_fill_gradient(
    low = "#E8F1FA",
    high = "#285A9F",
    limits = c(
      0.35,
      0.82
    ),
    name = "GS"
  )+
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    
    # --------------------------------------------------------
    # Axis text
    # --------------------------------------------------------
    
    axis.text.x = element_text(
      family = "Arial",
      size = 9,
      color = "black"
    ),
    
    axis.text.y = element_text(
      family = "Arial",
      size = 10,
      color = "black",
      face = "bold"
    ),
    
    # --------------------------------------------------------
    # Axis
    # --------------------------------------------------------
    
    axis.line = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      color = "black"
    ),
    
    # --------------------------------------------------------
    # Grid
    # --------------------------------------------------------
    
    panel.grid = element_blank(),
    
    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    legend.title = element_text(
      family = "Arial",
      size = 9,
      colour = "black"
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 8.5,
      colour = "black"
    ),
    
    legend.key.height = grid::unit(
      0.55,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.35,
      "cm"
    ),
    
    legend.spacing.y = grid::unit(
      0.08,
      "cm"
    ),
    
    legend.box = "vertical",
    
    plot.margin = margin(
      t = 5,
      r = 8,
      b = 5,
      l = 5
    )
  )




# ------------------------------------------------------------
# 5. Display plot in RStudio
# ------------------------------------------------------------

print(
  p_sf4d
)


# ------------------------------------------------------------
# 6. Output directory
# ------------------------------------------------------------

sf4d_plot_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "SF4D"
)

dir.create(
  sf4d_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 7. Final Cairo PDF
# ------------------------------------------------------------

sf4d_pdf <- file.path(
  sf4d_plot_dir,
  "Paper1_Supplementary_Figure4D_FINAL_Cairo.pdf"
)


if (file.exists(sf4d_pdf)) {
  file.remove(sf4d_pdf)
}


CairoPDF(
  sf4d_pdf,
  width = 6.6,
  height = 4.6,
  family = "Arial"
)


print(
  p_sf4d
)


dev.off()


# ------------------------------------------------------------
# 8. File integrity check
# ------------------------------------------------------------

stopifnot(
  file.exists(sf4d_pdf)
)


sf4d_file_info <- file.info(
  sf4d_pdf
)


stopifnot(
  !is.na(sf4d_file_info$size),
  sf4d_file_info$size > 10000
)


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "STEP 19D-3: SF4D Cairo PDF export\n",
  "============================================================\n",
  "\nOutput file:\n",
  sf4d_pdf,
  "\n\nFile size:",
  round(
    sf4d_file_info$size / 1024,
    1
  ),
  "KB\n",
  "\nSTATUS: PASS / LOCKED\n",
  "============================================================\n"
)


# ============================================================
# STEP 19D-1
# SUPPLEMENTARY TABLE X
# Tumor-only WGCNA characteristics of prioritized study genes
# Correct ENSEMBL-to-SYMBOL mapping
# ============================================================

cat("\n============================================================\n")
cat("STEP 19D-1: SUPPLEMENTARY TABLE X — CORE GENE WGCNA CHARACTERISTICS\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define LOCKED core-gene mapping
# ------------------------------------------------------------

core_gene_map <- data.frame(
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  ENSEMBL = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Extract LOCKED gene-level TME statistics
#    GS_all$Gene contains ENSEMBL IDs
# ------------------------------------------------------------

core_tme_table <- GS_all %>%
  filter(
    Gene %in% core_gene_map$ENSEMBL
  ) %>%
  select(
    Gene,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  ) %>%
  left_join(
    core_gene_map,
    by = c("Gene" = "ENSEMBL")
  ) %>%
  select(
    Gene,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  )


# ------------------------------------------------------------
# 3. Rename ENSEMBL-derived Gene column explicitly
# ------------------------------------------------------------

colnames(core_tme_table)[1] <- "ENSEMBL"


# ------------------------------------------------------------
# 4. Reorder to AOC3 / MS4A7 / IL7R
# ------------------------------------------------------------

core_tme_table <- core_gene_map %>%
  left_join(
    core_tme_table,
    by = "ENSEMBL"
  ) %>%
  select(
    Gene,
    ENSEMBL,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  )


# ------------------------------------------------------------
# 5. Extract LOCKED core-gene module positioning
# ------------------------------------------------------------

core_position_table <- sf4b_plot_df %>%
  select(
    Gene,
    Module,
    Module_Size,
    MM_Rank,
    MM_Percentile
  )


# ------------------------------------------------------------
# 6. Merge positioning + TME statistics
# ------------------------------------------------------------

supp_table_x <- core_position_table %>%
  left_join(
    core_tme_table,
    by = "Gene"
  )


# ------------------------------------------------------------
# 7. Reorder columns for Supplementary Table X
# ------------------------------------------------------------

supp_table_x <- supp_table_x %>%
  select(
    Gene,
    ENSEMBL,
    Module,
    Module_Size,
    MM_Rank,
    MM_Percentile,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  )


# ------------------------------------------------------------
# 8. Validation
# ------------------------------------------------------------

cat("\nSupplementary Table X:\n")
print(supp_table_x)


cat("\n------------------------------------------------------------\n")
cat("VALIDATION\n")
cat("------------------------------------------------------------\n")

cat("Rows:", nrow(supp_table_x), "\n")
cat("Columns:", ncol(supp_table_x), "\n")

cat("\nCore genes:\n")
print(supp_table_x$Gene)

cat("\nMissing GS values:\n")
print(
  supp_table_x %>%
    select(
      Gene,
      GS_Immune,
      GS_Stromal,
      GS_ESTIMATE
    )
)

cat("\nMissing-value check:\n")
print(
  colSums(
    is.na(
      supp_table_x
    )
  )
)


# ------------------------------------------------------------
# 9. LOCKED-value checks
# ------------------------------------------------------------

stopifnot(
  nrow(supp_table_x) == 3
)

stopifnot(
  identical(
    supp_table_x$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )
)

stopifnot(
  all(!is.na(supp_table_x$GS_Immune))
)

stopifnot(
  all(!is.na(supp_table_x$GS_Stromal))
)

stopifnot(
  all(!is.na(supp_table_x$GS_ESTIMATE))
)


# ------------------------------------------------------------
# 10. Explicit numerical audit against LOCKED values
# ------------------------------------------------------------

expected_gs <- data.frame(
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  GS_Immune = c(
    0.3628795,
    0.6866628,
    0.7688851
  ),
  GS_Stromal = c(
    0.6555792,
    0.6460729,
    0.7067829
  ),
  GS_ESTIMATE = c(
    0.5399233,
    0.7230010,
    0.8045778
  )
)


audit_gs <- supp_table_x %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  ) %>%
  left_join(
    expected_gs,
    by = "Gene",
    suffix = c("_actual", "_expected")
  ) %>%
  mutate(
    GS_Immune_match =
      abs(GS_Immune_actual - GS_Immune_expected) < 1e-6,
    
    GS_Stromal_match =
      abs(GS_Stromal_actual - GS_Stromal_expected) < 1e-6,
    
    GS_ESTIMATE_match =
      abs(GS_ESTIMATE_actual - GS_ESTIMATE_expected) < 1e-6
  )


print(audit_gs)


stopifnot(
  all(audit_gs$GS_Immune_match),
  all(audit_gs$GS_Stromal_match),
  all(audit_gs$GS_ESTIMATE_match)
)


cat("\n============================================================\n")
cat("STEP 19D-1 STATUS: PASS / LOCKED\n")
cat("Supplementary Table X uses LOCKED GS and module-position values.\n")
cat("No GS, MM, rank, percentile, or module assignment was recalculated.\n")
cat("============================================================\n")


# ============================================================
# STEP 19D-2
# Supplementary Table X
# Formatting and export
# ============================================================

cat("\n============================================================\n")
cat("STEP 19D-2: SUPPLEMENTARY TABLE X — FORMATTING AND EXPORT\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Output directory
# ------------------------------------------------------------

sf4_dir <- file.path(
  "Paper1_STEP19_PREP_20260914",
  "audit",
  "Supplementary_Table_X"
)

if (!dir.exists(sf4_dir)) {
  dir.create(
    sf4_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Confirm STEP 19D-1 object exists
# ------------------------------------------------------------

stopifnot(
  exists("supp_table_x")
)

stopifnot(
  nrow(supp_table_x) == 3
)

stopifnot(
  identical(
    supp_table_x$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )
)


# ------------------------------------------------------------
# 2. Create formal module labels
# ------------------------------------------------------------

supp_table_x_formal <- supp_table_x %>%
  mutate(
    Module = case_when(
      Module == "blue" ~ "ME2/blue",
      Module == "magenta" ~ "ME9/magenta",
      TRUE ~ Module
    )
  )


# ------------------------------------------------------------
# 3. Construct FINAL presentation table
#
#    Only formatting:
#    - no new statistical calculation
#    - no recalculation of GS/MM
#    - no recalculation of FDR
# ------------------------------------------------------------

supp_table_x_final <- supp_table_x_formal %>%
  select(
    Gene,
    ENSEMBL,
    Module,
    Module_Size,
    MM_Rank,
    MM_Percentile,
    GS_Immune,
    FDR_Immune,
    GS_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    FDR_ESTIMATE
  )


# ------------------------------------------------------------
# 4. Apply display formatting
# ------------------------------------------------------------

supp_table_x_display <- supp_table_x_final %>%
  mutate(
    
    # Module size / rank
    Module_Size = as.integer(Module_Size),
    MM_Rank = as.integer(MM_Rank),
    
    # MM percentile: 2 decimal places
    MM_Percentile = sprintf(
      "%.2f",
      MM_Percentile
    ),
    
    # GS: 3 decimal places
    GS_Immune = sprintf(
      "%.3f",
      GS_Immune
    ),
    
    GS_Stromal = sprintf(
      "%.3f",
      GS_Stromal
    ),
    
    GS_ESTIMATE = sprintf(
      "%.3f",
      GS_ESTIMATE
    ),
    
    # FDR: scientific notation
    FDR_Immune = formatC(
      FDR_Immune,
      format = "e",
      digits = 3
    ),
    
    FDR_Stromal = formatC(
      FDR_Stromal,
      format = "e",
      digits = 3
    ),
    
    FDR_ESTIMATE = formatC(
      FDR_ESTIMATE,
      format = "e",
      digits = 3
    )
  )


# ------------------------------------------------------------
# 5. Rename columns for formal Supplementary Table
# ------------------------------------------------------------

colnames(supp_table_x_display) <- c(
  "Gene",
  "ENSEMBL ID",
  "Module",
  "Module size",
  "MM rank",
  "MM percentile",
  "GS–Immune",
  "FDR–Immune",
  "GS–Stromal",
  "FDR–Stromal",
  "GS–ESTIMATE",
  "FDR–ESTIMATE"
)


# ------------------------------------------------------------
# 6. Print final table
# ------------------------------------------------------------

cat("\nFINAL SUPPLEMENTARY TABLE X:\n")
print(
  supp_table_x_display,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Export CSV
# ------------------------------------------------------------

csv_file <- file.path(
  sf4_dir,
  "Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes.csv"
)

write.csv(
  supp_table_x_display,
  csv_file,
  row.names = FALSE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 8. Export XLSX
# ------------------------------------------------------------

xlsx_file <- file.path(
  sf4_dir,
  "Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes.xlsx"
)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "Package 'openxlsx' is required for XLSX export."
  )
}

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  wb,
  "Supplementary Table X"
)

openxlsx::writeData(
  wb,
  sheet = "Supplementary Table X",
  x = supp_table_x_display,
  startRow = 1,
  startCol = 1,
  withFilter = FALSE
)


# ------------------------------------------------------------
# 9. XLSX formatting
# ------------------------------------------------------------

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

body_style <- openxlsx::createStyle(
  valign = "center"
)

openxlsx::addStyle(
  wb,
  sheet = "Supplementary Table X",
  style = header_style,
  rows = 1,
  cols = 1:ncol(supp_table_x_display),
  gridExpand = TRUE
)

openxlsx::addStyle(
  wb,
  sheet = "Supplementary Table X",
  style = body_style,
  rows = 2:(nrow(supp_table_x_display) + 1),
  cols = 1:ncol(supp_table_x_display),
  gridExpand = TRUE,
  stack = TRUE
)


# ------------------------------------------------------------
# 10. Set column widths
# ------------------------------------------------------------

openxlsx::setColWidths(
  wb,
  sheet = "Supplementary Table X",
  cols = 1,
  widths = 12
)

openxlsx::setColWidths(
  wb,
  sheet = "Supplementary Table X",
  cols = 2,
  widths = 20
)

openxlsx::setColWidths(
  wb,
  sheet = "Supplementary Table X",
  cols = 3,
  widths = 15
)

openxlsx::setColWidths(
  wb,
  sheet = "Supplementary Table X",
  cols = 4:6,
  widths = 15
)

openxlsx::setColWidths(
  wb,
  sheet = "Supplementary Table X",
  cols = 7:12,
  widths = 15
)


# ------------------------------------------------------------
# 11. Freeze header row
# ------------------------------------------------------------

openxlsx::freezePane(
  wb,
  sheet = "Supplementary Table X",
  firstRow = TRUE
)


# ------------------------------------------------------------
# 12. Save XLSX
# ------------------------------------------------------------

openxlsx::saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 13. Save final RDS
#
#     This preserves the formatted table object.
# ------------------------------------------------------------

rds_file <- file.path(
  sf4_dir,
  "Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes_FINAL.rds"
)

saveRDS(
  supp_table_x_display,
  rds_file
)


# ------------------------------------------------------------
# 14. Create table description / audit note
# ------------------------------------------------------------

txt_file <- file.path(
  sf4_dir,
  "Supplementary_Table_X_README.txt"
)

description_lines <- c(
  "Supplementary Table X",
  "Tumor-only WGCNA characteristics of prioritized study genes",
  "",
  "Study genes:",
  "AOC3, MS4A7, IL7R",
  "",
  "WGCNA cohort:",
  "513 TCGA-LUAD tumor samples",
  "",
  "Module labels:",
  "ME2/blue and ME9/magenta",
  "",
  "MM percentile:",
  "Displayed to two decimal places.",
  "",
  "Gene significance (GS):",
  "Displayed to three decimal places.",
  "",
  "P values / FDR:",
  "FDR values are displayed in scientific notation.",
  "",
  "Data provenance:",
  "All GS, FDR, module assignment, module size, MM rank, and MM percentile",
  "values were derived from locked STEP 19 results.",
  "",
  "No new WGCNA calculation was performed in STEP 19D-2.",
  "No GS, MM, rank, percentile, or module assignment was recalculated.",
  "",
  "Corresponding supplementary analyses:",
  "Supplementary Figure 4A–4D."
)

writeLines(
  description_lines,
  txt_file,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 15. File existence and size check
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("EXPORT CHECK\n")
cat("------------------------------------------------------------\n")

export_files <- c(
  csv_file,
  xlsx_file,
  rds_file,
  txt_file
)

export_check <- data.frame(
  File = basename(export_files),
  Exists = file.exists(export_files),
  Size_kB = round(
    file.size(export_files) / 1024,
    2
  ),
  stringsAsFactors = FALSE
)

print(export_check)


# ------------------------------------------------------------
# 16. Final content check
# ------------------------------------------------------------

stopifnot(
  all(file.exists(export_files))
)

stopifnot(
  nrow(supp_table_x_display) == 3
)

stopifnot(
  all(
    supp_table_x_display$Gene ==
      c("AOC3", "MS4A7", "IL7R")
  )
)


# ------------------------------------------------------------
# 17. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 19D-2 STATUS: PASS / LOCKED\n")
cat("============================================================\n")

cat(
  "Supplementary Table X formatted and exported successfully.\n"
)

cat(
  "CSV:", csv_file, "\n"
)

cat(
  "XLSX:", xlsx_file, "\n"
)

cat(
  "RDS:", rds_file, "\n"
)

cat(
  "README:", txt_file, "\n"
)

cat(
  "No statistical values were recalculated.\n"
)

cat(
  "============================================================\n"
)



# ============================================================
# STEP 19D-3
# END-OF-WORK SESSION WORKSPACE CHECKPOINT
# Paper 1 — Supplementary Figure 4 + Supplementary Table X
# ============================================================
# ============================================================
# STEP 19D-3 RECOVERY CHECK
# Check existing STEP 19 archived files
# ============================================================

cat("\n============================================================\n")
cat("STEP 19D-3 RECOVERY CHECK\n")
cat("============================================================\n")

prep_dir <- "Paper1_STEP19_PREP_20260914"

cat("\nPreparation directory:\n")
cat(normalizePath(prep_dir, winslash = "/", mustWork = FALSE), "\n")


# ------------------------------------------------------------
# 1. Check directory
# ------------------------------------------------------------

cat("\nDirectory exists:\n")
print(dir.exists(prep_dir))


# ------------------------------------------------------------
# 2. List RDS files
# ------------------------------------------------------------

rds_files <- list.files(
  prep_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

cat("\nRDS files found:", length(rds_files), "\n")

if (length(rds_files) > 0) {
  print(
    data.frame(
      File = normalizePath(
        rds_files,
        winslash = "/",
        mustWork = FALSE
      ),
      Size_kB = round(
        file.size(rds_files) / 1024,
        2
      ),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 3. List RData files
# ------------------------------------------------------------

rdata_files <- list.files(
  prep_dir,
  pattern = "\\.(RData|rda)$",
  recursive = TRUE,
  full.names = TRUE
)

cat("\nRData files found:", length(rdata_files), "\n")

if (length(rdata_files) > 0) {
  print(
    data.frame(
      File = normalizePath(
        rdata_files,
        winslash = "/",
        mustWork = FALSE
      ),
      Size_MB = round(
        file.size(rdata_files) / (1024^2),
        2
      ),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 4. Search specifically for STEP 19 / SF4 files
# ------------------------------------------------------------

step19_files <- list.files(
  prep_dir,
  pattern = "SF4|STEP19|Supplementary_Table_X",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\nSTEP 19 / SF4 / Supplementary Table files:\n")

if (length(step19_files) > 0) {
  print(
    normalizePath(
      step19_files,
      winslash = "/",
      mustWork = FALSE
    )
  )
} else {
  cat("NONE FOUND\n")
}


cat("\n============================================================\n")
cat("RECOVERY CHECK COMPLETED\n")
cat("============================================================\n")

# ============================================================
# STEP 19D-3
# RECOVER LOCKED WORKSPACE + STEP 19 FINAL OBJECTS
# ============================================================

cat("\n============================================================\n")
cat("STEP 19D-3: RECOVER LOCKED WORKSPACE + FINAL STEP 19 OBJECTS\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

prep_dir <- "Paper1_STEP19_PREP_20260914"

pre_step19_workspace <- file.path(
  prep_dir,
  "R_workspace",
  "Paper1_PRE_STEP19_WORKSPACE_20260914.RData"
)

sf4c_rds <- file.path(
  prep_dir,
  "audit",
  "SF4C",
  "SF4C_statistics_LOCKED_20260914.rds"
)

supp_table_rds <- file.path(
  prep_dir,
  "audit",
  "Supplementary_Table_X",
  "Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes_FINAL.rds"
)


# ------------------------------------------------------------
# 2. Confirm archive files exist
# ------------------------------------------------------------

stopifnot(
  file.exists(pre_step19_workspace)
)

stopifnot(
  file.exists(sf4c_rds)
)

stopifnot(
  file.exists(supp_table_rds)
)


cat("\nArchive files confirmed.\n")


# ------------------------------------------------------------
# 3. Load PRE-STEP19 complete workspace
# ------------------------------------------------------------

cat("\nLoading PRE-STEP19 workspace...\n")

load(
  pre_step19_workspace,
  envir = .GlobalEnv
)

cat("PRE-STEP19 workspace loaded.\n")


# ------------------------------------------------------------
# 4. Validate critical locked objects
# ------------------------------------------------------------

locked_objects <- c(
  "module_tme_results",
  "GS_all",
  "MM_all",
  "moduleColors",
  "ME_color_formal_raw",
  "datExpr_MAD",
  "estimate_scores_wgcna",
  "module_genes",
  "core_visualization_table",
  "core_position_final",
  "core_descriptive_profile",
  "candidate_master_18D"
)

locked_status <- data.frame(
  Object = locked_objects,
  Exists = vapply(
    locked_objects,
    exists,
    logical(1),
    inherits = FALSE
  ),
  stringsAsFactors = FALSE
)

cat("\nLocked object recovery status:\n")
print(
  locked_status,
  row.names = FALSE
)

stopifnot(
  all(locked_status$Exists)
)


# ------------------------------------------------------------
# 5. Recover SF4C locked statistics
# ------------------------------------------------------------

cat("\nRecovering SF4C locked statistics...\n")

sf4c_stats <- readRDS(
  sf4c_rds
)

cat(
  "SF4C object class:",
  paste(class(sf4c_stats), collapse = ";"),
  "\n"
)


# ------------------------------------------------------------
# 6. Recover Supplementary Table X
# ------------------------------------------------------------

cat("\nRecovering Supplementary Table X...\n")

supp_table_x_display <- readRDS(
  supp_table_rds
)

cat(
  "Supplementary Table X dimensions:",
  nrow(supp_table_x_display),
  "x",
  ncol(supp_table_x_display),
  "\n"
)


# ------------------------------------------------------------
# 7. Restore the underlying Supplementary Table X object
#    from the locked GS/module data
#
#    This is object reconstruction only.
#    No statistical calculation is performed.
# ------------------------------------------------------------

core_gene_map <- data.frame(
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  ENSEMBL = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  stringsAsFactors = FALSE
)


core_tme_table <- GS_all %>%
  filter(
    Gene %in% core_gene_map$ENSEMBL
  ) %>%
  select(
    Gene,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  ) %>%
  left_join(
    core_gene_map,
    by = c("Gene" = "ENSEMBL")
  ) %>%
  select(
    Gene,
    ENSEMBL,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  )


core_position_table <- core_position_final %>%
  filter(
    Symbol %in% c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  ) %>%
  transmute(
    Gene = Symbol,
    Module,
    Module_Size,
    MM_Rank = MM_rank,
    MM_Percentile = MM_percentile
  )


supp_table_x <- core_position_table %>%
  left_join(
    core_tme_table,
    by = "Gene"
  )


supp_table_x_final <- supp_table_x %>%
  select(
    Gene,
    ENSEMBL,
    Module,
    Module_Size,
    MM_Rank,
    MM_Percentile,
    GS_Immune,
    P_Immune,
    FDR_Immune,
    GS_Stromal,
    P_Stromal,
    FDR_Stromal,
    GS_ESTIMATE,
    P_ESTIMATE,
    FDR_ESTIMATE
  )


# ------------------------------------------------------------
# 8. Restore SF4B plotting data
# ------------------------------------------------------------

sf4b_plot_df <- core_position_final %>%
  filter(
    Symbol %in% c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  ) %>%
  transmute(
    Gene = Symbol,
    Module,
    ME,
    Module_Size,
    MM,
    MM_Rank = MM_rank,
    MM_Percentile = MM_percentile
  )


# ------------------------------------------------------------
# 9. Restore SF4D source data
# ------------------------------------------------------------

sf4d_df <- GS_all %>%
  filter(
    Gene %in% core_gene_map$ENSEMBL
  ) %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  ) %>%
  left_join(
    core_gene_map,
    by = c("Gene" = "ENSEMBL")
  ) %>%
  select(
    Gene,
    ENSEMBL,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  )


# ------------------------------------------------------------
# 10. Restore SF4D long-format plotting data
# ------------------------------------------------------------

sf4d_long <- sf4d_df %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      GS_Immune,
      GS_Stromal,
      GS_ESTIMATE
    ),
    names_to = "Trait",
    values_to = "GS"
  ) %>%
  mutate(
    Trait = dplyr::recode(
      Trait,
      GS_Immune = "ImmuneScore",
      GS_Stromal = "StromalScore",
      GS_ESTIMATE = "ESTIMATEScore"
    ),
    Abs_GS = abs(GS)
  )


# ------------------------------------------------------------
# 11. Validate recovered SF4D values
# ------------------------------------------------------------

expected_sf4d <- data.frame(
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  GS_Immune = c(
    0.3628795,
    0.6866628,
    0.7688851
  ),
  GS_Stromal = c(
    0.6555792,
    0.6460729,
    0.7067829
  ),
  GS_ESTIMATE = c(
    0.5399233,
    0.7230010,
    0.8045778
  )
)

sf4d_audit <- sf4d_df %>%
  select(
    Gene,
    GS_Immune,
    GS_Stromal,
    GS_ESTIMATE
  ) %>%
  left_join(
    expected_sf4d,
    by = "Gene",
    suffix = c("_actual", "_expected")
  )

stopifnot(
  all(
    abs(
      sf4d_audit$GS_Immune_actual -
        sf4d_audit$GS_Immune_expected
    ) < 1e-6
  )
)

stopifnot(
  all(
    abs(
      sf4d_audit$GS_Stromal_actual -
        sf4d_audit$GS_Stromal_expected
    ) < 1e-6
  )
)

stopifnot(
  all(
    abs(
      sf4d_audit$GS_ESTIMATE_actual -
        sf4d_audit$GS_ESTIMATE_expected
    ) < 1e-6
  )
)


# ------------------------------------------------------------
# 12. Final Supplementary Table X validation
# ------------------------------------------------------------

stopifnot(
  nrow(supp_table_x_display) == 3
)

stopifnot(
  identical(
    supp_table_x_display$Gene,
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )
)


# ------------------------------------------------------------
# 13. Define final checkpoint
# ------------------------------------------------------------

checkpoint_dir <- file.path(
  prep_dir,
  "R_workspace"
)

checkpoint_file <- file.path(
  checkpoint_dir,
  "Paper1_STEP19D3_END_WORKSPACE_20260914.RData"
)


# ------------------------------------------------------------
# 14. Save complete recovered workspace
# ------------------------------------------------------------

cat("\nSaving final recovered workspace...\n")

save.image(
  file = checkpoint_file
)


# ------------------------------------------------------------
# 15. Verify checkpoint
# ------------------------------------------------------------

stopifnot(
  file.exists(checkpoint_file)
)

stopifnot(
  file.size(checkpoint_file) > 0
)


cat("\nCheckpoint file:\n")
cat(
  normalizePath(
    checkpoint_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "Size:",
  round(
    file.size(checkpoint_file) /
      (1024^2),
    2
  ),
  "MB\n"
)


# ------------------------------------------------------------
# 16. Final checkpoint status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 19D-3 STATUS: PASS / LOCKED\n")
cat("============================================================\n")

cat(
  "PRE-STEP19 locked workspace recovered.\n"
)

cat(
  "SF4C locked statistics recovered from archived RDS.\n"
)

cat(
  "Supplementary Table X recovered from archived RDS.\n"
)

cat(
  "SF4B/SF4D source objects restored from locked data.\n"
)

cat(
  "No WGCNA analysis was rerun.\n"
)

cat(
  "No statistical result was recalculated.\n"
)

cat(
  "Final workspace checkpoint saved successfully.\n"
)

cat(
  "============================================================\n"
)

