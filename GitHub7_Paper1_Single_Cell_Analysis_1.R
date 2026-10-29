##GitHub7_Paper1_Single_Cell_Analysis_1

###Redraw Figure 9A–9B (IL7R, AOC3, MS4A7)
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

setwd("input_dir/Figure 9A_9B_umap")
getwd()

###############################################################
##
## Figure 9A
## NSCLC_GSE131907
## Celltype UMAP
##
## Final Version V3
##
###############################################################

options(stringsAsFactors = FALSE)

set.seed(20250724)

library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(Matrix)

theme_set(theme_bw())

###############################################################
## Output directory
###############################################################

dir.create("Figure9",showWarnings = FALSE)

###############################################################
## Read expression matrix
###############################################################

expr <- Read10X_h5(
  "input_dir/Figure 9C Dot Plot (Publication-ready)/NSCLC_GSE131907_expression.h5"
)

dim(expr)

head(colnames(expr))


cat("Genes :",nrow(expr),"\n")

cat("Cells :",ncol(expr),"\n")

###############################################################
## Create Seurat Object
###############################################################

sce <- CreateSeuratObject(
  
  counts = expr,
  
  project = "NSCLC_GSE131907",
  
  min.cells = 0,
  
  min.features = 0
  
)

sce

###############################################################
## Read metadata
###############################################################

meta <- fread(
  
  "input_dir/Figure 9C Dot Plot (Publication-ready)/NSCLC_GSE131907_CellMetainfo_table.tsv",
  
  data.table = FALSE
  
)

dim(meta)

head(meta)

rownames(meta) <- meta$Cell

###############################################################
## Match barcodes
###############################################################

all(
  colnames(sce) %in%
    rownames(meta)
)
#then：
meta <- meta[
  colnames(sce),
]
#check again：
sce@meta.data <- meta

###############################################################
## Add UMAP coordinates
###############################################################

umap <- as.matrix(
  
  meta[
    ,
    c(
      "UMAP_1",
      "UMAP_2"
    )
  ]
  
)

rownames(umap) <- rownames(meta)

sce[["umap"]] <-
  
  CreateDimReducObject(
    
    embeddings = umap,
    
    key = "UMAP_",
    
    assay = "RNA"
    
  )

#then check：
Embeddings(
  sce,
  "umap"
)[1:5,]

###############################################################
## Step 7
## Set Celltype
## Set Color
## Optimize plotting order
###############################################################

##------------------------------------------------------------
## 7.1  Set Celltype (fixed order)
##------------------------------------------------------------

cell.order <- c(
  "B",
  "CD4Tconv",
  "CD8T",
  "CD8Tex",
  "DC",
  "Endothelial",
  "Epithelial",
  "Fibroblasts",
  "Mast",
  "Mono/Macro",
  "Oligodendrocyte",
  "Plasma"
)

sce$Celltype <- factor(
  sce$`Celltype (major-lineage)`,
  levels = cell.order
)

meta$Celltype <- factor(
  meta$`Celltype (major-lineage)`,
  levels = cell.order
)

Idents(sce) <- "Celltype"

##------------------------------------------------------------
## 7.2  Check Celltype
##------------------------------------------------------------

table(Idents(sce))

levels(Idents(sce))

##------------------------------------------------------------
## 7.3  Publication colour palette
## (close to TISCH2)
##------------------------------------------------------------
cell.cols <- c(
  
  "B"               ="#E68613",
  
  "CD4Tconv"        ="#5968C5",
  
  "CD8T"            ="#43B8A6",
  
  "CD8Tex"          ="#97CC3C",
  
  "DC"              ="#C65CB3",
  
  "Endothelial"     ="#2D8A74",
  
  "Epithelial"      ="#E2AE14",
  
  "Fibroblasts"     ="#3D8DCB",
  
  "Mast"            ="#7A55AE",
  
  "Mono/Macro"      ="#F2665A",
  
  "Oligodendrocyte" ="#A5AAA2",
  
  "Plasma"          ="#C2C431"
  
)

##------------------------------------------------------------
## 7.4  Check colour names
##------------------------------------------------------------

stopifnot(
  identical(
    names(cell.cols),
    levels(Idents(sce))
  )
)

##------------------------------------------------------------
## 7.5  Optimize drawing order
##
## Small clusters are drawn LAST,
## preventing them from being hidden
## by large clusters.
##------------------------------------------------------------

cluster.size <- table(meta$Celltype)

cluster.size

##------------------------------------------------------------
## Order by cluster size
##------------------------------------------------------------

draw.order <- names(
  sort(
    cluster.size,
    decreasing = TRUE
  )
)

draw.order

##------------------------------------------------------------
## Reorder cells
## Large clusters first
## Small clusters last
##------------------------------------------------------------

meta$DrawOrder <- factor(
  meta$Celltype,
  levels = draw.order
)

meta <- meta[
  order(meta$DrawOrder),
]

## Keep Seurat metadata consistent
sce <- subset(
  sce,
  cells = rownames(meta)
)

###############################################################
## Check
###############################################################

head(meta[, c(
  "Celltype",
  "DrawOrder"
)])

table(meta$Celltype)


###############################################################
#Step 8.1 Font initialization

###############################################################
##
## Figure 9A Final V4 (Continuation Version)
##
## Step 8.1
## Font initialization
##
###############################################################

library(ggplot2)
library(dplyr)
library(ggrepel)

## Windows Arial
windowsFonts(
  Arial = windowsFont("Arial")
)

## Publication font
base_family <- "sans"

message("✓ Windows font initialized.")


#Figure 9A Final V4
#Step 8.2（Final）
###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.2
## Extract UMAP & Prepare metadata
##
###############################################################

##------------------------------------------------------------
## Check Seurat object
##------------------------------------------------------------

if (!exists("sce")) {
  
  stop(
    "Seurat object 'sce' was not found.\nPlease finish Step 1–7 first."
  )
  
}

##------------------------------------------------------------
## Check metadata object
##------------------------------------------------------------

if (!exists("meta")) {
  
  stop(
    "Metadata object 'meta' was not found."
  )
  
}

if (nrow(meta) == 0) {
  
  stop(
    "Metadata is empty."
  )
  
}

##------------------------------------------------------------
## Check UMAP reduction
##------------------------------------------------------------

if (!("umap" %in% names(sce@reductions))) {
  
  stop(
    "UMAP reduction was not found."
  )
  
}

##------------------------------------------------------------
## Extract UMAP coordinates
##------------------------------------------------------------

umap_df <-
  
  as.data.frame(
    
    Embeddings(
      sce,
      reduction = "umap"
    )
    
  )

##------------------------------------------------------------
## Add Cell ID
##------------------------------------------------------------

umap_df$Cell <- rownames(umap_df)

##------------------------------------------------------------
## Check Cell names
##------------------------------------------------------------

missing_cells <-
  
  setdiff(
    
    rownames(umap_df),
    
    rownames(meta)
    
  )

if(length(missing_cells) > 0){
  
  stop(
    
    paste(
      
      length(missing_cells),
      
      "cells are missing in metadata."
      
    )
    
  )
  
}

##------------------------------------------------------------
## Merge Celltype
##------------------------------------------------------------

if (!("Celltype" %in% colnames(meta))) {
  
  stop(
    "Column 'Celltype' was not found in metadata."
  )
  
}

umap_df$Celltype <-
  
  meta[
    rownames(umap_df),
    "Celltype"
  ]

##------------------------------------------------------------
## Merge Sample (optional)
##------------------------------------------------------------

if ("Sample" %in% colnames(meta)) {
  
  umap_df$Sample <-
    
    meta[
      rownames(umap_df),
      "Sample"
    ]
  
}

##------------------------------------------------------------
## Merge Patient (optional)
##------------------------------------------------------------

if ("Patient" %in% colnames(meta)) {
  
  umap_df$Patient <-
    
    meta[
      rownames(umap_df),
      "Patient"
    ]
  
}

##------------------------------------------------------------
## Remove incomplete observations
##------------------------------------------------------------

n_before <- nrow(umap_df)

umap_df <-
  
  umap_df %>%
  
  filter(
    
    !is.na(UMAP_1),
    
    !is.na(UMAP_2),
    
    !is.na(Celltype)
    
  )

n_after <- nrow(umap_df)

message(
  
  "Cells removed: ",
  
  n_before - n_after
  
)

##------------------------------------------------------------
## QC 1
## Check duplicated Cell IDs
##------------------------------------------------------------

if(anyDuplicated(umap_df$Cell) > 0){
  
  stop(
    
    "Duplicated Cell IDs detected."
    
  )
  
}

##------------------------------------------------------------
## QC 2
## Check Celltype
##------------------------------------------------------------

cat("\n")

cat("----------------------------------------\n")

cat("Celltype summary\n")

cat("----------------------------------------\n")

print(
  
  table(
    
    umap_df$Celltype
    
  )
  
)

cat("----------------------------------------\n")

##------------------------------------------------------------
## QC 3
## Check dimensions
##------------------------------------------------------------

cat(
  
  "Number of cells :",
  
  nrow(umap_df),
  
  "\n"
  
)

cat(
  
  "Number of Celltypes :",
  
  length(
    
    unique(
      
      umap_df$Celltype
      
    )
    
  ),
  
  "\n"
  
)

##------------------------------------------------------------
## Preview
##------------------------------------------------------------

cat("\n")

cat("Preview of plotting data:\n")

head(umap_df)

#Step 8.3（Final V4）
###############################################################
##
## Step 8.3
## Publication cell order
##
###############################################################
celltype_levels <-c(
  "B",
  "CD4Tconv",
  "CD8T",
  "CD8Tex",
  "DC",
  "Endothelial",
  "Epithelial",
  "Fibroblasts",
  "Mast",
  "Mono/Macro",
  "Oligodendrocyte",
  "Plasma"
)

###############################################################
## QC
###############################################################

missing_levels <-
  
  setdiff(
    
    unique(umap_df$Celltype),
    
    celltype_levels
    
  )

if(length(missing_levels) > 0){
  
  stop(
    
    paste(
      
      "Undefined Celltype:",
      
      paste(missing_levels, collapse=", ")
      
    )
    
  )
  
}

###############################################################
## Fixed factor levels
###############################################################

umap_df$Celltype <-
  
  factor(
    
    umap_df$Celltype,
    
    levels = celltype_levels
    
  )

###############################################################
## Automatic DrawOrder
###############################################################

umap_df$DrawOrder <-
  
  as.numeric(
    
    umap_df$Celltype
    
  )

###############################################################
## Sort
###############################################################

umap_df <-
  
  umap_df %>%
  
  arrange(
    
    DrawOrder
    
  )

###############################################################
## QC
###############################################################

cat("\n")

cat("--------------------------------\n")

cat("Publication drawing order\n")

cat("--------------------------------\n")

print(
  
  levels(
    
    umap_df$Celltype
    
  )
  
)

cat("--------------------------------\n")

#Figure 9A Final V4
#Step 8.4 Fixed color palette + Color QC
###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.4
## Fixed color palette + Color QC
##
###############################################################

##------------------------------------------------------------
## Fixed publication color palette
##
## Color-blind friendly
## Adobe Illustrator compatible
## Consistent across all figures
##------------------------------------------------------------
celltype_colors <- c(
  
  "B"               ="#E68613",
  
  "CD4Tconv"        ="#5968C5",
  
  "CD8T"            ="#43B8A6",
  
  "CD8Tex"          ="#97CC3C",
  
  "DC"              ="#C65CB3",
  
  "Endothelial"     ="#2D8A74",
  
  "Epithelial"      ="#E2AE14",
  
  "Fibroblasts"     ="#3D8DCB",
  
  "Mast"            ="#7A55AE",
  
  "Mono/Macro"      ="#F2665A",
  
  "Oligodendrocyte" ="#A5AAA2",
  
  "Plasma"          ="#C2C431"
  
)

##------------------------------------------------------------
## QC 1
## Duplicate color names
##------------------------------------------------------------

if(anyDuplicated(names(celltype_colors)) > 0){
  
  stop("Duplicate Celltype names detected in celltype_colors.")
  
}

##------------------------------------------------------------
## QC 2
## Missing Celltype
##------------------------------------------------------------

missing_color <-
  
  setdiff(
    
    levels(umap_df$Celltype),
    
    names(celltype_colors)
    
  )

if(length(missing_color) > 0){
  
  stop(
    
    paste(
      
      "Missing colors for:",
      
      paste(missing_color, collapse = ", ")
      
    )
    
  )
  
}

##------------------------------------------------------------
## QC 3
## Extra Celltype
##------------------------------------------------------------

unused_color <-
  
  setdiff(
    
    names(celltype_colors),
    
    levels(umap_df$Celltype)
    
  )

if(length(unused_color) > 0){
  
  warning(
    
    paste(
      
      "Unused colors:",
      
      paste(unused_color, collapse = ", ")
      
    )
    
  )
  
}

##------------------------------------------------------------
## Reorder colors
##------------------------------------------------------------

celltype_colors <-
  
  celltype_colors[
    
    levels(
      
      umap_df$Celltype
      
    )
    
  ]

##------------------------------------------------------------
## QC 4
## Hex color validation
##------------------------------------------------------------

valid_hex <-
  
  grepl(
    
    "^#[A-Fa-f0-9]{6}$",
    
    celltype_colors
    
  )

if(!all(valid_hex)){
  
  stop("Invalid HEX color detected.")
  
}

##------------------------------------------------------------
## QC 5
## Color uniqueness
##------------------------------------------------------------

if(anyDuplicated(celltype_colors) > 0){
  
  warning(
    
    "Some Celltypes share identical colors."
    
  )
  
}

##------------------------------------------------------------
## Print color table
##------------------------------------------------------------

cat("\n")

cat("----------------------------------------\n")

cat("Publication color palette\n")

cat("----------------------------------------\n")

print(
  
  data.frame(
    
    Celltype = names(celltype_colors),
    
    Color = unname(celltype_colors),
    
    row.names = NULL
    
  )
  
)

cat("----------------------------------------\n")

##------------------------------------------------------------
## Preview palette
##------------------------------------------------------------

palette_df <-
  
  data.frame(
    
    Celltype = factor(
      
      names(celltype_colors),
      
      levels = names(celltype_colors)
      
    ),
    
    x = seq_along(celltype_colors),
    
    y = 1
    
  )

palette_plot <-
  
  ggplot(
    
    palette_df,
    
    aes(
      
      x,
      
      y,
      
      fill = Celltype
      
    )
    
  ) +
  
  geom_tile(
    
    height = 0.8
    
  ) +
  
  scale_fill_manual(
    
    values = celltype_colors
    
  ) +
  
  theme_void(
    
    base_family = base_family
    
  ) +
  
  theme(
    
    legend.position = "none",
    
    axis.text = element_blank(),
    
    axis.ticks = element_blank()
    
  )

print(palette_plot)

rm(palette_df)

#Step 8.5 Calculate label coordinates
###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.5
## Calculate label coordinates
##
###############################################################

library(dplyr)

##------------------------------------------------------------
## Calculate label positions
##------------------------------------------------------------

label_df <-
  
  umap_df %>%
  
  group_by(Celltype) %>%
  
  summarise(
    
    UMAP_1 = median(UMAP_1),
    
    UMAP_2 = median(UMAP_2),
    
    nCells = n(),
    
    xmin = min(UMAP_1),
    
    xmax = max(UMAP_1),
    
    ymin = min(UMAP_2),
    
    ymax = max(UMAP_2),
    
    .groups = "drop"
    
  )

##------------------------------------------------------------
## Keep Celltype order
##------------------------------------------------------------

label_df$Celltype <-
  
  factor(
    
    label_df$Celltype,
    
    levels = levels(umap_df$Celltype)
    
  )

label_df <-
  
  label_df %>%
  
  arrange(Celltype)

##------------------------------------------------------------
## QC 1
## Missing labels
##------------------------------------------------------------

missing_label <-
  
  setdiff(
    
    levels(umap_df$Celltype),
    
    label_df$Celltype
    
  )

if(length(missing_label) > 0){
  
  stop(
    
    paste(
      
      "Missing label coordinates for:",
      
      paste(missing_label, collapse = ", ")
      
    )
    
  )
  
}

##------------------------------------------------------------
## QC 2
## Median inside UMAP range
##------------------------------------------------------------

inside_x <-
  
  label_df$UMAP_1 >= label_df$xmin &
  label_df$UMAP_1 <= label_df$xmax

inside_y <-
  
  label_df$UMAP_2 >= label_df$ymin &
  label_df$UMAP_2 <= label_df$ymax

if(!all(inside_x & inside_y)){
  
  warning(
    
    "Some label positions are outside the cluster range."
    
  )
  
}

##------------------------------------------------------------
## QC 3
## Cell numbers
##------------------------------------------------------------

if(any(label_df$nCells == 0)){
  
  stop(
    
    "Empty Celltype detected."
    
  )
  
}

##------------------------------------------------------------
## Print summary
##------------------------------------------------------------

cat("\n")

cat("----------------------------------------\n")

cat("Label coordinate summary\n")

cat("----------------------------------------\n")

print(
  
  label_df[, c(
    
    "Celltype",
    
    "nCells",
    
    "UMAP_1",
    
    "UMAP_2"
    
  )]
  
)

cat("----------------------------------------\n")

##------------------------------------------------------------
## Preview
##------------------------------------------------------------

head(label_df)

###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.6
## Publication Theme
##
###############################################################

##------------------------------------------------------------
## Base publication theme
##------------------------------------------------------------

publication_theme <-
  
  theme_classic(
    
    base_family = base_family
    
  ) +
  
  theme(
    
    #######################################################
    ## Text
    #######################################################
    
    text = element_text(
      
      family = base_family,
      
      colour = "black"
      
    ),
    
    #######################################################
    ## Axis title
    #######################################################
    
    axis.title = element_blank(),
    
    #######################################################
    ## Axis text
    #######################################################
    
    axis.text.x = element_text(
      
      size = 10,
      
      colour = "black"
      
    ),
    
    axis.text.y = element_text(
      
      size = 10,
      
      colour = "black"
      
    ),
    
    #######################################################
    ## Axis ticks
    #######################################################
    
    axis.ticks = element_line(
      
      colour = "black",
      
      linewidth = 0.4
      
    ),
    
    axis.ticks.length =
      
      unit(
        
        1,
        
        "mm"
        
      ),
    
    #######################################################
    ## Axis line
    #######################################################
    
    axis.line = element_line(
      
      colour = "black",
      
      linewidth = 0.35
      
    ),
    
    #######################################################
    ## Panel
    #######################################################
    
    panel.background =
      
      element_blank(),
    
    panel.border =
      
      element_blank(),
    
    panel.grid.major =
      
      element_blank(),
    
    panel.grid.minor =
      
      element_blank(),
    
    #######################################################
    ## Legend
    #######################################################
    
    legend.position = "right",
    
    legend.direction = "vertical",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      
      size = 8,
      
      family = base_family,
      
      colour = "black"
      
    ),
    
    legend.key =
      
      element_blank(),
    
    legend.key.size =
      
      unit(
        
        2.6,
        
        "mm"
        
      ),
    
    legend.spacing.y =
      
      unit(
        
        1,
        
        "mm"
        
      ),
    
    #######################################################
    ## Plot
    #######################################################
    
    plot.title =
      
      element_blank(),
    
    plot.subtitle =
      
      element_blank(),
    
    plot.caption =
      
      element_blank(),
    
    #######################################################
    ## Margin
    #######################################################
    
    plot.margin =
      
      margin(
        
        3,
        
        3,
        
        3,
        
        3,
        
        unit = "mm"
        
      )
    
  )

###############################################################
## QC
###############################################################

message("✓ Publication theme created.")

###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.7
## Draw UMAP
##
###############################################################

library(ggplot2)
library(ggrepel)

##------------------------------------------------------------
## Build Figure
##------------------------------------------------------------

p <-
  
  ggplot(
    
    umap_df,
    
    aes(
      
      x = UMAP_1,
      
      y = UMAP_2,
      
      colour = Celltype
      
    )
    
  ) +
  
  ###############################################################
## UMAP points
###############################################################
geom_point(
  size = 0.22,
  shape = 16,
  stroke = 0,
  alpha = 1
)+
  
  ###############################################################
## Celltype labels
###############################################################

geom_text_repel(
  
  data = label_df,
  
  aes(
    
    x = UMAP_1,
    
    y = UMAP_2,
    
    label = Celltype
    
  ),
  
  family = base_family,
  
  fontface = "plain",
  
  colour = "black",
  
  size = 2.8,
  
  seed = 1234,
  
  box.padding = 0.35,
  
  point.padding = 0.15,
  
  force = 1.2,
  
  max.overlaps = Inf,
  
  min.segment.length = 0,
  
  segment.color = NA,
  
  segment.size = 0,
  
  show.legend = FALSE
  
) +
  
  ###############################################################
## Fixed colours
###############################################################
scale_color_manual(
  values = celltype_colors,
  drop = FALSE
) +
  
  ###############################################################
## Keep aspect ratio
###############################################################

coord_equal() +
  
  ###############################################################
## Publication theme
###############################################################

publication_theme +
  
  ###############################################################
## Legend
###############################################################

guides(
  
  colour = guide_legend(
    
    override.aes = list(
      
      size = 2
      
    )
    
  )
  
)

###############################################################
## QC
###############################################################

print(p)

###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.8
## Export publication figures
##
###############################################################

##------------------------------------------------------------
## Output directory
##------------------------------------------------------------

outdir <- "figures"

if (!dir.exists(outdir)) {
  
  dir.create(
    outdir,
    recursive = TRUE
  )
  
}

##------------------------------------------------------------
## Figure name
##------------------------------------------------------------

figure_name <- "Figure9A_Final_V4_3"

##------------------------------------------------------------
## Publication size
##
## Single-column figure
##------------------------------------------------------------

fig_width  <- 120      # mm
fig_height <-90      # mm

##------------------------------------------------------------
## Export PDF (vector)
##------------------------------------------------------------

ggsave(
  
  filename = file.path(
    outdir,
    paste0(figure_name, ".pdf")
  ),
  
  plot = p,
  
  device = cairo_pdf,
  
  width = fig_width,
  
  height = fig_height,
  
  units = "mm",
  
  bg = "white"
  
)

message("✓ PDF exported.")

##------------------------------------------------------------
## Export TIFF
##------------------------------------------------------------

ggsave(
  
  filename = file.path(
    outdir,
    paste0(figure_name, ".tiff")
  ),
  
  plot = p,
  
  width = fig_width,
  
  height = fig_height,
  
  units = "mm",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white"
  
)

message("✓ TIFF exported.")

##------------------------------------------------------------
## Export PNG
##------------------------------------------------------------

ggsave(
  
  filename = file.path(
    outdir,
    paste0(figure_name, ".png")
  ),
  
  plot = p,
  
  width = fig_width,
  
  height = fig_height,
  
  units = "mm",
  
  dpi = 600,
  
  bg = "white"
  
)

message("✓ PNG exported.")

##------------------------------------------------------------
## Export SVG
##------------------------------------------------------------

ggsave(
  
  filename = file.path(
    outdir,
    paste0(figure_name, ".svg")
  ),
  
  plot = p,
  
  width = fig_width,
  
  height = fig_height,
  
  units = "mm",
  
  bg = "white"
  
)

message("✓ SVG exported.")

###############################################################
## Export summary
###############################################################

cat("\n")

cat("----------------------------------------\n")

cat("Figure exported successfully\n")

cat("----------------------------------------\n")

cat("Output directory :", normalizePath(outdir), "\n")

cat("Figure size      :", fig_width, "×", fig_height, "mm\n")

cat("Formats          : PDF / TIFF / PNG / SVG\n")

cat("----------------------------------------\n")

###############################################################
##
## Figure 9A Final V4 (CSBJ Submission Version)
##
## Step 8.9
## Submission Quality Control
##
###############################################################

cat("\n")
cat("========================================\n")
cat(" Figure 9A Submission QC\n")
cat("========================================\n")

###############################################################
## QC 1
## Check plotting object
###############################################################

if (!exists("p")) {
  
  stop("Plot object 'p' was not found.")
  
}

cat("✓ Plot object exists.\n")

###############################################################
## QC 2
## Check Celltype
###############################################################
expected_celltypes <- c(
  "B",
  "CD4Tconv",
  "CD8T",
  "CD8Tex",
  "DC",
  "Endothelial",
  "Epithelial",
  "Fibroblasts",
  "Mast",
  "Mono/Macro",
  "Oligodendrocyte",
  "Plasma"
)

if (!identical(
  
  levels(umap_df$Celltype),
  
  expected_celltypes
  
)) {
  
  warning(
    
    "Celltype order differs from publication version."
    
  )
  
} else {
  
  cat("✓ Celltype order verified.\n")
  
}

###############################################################
## QC 3
## Check colour mapping
###############################################################

if (!all(
  
  expected_celltypes %in%
  
  names(celltype_colors)
  
)) {
  
  stop(
    
    "Color mapping is incomplete."
    
  )
  
}

cat("✓ Color mapping verified.\n")

###############################################################
## QC 4
## Check labels
###############################################################

if (
  
  nrow(label_df) != length(expected_celltypes)
  
) {
  
  warning(
    
    "Unexpected number of labels."
    
  )
  
} else {
  
  cat("✓ Label number verified.\n")
  
}

###############################################################
## QC 5
## Check exported files
###############################################################

files <- c(
  
  "Figure9A_Final_V4.pdf",
  
  "Figure9A_Final_V4.tiff",
  
  "Figure9A_Final_V4.png"
  
)

for(f in files){
  
  fullfile <-
    
    file.path(
      
      outdir,
      
      f
      
    )
  
  if(file.exists(fullfile)){
    
    cat(
      
      "✓",
      
      f,
      
      "\n"
      
    )
    
  }else{
    
    warning(
      
      paste(
        
        f,
        
        "was not exported."
        
      )
      
    )
    
  }
  
}

###############################################################
## QC 6
## File size
###############################################################

for(f in files){
  
  fullfile <-
    
    file.path(
      
      outdir,
      
      f
      
    )
  
  if(file.exists(fullfile)){
    
    info <- file.info(fullfile)
    
    cat(
      
      sprintf(
        
        "%-25s %.2f MB\n",
        
        f,
        
        info$size/1024^2
        
      )
      
    )
    
  }
  
}

###############################################################
## QC 7
## Session information
###############################################################

cat("\n")

cat("----------------------------------------\n")

cat("R version\n")

cat("----------------------------------------\n")

cat(
  
  R.version.string,
  
  "\n"
  
)

cat(
  
  "ggplot2:",
  
  as.character(
    
    packageVersion(
      
      "ggplot2"
      
    )
    
  ),
  
  "\n"
  
)

cat(
  
  "Seurat :", 
  
  as.character(
    
    packageVersion(
      
      "Seurat"
      
    )
    
  ),
  
  "\n"
  
)

cat(
  
  "showtext:",
  
  as.character(
    
    packageVersion(
      
      "showtext"
      
    )
    
  ),
  
  "\n"
  
)

cat("----------------------------------------\n")

###############################################################
## Final summary
###############################################################

cat("\n")

cat("========================================\n")

cat(" Figure 9A QC completed.\n")
cat(" Ready for Illustrator / Submission.\n")
cat("========================================\n")


# Figure 9B
## Final V4 (CSBJ Submission Version)
##
## Part 1
## build_plot_df()
##
###############################################################

build_plot_df <- function(
    sce,
    genes = c("AOC3", "IL7R", "MS4A7")
){
  
  ###########################################################
  ## Check Seurat object
  ###########################################################
  
  if (!inherits(sce, "Seurat")) {
    
    stop("Input must be a Seurat object.")
    
  }
  
  ###########################################################
  ## Check RNA assay
  ###########################################################
  
  if (!("RNA" %in% Assays(sce))) {
    
    stop("RNA assay was not found.")
    
  }
  
  DefaultAssay(sce) <- "RNA"
  
  ###########################################################
  ## Check UMAP
  ###########################################################
  
  if (!("umap" %in% names(sce@reductions))) {
    
    stop("UMAP reduction was not found.")
    
  }
  
  ###########################################################
  ## Check counts layer
  ###########################################################
  
  if (!("counts" %in% Layers(sce[["RNA"]]))) {
    
    stop("Counts layer was not found.")
    
  }
  
  ###########################################################
  ## Check genes
  ###########################################################
  
  gene_exist <- genes %in% rownames(sce)
  
  if (any(!gene_exist)) {
    
    stop(
      paste(
        "Gene(s) not found:",
        paste(genes[!gene_exist], collapse = ", ")
      )
    )
    
  }
  
  ###########################################################
  ## Read counts
  ###########################################################
  
  message("Reading counts layer ...")
  
  counts <- LayerData(
    
    object = sce,
    
    assay = "RNA",
    
    layer = "counts"
    
  )
  
  ###########################################################
  ## Library size
  ###########################################################
  
  message("Calculating library size ...")
  
  lib.size <- Matrix::colSums(counts)
  
  ###########################################################
  ## Extract genes only
  ###########################################################
  
  message("Extracting target genes ...")
  
  gene_counts <- counts[
    genes,
    ,
    drop = FALSE
  ]
  
  ###########################################################
  ## LogNormalize
  ##
  ## log1p(count / library *10000)
  ###########################################################
  
  message("Log-normalizing target genes ...")
  
  gene_expr <-
    
    t(
      
      t(gene_counts) /
        
        lib.size *
        
        10000
      
    )
  
  gene_expr <- log1p(gene_expr)
  
  gene_expr <- as.matrix(gene_expr)
  
  ###########################################################
  ## UMAP
  ###########################################################
  
  message("Extracting UMAP ...")
  
  umap_df <-
    
    as.data.frame(
      
      Embeddings(
        
        sce,
        
        reduction = "umap"
        
      )
      
    )
  
  ###########################################################
  ## Metadata
  ###########################################################
  
  meta <- sce@meta.data
  
  ###########################################################
  ## Build dataframe
  ###########################################################
  
  plot_df <- data.frame(
    
    Cell = rownames(meta),
    
    UMAP_1 = umap_df[,1],
    
    UMAP_2 = umap_df[,2],
    
    Celltype = meta$Celltype,
    
    Sample = meta$Sample,
    
    Patient = meta$Patient,
    
    stringsAsFactors = FALSE
    
  )
  
  ###########################################################
  ## Add gene expression
  ###########################################################
  
  for(g in genes){
    
    plot_df[[g]] <-
      
      as.numeric(
        
        gene_expr[g, ]
        
      )
    
  }
  
  ###########################################################
  ## Draw order
  ##
  ## low expression first
  ###########################################################
  
  plot_df$DrawOrder <- seq_len(
    
    nrow(plot_df)
    
  )
  
  ###########################################################
  ## Summary
  ###########################################################
  
  message("----------------------------------------")
  message("Figure 9B plot_df created successfully.")
  message("----------------------------------------")
  message("Cells : ", nrow(plot_df))
  message("Genes : ", paste(genes, collapse = ", "))
  message("----------------------------------------")
  
  return(plot_df)
  
}

###############################################################
##
## Figure 9B
## Final V5 (CSBJ Submission Version)
##
## Step 2.1
## Function framework + data check + adaptive cutoff
##
###############################################################

plot_feature_umap_v5 <- function(
    gene,
    cutoff.quantile = 0.05,
    clip.quantile   = 0.95
){
  
  ###########################################################
  ## Check plot_df
  ###########################################################
  
  if(!exists("plot_df", inherits = TRUE)){
    
    stop(
      "plot_df was not found.\n",
      "Please run:\n",
      "plot_df <- build_plot_df(sce)"
    )
    
  }
  
  ###########################################################
  ## Copy dataframe
  ###########################################################
  
  df <- plot_df
  
  ###########################################################
  ## Check gene
  ###########################################################
  
  if(!(gene %in% colnames(df))){
    
    stop(
      paste0(
        "Gene '",
        gene,
        "' was not found in plot_df."
      )
    )
    
  }
  
  ###########################################################
  ## Expression
  ###########################################################
  
  df$Expression <- df[[gene]]
  
  ###########################################################
  ## Positive cells
  ###########################################################
  
  positive <- df[df$Expression > 0.10, ]
  
  if(nrow(positive) == 0){
    
    stop(
      paste0(
        "No positive cells detected for ",
        gene,
        "."
      )
    )
    
  }
  
  ###########################################################
  ## Adaptive cutoff
  ###########################################################
  if(gene=="AOC3"){
    cutoff.quantile <-0.55
  }
  
  if(gene=="IL7R"){
    cutoff.quantile <-0.35
  }
  
  if(gene=="MS4A7"){
    cutoff.quantile <-0.45
  }
  
  cutoff <- quantile(
    positive$Expression,
    cutoff.quantile,
    na.rm=TRUE
  )
  
  ###########################################################
  ## Keep cells above cutoff
  ###########################################################
  
  positive <- positive[
    positive$Expression >= cutoff,
  ]
  
  ###########################################################
  ## Draw low first
  ###########################################################
  
  positive <- positive[
    order(positive$Expression),
  ]
  
  ###########################################################
  ## Clip upper values
  ###########################################################
  
  upper <- max(
    
    as.numeric(
      
      quantile(
        
        positive$Expression,
        
        probs = clip.quantile,
        
        na.rm = TRUE
        
      )
      
    ),
    
    cutoff
    
  )
  
  positive$Expression <-
    
    pmin(
      positive$Expression,
      upper
    )
  
  ###########################################################
  ## Legend setting
  ###########################################################
  
  if(gene == "AOC3"){
    
    legend.limit  <- c(0, 3.8)
    
    legend.breaks <- seq(0, 3.5, by = 0.5)
    
    legend.labels <- sprintf("%.1f", legend.breaks)
    
  }else{
    
    legend.limit  <- c(0, 4.3)
    
    legend.breaks <- seq(0, 4.0, by = 0.5)
    
    legend.labels <- sprintf("%.1f", legend.breaks)
    
  }
  
 
  ###########################################################
  ###########################################################
  ##
  ## Step 2.2
  ## Background + Expression + Colour
  ##
  ###########################################################
  
  ###########################################################
  ## TISCH2 colour palette
  ###########################################################
  expr.cols <- c(
    
    "#D9D9D9",
    
    "#D2D0F6",
    
    "#ACA2E9",
    
    "#7E6AD4",
    
    "#5941C0",
    
    "#2C1B8F"
    
  )
  
  ###########################################################
  ## Background dataframe
  ###########################################################
  
  background <- df
  
  ###########################################################
  ## Build ggplot
  ###########################################################
  
  p <-
    
    ggplot(
      
      background,
      
      aes(
        
        x = UMAP_1,
        
        y = UMAP_2
        
      )
      
    ) +
    
    ###########################################################
  ## Background cells
  ###########################################################
  geom_point(
    
    colour = "#E0E0E0",
    
    size = 0.10,
    
    shape = 16,
    
    stroke = 0,
    
    alpha = 1
    
  ) +
    
    ###########################################################
  ## Positive cells
  ###########################################################
  geom_point(
    
    data = positive,
    
    aes(
      
      colour = Expression
      
    ),
    
    size = 0.12,
    
    shape = 16,
    
    stroke = 0,
    
    alpha = 0.80
    
  ) +
    
    ###########################################################
  ## Expression colour scale
  ###########################################################
  
  scale_colour_gradientn(
    
    colours = expr.cols,
    
    limits = legend.limit,
    
    breaks = legend.breaks,
    
    labels = legend.labels,
    
    oob = scales::squish,
    
    name = NULL,
    
    guide = guide_colourbar(
      
      direction = "vertical",
      
      title = NULL,
      
      ticks = TRUE,
      
      draw.ulim = TRUE,
      
      draw.llim = TRUE,
      
      frame.colour = NA,
      
      ticks.colour = "grey60"
      
    )
    
  ) +
    
    ###########################################################
  ## Keep UMAP aspect ratio
  ###########################################################
  
  coord_equal(
    
    expand = FALSE
    
  )
  
  
  ###########################################################
  ###########################################################
  ##
  ## Step 2.3
  ## Theme + Axis + Legend + Return
  ##
  ###########################################################
  
  p <-
    
    p +
    
    ###########################################################
  ## Axis
  ###########################################################
  
  scale_x_continuous(
    
    expand = expansion(mult = 0.01)
    
  ) +
    
    scale_y_continuous(
      
      expand = expansion(mult = 0.01)
      
    ) +
    
    ###########################################################
  ## Labels
  ###########################################################
  
  labs(
    
    x = NULL,
    
    y = NULL,
    
    title = NULL
    
  ) +
    
    ###########################################################
  ## Theme
  ###########################################################
  
  theme_classic(
    
    base_family = "Arial"
    
  ) +
    
    theme(
      
      #######################################################
      ## Panel
      #######################################################
      
      panel.background = element_blank(),
      
      panel.border = element_blank(),
      
      #######################################################
      ## Grid
      #######################################################
      
      panel.grid.major = element_line(
        
        colour = "grey88",
        
        linewidth = 0.25
        
      ),
      
      panel.grid.minor = element_blank(),
      
      #######################################################
      ## Axis
      #######################################################
      
      axis.line = element_line(
        
        colour = "black",
        
        linewidth = 0.30
        
      ),
      
      axis.ticks = element_line(
        
        colour = "black",
        
        linewidth = 0.30
        
      ),
      
      axis.ticks.length = unit(
        
        1,
        
        "mm"
        
      ),
      
      axis.text = element_blank(),
      
      axis.title = element_blank(),
      
      #######################################################
      ## Legend
      #######################################################
      legend.spacing.x = unit(2, "mm"),
      
      legend.position = "right",
      
      legend.direction = "vertical",
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        
        family = "Arial",
        
        size = 7,
        
        margin=margin(
          
          l=6),
        
        colour = "black"
        
      ),
      
      legend.key.height = unit(
        
        18,
        
        "mm"
        
      ),
      
      legend.key.width = unit(
        
        2,
        
        "mm"
        
      ),
      
      #######################################################
      ## Margin
      #######################################################
      
      plot.margin = margin(
        
        2,
        
        2,
        
        2,
        
        2,
        
        unit = "mm"
        
      )
      
    ) +
    
    ###########################################################
  ## Legend polish
  ###########################################################
  guides(
    
    colour = guide_colourbar(
      
      direction = "vertical",
      
      title = NULL,
      
      title.position = "top",
      
      ticks = TRUE,
      
      draw.ulim = TRUE,
      
      draw.llim = TRUE,
      
      frame.colour = NA,
      
      ticks.colour="grey60",
      
      ticks.linewidth=0.15,
      
      barheight = unit(
        
        47,
        
        "mm"
        
      ),
      
      barwidth = unit(
        
        3.5,
        
        "mm"
        
      )
      
    )
    
  )
  
  ###########################################################
  ## Return
  ###########################################################
  
  return(p)
  
}


plot_df <- build_plot_df(sce)


if (!exists("plot_df")) {
  
  message("plot_df not found.")
  
  message("Building plot_df ...")
  
  plot_df <<- build_plot_df(sce)
  
  message("Done.")
  
}


## AOC3
p_aoc3 <- plot_feature_umap_v5("AOC3")
print(p_aoc3)

## IL7R
p_il7r <- plot_feature_umap_v5("IL7R")
print(p_il7r)

## MS4A7
p_ms4a7 <- plot_feature_umap_v5("MS4A7")
print(p_ms4a7)  

##Figure 9B Final V5
#Step 3（CSBJ Submission Export）
###############################################################
##
## Figure 9B
## Final V5
##
## Step 3
## Export function
##
###############################################################

export_feature_umap <- function(
    p,
    gene,
    outdir = "figures/Figure9",
    width = 65,
    height = 60,
    dpi = 600
){
  
  ###########################################################
  ## Create folder
  ###########################################################
  
  if(!dir.exists(outdir)){
    dir.create(
      outdir,
      recursive = TRUE
    )
  }
  
  ###########################################################
  ## File names
  ###########################################################
  
  pdf.file <- file.path(
    outdir,
    paste0(gene, ".pdf")
  )
  
  png.file <- file.path(
    outdir,
    paste0(gene, ".png")
  )
  
  tif.file <- file.path(
    outdir,
    paste0(gene, ".tiff")
  )
  
  ###########################################################
  ## PDF (Vector)
  ###########################################################
  
  ggsave(
    filename = pdf.file,
    plot = p,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "mm",
    bg = "white"
  )
  
  ###########################################################
  ## PNG
  ###########################################################
  
  ggsave(
    filename = png.file,
    plot = p,
    width = width,
    height = height,
    units = "mm",
    dpi = dpi,
    bg = "white"
  )
  
  ###########################################################
  ## TIFF
  ###########################################################
  
  ggsave(
    filename = tif.file,
    plot = p,
    width = width,
    height = height,
    units = "mm",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  
  ###########################################################
  ## Message
  ###########################################################
  
  message("-----------------------------------")
  message("Export finished:")
  message(pdf.file)
  message(png.file)
  message(tif.file)
  message("-----------------------------------")
  
}

###############################################################
##
## Export Figure 9B
##
###############################################################

export_feature_umap(
  p_aoc3,
  "AOC3"
)

export_feature_umap(
  p_il7r,
  "IL7R"
)

export_feature_umap(
  p_ms4a7,
  "MS4A7"
)

# ===============================
# Figure 9C Dot Plot (Publication-ready)
# ===============================

# 1. Load package
library(tidyverse)
library(ggplot2)

# 2. Read data

data <- read.delim("NSCLC_GSE131907_AllDiffGenes_table.tsv")

colnames(data)


# 1. Extract target data
plot_data <- data %>%
  filter(Gene %in% c("IL7R", "MS4A7", "AOC3")) %>%
  mutate(Celltype = Celltype..major.lineage.)

# 2. Summary (to avoid duplication)
plot_data <- plot_data %>%
  group_by(Celltype, Gene) %>%
  summarise(
    avg_log2FC = mean(log2FC, na.rm = TRUE),
    pct_exp = mean(Percentage...., na.rm = TRUE),
    .groups = "drop"
  )

# 3.Setup sequence
plot_data$Gene <- factor(plot_data$Gene, levels = c("IL7R", "MS4A7", "AOC3"))

plot_data$Celltype <- factor(plot_data$Celltype,
                             levels = c("CD4Tconv", "CD8T", "CD8Tex",
                                        "B", "Plasma",
                                        "Mono/Macro", "DC", "Mast",
                                        "Endothelial", "Fibroblasts", "Epithelial"))

# 4. Drawing
p <- ggplot(plot_data, aes(x = Celltype, y = Gene)) +
  geom_point(aes(size = pct_exp, color = avg_log2FC)) +
  
  scale_size(range = c(1, 8), name = "Percent Expressed") +
  
  scale_color_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    name = "Relative Expression"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      face = "plain",
      size = 13, 
      color = "black"),
    axis.text.y = element_text(
      size = 13, 
      face = "italic", 
      color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12, 
      color = "black"),
    legend.text = element_text(
      face = "plain",
      size = 11, 
      color = "black")
  )

print(p)

# PDF
ggsave("Figure9C-3.pdf", p, width = 8, height = 4.5)

# TIFF
ggsave("Figure9C-3.tiff", p, width = 8, height = 4.5, dpi = 1200, compression = "lzw")

# PNG
ggsave("Figure9C-3.png", p, width = 8, height = 4.5, dpi = 300)

unique(data$Gene) %>% grep("MS4A7", ., value = TRUE)


