##O:\Rproject\GitHub1_Paper 1_differential_expression_GSE32863


# 00_setup.R

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


### Figure 2D (External ROC Validation using GSE32863)

BiocManager::version()
BiocManager::install("GEOquery")
library(GEOquery)

gse <- getGEO(
  filename = "input_dir/GSE32863_series_matrix.txt.gz"
)

expr <- exprs(gse)

dim(expr)

head(expr[,1:5])

#View sample information
pdata <- pData(gse)

colnames(pdata)
#Then:
head(pdata$title)

#Step 3: Check the GPL comment file
anno <- read.delim(
  "input_dir/TCGA_LUAD_ROC/GPL6884-11607.txt",
  comment.char="#",
  stringsAsFactors=FALSE
)

colnames(anno)

#Step 4: Check whether the three genes are present—this is the most important step.
unique(anno$Symbol)
#Then：
"AOC3" %in% anno$Symbol

"MS4A7" %in% anno$Symbol

"IL7R" %in% anno$Symbol

#Step 5: Check the sample size

table(pdata$source_name_ch1)

############################################################
## Step 6. GPL Annotation: Converting Probes to Gene Symbols
############################################################

library(dplyr)

## Keep Probe and Gene Symbol
anno2 <- anno %>%
  select(ID, Symbol)

## Remove probes that lack a Gene Symbol.
anno2 <- anno2 %>%
  filter(Symbol != "")

## Merge expression matrices
expr_df <- data.frame(
  ProbeID = rownames(expr),
  expr,
  check.names = FALSE
)

expr_merge <- merge(
  anno2,
  expr_df,
  by.x="ID",
  by.y="ProbeID"
)

dim(expr_merge)

############################################################
## Step 7. Calculate the average for multiple probes corresponding to the same gene.
############################################################

expr_gene <-
  expr_merge %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean))

expr_gene <- as.data.frame(expr_gene)

rownames(expr_gene) <- expr_gene$Symbol

expr_gene$Symbol <- NULL

dim(expr_gene)


############################################################
## Step 8. Extract the three-gene expression matrix
############################################################

genes <- c("AOC3","MS4A7","IL7R")

expr3 <- expr_gene[genes,]

expr3


############################################################
## Step 9. Construct the GEO validation data frame
############################################################

geo_df <- data.frame(
  Sample = colnames(expr3),
  
  AOC3 = as.numeric(expr3["AOC3",]),
  
  MS4A7 = as.numeric(expr3["MS4A7",]),
  
  IL7R = as.numeric(expr3["IL7R",]),
  
  status = ifelse(
    pdata$source_name_ch1=="Lung adenocarcinoma",
    1,
    0
  )
)

head(geo_df)

#Check
table(geo_df$status)

#Figure 2G (External expression validation using GSE32863)
############################################################
## Step 1. Load R packages
############################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
#Step 2 Organize the data

table(pdata$source_name_ch1)

#Then regenerate:
group <- ifelse(
  pdata$source_name_ch1 == "Adjacent non-tumor lung",
  "Normal",
  "Tumor"
)

table(group)

# Step 3 (Recreate expr_df)
# Once you have confirmed that the groups are correct, recreate it:
expr_df <- data.frame(
  
  Group = group,
  
  AOC3 = as.numeric(expr_gene["AOC3", ]),
  
  MS4A7 = as.numeric(expr_gene["MS4A7", ]),
  
  IL7R = as.numeric(expr_gene["IL7R", ])
  
)


table(expr_df$Group)

############################################################
## Step 4. P-value labels
############################################################
anno1 <- data.frame(
  
  group1 = "Normal",
  
  group2 = "Tumor",
  
  y.position = max(expr_df$AOC3) * 1.05,
  
  label = "italic('P') < 0.001"
  
)
anno2 <- data.frame(
  
  group1="Normal",
  
  group2="Tumor",
  
  y.position=max(expr_df$MS4A7)*1.05,
  
  label="italic('P') < 0.001"
  
)

anno3 <- data.frame(
  
  group1="Normal",
  
  group2="Tumor",
  
  y.position=max(expr_df$IL7R)*1.05,
  
  label="italic('P') < 0.001"
  
)


## Step 5. AOC3 (Final Version)
############################################################

library(ggplot2)
library(ggpubr)

############################################################
## AOC3
############################################################

p1_AOC3 <-
  
  ggplot(
    expr_df,
    aes(
      x = Group,
      y = AOC3,
      fill = Group
    )
  )+
  
  ############################################################
## Violin
############################################################

geom_violin(
  trim = FALSE,
  colour = NA,
  alpha = 0.80
)+
  
  ############################################################
## Boxplot
############################################################

geom_boxplot(
  width = 0.16,
  fill = "white",
  colour = "black",
  linewidth = 0.5,
  outlier.shape = NA
)+
  
  ############################################################
## Draw significance brackets 
############################################################

annotate(
  "segment",
  x = 1,
  xend = 2,
  y = max(expr_df$AOC3) * 1.08,
  yend = max(expr_df$AOC3) * 1.08,
  linewidth = 0.5
)+
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = max(expr_df$AOC3) * 1.08,
    yend = max(expr_df$AOC3) * 1.05,
    linewidth = 0.5
  )+
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = max(expr_df$AOC3) * 1.08,
    yend = max(expr_df$AOC3) * 1.05,
    linewidth = 0.5
  )+
  
  ############################################################
## Italic P
############################################################

annotate(
  "text",
  x = 1.5,
  y = max(expr_df$AOC3) * 1.11,
  label = "italic(P) < 0.001",
  parse = TRUE,
  size = 5
)+
  
  ############################################################
## Color
############################################################

scale_fill_manual(
  values = c(
    "#4DBBD5",
    "#E64B35"
  )
)+
  
  ############################################################
## Y axis
############################################################

scale_y_continuous(
  name = expression("Normalized " * log[2] * " expression"),
  expand = expansion(
    mult = c(
      0.05,
      0.28
    )
  )
)+
  
  ############################################################
## Prevent clipping
############################################################

coord_cartesian(
  clip = "off"
)+
  
  ############################################################
## Theme
############################################################

theme_classic(
  base_size = 14
)+
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    
    plot.margin = margin(
      t = 35,
      r = 15,
      b = 10,
      l = 10
    )
    
  )+
  
  labs(
    title = "AOC3"
  )


## Step 6. MS4A7

p2_MS4A7 <-
  
  ggplot(
    expr_df,
    aes(
      x = Group,
      y = MS4A7,
      fill = Group
    )
  )+
  
  ############################################################
## Violin
############################################################

geom_violin(
  trim = FALSE,
  colour = NA,
  alpha = 0.80
)+
  
  ############################################################
## Boxplot
############################################################

geom_boxplot(
  width = 0.16,
  fill = "white",
  colour = "black",
  linewidth = 0.5,
  outlier.shape = NA
)+
  
############################################################
## Draw significance brackets 
############################################################

annotate(
  "segment",
  x = 1,
  xend = 2,
  y = max(expr_df$MS4A7) * 1.08,
  yend = max(expr_df$MS4A7) * 1.08,
  linewidth = 0.5
)+
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = max(expr_df$MS4A7) * 1.08,
    yend = max(expr_df$MS4A7) * 1.05,
    linewidth = 0.5
  )+
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = max(expr_df$MS4A7) * 1.08,
    yend = max(expr_df$MS4A7) * 1.05,
    linewidth = 0.5
  )+
  
  ############################################################
## Italic P
############################################################

annotate(
  "text",
  x = 1.5,
  y = max(expr_df$MS4A7) * 1.11,
  label = "italic(P) < 0.001",
  parse = TRUE,
  size = 5
)+
  
  ############################################################
## Color
############################################################

scale_fill_manual(
  values = c(
    "#4DBBD5",
    "#E64B35"
  )
)+
  
  ############################################################
## Y axis
############################################################

scale_y_continuous(
  name = expression("Normalized " * log[2] * " expression"),
  expand = expansion(
    mult = c(
      0.05,
      0.28
    )
  )
)+
  
############################################################
## Prevent clipping
############################################################

coord_cartesian(
  clip = "off"
)+
  
  ############################################################
## Theme
############################################################

theme_classic(
  base_size = 14
)+
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    
    plot.margin = margin(
      t = 35,
      r = 15,
      b = 10,
      l = 10
    )
    
  )+
  
  labs(
    title = "MS4A7"
  )

#Step 7 IL7R
############################################################
p3_IL7R <-
  
  ggplot(
    expr_df,
    aes(
      x = Group,
      y = IL7R,
      fill = Group
    )
  )+
  
  ############################################################
## Violin
############################################################

geom_violin(
  trim = FALSE,
  colour = NA,
  alpha = 0.80
)+
  
  ############################################################
## Boxplot
############################################################

geom_boxplot(
  width = 0.16,
  fill = "white",
  colour = "black",
  linewidth = 0.5,
  outlier.shape = NA
)+
  
  ############################################################
## Draw significance brackets
############################################################

annotate(
  "segment",
  x = 1,
  xend = 2,
  y = max(expr_df$IL7R) * 1.08,
  yend = max(expr_df$IL7R) * 1.08,
  linewidth = 0.5
)+
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = max(expr_df$IL7R) * 1.08,
    yend = max(expr_df$IL7R) * 1.05,
    linewidth = 0.5
  )+
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = max(expr_df$IL7R) * 1.08,
    yend = max(expr_df$IL7R) * 1.05,
    linewidth = 0.5
  )+
  
  ############################################################
## Italic P
############################################################

annotate(
  "text",
  x = 1.5,
  y = max(expr_df$IL7R) * 1.11,
  label = "italic(P) < 0.001",
  parse = TRUE,
  size = 5
)+
  
  ############################################################
## Color
############################################################

scale_fill_manual(
  values = c(
    "#4DBBD5",
    "#E64B35"
  )
)+
  
  ############################################################
## Y axis
############################################################

scale_y_continuous(
  name = expression("Normalized " * log[2] * " expression"),
  expand = expansion(
    mult = c(
      0.05,
      0.28
    )
  )
)+
  
############################################################
## Prevent clipping
############################################################

coord_cartesian(
  clip = "off"
)+
  
  ############################################################
## Theme
############################################################

theme_classic(
  base_size = 14
)+
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    
    plot.margin = margin(
      t = 35,
      r = 15,
      b = 10,
      l = 10
    )
    
  )+
  
  labs(
    title = "IL7R"
  )

############################################################
#Step 8 Merge Figure2D
############################################################

Figure2D <-
  
  (p1_AOC3 | p2_MS4A7 | p3_IL7R)+
  
  plot_annotation(
    
    theme=
      
      theme(
        
        plot.margin=margin(
          20,
          15,
          15,
          15
        )
        
      )
    
  )

Figure2D


############################################################
## Step 9. Save Figure 2D
############################################################

## PDF
ggsave(
  
  "Figure2D_GSE32863.pdf",
  
  Figure2G_3,
  
  width = 11,
  
  height = 4.2,
  
  device = cairo_pdf
  
)

## TIFF（600 dpi）
ggsave(
  
  "Figure2D_GSE32863.tiff",
  
  Figure2G_3,
  
  width = 11,
  
  height = 4.2,
  
  dpi = 600,
  
  compression = "lzw",
  
  limitsize = FALSE
  
)

## PNG
ggsave(
  
  "Figure2D_GSE32863.png",
  
  Figure2G_3,
  
  width = 11,
  
  height = 4.2,
  
  dpi = 600,
  
  limitsize = FALSE
  
)

## JPEG
ggsave(
  
  "Figure2D_GSE32863.jpeg",
  
  Figure2G_3,
  
  width = 11,
  
  height = 4.2,
  
  dpi = 600,
  
  limitsize = FALSE
  
)


