# ml Anaconda3/2020.02
# source activate R-4.0.2-BABS
# R

library(tidyverse)
library(DESeq2)
library(emmeans)
library(plotly)
library(ggrepel)
library(biomaRt)
library("RColorBrewer")
library("gplots")
library(openxlsx)

library(ComplexHeatmap)
library(circlize)

load( "/camp/stp/babs/working/sopenam/projects/tybulewiczv/eva.lana-elola/12_RNASeq_HumanEmbryonicHearts_DownSyndrome/results/DESeq/DESeq_analysis.RData" ) 


projectDir <- "/camp/stp/babs/working/sopenam/projects/"
### Change the following entries for each project accordingly

PI <- "tybulewiczv/"
scientist <- "eva.lana-elola/"
project <- "12_RNASeq_HumanEmbryonicHearts_DownSyndrome/"

#directory where to export results
DESeq.Dir <- paste0(projectDir,PI, scientist, project, "/results/DESeq/")


coldat <- colData(rld1)

Oxphos <- c( "ABCB7","COX10","NDUFB2", "SDHA", "UQCRC2")
Oxphos <- merge(Oxphos, lookup, by.x=1, by.y=2, sort=FALSE)
cellcycle <- c( "CHEK1", "KIF2C", "RAD21", "RBL1",  "SMC4")
cellcycle <- merge(cellcycle, lookup, by.x=1, by.y=2, sort=FALSE)

## Plots for Dp1tyb/dyrk1a & Cell Cycle

mat1 <- assay(rld1)[ rownames(assay(rld1)) %in% cellcycle$ensembl_gene_id, ]
colnames(mat1) <- rld1$treatment

mat2 <- t(scale(t(mat1)))
# define colours for scale
col_fun = colorRamp2(c(min(mat2), 0, max(mat2)), c("blue", "white", "red"))
col.pal <-  brewer.pal(10, "Paired")

column_ha = HeatmapAnnotation(treatment = coldat$treatment, 
                              col = list(treatment = c("control" = col.pal[2], 
                                                    "DS" = col.pal[4])
                              ))


jpeg( paste0(DESeq.Dir, "Heatmap_DS_control_CellCycle_20210518.jpeg"), width = 768, type= "cairo")
#par( mar=c(5, 10, 5, 2))
#par(oma = c(3,3,3,3))
 Heatmap(mat2, col= col_fun, 
        name = " ",
        cluster_rows = FALSE,
        #clustering_distance_columns = dst,
        #clustering_method_columns = mtd,
        top_annotation= column_ha,
        column_labels = coldat$treatment,
        #show_column_names= FALSE
        row_labels = cellcycle[,1],
        show_row_names = TRUE,
        column_order = c( which(colData$treatment == "control"),
                          which(colData$treatment == "DS")),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        width=unit(2,'in'),
        height= unit(2, 'in'))
dev.off()


## Plots for D1`tyb/dyrk1a & Cell Oxphos
mat1 <- assay(rld1)[ rownames(assay(rld1)) %in% Oxphos$ensembl_gene_id, ]
colnames(mat1) <- rld1$treatment

mat2 <- t(scale(t(mat1)))
# define colours for scale
col_fun = colorRamp2(c(min(mat2), 0, max(mat2)), c("blue", "white", "red"))
col.pal <-  brewer.pal(10, "Paired")

column_ha = HeatmapAnnotation(treatment = coldat$treatment, 
                              col = list(treatment = c("control" = col.pal[2], 
                                                      "DS" = col.pal[4])
                              ))


jpeg( paste0(DESeq.Dir, "Heatmap_DS_control_Oxphos_20210518.jpeg"), width = 768, type= "cairo")
#par( mar=c(5, 10, 5, 2))
#par(oma = c(3,3,3,3))
Heatmap(mat2, col= col_fun, 
        name = " ",
        cluster_rows = FALSE,
        #clustering_distance_columns = dst,
        #clustering_method_columns = mtd,
        top_annotation= column_ha,
        column_labels = coldat$treatment,
        #show_column_names= FALSE
        row_labels = Oxphos[,1],
        show_row_names = TRUE,
        column_order = c( which(colData$treatment == "control"),
                          which(colData$treatment == "DS")),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        width=unit(2,'in'),
        height= unit(2, 'in'))
dev.off()



