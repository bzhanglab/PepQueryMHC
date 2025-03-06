library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(SeuratData)
library(gridExtra)
library(networkD3)
library(patchwork)
library(htmlwidgets)
source("Color.R")
options(future.globals.maxSize = 1e9)
w_data_frame <- LoadSeuratRds(file = "meta/MSK_LUAD_Primary_seurat_object.rds")

########## ExtFig. 9: UMAP by samples, clusters, and cell annotations ##########

#### Naming convention ###
w_data_frame$sample <- gsub("PRIMARY_TUMOR", "primary tumor", w_data_frame$sample, fixed = T)
w_data_frame$sample <- gsub("_", "-", w_data_frame$sample, fixed = T)
w_data_frame$type <- gsub("PRIMARY_TUMOR", "Primary tumor", w_data_frame$type, fixed = T)
w_data_frame$sample <- gsub("-primary tumor", "", w_data_frame$sample, fixed=T)

unique(w_data_frame$sample)

umap_data <- w_data_frame@reductions[['umap.cca']]
umap_data <- umap_data@cell.embeddings %>% as.data.frame()
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
umap_data$cluster <- w_data_frame$cca_clusters
umap_data$sample <- w_data_frame$sample
umap_data$type <-  w_data_frame$type
umap_data$cell_annotation <- w_data_frame$cell_annotation

ids <- w_data_frame@active.ident %>% as.data.frame()
umap_data$cell_id <- rownames(ids)

umap_data$sample <- factor(x = umap_data$sample, levels = c("LX653",
                                                            "LX661",
                                                            "LX675",
                                                            "LX676",
                                                            "LX679",
                                                            "LX680",
                                                            "LX682",
                                                            "LX684"))

umap_data$type <- factor(x = umap_data$type, levels = c("Primary tumor"))


##### ExtFig. 9a: By samples #####
plot <- ggplot(umap_data, 
               aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(size = 0.3, alpha = 0.5, shape = 19, stroke = 0, aes(color = sample)) + 
  labs(title = "", x = "UMAP-1", y = "UMAP-2") + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5), title = "Sample")) +
  scale_color_manual(values = scRNA_sample_color)
plot
ggsave("ExtFig9a.png", plot, width = 7, height = 6, dpi = 600)

##### ExtFig. 9b: By clusters #####
plot <- DimPlot(object = w_data_frame, reduction = "umap.cca", group.by = c("seurat_clusters"), label = T, 
                repel= T, alpha = 1, pt.size = 1, raster = T, raster.dpi = c(600,600)) +
  labs(title = "", x = "UMAP-1", y = "UMAP-2") + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5), title = "Cluster")) + 
  scale_color_manual(values = scRNA_cluster_color)

plot
ggsave("ExtFig9b.png", plot, width = 7, height = 6, dpi = 600)

##### ExtFig. 9c: By cell annotations #####
umap_data$unified_cell_annotation <- umap_data$cell_annotation
umap_data[umap_data$unified_cell_annotation == "Naive CD8+ T cell", ]$unified_cell_annotation <- "T cell"
umap_data[umap_data$unified_cell_annotation == "Conventional T(Tconv) cell", ]$unified_cell_annotation <- "T cell"
umap_data[umap_data$unified_cell_annotation == "Regulatory T(Treg) cell", ]$unified_cell_annotation <- "T cell"
umap_data[umap_data$unified_cell_annotation == "Germinal center B cell", ]$unified_cell_annotation <- "B cell"
umap_data[umap_data$unified_cell_annotation == "Plasmacytoid dendritic cell(pDC)", ]$unified_cell_annotation <- "Dendritic cell"

plot <- ggplot(umap_data, 
               aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(size = 0.3, alpha = 1, shape = 19, stroke = 0, aes(color = unified_cell_annotation)) + 
  labs(title = "", x = "UMAP-1", y = "UMAP-2") + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5), title = "Cell annotation")) +
  scale_color_manual(values = scRNA_cell_annotation_color)
plot
ggsave("ExtFig9c.png", plot, width = 7, height = 6, dpi = 600)
