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
library(reshape2)
library(readxl)
source("Color.R")
options(future.globals.maxSize = 1e9)

########## Fig. 2b: Highly-variable peptides ##########
w_data_frame <- LoadSeuratRds(file = "meta/MSK_LUAD_Primary_seurat_object.rds")
length(unique(w_data_frame$cca_clusters)) ## 21 clusters

#### Naming convention ###
w_data_frame$sample <- gsub("PRIMARY_TUMOR", "primary tumor", w_data_frame$sample, fixed = T)
w_data_frame$sample <- gsub("_", "-", w_data_frame$sample, fixed = T)
w_data_frame$type <- gsub("PRIMARY_TUMOR", "Primary tumor", w_data_frame$type, fixed = T)

unique(w_data_frame$sample) ## samples



###### Load the top 10% most abundant peptides #####
RPHT_data <- read.csv("meta/PepQueryMHC_top10p_result.tsv", header = T, sep = "\t")
# drop Other and null
RPHT_data <- RPHT_data %>% filter(!grepl("(Other)|(Null)", Cell_id)) %>% as.data.frame()

# add cell anntations
cell_annotations <- w_data_frame$cell_annotation %>% as.data.frame()
cell_annotations$Cell_id <- rownames(cell_annotations)
split_ids <- sapply(strsplit(cell_annotations$Cell_id,"_"), `[`, c(1,2))

cell_annotations$Sample <- paste(split_ids[1,], split_ids[2,], sep = "_")
cell_annotations$Type <- sapply(strsplit(cell_annotations$Sample,"_"), `[`, 2)
cell_annotations[cell_annotations$Type == "PRIMARY", ]$Type <- "Primary tumor"
colnames(cell_annotations) <- c("Annotation", "Cell_id", "Sample", "Type")

# merge detailed cell types into broader categories
cell_annotations[cell_annotations$Annotation == "Naive CD8+ T cell", ]$Annotation <- "T cell"
cell_annotations[cell_annotations$Annotation == "Conventional T(Tconv) cell", ]$Annotation <- "T cell"
cell_annotations[cell_annotations$Annotation == "Regulatory T(Treg) cell", ]$Annotation <- "T cell"
cell_annotations[cell_annotations$Annotation == "Germinal center B cell", ]$Annotation <- "B cell"
cell_annotations[cell_annotations$Annotation == "Plasmacytoid dendritic cell(pDC)", ]$Annotation <- "Dendritic cell"

total_annotation_count_per_sample <- cell_annotations %>% count(Annotation, Sample)
colnames(total_annotation_count_per_sample) <- c("Annotation", "Sample", "Total")

## load MHC_I/_II annotations
MHC_I_annotations <- read_excel("SupplementaryTable1.xlsx", sheet = "LUAD_MHC-I")
MHC_II_annotations <- read_excel("SupplementaryTable1.xlsx", sheet = "LUAD_MHC-II")

## refine the fields
MHC_I_annotations <- MHC_I_annotations[, c("Sequence", "Gene", "Type")]
colnames(MHC_I_annotations) <- c("Peptide", "Gene-MHC-I", "Type-MHC-I")
MHC_II_annotations <- MHC_II_annotations[, c("Sequence", "Gene", "Type")]
colnames(MHC_II_annotations) <- c("Peptide", "Gene-MHC-II", "Type-MHC-II")

RPHT_cells <- RPHT_data
RPHT_cells <- left_join(cell_annotations, RPHT_cells, by = "Cell_id")

## Top 50 CVmad
sample_types <- unique(RPHT_cells$Type)
cell_types <- unique(RPHT_cells$Annotation)

prop_list <- list()
count_list <- list()
prop_cutoff <- 0.5
for(sample_type in sample_types) {
  prop_frame <- data.frame()
  count_frame <- data.frame()
  tmp <- RPHT_cells[RPHT_cells$Type == sample_type, ]
  
  prop_cut <- data.frame()
  for(cell_type in cell_types) {
    sub_data <- tmp[tmp$Annotation == cell_type, ]
    n <- nrow(sub_data)
    row_count <- apply(sub_data %>% subset(select = -c(Annotation, Cell_id, Sample, Type)),2,function(x){sum(x > 0)}) %>% as.data.frame()
    
    row_prop <- row_count / n
    colnames(row_prop) <- cell_type
    colnames(row_count) <- cell_type
    
    row_prop$Peptide <- rownames(row_prop)
    row_count$Peptide <- rownames(row_count)
    row_prop <- row_prop[, c("Peptide", cell_type)]
    row_count <- row_count[, c("Peptide", cell_type)]
    
    if(nrow(prop_frame) == 0) {
      prop_frame <- row_prop
      count_frame <- row_count
      prop_cut <- row_prop[, cell_type] > prop_cutoff
    } else {
      prop_frame <- left_join(prop_frame, row_prop, by = "Peptide")
      count_frame <- left_join(count_frame, row_count, by = "Peptide")
      
      prop_cut <- prop_cut | (row_prop[, cell_type] > prop_cutoff)
    }
  }
  # upper and lower 0.25 => resulting in too small sd => bigger values (over than 20)
  # just note...
  
  prop_frame$MAD <- apply(prop_frame %>% subset(select = -c(Peptide)),1,function(x){mad(x, na.rm = T)})
  prop_frame$Median <- apply(prop_frame %>% subset(select = -c(Peptide, MAD)),1,function(x){median(x, na.rm = T)})
  
  prop_frame$CVmad <- prop_frame$MAD / prop_frame$Median
  prop_frame <- left_join(prop_frame, MHC_I_annotations, by = "Peptide")
  prop_frame <- left_join(prop_frame, MHC_II_annotations, by = "Peptide")
  
  prop_frame[!prop_cut, ]$CVmad <- -abs(prop_frame[!prop_cut, ]$CVmad)
  prop_list[[sample_type]] <- prop_frame
  count_list[[sample_type]] <- count_frame
}


## select top 50 peptides
top_peptides <- c()
for(sample_type in sample_types) {
  sub_data <- prop_list[[sample_type]]
  sub_data <- sub_data[order(sub_data$CVmad, decreasing = T), ]
  
  sub_data$GeneName <- sub_data$`Gene-MHC-I`
  sub_data[is.na(sub_data$`Gene-MHC-I`), ]$GeneName <- sub_data[is.na(sub_data$`Gene-MHC-I`), ]$`Gene-MHC-II`
  
  prop_list[[sample_type]] <- sub_data
  
  sub_data <- sub_data[!duplicated(sub_data$GeneName), ]
  top_peptides <- append(top_peptides, sub_data[seq(1,50,by=1), ]$Peptide)
}

top_peptides <- top_peptides[!duplicated(top_peptides)]

## select list containing top X peptides
top_frac_list <- list()
top_count_list <- list()
for(sample_type in sample_types) {
  sub_data <- prop_list[[sample_type]]
  top_frac_list[[sample_type]] <- sub_data[sub_data$Peptide %in% top_peptides, ]
  
  sub_data <- count_list[[sample_type]]
  top_count_list[[sample_type]] <- sub_data[sub_data$Peptide %in% top_peptides, ]
}


# Load libraries
unique(sample_types)
scPercent_primary_tumor <- top_frac_list[["Primary tumor"]] %>% subset(select = -c(MAD, Median, CVmad, `Gene-MHC-I`, `Type-MHC-I`, `Gene-MHC-II`, `Type-MHC-II`, GeneName))
count_primary_tumor <- top_count_list[["Primary tumor"]]


# Reshape data to long format
scPercent_primary_tumor_long <- melt(scPercent_primary_tumor, id.vars = "Peptide", value.name = "value")
count_primary_tumor_long <- melt(count_primary_tumor, id.vars = "Peptide", value.name = "count")

# Create dot plot
data_long <- rbind(scPercent_primary_tumor_long)
count_long <- rbind(count_primary_tumor_long)

data_long <- left_join(data_long, count_long, by = c("Peptide", "variable"))
data_long <- left_join(data_long, top_frac_list[['Primary tumor']][, c("Peptide", "GeneName")], by = "Peptide")
data_long <- data_long[order(data_long$GeneName), ]

## assign peptide id
peptide_ids <- data.frame()
peptides <- unique(data_long$Peptide)
for(id in c(1:length(peptides))) {
  
  peptide_id <- paste("PEPT", id, sep = "")
  
  if(nrow(peptide_ids) == 0) {
    peptide_ids <- as.data.frame(t(c(peptides[id], peptide_id)))
  } else {
    peptide_ids <- rbind(peptide_ids, c(peptides[id], peptide_id))
  }
}
colnames(peptide_ids) <- c("Peptide", "PeptideId")

data_long <- left_join(data_long, peptide_ids, by = "Peptide")
data_long$Present <- paste(data_long$PeptideId, " (",data_long$GeneName,")", sep = "")
data_long$Present <- factor(data_long$Present, levels = unique(data_long$Present))
data_long$variable <- factor(data_long$variable, 
                             levels = c("Epithelial cell", "Endothelial cell", "Mesenchymal cell", 
                                        "Plasma cell", "Fibroblast", "Mast cell", 
                                        "Monocyte", "Dendritic cell", 
                                        "T cell", "B cell", "Natural killer cell"))

plot <- ggplot(data_long[data_long$count > 0, ], aes(x = Present, y = variable)) +
  geom_point(aes(color = log10(count), size = value)) +
  scale_color_gradient(low = mypal[2], high = mypal[1]) +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(color = "black", size = 11, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0), 
        legend.text = element_text(size = 10), 
        title = element_text(color = "black", size = 12),
        legend.title = element_text(size = 10), legend.position = "right",
        panel.spacing.y = unit(1, "lines")) +
  labs(x = "Peptides", y = "Cell types", color = expression("Log"[10]("scSize")), size = "scProportion")
plot

## for size 50
ggsave(filename = "Fig2b.png", plot = plot, width = 10, height = 5, dpi = 600)

########## Fig. 2c: UMAP examples ##########
## umap data
umap_data <- w_data_frame@reductions[['umap.cca']]
umap_data <- umap_data@cell.embeddings %>% as.data.frame()
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
umap_data$cluster <- w_data_frame$cca_clusters
umap_data$sample <- w_data_frame$sample
umap_data$type <-  w_data_frame$type
umap_data$cell_annotation <- w_data_frame$cell_annotation

ids <- w_data_frame@active.ident %>% as.data.frame()
umap_data$cell_id <- rownames(ids)

umap_data[umap_data$cell_annotation == "Naive CD8+ T cell", ]$cell_annotation <- "T cell"
umap_data[umap_data$cell_annotation == "Conventional T(Tconv) cell", ]$cell_annotation <- "T cell"
umap_data[umap_data$cell_annotation == "Regulatory T(Treg) cell", ]$cell_annotation <- "T cell"
umap_data[umap_data$cell_annotation == "Germinal center B cell", ]$cell_annotation <- "B cell"
umap_data[umap_data$cell_annotation == "Plasmacytoid dendritic cell(pDC)", ]$cell_annotation <- "Dendritic cell"
umap_data$cell_annotation <- factor(umap_data$cell_annotation, 
                                    levels = rev(c("Epithelial cell", "Endothelial cell", "Mesenchymal cell", 
                                                   "Plasma cell", "Fibroblast", "Mast cell", 
                                                   "Monocyte", "Dendritic cell", 
                                                   "T cell", "B cell", "Natural killer cell")))


example_peptides <- c("YAYEPADTA", "QTPGFLKLL", "KASYSGVSLF", "REIKLLSKEL")

plots <- list()
for(idx in c(1:length(example_peptides))) {
  tmp <- umap_data
  peptide <- example_peptides[idx]
  gene <- data_long[data_long$Peptide == peptide, ]$GeneName
  
  tmp$size <- 0.5
  tmp$stroke <- 0
  tmp$alpha <- "N"
  
  ## select peptide
  selected_cell_ids <- RPHT_data[RPHT_data[, peptide] > 0, ]$Cell_id
  
  tmp[tmp$cell_id %in% selected_cell_ids, ]$size <- 0.5
  tmp[tmp$cell_id %in% selected_cell_ids, ]$stroke <- 0
  tmp[tmp$cell_id %in% selected_cell_ids, ]$alpha <- "Y"
  
  plot <- ggplot(tmp, 
                 aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(size = tmp$size, stroke = tmp$stroke, 
               aes(fill = cell_annotation, color = cell_annotation, alpha = alpha)) + 
    labs(title = paste(peptide," (",gene,")", sep = ""), x = "UMAP-1", y = "UMAP-2") + 
    theme_classic() + 
    scale_x_continuous(limits = c(-16, 16)) +
    scale_y_continuous(limits = c(-16, 16)) +
    theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12), 
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.ticks.x = element_line(color="black"),
          axis.ticks.y = element_line(color="black"),
          legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15), 
          legend.position = "none",
          strip.text = element_text(size = 12, color = "black")) +
    
    scale_color_manual(values = scRNA_cell_annotation_color) +
    scale_fill_manual(values = scRNA_cell_annotation_color) + 
    scale_alpha_manual(values = c("Y" = 1, "N" = 0.05))
  
  plots[[idx]] <- plot
}


combine_plots <- wrap_plots(plots)
plot <- combine_plots + plot_layout(nrow = 2, byrow = F, guides = "collect", axis_titles = "collect") + 
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size=3), title = "Cell annotation")) +
  guides(fill = "none") + guides(alpha = "none")
plot
ggsave(filename = "Fig2c.png", plot = plot, width = 10, height = 7, dpi = 600)
