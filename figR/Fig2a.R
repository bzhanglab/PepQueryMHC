library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
source("GTExMapper.R")
source("Color.R")

pepquerymhc_tumor_antigens <- read_excel(path = "SupplementaryTable3.xlsx", sheet = "final tumor antigens (PQM)")
bamquery_tumor_antigens <- read_excel(path = "SupplementaryTable3.xlsx", sheet = "final tumor antigens (BQ)")
cutoff_1p <-read_excel(path = "SupplementaryTable3.xlsx", sheet = "p-value cutoff (log2(RPHM+1))")

########## Fig. 2a: Tumor antigen heatmap ##########
dev.off()
png("Fig2a.png", width = 6, height =3, units = 'in', res = 600)

bamquery <- bamquery_tumor_antigens %>% subset(select = c("LUAD_Tumor_Avg", "LUAD_Normal_Avg"))
bamquery <- bamquery %>% mutate(across(grep("Avg", colnames(bamquery)), ~ if_else(. == 0, NA_real_, .)))
colnames(bamquery) <- c("LUAD (10) ", "NAT (10) ")

heat_data <- pepquerymhc_tumor_antigens

heat_data$MHC <- factor(heat_data$MHC, 
                        levels = c("MHC-I", "MHC-II"))
heat_data$Label <- factor(heat_data$Label, levels = c("TSA", "TAA"))

category_count <- heat_data %>% dplyr::count(MHC)
group <- rep(unique(category_count$MHC), times = category_count$n)
tools <- rep(c("PepQueryMHC", "PepQueryMHC", "PepQueryMHC", "PepQueryMHC", "PepQueryMHC", "PepQueryMHC",
               "BamQuery", "BamQuery"))
tools <- factor(tools, levels = c("PepQueryMHC", "BamQuery"))

#### below to NA
heat_data <- 
  heat_data %>% mutate(across(grep("(GTEX)|(LUAD)|(Avg)|(mTEC)|(Max)|(Count)", colnames(heat_data)), ~ if_else(. == 0, NA_real_, .)))

heat_data <- heat_data %>% mutate(across(Label, ~ if_else(. == "TSA", 1, 0)))
heat_data <- heat_data %>% mutate(across(IEDB_MHC, ~ if_else(. == "Yes", 1, 0)))
heat_data <- heat_data %>% mutate(across(HLA_LIGAND_ATLAS, ~ if_else(. == "Yes", 1, 0)))
heat_data <- heat_data %>% mutate(across(IEAtlas_normal, ~ if_else(. == "Yes", 1, 0)))
heat_data <- heat_data %>% mutate(across(IEAtlas_cancer, ~ if_else(. == "Yes", 1, 0)))
heat_data <- heat_data %>% mutate(across(caAtlas_CT_antigen, ~ if_else(. == "Yes", 1, 0)))
heat_data <- heat_data %>% mutate(across(caAtlas_CA_antigen, ~ if_else(. == "Yes", 1, 0)))
heat_data$Report <- rowSums(heat_data[, grepl("(IEDB_MHC)|(HLA_LIGAND)|(IEAtlas)|(caAtlas)", colnames(heat_data))])
heat_data <- heat_data %>% mutate(across(Report, ~ if_else(. > 0, 1, 0)))

column_ha = HeatmapAnnotation(Label = anno_simple(heat_data$Label, col = c(`1` = "black", `0` = "white"), 
                                                  height = unit(3, "mm"), border = T),
                              Report = anno_simple(heat_data$Report, col = c(`1` = "black", `0` = "white"), 
                                                   height = unit(3, "mm"), border = T),
                              show_legend = F, 
                              show_annotation_name = T, 
                              annotation_name_gp = gpar(fontsize= 8)
)

heat_data[heat_data$Comment == ".", ]$Comment <- ""
heat_data[heat_data$Comment != "", ]$Comment <- paste(", ", heat_data[heat_data$Comment != "", ]$Comment, sep = "")

row_name <- heat_data$Peptide_ID
heat_data <- as.data.frame(heat_data[, grep("(GTEx_All_Avg)|(Testis_Tissue_Avg)|(Ovary_Tissue_Avg)|(Tumor_Avg)|(Normal_Avg)|(mTEC_Avg)", colnames(heat_data))])
heat_data <- heat_data[, order(colnames(heat_data))]
rownames(heat_data) <- row_name

colnames(heat_data)
colnames(heat_data) <- gsub("Testis_Tissue_Avg", "Testis (50)", colnames(heat_data))
colnames(heat_data) <- gsub("Ovary_Tissue_Avg", "Ovary (50)", colnames(heat_data))
colnames(heat_data) <- gsub("GTEx_All_Avg", "oGTEx (2303)", colnames(heat_data))
colnames(heat_data) <- gsub("LUAD_Tumor_Avg", "LUAD (10)", colnames(heat_data))
colnames(heat_data) <- gsub("LUAD_Normal_Avg", "NAT (10)", colnames(heat_data))
colnames(heat_data) <- gsub("mTEC_Avg", "mTEC (8)", colnames(heat_data))

## append bamquery
heat_data <- cbind(heat_data, bamquery)

heat_data <- as.data.frame(t(heat_data))
heat_data$Class <- rownames(heat_data)
heat_data$Class <- factor(heat_data$Class, levels = c("LUAD (10)", "NAT (10)", "Testis (50)", "Ovary (50)", "oGTEx (2303)", "mTEC (8)",
                                                      "LUAD (10) ", "NAT (10) "))
heat_data <- heat_data[order(heat_data$Class), ]
heat_data <- heat_data %>% subset(select = -c(Class))

mycols <- colorRamp2(breaks = c(0, 2, 4), 
                     colors = c(mypal[2],"white", mypal[1]))

plot <- Heatmap(heat_data, column_split = group,row_split = tools,
                name = "Log2(RPHM+1)", #title of legend
                column_title = "", row_title = "",
                row_names_gp = gpar(fontsize = 8, color = "black"), # Text size for row names
                
                row_title_gp = gpar(fontsize = 8, color = "black"),
                column_title_gp = gpar(fontsize = 8, color = "black"), column_title_side = "bottom",
                
                column_names_gp = gpar(fontsize = 8, color = "black"), column_names_side = "top",
                show_heatmap_legend = F,
                show_row_names = T,
                show_column_names = T,
                
                cluster_columns = F,  
                show_column_dend = F,
                cluster_rows = F,
                col = mycols, 
                
                heatmap_width = 4.5 * unit(1, "in"),
                heatmap_height = 2.6 * unit(1, "in"),
                
                top_annotation = column_ha,
                
                
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(is.na(heat_data[i, j]) | heat_data[i, j] == 0) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(fill = "#D3D1DF", lwd = 1, col = "black"))
                  }else if(heat_data[i, j] > 2) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "black", lwd = 1, fill = fill))
                    #grid.text("+", x, y, 
                    #          gp = gpar(col = "black", lwd = 1, fill = fill, fontsize = 25))
                  }else {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "black", lwd = 1, fill = fill))
                  }
                })

draw(plot, 
     legend_grouping = "original", 
     padding = unit(c(1, 0, 1, 0), "in"))

lgd <- Legend(legend_width = unit(0.8, "in"), grid_height = unit(0.1, "in"),
              direction = "vertical",
              title_position = "leftcenter-rot",
              at = c(0, 2, 4), 
              labels = c(0, 2, 4),
              col_fun = mycols,
              labels_gp = gpar(fontsize = 8, color = "black"),
              title_gp = gpar(fontsize = 8, color = "black", fontface = "bold"),
              title = expression("Mean"("Log"[2]("RPHM+1"))))

draw(lgd, x = unit(0.95, "npc"), y = unit(0.4, "npc"), just = c("right", "top"))


lgd_label <- Legend(legend_width = unit(0.8, "in"), grid_height = unit(0.1, "in"),
                    direction = "horizontal",
                    title_position = "topcenter",
                    labels = c("TSA", "TAA"),
                    legend_gp = gpar(fill = c("TSA" = "black", "TAA" = "white")), border = "black",
                    labels_gp = gpar(fontsize = 8, color = "black"),
                    title_gp = gpar(fontsize = 8, color = "black", fontface = "bold"),
                    title = "Label")

report_label <- Legend(legend_width = unit(0.8, "in"), grid_height = unit(0.1, "in"),
                       direction = "horizontal",
                       title_position = "topcenter",
                       labels = c("Yes", "No"),
                       legend_gp = gpar(fill = c("Yes" = "black", "No" = "white")), border = "black",
                       labels_gp = gpar(fontsize = 8, color = "black"),
                       title_gp = gpar(fontsize = 8, color = "black", fontface = "bold"),
                       title = "Report")

draw(lgd_label, x = unit(0.945, "npc"), y = unit(0.8, "npc"), just = c("right", "top"))
draw(report_label, x = unit(0.94, "npc"), y = unit(0.62, "npc"), just = c("right", "top"))

dev.off()
