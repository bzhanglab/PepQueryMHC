library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(tibble)
library(tidyr)
source("Color.R")

########## Fig. 2d: Proportion of type ##########
annotations <- read_excel("TableS5.xlsx", sheet = "C1R_Top_annotation")

annotations$Category <- "Reference"
annotations[annotations$`PepQueryMHC:Unique_class_code` != "PC", ]$Category <- "Non-reference"
annotations[annotations$`PepQueryMHC:Unique_class_code` == ".", ]$Category <- "Not found"

category_count <- annotations %>% dplyr::count(Type,Category)
category_count$prop <- as.double(category_count$n)

category_count[category_count$Type == "Cis-spliced", ]$prop <- category_count[category_count$Type == "Cis-spliced", ]$prop / 5
category_count[category_count$Type != "Cis-spliced", ]$prop <- category_count[category_count$Type != "Cis-spliced", ]$prop / 23
category_count$Category <- factor(category_count$Category, levels = c("Reference", "Non-reference", "Not found"))

plot_time_top <- ggplot(category_count, 
                        aes(x = Type, y = prop, fill = as.factor(Category))) + 
  geom_bar(stat = "identity", position = "stack", size = 0.1, color = "black") + 
  labs(title = "", x = "", y = "Proportion") + 
  guides(fill = guide_legend(title = "", nrow = 3)) +
  theme_minimal() + 
  scale_y_continuous(n.breaks = 5, expand = c(0,0)) +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 12, hjust = 0.5), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"), 
        legend.margin = margin(b = 20, r = 20),
        legend.text = element_text(size = 12), 
        legend.title = element_blank(),
        legend.position= "top", 
        strip.placement = "outside", 
        strip.text.y = element_text(color = "black", size = 12), 
        strip.text.x = element_text(color = "black", size = 12),
        plot.background = element_rect(fill = "white", colour = "white"), panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = c("Reference" = "#00A087FF", 
                               "Non-reference" = "#F39B7FFF",
                               "Not found" = "#DC0000FF"))

plot_time_top
ggsave(filename = "Fig2d.png", plot = plot_time_top, height = 7, width = 1.8, dpi = 600)



########## Fig.2e: RPHM value per each peptide ##########
pepquerymhc <- read_excel("TableS5.xlsx", sheet = "PepQueryMHC")
bamquery <- read_excel("TableS5.xlsx", sheet = "BamQuery")
cutoffs <- read_excel("TableS3.xlsx", sheet = "p-value cutoff (log2(RPHM+1))")
input_data <- read_excel("TableS1.xlsx", sheet = "Spliced-MHC-I")
names(input_data)[names(input_data) == 'Sequence'] <- "DeNovoSequence"

## set length
bamquery$Length <- nchar(bamquery$Sequence)
pepquerymhc$Length <- nchar(pepquerymhc$Sequence)

## IL sequence setting
bamquery$DeNovoSequence <- gsub("I", "L", bamquery$Sequence)
pepquerymhc$DeNovoSequence <- gsub("I", "L", pepquerymhc$Sequence)

## load cutoff
cutoffs <- cutoffs[cutoffs$MHC == "MHC-I", c("Length", "1% cutoff")]
bamquery <- left_join(bamquery, cutoffs, by = "Length")
bamquery$`C1R-Max` <- apply(bamquery[, c("C1R", "C1R-B5701")],1,max)
pepquerymhc <- left_join(pepquerymhc, cutoffs, by = "Length")
pepquerymhc$`C1R-Max` <- apply(pepquerymhc[, c("C1R", "C1R-B5701")],1,max)

## filter below cutoff
nrow(bamquery) ## 106
bamquery_filter <- bamquery[bamquery$`C1R-Max` > bamquery$`1% cutoff`, ]
nrow(bamquery_filter) ## 45

nrow(pepquerymhc) ## 364
pepquerymhc_filter <- pepquerymhc[pepquerymhc$`C1R-Max` > pepquerymhc$`1% cutoff`, ]
nrow(pepquerymhc_filter) ## 101

## select DeNovoSequence with the max
bamquery_filter <- bamquery_filter[order(bamquery_filter$`C1R-Max`, decreasing = T), ]
bamquery_filter <- bamquery_filter[!duplicated(bamquery_filter$DeNovoSequence), ]
pepquerymhc_filter <- pepquerymhc_filter[order(pepquerymhc_filter$`C1R-Max`, decreasing = T), ]
pepquerymhc_filter <- pepquerymhc_filter[!duplicated(pepquerymhc_filter$DeNovoSequence), ]

## log transform
bamquery_filter$`Log2(C1R-Max+1)` <- log2(bamquery_filter$`C1R-Max`+1)
pepquerymhc_filter$`Log2(C1R-Max+1)` <- log2(pepquerymhc_filter$`C1R-Max`+1)

## assign id
bamquery_filter$Id <- paste(bamquery_filter$Sequence, bamquery_filter$Location, bamquery_filter$Strand, sep ="_")
pepquerymhc_filter$Id <- paste(pepquerymhc_filter$Sequence, pepquerymhc_filter$Location, pepquerymhc_filter$Strand, sep ="_")

## assign type
bamquery_filter <- left_join(input_data[, c("DeNovoSequence", "Reported fusion", "Type")], bamquery_filter, by = "DeNovoSequence")
pepquerymhc_filter <- left_join(input_data[, c("DeNovoSequence", "Reported fusion", "Type")], pepquerymhc_filter, by = "DeNovoSequence")


## expression data
plot_bamquery <- bamquery_filter[, c("Type", "DeNovoSequence", "Sequence", "Location", "Strand", "Log2(C1R-Max+1)")]
plot_pepqueryMHC <- pepquerymhc_filter[, c("Type","DeNovoSequence", "Sequence", "Location", "Strand", "Log2(C1R-Max+1)")]

plot_bamquery$Tool <- "BamQuery"
plot_pepqueryMHC$Tool <- "PepQueryMHC"

plot_data2 <- rbind(plot_pepqueryMHC, plot_bamquery)
# NA to zero
plot_data2[is.na(plot_data2$`Log2(C1R-Max+1)`), ]$`Log2(C1R-Max+1)` <- 0

plot_data2$`Log2(C1R-Max+1)` <- as.numeric(plot_data2$`Log2(C1R-Max+1)`)


plot_data2 <- plot_data2[order(plot_data2$DeNovoSequence, decreasing = F), ]
plot_data2 <- plot_data2[order(plot_data2$`Log2(C1R-Max+1)`, decreasing = F), ]
plot_data2 <- plot_data2[order(plot_data2$Type, decreasing = T), ]
plot_data2$DeNovoSequence <- factor(plot_data2$DeNovoSequence, levels = unique(plot_data2$DeNovoSequence))

label_colors <- c("Cis-spliced" = "#4DBBD5FF", "Trans-spliced" = "#E64B35FF")
text_colors <- c(rep(label_colors[2], 23), rep(label_colors[1], 5))

## nullify
plot_data2[plot_data2$`Log2(C1R-Max+1)` == 0, ]$`Log2(C1R-Max+1)` <- NA

plot <- ggplot(plot_data2, aes(x = Tool, y = DeNovoSequence)) +
  geom_point(aes(color = Type, size = `Log2(C1R-Max+1)`)) +
  scale_color_manual(values=label_colors) + 
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = text_colors, size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 0), 
        legend.text = element_text(size = 12), 
        legend.margin = margin(b = 20, r = 20),
        plot.margin = margin(t = 20),
        title = element_text(color = "black", size = 12), 
        legend.title = element_text(size = 12), 
        legend.position = "right",
        panel.spacing.y = unit(1, "lines")) +
  labs(x = "", y = "Cis-/Trans-spliced peptides", size = expression("Log"[2]("RPHM+1")), size = "Type") +
  guides(color = guide_legend(title = "Type", override.aes = list(size=5)))
plot
ggsave(filename = "Fig2e.png", plot = plot, height = 7, width = 4, dpi = 600)
