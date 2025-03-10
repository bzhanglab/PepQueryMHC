library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
source("GTExMapper.R")
source("Color.R")


########## ExtFig. 7,8: Overall expression across the tumor antigens ##########
tumor_antigens <- read_excel(path = "SupplementaryTable3.xlsx", sheet = "final tumor antigens (PQM)")
cutoff_1p <-read_excel(path = "SupplementaryTable3.xlsx", sheet = "cutoff (log2(RPHM+1))")


### Extend version
target_sequences <- tumor_antigens$Peptide
target_data <- tumor_antigens[, grep("GTEX-|LUAD_|mTEC", colnames(tumor_antigens), fixed = F)]
target_data <- as.data.frame(target_data[, -grep("Avg", colnames(target_data), fixed = F)])
target_data <- as.data.frame(target_data[, -grep("Max", colnames(target_data), fixed = F)])
## the length of target_data must be 2431!


target_data$Sequence <- tumor_antigens$Peptide
target_data$ID <- tumor_antigens$Peptide_ID
target_data$ID <- factor(target_data$ID, levels = target_data$ID)
target_data$Label <- tumor_antigens$Label
target_data$MHC <- tumor_antigens$MHC
target_data$Gene <- tumor_antigens$Gene
target_data$Comment <- tumor_antigens$Category
target_data[target_data$Comment == "Neoantigen", ]$Comment <- tumor_antigens[tumor_antigens$Category == "Neoantigen", ]$Category

target_data$Mark <- "Not selected"
target_data[target_data$Sequence %in% target_sequences, ]$Mark <- "Selected"
target_data$Label <- factor(target_data$Label, levels = c("TSA", "TAA"))

target_data <- pivot_longer(target_data, cols = grep("GTEX-|LUAD_|mTEC", colnames(target_data)))

target_data$Sample <- paste(sapply(strsplit(target_data$name,".", fixed =T), `[`, 2), "(50)", sep =" ")
target_data[grep("-T$", target_data$name), ]$Sample <- "Lung adenocarcinoma (10)"
target_data[grep("-A$", target_data$name), ]$Sample <- "Normal adjacent tumor (10)"
target_data[grep("mTEC", target_data$name), ]$Sample <- "mTEC (8)"
## manually replace sample numbers
target_data[grep("Bladder", target_data$Sample, fixed = T), ]$Sample <- "Bladder (21)"
target_data[grep("Cervix ectocervix", target_data$Sample, fixed = T), ]$Sample <- "Cervix ectocervix (9)"
target_data[grep("Cervix endocervix", target_data$Sample, fixed = T), ]$Sample <- "Cervix endocervix (10)"
target_data[grep("Fallopian tube", target_data$Sample, fixed = T), ]$Sample <- "Fallopian tube (9)"
target_data[grep("Kidney medulla", target_data$Sample, fixed = T), ]$Sample <- "Kidney medulla (4)"
## 53 + 45*50 = 2303: good.

target_data[target_data$Sample == "Testis", ]$Sample <- "Testis (50)"
target_data[target_data$Sample == "Ovary", ]$Sample <- "Ovary (50)"

oGTEx_samples <- unique(target_data[setdiff(grep("GTEX", target_data$name, fixed = T), grep("(.Testis)|(.Ovary)", target_data$name, fixed = F)), ]$Sample)
target_data$Sample <- factor(target_data$Sample, levels = c("Lung adenocarcinoma (10)", "Normal adjacent tumor (10)", "Testis (50)", "Ovary (50)", oGTEx_samples, "mTEC (8)"))

target_data$Length <- nchar(target_data$Sequence)
target_data$Pass <- 1.99
target_data$IsOutlier <- F
target_data[target_data$value == 0, ]$value <- NA

lengths <- unique(target_data[target_data$MHC == "MHC-I", ]$Length)
for(len in lengths) {
  cutoff <- cutoff_1p[cutoff_1p$MHC == "MHC-I" & cutoff_1p$Length == len, ]$`1% cutoff`
  tmp <- target_data[target_data$MHC == "MHC-I" & target_data$Length == len, ] %>% mutate(across(value, ~ if_else(. < cutoff, 1.99, 2)))
  
  outliers <- target_data[target_data$MHC == "MHC-I" & target_data$Length == len, ] %>%
    group_by(Sample, ID) %>%
    mutate(IsOutlier = value < quantile(value, 0.25, na.rm = T) - 1.5 * IQR(value, na.rm = T) | 
             value > quantile(value, 0.75, na.rm = T) + 1.5 * IQR(value, na.rm = T))
  
  target_data[target_data$MHC == "MHC-I" & target_data$Length == len, ]$Pass <- tmp$value
  target_data[target_data$MHC == "MHC-I" & target_data$Length == len, ]$IsOutlier <- outliers$IsOutlier
  
}

lengths <- unique(target_data[target_data$MHC == "MHC-II", ]$Length)
for(len in lengths) {
  cutoff <- cutoff_1p[cutoff_1p$MHC == "MHC-II" & cutoff_1p$Length == len, ]$`1% cutoff`
  tmp <- target_data[target_data$MHC == "MHC-II" & target_data$Length == len, ] %>% mutate(across(value, ~ if_else(. < cutoff, 1.99, 2)))
  
  outliers <- target_data[target_data$MHC == "MHC-II" & target_data$Length == len, ] %>%
    group_by(Sample, ID) %>%
    mutate(IsOutlier = value < quantile(value, 0.25, na.rm = T) - 1.5 * IQR(value, na.rm = T) | 
             value > quantile(value, 0.75, na.rm = T) + 1.5 * IQR(value, na.rm = T))
  
  target_data[target_data$MHC == "MHC-II" & target_data$Length == len, ]$Pass <- tmp$value
  target_data[target_data$MHC == "MHC-II" & target_data$Length == len, ]$IsOutlier <- outliers$IsOutlier
}


colnames(target_data) <- c("Sequence", "ID","Label", "MHC","Gene","Comment", "Mark","Data", "RPHM", "Sample", "Length", "Pass", "IsOutlier")

target_data$Category <- "GTEx"
target_data[grep("LUAD_", target_data$Data, fixed = T), ]$Category <- sapply(strsplit(target_data[grep("LUAD_", target_data$Data), ]$Data,"_", fixed =T), `[`, 2)
target_data[grep("LUAD_", target_data$Data, fixed = T), ]$Category <- gsub("(-A)|(-T)", "",target_data[grep("LUAD_", target_data$Data, fixed = T), ]$Category)
target_data[target_data$Sample == "mTEC (8)", ]$Category <- "mTEC"


target_data <- target_data[!is.na(target_data$RPHM), ]
target_data$Category <- factor(target_data$Category, levels =  c("GTEx",
                                                                 "mTEC",
                                                                 "C3L-01632",
                                                                 "C3L-02549",
                                                                 "C3N-00169",
                                                                 "C3N-00199",
                                                                 "C3N-00547",
                                                                 "C3N-00579",
                                                                 "C3N-01016",
                                                                 "C3N-01024",
                                                                 "C3N-01416",
                                                                 "C3N-02145"))

target_data$Details <- paste(target_data$ID, ": ", target_data$Sequence, " (", target_data$Gene, ", ",target_data$Comment, ")", sep ="")
target_data <- target_data[order(target_data$ID), ]
target_data$Details <- factor(target_data$Details, levels = unique(target_data$Details))

## shape
target_data$Shape <- 19
target_data[target_data$IsOutlier, ]$Shape <- 17

## left/right mark
target_data[target_data$ID %in% unique(target_data$ID)[1:9], ]$Mark = "1"
target_data[target_data$ID %in% unique(target_data$ID)[10:21], ]$Mark = "2"

## plot ExtFig. 7
plot <- ggplot(target_data[target_data$Mark == "1", ], aes(x = Sample, y= RPHM)) +
  geom_boxplot(notch=FALSE, outliers = F) +
  labs(title = "", x = "", y = expression("Log"[2]("RPHM+1"))) + 
  geom_jitter(width = 0.2, height = 0, aes(color = Category), 
              size = 2.5,
              alpha = 0.7, 
              shape = target_data[target_data$Mark == "1", ]$Shape) +
  facet_wrap(Details ~ ., scales = "free_y",  strip.position = "top", ncol = 1) + 
  theme_minimal() + 
  scale_y_continuous(limits = c(0, max(target_data$RPHM+2)), expand = c(0,0)) +
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        title = element_text(color = "black", size = 12),
        legend.title = element_text(size = 15), legend.position = "right",
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.margin = margin(0, 20, 0, 100),
        panel.spacing.y = unit(1, "lines")) +
  geom_text(aes(label=..count..), y=max(target_data$RPHM+1.1), stat='count', colour="black", size=5) +
  guides(color = guide_legend(title = "Sample", ncol = 1), fill = "none") +
  scale_color_manual(values = sample_color)

plot
ggsave("ExtFig7.png", plot, width = 25, height = 18, dpi = 600)

## plot ExtFig. 8
plot <- ggplot(target_data[target_data$Mark == "2", ], aes(x = Sample, y= RPHM)) +
  geom_boxplot(notch=FALSE, outliers = F) +
  labs(title = "", x = "", y = expression("Log"[2]("RPHM+1"))) + 
  geom_jitter(width = 0.2, height = 0, aes(color = Category), 
              size = 2.5,
              alpha = 0.7,
              shape = target_data[target_data$Mark == "2", ]$Shape) +
  facet_wrap(Details ~ ., scales = "free_y",  strip.position = "top", ncol = 1) + 
  theme_minimal() + 
  scale_y_continuous(limits = c(0, max(target_data$RPHM+2)), expand = c(0,0)) +
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        title = element_text(color = "black", size = 12),
        legend.title = element_text(size = 15), legend.position = "right",
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.margin = margin(0, 20, 0, 100),
        panel.spacing.y = unit(1, "lines")) +
  geom_text(aes(label=..count..), y=max(target_data$RPHM+1.1), stat='count', colour="black", size=5) +
  guides(color = guide_legend(title = "Sample", ncol = 1), fill = "none") +
  scale_color_manual(values = sample_color)

plot
ggsave("ExtFig8.png", plot, width = 25, height = 18, dpi = 600)
