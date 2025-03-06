library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(tidyr)
library(readxl)
source("Color.R")


########## ExtFig. 2a: The number of reads across 10 LUAD and 10 NAT samples ##########
libSize <- read_excel(path ="SupplementaryTable1.xlsx", sheet = "Reads")
## reads indicate the number of both primary and unmapped reads
libSize <- libSize[order(libSize$Reads), ]
libSize$Sample <- factor(libSize$Sample, levels = libSize$Sample)
libSize$Reads <- libSize$Reads / 1000000

plot <- ggplot(libSize, aes(x = as.factor(Sample), y = Reads, color = Type)) + 
  geom_point(size = 3) +
  labs(title = "", x = "Sample", y = expression("Number of reads (million)")) + 
  guides(color = guide_legend(title = "Type")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(hjust = 1, color = "black", size = 15, angle = 45),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.background = element_rect(colour="black", fill="white"),
        legend.position=c(0.9, 0.25)) +
  scale_color_manual(values = c("Tumor" = tumor_color, 
                                "Normal" = normal_color))
plot

ggsave(filename = "ExtFig2a.png", plot = plot, height = 4, width = 7, dpi = 600)

########## ExtFig. 2b: The number of MHC-I peptides according to lengths ##########
func_count <- function(data, label) {
  data <- data[!duplicated(data$Sequence), ]
  data$Length <- nchar(data$Sequence)
  category_count <- data %>% dplyr::count(Length)
  category_count$p <- category_count$n / sum(category_count$n)
  category_count$Label <- label
  return(category_count)
}

data_MHC_I <- read_excel("SupplementaryTable1.xlsx", sheet = "LUAD_MHC-I")
MHC_I <- func_count(data_MHC_I, "LUAD pMHC-I")

plot <- ggplot(MHC_I, aes(x = as.factor(Length), y = n, fill = n)) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), size = 0.5, color = "black") +
  labs(title = "", x = "Length", y = "Number of pMHC-Is") + 
  guides(color = guide_legend(title = "Type")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_gradient(low = "#4DBBD5FF",
                      high = "#E64B35FF")
plot

ggsave(filename = "ExtFig2b.png", plot = plot, height = 4, width = 3, dpi = 600)

########## ExtFig. 2c: The number of MHC-II peptides according to lengths ##########
data_MHC_II <- read_excel("SupplementaryTable1.xlsx", sheet = "LUAD_MHC-II")
MHC_II <- func_count(data_MHC_II, "LUAD pMHC-II")

plot <- ggplot(MHC_II, aes(x = Length, y = n, fill = n)) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), size = 0.5, color = "black") +
  labs(title = "", x = "Length", y = "Number of pMHC-IIs") + 
  guides(color = guide_legend(title = "Type")) + 
  scale_x_continuous(breaks = seq(9, max(MHC_II$Length), 3)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_gradient(low = "#4DBBD5FF",
                      high = "#E64B35FF")
plot

ggsave(filename = "ExtFig2c.png", plot = plot, height = 4, width = 5, dpi = 600)


