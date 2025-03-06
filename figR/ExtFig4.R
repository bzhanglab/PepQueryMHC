library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source("Color.R")

########## ExtFig. 4: Correlation plots ##########
reads_9mers <- read_excel("TableS2.xlsx", sheet = "40k 9mers (log2_RPHM+1)")
## calculate correlation for all samples
sample_list <- unique(reads_9mers$Sample)
plot_list <- list()
for(target_sample in sample_list) {
  sub_data <- reads_9mers[reads_9mers$Sample == target_sample, ]
  
  plot <- ggplot(sub_data, aes(x = BamQuery, y = PepQueryMHC)) + 
    geom_point(size = 1.5, color = mypal[4], alpha = 0.5) +
    labs(title = "", x="", y="") + 
    guides(color = guide_legend(title = "Type")) +
    theme_minimal() + 
    theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
          axis.text.y = element_text(color = "black", size = 15), 
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.ticks.x = element_line(color="black"),
          axis.ticks.y = element_line(color="black"),
          legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15)) +
    # linear fitting
    geom_smooth(method = 'lm', color= mypal[1], fill=  mypal[2]) +
    stat_cor(alternative = "two.sided",show.legend = FALSE, color = "black", 
             label.x.npc = "left", label.y.npc = "top", r.digits = 3, p.digits = 3, label.sep = "\n", size = 4)
  
  plot_list <- append(plot_list, list(plot))
}

plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 4, labels = sample_list)
plot
ggsave("ExtFig4.png", plot, width = 20, height = 16, dpi = 600)