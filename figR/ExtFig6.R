library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
source("Color.R")
source("GTExMapper.R")


## calculate mean
cal_mean <- function(data) {
  
  ## calculate tissue-wise average
  tissues <- colnames(data[, grep("GTEX", colnames(data), fixed= T)])
  tissues <- sapply(strsplit(tissues,".", fixed = T), `[`, 2)
  tissues <- unique(tissues)
  print(paste("Number of tissues: ", length(tissues), sep = ""))
  data[is.na(data)] <- 0
  
  for(tissue in tissues) {
    data[paste(tissue,"_Tissue_Avg", sep ="")] <- 
      apply(data[, grep(tissue, colnames(data))],1,function(x){mean(x[x>0])})
  }
  
  data[is.na(data)] <- 0
  ## Except for Testis, Ovary
  # cancer-germline antigens
  tmp <- data[, -grep("Testis", colnames(data))]
  tmp <- tmp[, -grep("Ovary", colnames(tmp))]
  print(length(data)-length(tmp))
  
  tmp <- tmp[, grep("GTEX", colnames(tmp), fixed= T)]
  print(length(tmp))
  data$GTEx_All_Avg <- apply(tmp,1,function(x){mean(x[x>0])})
  
  data$LUAD_Normal_Avg <- apply(data[, grep("-A$", colnames(data))],1,function(x){mean(x[x>0])})
  data$LUAD_Tumor_Avg <- apply(data[, grep("-T$", colnames(data))],1,function(x){mean(x[x>0])})
  
  data[is.na(data)] <- 0
  
  return (data)
}

draw_cutoff_figure <- function(data, y_label = "y_label") {
  
  
  length_dist_gtex <- data[data$GTEx_All_Avg > 0, c("Length", "GTEx_All_Avg")]
  length_dist_gtex$Category <- "GTEx"
  length_dist_luad_tumor <- data[data$LUAD_Tumor_Avg > 0, c("Length", "LUAD_Tumor_Avg")]
  length_dist_luad_tumor$Category <- "LUAD"
  length_dist_luad_normal <- data[data$LUAD_Normal_Avg > 0, c("Length", "LUAD_Normal_Avg")]
  length_dist_luad_normal$Category <- "NAT"
  
  colnames(length_dist_gtex) <- c("Length", "RPHM", "Category")
  colnames(length_dist_luad_tumor) <- c("Length", "RPHM", "Category")
  colnames(length_dist_luad_normal) <- c("Length", "RPHM", "Category")
  
  length_dist <- rbind(length_dist_gtex, length_dist_luad_normal, length_dist_luad_tumor)
  
  lengths <- unique(length_dist$Length)
  categories <- unique(length_dist$Category)
  cutoff_per_length <- data.frame()
  
  for(category in categories) {
    for(len in lengths) {
      p1 <- quantile(length_dist[length_dist$Length == len & length_dist$Category == category, c("RPHM")], 0.01)
      p5 <- quantile(length_dist[length_dist$Length == len & length_dist$Category == category, c("RPHM")], 0.05)
      len_data <- as.data.frame(t(c(category, len, p1, p5)))
      if(nrow(cutoff_per_length) == 0) {
        cutoff_per_length <- len_data
      } else {
        cutoff_per_length <- rbind(cutoff_per_length, len_data)
      }
    }
  }
  
  
  colnames(cutoff_per_length) <- c("Category", "Length", "1%", "5%")
  cutoff_per_length$Length <- as.numeric(cutoff_per_length$Length)
  cutoff_per_length$`1%` <- as.numeric(cutoff_per_length$`1%`)
  cutoff_per_length$`5%` <- as.numeric(cutoff_per_length$`5%`)
  
  length_dist <- left_join(length_dist, cutoff_per_length, by = c("Length", "Category"))
  
  ## present only GTEx except Testis and Ovary
  cutoff_per_length <- cutoff_per_length[cutoff_per_length$Category == "GTEx", ]
  length_dist <- length_dist[length_dist$Category == "GTEx", ]
  
  plot <- ggplot(length_dist, aes(x = RPHM)) +
    geom_histogram(aes(y=after_stat(count)), colour=mypal[2], fill="white", bins = 300)+
    labs(title = "", x = expression("Mean"("log"[2]("RPHM+1"))), y = y_label) + 
    theme_classic() + 
    geom_vline(aes(xintercept = `1%`), color = mypal[1], linewidth = 1) +
    geom_vline(aes(xintercept = `5%`), color = mypal[4], linewidth = 1) +
    scale_y_continuous(expand = c(0,0), n.breaks = 2) +
    scale_x_continuous(breaks = seq(0, 16), expand=c(0,0.2)) +
    facet_grid(Length ~ ., scales = "free_y") +
    theme(axis.text.x = element_text(color = "black", size = 15, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 15), 
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          legend.text = element_text(size = 15), 
          title = element_text(color = "black", size = 12),
          legend.title = element_text(size = 15), 
          strip.text.y = element_text(color = "black", angle = 0),
          strip.text = element_text(size = 13, color = "black"), panel.spacing = unit(1, "lines")) +
    geom_text(aes(label=format(round(`1%`, 3), nsmall = 3)), y=0.5, x=length_dist$`1%`-0.05 , colour="black", 
              size=4, angle = 90, hjust = 0, vjust=0) +
    geom_text(aes(label=format(round(`5%`, 3), nsmall = 3)), y=0.5, x=length_dist$`5%`+0.25 , colour="black", 
              size=4, angle = 90, hjust = 0, vjust=0)
  return(plot)
}

########## ExtFig. 6b: RPHM thresholds for class I ##########
## load data
MHC_I <- read.csv(
  "meta/LUAD_MHC_I_GTEX.tsv",
  header = T,
  check.names = FALSE,
  sep="\t"
)

MHC_I_mean <- cal_mean(MHC_I)
MHC_I_plot <- draw_cutoff_figure(MHC_I_mean, "Number of MHC-I-bound peptides")
ggsave("ExtFig6b.png", MHC_I_plot, width = 12, height = 5, dpi = 600)

########## ExtFig. 6c: RPHM thresholds for class II ##########
MHC_II <- read.csv(
  "meta/LUAD_MHC_II_GTEX.tsv",
  header = T,
  check.names = FALSE,
  sep="\t"
)

MHC_II_mean <- cal_mean(MHC_II)
MHC_II_plot <- draw_cutoff_figure(MHC_II_mean, "Number of MHC-II-bound peptides")
ggsave("ExtFig6c.png", MHC_II_plot, width = 12, height = 19, dpi = 600)
