library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source("Color.R")

reads_9mers <- read_excel("TableS2.xlsx", sheet = "40k 9mers (log2_RPHM+1)")
########## Fig. 1b: Performance measurement ##########
data <- read_excel("TableS2.xlsx", sheet = "Performance")
data_mem <- data %>% subset(select = c(Tool, Size, `Memory (GB)`))
data_run <- data %>% subset(select = c(Tool, Size, `Time (hour)`))
data_mem$Type <- "Memory usage"
data_run$Type <- "Runtime"

colnames(data_mem) <- c("Tool", "Size", "Estimate", "Type")
colnames(data_run) <- c("Tool", "Size", "Estimate", "Type")

data <- rbind(data_mem, data_run)

data$Size <- factor(data$Size, levels = c("1k", "5k", "10k", "20k", "30k", "40k"))
data$Type <- factor(data$Type, levels = c("Runtime", "Memory usage"))
plot_time_top <- ggplot(data, 
                        aes(x = Size, y = Estimate, color = Tool)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=0.5, notch=FALSE, outliers = F) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5) +
  facet_wrap(~Type, scales = "free_y", ncol = 2, 
             strip.position = "left", 
             labeller = as_labeller(c(`Memory usage` = "Memory (GB)", `Runtime` = "Runtime (hour)") ) ) +
  labs(title = "", x = "Number of 9-mers", y = "Runtime (hour)") + 
  guides(color = guide_legend(title = "Tool")) +
  theme_minimal() + 
  scale_y_continuous(n.breaks = 10) +
  theme(axis.text.x = element_text(hjust = 0.5, color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),
        legend.text = element_text(size = 15), 
        legend.title = element_blank(),
        legend.position= "top", 
        strip.placement = "outside", strip.text.y = element_text(color = "black", size = 15),
        plot.background = element_rect(fill = "white", colour = "white"), panel.spacing = unit(1, "lines")) +
  scale_color_manual(values = c("#4DBBD5FF", "#E64B35FF"))

plot_time_top
ggsave("Fig1b.png", plot_time_top, width = 8, height = 6, dpi = 600)



########## Fig. 1c: Estimated performance ##########
data <- read_excel("TableS2.xlsx", sheet = "Estimation")
data <- data %>% subset(select = c(Tool, Size, `Time (min)`))
colnames(data) <- c("Tool", "Size", "Estimate")
data[data$Size == "MHC-I (8-12-mers)", ]$Size <- "MHC-I\n(8-12-mers)"
data[data$Size == "MHC-II (7-25-mers)", ]$Size <- "MHC-II\n(7-25-mers)"

data$Size <- factor(data$Size, levels = c("MHC-I\n(8-12-mers)", "MHC-II\n(7-25-mers)"))

data_max <- data
data_max <- data_max[!duplicated(paste(data_max$Tool, data_max$Size)), ]
data_max[data_max$Tool == "PepQueryMHC" & data_max$Size == "MHC-I\n(8-12-mers)", ]$Estimate <- 
  max(data[data$Tool == "PepQueryMHC" & data$Size == "MHC-I\n(8-12-mers)", ]$Estimate)

data_max[data_max$Tool == "PepQueryMHC" & data_max$Size == "MHC-II\n(7-25-mers)", ]$Estimate <- 
  max(data[data$Tool == "PepQueryMHC" & data$Size == "MHC-II\n(7-25-mers)", ]$Estimate)


data_bamquery <- data_max[data_max$Tool == "BamQuery (estimated)", c("Size", "Estimate")]
data_pepqueryMHC <- data_max[data_max$Tool == "PepQueryMHC", c("Size", "Estimate")]

data_bamquery$PepQueryMHC <- data_pepqueryMHC$Estimate
colnames(data_bamquery) <- c("Size", "BamQuery (estimated)", "PepQueryMHC")
data <- data_bamquery
data$FoldChange <- data$`BamQuery (estimated)`/data$PepQueryMHC
data_max$Estimate <- log10(data_max$Estimate) ## log10

plot_time_top <- ggplot(data_max, 
                        aes(x = Tool, y = Estimate, fill = as.factor(Tool))) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "total"), size = 0.1, color = "black") + 
  labs(title = "", x = "", y = expression("Log"[10](Runtime(min)))) + 
  facet_wrap(Size ~ ., scales = 'free_y', strip.position = "bottom") + 
  guides(fill = guide_legend(title = "Tool", nrow = 2)) +
  theme_minimal() + 
  scale_y_continuous(n.breaks = 6, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 15, angle = 90, hjust = 0.5), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"), 
        legend.margin = margin(b = 20, r = 20),
        legend.text = element_text(size = 15), 
        legend.title = element_blank(),
        legend.position= "top", 
        strip.placement = "outside", strip.text.y = element_text(color = "black", size = 15), 
        strip.text.x = element_text(color = "black", size = 15),
        plot.background = element_rect(fill = "white", colour = "white"), panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))

plot_time_top
ggsave("Fig1c.png", plot_time_top, width = 3.4, height = 6, dpi = 600)



########## Fig. 1d ##########
## calculate single sample correlation (C3N-00547-T)
plot <- ggplot(reads_9mers[reads_9mers$Sample == "C3N-00547-T",], aes(x = BamQuery, y = PepQueryMHC)) + 
  geom_point(size = 1.5, color = mypal[4], alpha = 0.5) +
  labs(title = "C3N-00547-T (40k 9-mers)", x=expression("BamQuery: Log"[2](RPHM+1)), y=expression("PepQueryMHC: Log"[2](RPHM+1))) + 
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
plot
ggsave("Fig1d.png", plot, width = 5, height = 5, dpi = 600)



########## Fig. 1e ##########
## !! this figure needs bam file.
library(Gviz)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

## bam files
# bai file should be in the same folder
bam_file1 <- "C3N-00547-T.markdup.sorted.bam"
bam_file2 <- "C3N-00547-A.markdup.sorted.bam"
## gene model
gene_models <- import("gencode.v38.annotation.gtf", format = "GTF")
gene_models <- as.data.frame(gene_models)

names(gene_models)[names(gene_models) == 'seqnames'] <- 'chromosome'
names(gene_models)[names(gene_models) == 'gene_id'] <- 'gene'
names(gene_models)[names(gene_models) == 'exon_number'] <- 'exon'
names(gene_models)[names(gene_models) == 'transcript_id'] <- 'transcript'
names(gene_models)[names(gene_models) == 'gene_name'] <- 'symbol'

scheme <- getScheme()
scheme$GeneRegionTrack$fill <- mypal[2]
scheme$GeneRegionTrack$col <- NA
scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
scheme$GeneRegionTrack$collapseTranscripts <- FALSE
scheme$GeneRegionTrack$thinBoxFeature <- c("UTR")
scheme$GeneRegionTrack$col.line="black"
scheme$GeneRegionTrack$fontcolor.item="black"
scheme$GeneRegionTrack$background.panel <- "white"
scheme$GeneRegionTrack$background.title <- "white"
scheme$DataTrack$col.sampleNames <- "black"
scheme$DataTrack$background.title <- "white"

addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")
options(ucscChromosomeNames=FALSE)
getOption("Gviz.scheme")

## draw HUEW1 somatic mutations
from <- 53559493 - 5
to <-   53559519 + 18
maTrack <- GeneRegionTrack(rstarts = c(53559493), 
                           rends = c(53559519), 
                           strand = c("-"),
                           id = c("1"), group = c("1"),
                           genome="hg38", name = "SAAADILLL", chromosome = "chrX", transcript = " ")
maTrack@dp@pars$fill <- mypal[1]
maTrack@dp@pars$fontsize <- 0

aTrack1 <- AlignmentsTrack(range = bam_file1, genome = "hg38",
                           name = " ",
                           chromosome = "chrX", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack1@dp@pars$yTicksAt <- c(50, 100, 150, 200)
aTrack1@dp@pars$fontsize.title <- 15
aTrack1@dp@pars$col.axis <- "black"
aTrack1@dp@pars$background.panel <- "white"
aTrack1@dp@pars$background.title <- "white"
aTrack1@dp@pars$col <- mypal[4]

aTrack2 <- AlignmentsTrack(range = bam_file2, genome = "hg38",
                           name = " ",
                           chromosome = "chrX", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack2@dp@pars$yTicksAt <- c(50, 100, 150, 200)
aTrack2@dp@pars$fontsize.title <- 15
aTrack2@dp@pars$col.axis <- "black"
aTrack2@dp@pars$background.panel <- "white"
aTrack2@dp@pars$background.title <- "white"
aTrack2@dp@pars$col <- mypal[4]

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chrX")
idxTrack@dp@pars$fontcolor <- "black"
idxTrack@dp@pars$fontsize <- 15

axTrack <- GenomeAxisTrack()
axTrack@dp@pars$fontcolor <- "black"
axTrack@dp@pars$fontsize <- 15


huwe1 <- gene_models[gene_models$symbol == "HUWE1" & gene_models$type == "exon", ]
gene_track <- GeneRegionTrack(range = huwe1, 
                              genome = "hg38", chromosome = "chrX", 
                              from = from, to = to,
                              transcriptAnnotation = "symbol", collapseTranscripts = "longest")
gene_track@dp@pars$fontcolor.group <- "white"


sTrack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg38, chromosome = "chrX")
png("Fig1e.png", width = 7.5, height = 6, res = 600, units = "in")
plotTracks(list(idxTrack,axTrack, aTrack1, aTrack2, sTrack, maTrack, gene_track), sizes = c(1,2,8,8,1,1,1), 
           from = from, to = to, showTitle = T, stacking = "squish", title.width = 1, 
           type = c("coverage", "pileup"))
dev.off()


