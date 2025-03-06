library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(Gviz)
library(rtracklayer)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
source("Color.R")

## !! this figure needs bam file.
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

########## ExtFig.5a: Alternative mapping ##########
from <- 20125056 - 10
to <-   20125082 + 10

maTrack <- GeneRegionTrack(rstarts = c(20125056), 
                           rends = c(20125082), 
                           strand = c("-"),
                           id = c("1"), group = c("1"),
                           genome="hg38", name = "GPLQSVQVF", chromosome = "chr19", transcript = " ")
maTrack@dp@pars$fill <- mypal[1]
maTrack@dp@pars$fontsize <- 0

aTrack1 <- AlignmentsTrack(range = bam_file1, genome = "hg38",
                           name = " ",
                           chromosome = "chr19", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack1@dp@pars$yTicksAt <- c(1, 2, 3)
aTrack1@dp@pars$fontsize.title <- 15
aTrack1@dp@pars$col.axis <- "black"
aTrack1@dp@pars$background.panel <- "white"
aTrack1@dp@pars$background.title <- "white"
aTrack1@dp@pars$col <- mypal[4]

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chr19")
idxTrack@dp@pars$fontcolor <- "black"
idxTrack@dp@pars$fontsize <- 15

axTrack <- GenomeAxisTrack()
axTrack@dp@pars$fontcolor <- "black"
axTrack@dp@pars$fontsize <- 15

znf90 <- gene_models[gene_models$symbol == "ZNF90" & gene_models$type == "exon", ]
gene_track <- GeneRegionTrack(range = znf90, 
                              genome = "hg38", chromosome = "chr19", 
                              from = from, to = to,
                              transcriptAnnotation = "symbol", 
                              collapseTranscripts = "longest")
gene_track@dp@pars$fontcolor.group <- "white"
sTrack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg38, chromosome = "chr19")

png("ExtFig5a.png", width = 7.5, height = 6, res = 600, units = "in")
plotTracks(list(idxTrack, axTrack, aTrack1, sTrack, maTrack, gene_track), sizes = c(1,2,6,1,1,1), 
           from = from, to = to, showTitle = T, stacking = "squish", title.width = 1, 
           type = c("coverage", "pileup"))
dev.off()



########## ExtFig. 5b: Synonymous mutations ##########
from <- 39435857 - 10
to <-   39435883 + 10

maTrack <- GeneRegionTrack(rstarts = c(39435857), 
                           rends = c(39435883), 
                           strand = c("-"),
                           id = c("1"), group = c("1"),
                           genome="hg38", name = "GPLQSVQVF", chromosome = "chr19", transcript = " ")
maTrack@dp@pars$fill <- mypal[1]
maTrack@dp@pars$fontsize <- 0

aTrack1 <- AlignmentsTrack(range = bam_file1, genome = "hg38",
                           name = " ",
                           chromosome = "chr19", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack1@dp@pars$yTicksAt <- c(1, 2, 3)
aTrack1@dp@pars$fontsize.title <- 15
aTrack1@dp@pars$col.axis <- "black"
aTrack1@dp@pars$background.panel <- "white"
aTrack1@dp@pars$background.title <- "white"
aTrack1@dp@pars$col <- mypal[4]


aTrack2 <- AlignmentsTrack(range = bam_file2, genome = "hg38",
                           name = " ",
                           chromosome = "chr19", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack2@dp@pars$yTicksAt <- c(1, 2, 3)
aTrack2@dp@pars$fontsize.title <- 15
aTrack2@dp@pars$col.axis <- "black"
aTrack2@dp@pars$background.panel <- "white"
aTrack2@dp@pars$background.title <- "white"
aTrack2@dp@pars$col <- mypal[4]

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chr19")
idxTrack@dp@pars$fontcolor <- "black"
idxTrack@dp@pars$fontsize <- 15

axTrack <- GenomeAxisTrack()
axTrack@dp@pars$fontcolor <- "black"
axTrack@dp@pars$fontsize <- 15

rps16 <- gene_models[gene_models$symbol == "RPS16" & gene_models$type == "exon", ]
gene_track <- GeneRegionTrack(range = rps16, 
                              genome = "hg38", chromosome = "chr19", 
                              from = from, to = to,
                              transcriptAnnotation = "symbol", collapseTranscripts = "longest")
gene_track@dp@pars$fontcolor.group <- "white"
sTrack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg38, chromosome = "chr19")

png("ExtFig5b.png", width = 7.5, height = 6, res = 600, units = "in")
plotTracks(list(idxTrack, axTrack, aTrack1, aTrack2, sTrack, maTrack, gene_track), sizes = c(1,2,6,6,1,1,1), 
           from = from, to = to, showTitle = T, stacking = "squish", title.width = 1, 
           type = c("coverage", "pileup"))
dev.off()



########## ExtFig. 5c: Matched read types ##########
percentMapping <- read_excel("TableS2.xlsx", sheet = "40k 9mers read type")

long_data <- pivot_longer(percentMapping, cols = c("Wildtype", "Variant", "Unmapped"))
long_data$value <- as.numeric(long_data$value)
long_data <- long_data[order(long_data$value, decreasing = T), ]
long_data$Sample <- factor(long_data$Sample, levels = unique(long_data$Sample))
long_data$name <- factor(long_data$name, levels = rev(c("Wildtype", "Variant", "Unmapped")))
plot <- ggplot(long_data, 
               aes(x = Sample, y = value, fill = as.factor(name))) + 
  geom_bar(stat = "identity", position = "stack", size = 0.1, color = "black") + 
  labs(title = "", x = "", y = "Proportion of matched reads") + 
  guides(fill = guide_legend(title = "", nrow = 1, reverse = T)) +
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
  scale_fill_manual(values = c("Wildtype" = "#00A087FF", 
                               "Variant" = "#F39B7FFF",
                               "Unmapped" = "#DC0000FF"))

plot

ggsave("ExtFig5c.png", plot, width = 5, height = 6, dpi = 600)



########## ExtFig. 5d: Poor reads ##########
## We manually put the sequencing qualities.
from <- 71681823 - 10
to <-   71681849 + 10

maTrack <- GeneRegionTrack(rstarts = c(71681823), 
                           rends = c(71681849), 
                           strand = c("-"),
                           id = c("1"), group = c("1"),
                           genome="hg38", name = "IMPHSVPGL", chromosome = "chr13", transcript = " ")
maTrack@dp@pars$fill <- mypal[1]
maTrack@dp@pars$fontsize <- 0

aTrack1 <- AlignmentsTrack(range = bam_file1, genome = "hg38",
                           name = " ",
                           chromosome = "chr13", start = from, end = to, 
                           showMismatches = T, mismatchSummary=T)

aTrack1@dp@pars$yTicksAt <- c(1, 2, 3)
aTrack1@dp@pars$fontsize.title <- 15
aTrack1@dp@pars$col.axis <- "black"
aTrack1@dp@pars$background.panel <- "white"
aTrack1@dp@pars$background.title <- "white"
aTrack1@dp@pars$col <- mypal[4]

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chr13")
idxTrack@dp@pars$fontcolor <- "black"
idxTrack@dp@pars$fontsize <- 15

axTrack <- GenomeAxisTrack()
axTrack@dp@pars$fontcolor <- "black"
axTrack@dp@pars$fontsize <- 15

dach1 <- gene_models[gene_models$symbol == "DACH1" & gene_models$type == "exon", ]
gene_track <- GeneRegionTrack(range = dach1, 
                              genome = "hg38", chromosome = "chr13", 
                              from = from, to = to,
                              transcriptAnnotation = "symbol", collapseTranscripts = "longest")
gene_track@dp@pars$fontcolor.group <- "white"
sTrack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg38, chromosome = "chr13")

png("ExtFig5d.png", width = 7.5, height = 6, res = 600, units = "in")
plotTracks(list(idxTrack, axTrack, aTrack1, sTrack, maTrack, gene_track), sizes = c(1,2,6,1,1,1), 
           from = from, to = to, showTitle = T, stacking = "squish", title.width = 1, 
           type = c("coverage", "pileup"))
dev.off()
