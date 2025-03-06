library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(Gviz)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
source("Color.R")


########## Fig. 2f: Gviz for cis-/trans-spliced peptide ##########
## !! this figure needs bam file.
## bam files
# bai file should be in the same folder
bam_file1 <- "C1R.markdup.sorted.bam"
bam_file2 <- "C1R-B5701.markdup.sorted.bam"

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
from <- 32157895 - 5000
to <-   32164766 + 5000

maTrack <- GeneRegionTrack(rstarts = c(32157895, 32164760), 
                           rends = c(32157911, 32164766), 
                           strand = c("-", "-"),
                           id = c("1", "1"), group = c("1"),
                           genome="hg38", name = "LTALSKLW", chromosome = "chr15", transcript = " ")
maTrack@dp@pars$fill <- mypal[1]
maTrack@dp@pars$fontsize <- 0

aTrack1 <- DataTrack(range = bam_file1, genome = "hg38", type = "l", 
                     name = " ",#C1R
                     chromosome = "chr15", start = from, end = to)
aTrack1@dp@pars$yTicksAt <- c(500, 1000, 1500, 2000)
aTrack1@dp@pars$fontsize.title <- 20
aTrack1@dp@pars$col.axis <- "black"
aTrack1@dp@pars$col <- mypal[4]

aTrack2 <- DataTrack(range = bam_file2, genome = "hg38", type = "l", 
                     name = " ", #C1R-B57:01
                     chromosome = "chr15", start = from, end = to)
aTrack2@dp@pars$yTicksAt <- c(500, 1000, 1500, 2000)
aTrack2@dp@pars$fontsize.title <- 20
aTrack2@dp@pars$col.axis <- "black"
aTrack2@dp@pars$col <- mypal[4]

idxTrack <- IdeogramTrack(genome="hg38", chromosome="chr15")
idxTrack@dp@pars$fontcolor <- "black"
idxTrack@dp@pars$fontsize <- 15

axTrack <- GenomeAxisTrack()
axTrack@dp@pars$fontcolor <- "black"
axTrack@dp@pars$fontsize <- 15

chrna7 <- gene_models[gene_models$symbol == "CHRNA7" & gene_models$type == "exon", ]
gene_track <- GeneRegionTrack(range = chrna7, 
                              genome = "hg38", chromosome = "chr15", 
                              from = from, to = to,
                              transcriptAnnotation = "symbol")
gene_track@dp@pars$fontcolor.group <- "white"

plot_list <- append(list(idxTrack,axTrack,aTrack1, aTrack2, maTrack), gene_track)

png("Fig2f.png", width = 8, height = 7, res = 600, units = "in")
plotTracks(plot_list, sizes = c(1,2,4,4,1,10), 
           from = from, to = to, showTitle = T, stacking = "squish", title.width = 1)

dev.off()

