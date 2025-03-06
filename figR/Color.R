library("ggsci")
library("ggplot2")
library("gridExtra")
library("scales")
library("pals")

mypal <- pal_npg("nrc", alpha = 1)(10)

mypal
show_col(mypal)


null_color <- mypal[6]
tumor_color <- mypal[1]
normal_color <- mypal[3]


sample_color <- c("C3L-01632" = mypal[1],
                  "C3L-02549" = mypal[2],
                  "C3N-00169" = mypal[3],
                  "C3N-00199" = mypal[4],
                  "C3N-00547" = mypal[5],
                  "C3N-00579" = mypal[6],
                  "C3N-01016" = mypal[7],
                  "C3N-01024" = "#AC0000FF",
                  "C3N-01416" = mypal[9],
                  "C3N-02145" = mypal[10],
                  "GTEx" = "grey",
                  "mTEC" = "black")

category_color <- c("Canonical" = mypal[1], "Non-canonical" = mypal[3])

scRNA_sample_color <- c("LX653" = mypal[1],
                        "LX661" = mypal[2],
                        "LX675" = mypal[3],
                        "LX676" = mypal[4],
                        "LX679" = mypal[5],
                        "LX680" = mypal[6],
                        "LX682" = mypal[7],
                        "LX684" = mypal[8])

scRNA_type_color <- c("Primary tumor" = mypal[1])

scRNA_cluster_color <- c("0" = cols25(n=25)[1],
                         "1" = cols25(n=25)[2],
                         "2" = cols25(n=25)[3],
                         "3" = cols25(n=25)[4],
                         "4" = cols25(n=25)[5],
                         "5" = cols25(n=25)[6],
                         "6" = cols25(n=25)[7],
                         "7" = cols25(n=25)[8],
                         "8" = cols25(n=25)[9],
                         "9" = cols25(n=25)[10],
                         "10" = cols25(n=25)[1],
                         "11" = cols25(n=25)[2],
                         "12" = cols25(n=25)[3],
                         "13" = cols25(n=25)[4],
                         "14" = cols25(n=25)[5],
                         "15" = cols25(n=25)[6],
                         "16" = cols25(n=25)[7],
                         "17" = cols25(n=25)[8],
                         "18" = cols25(n=25)[9],
                         "19" = cols25(n=25)[10],
                         "20" = cols25(n=25)[1],
                         "21" = cols25(n=25)[2],
                         "22" = cols25(n=25)[3],
                         "23" = cols25(n=25)[4],
                         "24" = cols25(n=25)[5])

scRNA_cell_annotation_color <- c("Epithelial cell" = cols25(n=25)[2],
                                 "Endothelial cell" = cols25(n=25)[1],
                                 "Fibroblast" = cols25(n=25)[3],
                                 "T cell" = cols25(n=25)[9],
                                 "B cell" = cols25(n=25)[8],
                                 "Dendritic cell" = cols25(n=25)[6],
                                 "Macrophage" = cols25(n=25)[7],
                                 "Mast cell" = cols25(n=25)[25],
                                 "Natural killer cell" = cols25(n=25)[12],
                                 "Plasma cell" = cols25(n=25)[4],
                                 "Mesenchymal cell" = cols25(n=25)[13],
                                 "Monocyte" = cols25(n=25)[14])

scRNA_cell_simple_annotation_color <- c("Epithelial cell" = mypal[1],
                                        "Other" = mypal[3])


