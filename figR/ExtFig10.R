library(networkD3)
library(htmlwidgets)

########## ExtFig. 10: Reannotation of previously reported cis- and trans-spliced pMHC-I peptides ##########
## deconvolution of locations
# after removing values below cutoff
# run annotate function in PepQueryMHC
annotations <- read_excel("SupplementaryTable5.xlsx", sheet = "C1R_Top_annotation")

# drop null
annotations <- annotations[annotations$Comparison != "Not found", ]

type_to_peptide <- annotations
type_to_peptide <- type_to_peptide %>% dplyr::count(Type, `PepQueryMHC:Sequence`)

peptide_to_locus <- annotations
peptide_to_locus <- peptide_to_locus %>% dplyr::count(`PepQueryMHC:Sequence`, `PepQueryMHC:Location`)

locus_to_gene <- annotations
locus_to_gene <- locus_to_gene %>% dplyr::count(`PepQueryMHC:Location`, `PepQueryMHC:Gene_name`)

gene_to_class <- annotations
gene_to_class <- gene_to_class %>% dplyr::count(`PepQueryMHC:Gene_name`, `PepQueryMHC:Unique_class_code`)

colnames(type_to_peptide) <- c("source", "target", "value", "group")
colnames(peptide_to_locus) <- c("source", "target", "value", "group")
colnames(locus_to_gene) <- c("source", "target", "value", "group")
colnames(gene_to_class) <- c("source", "target", "value", "group")

sankey_data <- rbind(type_to_peptide, peptide_to_locus, locus_to_gene, gene_to_class)

sankey_data <- sankey_data[order(sankey_data$target), ]

links <- data.frame(
  source = sankey_data$source, # e.g., c(1,1,2,2,1,2)
  target = sankey_data$target, # e.g., c(A,A,C,D,B,C)
  value = sankey_data$value # e.g., c(10,15,1,100,20,30)
)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 12)

saveWidget(p, file = "sankey_diagram.html")









