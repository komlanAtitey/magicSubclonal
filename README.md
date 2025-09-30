# magicSubclonal
magicSubclonal is a physics-informed framework for subclone discovery and gene identification from bulk transcriptomes. It models driver-gene dynamics via the Chemical Master Equation, pinpoints rare states, links them to non-drivers, and assigns clinical risk, yielding calibrated and reproducible results across cancers.

$~~$

## Usage <br>

**GEO-based assessment**<br>
source("magicClonal_GEO.R")<br>
GEO_number <- "GSE9891" <br>  
***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>
***Number of coloumn for visualization***<br>
gene_column_number <- 3<br>
***Number of samples***<br>
number_sample <- 75<br>
***Run***<br>
run <- magicClonal_GEO(GEO_number, genes_of_interest, number_sample, gene_column_number)<br>
***List of subclone genes***<br>
subclone_gene <- run$high_risk_genes_in_order_unique<br>
subclone_gene_GSE171415 <- subclone_gene<br>
***Example of visualization***<br>
Barplot_tumorHeterogeneity <- run$Barplot_tumorHeterogeneity<br>
print(Barplot_tumorHeterogeneity)<br>
viz_subclone <- run$viz_subclone<br>
print(viz_subclone)<br>

