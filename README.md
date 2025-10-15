# magicSubclonal: Physics-Informed Subclone Discovery from Bulk Transcriptomes
magicSubclonal (Model-Aware, Gene-driven Inference of Clones and Subclones) is a physics-informed computational framework for identifying subclonal states and influential genes from bulk transcriptomic data. It models the stochastic dynamics of driver genes using the Chemical Master Equation (CME) to capture rare expression states, link them to non-driver genes, and quantify clinical risk. The framework delivers calibrated, reproducible, and biologically interpretable insights across diverse cancer types.<br>


$~~$

## Input Data<br>
magicSubclonal accepts two forms of input: <br>
1.	GEO Datasets: The loader handles GEO accessions in multiple formats and standardizes them into a consistent gene × sample expression matrix, where rows are HGNC gene symbols. When a processed GEO ExpressionSet (Series Matrix) is available, the rows may represent probes or gene IDs, and columns represent samples with numeric expression values. <br>

2.	Custom Gene Expression Matrices: Users can directly provide a gene expression matrix where rows correspond to gene names and columns to samples, containing normalized numeric expression values. <br>

$~~$

## Usage <br>
**GEO-based assessment**<br>
source("magicSubclonal_GEO.R")<br>
GEO_number <- "GSE9891" <br>  
***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>
***Number of coloumn for visualization***<br>
gene_column_number <- 3<br>
***Number of samples***<br>
number_sample <- 75<br>
***Run***<br>
run <- magicSubclonal_GEO(GEO_number, genes_of_interest, number_sample, gene_column_number)<br>
***List of subclone genes***<br>
subclone_gene <- run$high_risk_genes_in_order_unique<br>
subclone_gene_GSE171415 <- subclone_gene<br>
***Example of visualization***<br>
Barplot_tumorHeterogeneity <- run$Barplot_tumorHeterogeneity<br>
print(Barplot_tumorHeterogeneity)<br>
viz_subclone <- run$viz_subclone<br>
print(viz_subclone)<br>

$~~$

**Data matrix assessment**<br>
source("magicSubclonal_matrix.R")<br>
***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>
***Number of coloumn for visualization***<br>  
gene_column_number <- 3<br>
***Number of samples***<br>
number_sample <- 75<br>
run <- magicClonal_matrix(input_data_matrix, genes_of_interest, number_sample, gene_column_number)<br>





