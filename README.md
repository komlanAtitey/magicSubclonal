# magicSubclonal: Physics-Informed Subclone Discovery from Bulk Transcriptomes
magicSubclonal (Model-Aware, Gene-driven Inference of Clones and Subclones) is an AI physics-informed computational framework for identifying subclonal states and influential genes from bulk transcriptomic data. It models the stochastic dynamics of driver genes using the Chemical Master Equation (CME) to capture rare expression states, link them to non-driver genes, and quantify clinical risk. The framework delivers calibrated, reproducible, and biologically interpretable insights across diverse cancer types.<br>
[Watch the video](Figure/clonal_expansion_qt.mp4)
![Clonal Expansion](Figure/clonal_expansion_qt.gif)


$~~$
## Input Data<br>
magicSubclonal accepts two forms of input: <br>
1.	GEO Datasets: The loader handles GEO accessions in multiple formats and standardizes them into a consistent gene × sample expression matrix, where rows are HGNC gene symbols. When a processed GEO ExpressionSet (Series Matrix) is available, the rows may represent probes or gene IDs, and columns represent samples with numeric expression values. <br>

2.	Custom Gene Expression Matrices
Users may directly provide a gene expression matrix as input.<br>

For **R scripts**: the matrix should have genes as rows and samples as columns, containing normalized numeric expression values.<br>

For **Python scripts**: the input should be a CSV file where the first column is labeled “GeneSymbol”, listing the gene names, followed by columns corresponding to the samples.<br>

## Runtime<br>
Runtime for magicSubclonal varies with dataset size, model settings (CME fitting, bootstrapping, resampling), and hardware configuration. On a MacBook Pro M2 Max (12-core CPU, 30/38-core GPU, up to 96 GB unified memory), processing took approximately 17, 8, 7, and 5 minutes for the ovarian (285 samples), lung (67 samples), DCIS (67 samples), and breast (19 samples) cohorts, respectively. <br>

$~~$

## Usage: Python script <br>
***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>

***Number of samples***<br>
number_sample <- 200<br>

***Number of coloumn for subplot grid layout***<br>
ncol <- 3<br>

***Execute in the terminal***<br>
python3 run_magicsubclonal.py --csv ovarian_data.csv --genes "TP53,BRCA1,BRCA2,ARID1A,PIK3CA,PTEN" --samples 200 --ncol 3

## Usage: R script <br>
#--------------<br>
**GEO-based assessment**<br>
#--------------<br>
source("magicSubclonal_GEO.R")<br>
GEO_number <- "GSE9891" <br>  

***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>

***Number of coloumn for subplot grid layout***<br>
gene_column_number <- 3<br>

***Number of samples***<br>
number_sample <- 75<br>

***Run***<br>
run <- magicSubclonal_GEO(GEO_number, genes_of_interest, number_sample, gene_column_number)<br>

***List of subclone genes***<br>
subclone_gene <- run$high_risk_genes_in_order_unique<br>
subclone_gene_GSE171415 <- subclone_gene<br>

***Example of visualization***<br>

Heatmap_tumorHeterogeneity <- run$Heatmap_tumorHeterogeneity<br>
print(Heatmap_tumorHeterogeneity)
![](Figure/magicSubclonal_hetero.png)

Barplot_tumorHeterogeneity <- run$Barplot_tumorHeterogeneity<br>
print(Barplot_tumorHeterogeneity)<br>
![](Figure/magicSubclonal_fraction.png)

viz_subclone <- run$viz_subclone<br>
print(viz_subclone)<br>
![](Figure/magicSubclonal_genes.png)

$~~$

#--------------<br>
**Data matrix assessment**<br>
#--------------<br>
source("magicSubclonal_matrix.R")<br>
***List of driver genes***<br>
genes_of_interest <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")<br>

***Number of coloumn for visualization***<br>  
gene_column_number <- 3<br>

***Number of samples***<br>
number_sample <- 75<br>
run <- magicClonal_matrix(input_data_matrix, genes_of_interest, number_sample, gene_column_number)<br>






