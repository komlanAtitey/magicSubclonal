setwd("/Users/atiteyk2/Documents/Master_Equation")
getwd()


## ============================================================
## Cancer-focused enrichment for a subclone gene set (ORA + viz)
## - Collect ALL results (p/q cutoffs = 1), then derive significance
## - Plot objects are created and saved; you can print them later
##   with: print(p_dot); print(p_dot_fallback); print(p_emap); print(p_go)
## ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("dplyr","tibble","readr","ggplot2","stringr","purrr")
bioc <- c("clusterProfiler","ReactomePA","org.Hs.eg.db","msigdbr","enrichplot","AnnotationDbi")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos="https://cloud.r-project.org")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(ggplot2)
  library(stringr); library(purrr)
  library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db)
  library(msigdbr); library(enrichplot); library(AnnotationDbi)
})

set.seed(42)
outdir <- "example_enrichment_out"; dir.create(outdir, showWarnings = FALSE)

## ---- 1) Inputs ----
# Use ONE of your lists (SYMBOLs)
#subclone_genes <- c("PRC1","EFNA4","IL20","PCSK1","ALKAL2","ZB16B","ALKAL2","CLSTN2")
#subclone_genes <- c("SST","LALBA","PCSK1","DGAT1","IGKV1D-13","H19","SOX11","NTS")
subclone_genes <- c("RN7SL1", "RN7SL2", "RN7SK", "RMRP", "EFS", "RPPH1", "CCDC3", "RPS16")
subclone_genes <- c("IGLV3-10", "S100A7", "MIA", "LTF", "AGR2", "HBB", "LTF", "GABRP")
  
subclone_genes <- unique(subclone_genes)

# Optional: provide expressed-gene universe (SYMBOLs) for better calibration
# expressed_genes <- rownames(expr_matrix)
expressed_genes <- NULL

## ---- 2) Helpers ----
map_symbols_to_entrez <- function(sym) {
  sym <- unique(na.omit(sym))
  clusterProfiler::bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

build_cancer_T2G <- function() {
  # Hallmarks
  H <- msigdbr(species = "Homo sapiens", category = "H") |>
    dplyr::transmute(term = gs_name, entrez_gene)
  
  # Full C2, then filter by subcat and name patterns (robust across msigdbr versions)
  C2 <- msigdbr(species = "Homo sapiens", category = "C2")
  
  KEGG <- C2 |>
    dplyr::filter(grepl("KEGG", gs_subcat, ignore.case = TRUE)) |>
    dplyr::filter(grepl("(cancer|carcinoma|leukemia|lymphoma|melanoma|glioma|sarcoma|MAPK|PI3K|AKT|MTOR|WNT|NOTCH|HEDGEHOG|TGFB|JAK|STAT|APOPTOSIS|CELL[_ ]?CYCLE)",
                        gs_name, ignore.case = TRUE)) |>
    dplyr::transmute(term = gs_name, entrez_gene)
  
  REACT <- C2 |>
    dplyr::filter(grepl("REACTOME", gs_subcat, ignore.case = TRUE)) |>
    dplyr::filter(grepl("(cancer|oncogenic|neoplasm|dna_repair|cell_cycle|apoptosis|metastasis|EMT|RAS|RAF|ERBB|EGFR|VEGF|MYC|TP53|RB1|NEUROENDOCRINE)",
                        gs_name, ignore.case = TRUE)) |>
    dplyr::transmute(term = gs_name, entrez_gene)
  
  # Oncogenic signatures (C6)
  C6 <- msigdbr(species = "Homo sapiens", category = "C6") |>
    dplyr::transmute(term = gs_name, entrez_gene)
  
  dplyr::bind_rows(H, KEGG, REACT, C6) |>
    dplyr::distinct(term, entrez_gene)
}

## ---- 3) Universe + mapping ----
if (is.null(expressed_genes)) {
  expressed_genes <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  message("NOTE: Using all human genes as the universe. For more power, pass your expressed genes.")
}
universe_df <- map_symbols_to_entrez(expressed_genes)
entrez_universe <- unique(universe_df$ENTREZID)

sub_df <- map_symbols_to_entrez(subclone_genes)
if (nrow(sub_df) == 0) stop("No subclone genes mapped. Check symbols/species.")
entrez_sub <- unique(sub_df$ENTREZID)

unmapped <- setdiff(subclone_genes, sub_df$SYMBOL)
if (length(unmapped)) message("Unmapped symbols (check aliases): ", paste(unmapped, collapse = ", "))

## ---- 4) Cancer-focused ORA (collect ALL results) ----
T2G_cancer <- build_cancer_T2G()
stopifnot(nrow(T2G_cancer) > 0)
T2N_cancer <- T2G_cancer |>
  dplyr::distinct(term) |>
  dplyr::transmute(term, term_name = term)

# Collect everything; do not filter here
ora_cancer_all <- enricher(
  gene          = entrez_sub,
  TERM2GENE     = T2G_cancer,
  TERM2NAME     = T2N_cancer,
  universe      = entrez_universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,          # <- collect all
  qvalueCutoff  = 1,          # <- collect all
  minGSSize     = 3,          # allow compact sets
  maxGSSize     = 2000
)

ora_tab <- as.data.frame(ora_cancer_all)
readr::write_tsv(ora_tab, file.path(outdir, "ORA_cancer_ALL.tsv"))

# Significant subset (q < 0.05)
sig_tab <- ora_tab |>
  dplyr::arrange(p.adjust) |>
  dplyr::filter(!is.na(p.adjust)) |>
  dplyr::filter(p.adjust < 0.05)

# Build an enrichResult containing only significant rows (if any)
ora_cancer_sig <- ora_cancer_all
if (nrow(sig_tab) > 0) {
  ora_cancer_sig@result <- dplyr::semi_join(ora_cancer_all@result, sig_tab, by = c("ID" = "ID"))
} else {
  ora_cancer_sig@result <- ora_cancer_all@result[0, ]
}

## ---- 5) GO BP (collect ALL results) ----
ego_bp_all <- enrichGO(
  gene          = entrez_sub,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  universe      = entrez_universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,         # collect all
  qvalueCutoff  = 1,         # collect all
  minGSSize     = 3,
  readable      = TRUE
)
write_tsv(as.data.frame(ego_bp_all), file.path(outdir, "ORA_GO_BP_ALL.tsv"))

ego_sig_tab <- as.data.frame(ego_bp_all) |>
  dplyr::filter(!is.na(p.adjust), p.adjust < 0.05)

## ---- 6) Build plot objects (no print in if-blocks) ----
p_dot <- NULL            # significant cancer terms
p_dot_fallback <- NULL   # preview (top 20 by p.adjust) if no significant
p_emap <- NULL           # emap of significant terms (if enough & has edges)
p_go <- NULL             # GO BP preview (top 20 by p.adjust)

# (a) Cancer ORA dotplot for significant terms (if any)
if (nrow(sig_tab) > 0) {
  p_dot <- tryCatch({
    enrichplot::dotplot(ora_cancer_sig, showCategory = min(20, nrow(sig_tab))) +
      ggplot2::ggtitle("Cancer-focused ORA (significant terms)")
  }, error = function(e) NULL)
}

# (b) Fallback/preview: top 20 by p.adjust from ALL cancer results
if (nrow(ora_tab) > 0 && is.null(p_dot)) {
  ora_top <- ora_cancer_all
  top_ids <- head(ora_tab$ID[order(ora_tab$p.adjust)], 20)
  ora_top@result <- dplyr::semi_join(ora_cancer_all@result,
                                     data.frame(ID = top_ids), by = "ID")
  p_dot_fallback <- tryCatch({
    enrichplot::dotplot(ora_top, showCategory = min(20, nrow(ora_top@result))) +
      ggplot2::ggtitle("Cancer-focused ORA (top terms, preview)")
  }, error = function(e) NULL)
}

# (c) emapplot only if enough significant terms AND non-zero similarity
if (nrow(sig_tab) >= 5) {
  ora_sim <- tryCatch(enrichplot::pairwise_termsim(ora_cancer_sig), error = function(e) NULL)
  if (!is.null(ora_sim)) {
    sim_mat <- tryCatch(attr(ora_sim, "termsim"), error = function(e) NULL)
    has_edges <- !is.null(sim_mat) && sum(sim_mat[upper.tri(sim_mat, diag = FALSE)] > 0, na.rm = TRUE) > 0
    if (has_edges) {
      p_emap <- tryCatch(
        enrichplot::emapplot(ora_sim, showCategory = min(30, nrow(sig_tab))),
        error = function(e) NULL
      )
    }
  }
}

# (d) GO BP preview plot (top 20 by p.adjust from ALL GO results)
ego_all_tab <- as.data.frame(ego_bp_all)
if (nrow(ego_all_tab) > 0) {
  ego_top <- ego_bp_all
  top_ids_go <- head(ego_all_tab$ID[order(ego_all_tab$p.adjust)], 20)
  ego_top@result <- dplyr::semi_join(ego_bp_all@result,
                                     data.frame(ID = top_ids_go), by = "ID")
  p_go <- tryCatch({
    enrichplot::dotplot(ego_top, showCategory = min(20, nrow(ego_top@result))) +
      ggplot2::ggtitle("GO BP ORA (top terms, preview)")
  }, error = function(e) NULL)
}

## ---- 7) Save plots if they exist (no printing here) ----
if (!is.null(p_dot)) {
  ggplot2::ggsave(file.path(outdir, "ORA_cancer_dotplot_SIG.png"),
                  plot = p_dot, width = 8, height = 5, dpi = 300)
}
if (!is.null(p_dot_fallback)) {
  ggplot2::ggsave(file.path(outdir, "ORA_cancer_dotplot_PREVIEW.png"),
                  plot = p_dot_fallback, width = 8, height = 5, dpi = 300)
}
if (!is.null(p_emap)) {
  ggplot2::ggsave(file.path(outdir, "ORA_cancer_emap_SIG.png"),
                  plot = p_emap, width = 9, height = 6, dpi = 300)
}
if (!is.null(p_go)) {
  ggplot2::ggsave(file.path(outdir, "ORA_GO_BP_dotplot_PREVIEW.png"),
                  plot = p_go, width = 8, height = 5, dpi = 300)
}

## ---- 8) Now you can print any plot independently in the console ----
## Example commands to run AFTER sourcing this script:
 print(p_dot)
 print(p_dot_fallback)
 print(p_emap)
 print(p_go)

cat("Done. Results saved to: ", normalizePath(outdir), "\n")






