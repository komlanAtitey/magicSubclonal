
magicSubclonal_GEO <- function(GEO_number, genes_of_interest, number_sample, gene_column_number) {
  
  
  # -----------------------------
  # Step 1: Load GEO dataset and Map driver genes 
  # ----------------------------
  sentence1 <- "Load GEO dataset and Map driver genes" # @@@@@@@@@@@@
  print(sentence1) # @@@@@@@@@@@@
  
  suppressPackageStartupMessages({
    library(GEOquery)
    library(AnnotationDbi)
  })
  
  # --- Detect whether this GSE is a probe-based ExpressionSet we can map with hgu133plus2.db
  use_stmt1 <- FALSE
  gse_try <- try(GEOquery::getGEO(GEO_number, GSEMatrix = TRUE, AnnotGPL = TRUE), silent = TRUE)
  
  if (!inherits(gse_try, "try-error") && !is.null(gse_try)) {
    esets <- if (is.list(gse_try)) gse_try else list(gse_try)
    dims  <- sapply(esets, function(es) nrow(Biobase::exprs(es)))
    if (any(dims > 0)) {
      gset <- esets[[ which.max(dims) ]]
      ids  <- rownames(Biobase::exprs(gset))
      plat <- annotation(gset)
      
      # Probe-like if rownames look like Affy probes (e.g., *_at) and platform is GPL570
      affy_like <- length(ids) > 0 && any(grepl("_at$", ids, ignore.case = TRUE))
      have_pkg  <- requireNamespace("hgu133plus2.db", quietly = TRUE)
      
      use_stmt1 <- affy_like && identical(plat, "GPL570") && have_pkg
    }
  }
  
  if (use_stmt1) {
    # =======================
    # ===== Statement 1 =====
    # =======================
    sentence0 <- "Load GEO dataset!" # @@@@@@@@@@@@
    print(sentence0)
    
    pheno_data <- Biobase::pData(gset)
    
    sentence1 <- "Map driver genes to probes!" # @@@@@@@@@@@@
    print(sentence1)
    
    suppressPackageStartupMessages(library(hgu133plus2.db))
    
    probes <- AnnotationDbi::select(
      hgu133plus2.db,
      keys    = genes_of_interest,
      columns = c("PROBEID", "SYMBOL"),
      keytype = "SYMBOL"
    )
    
    # ---- OUTPUT 1 (exact line preserved) ----
    probes <- probes[!is.na(probes$PROBEID), ]
    exprs_data <- Biobase::exprs(gset)
    
  } else {
    # =======================
    # ===== Statement 2 =====
    # =======================
    build_symbol_matrix_from_GEO <- function(GEO_number,
                                             genes_of_interest = character(0),
                                             prefer_supp = TRUE,
                                             keep_unmapped = TRUE,
                                             verbose = TRUE) {
      # deps (load quietly)
      suppressPackageStartupMessages({
        library(GEOquery)
        library(AnnotationDbi)
        library(org.Hs.eg.db)
        library(data.table)
        library(matrixStats)
        library(R.utils)
        library(Matrix)
      })
      
      # helpers (scoped inside so you don't clutter your workspace)
      collapse_by_symbol <- function(mat) {
        if (is.null(mat) || nrow(mat) == 0) {
          return(matrix(numeric(0), nrow = 0, ncol = ncol(mat),
                        dimnames = list(character(0), colnames(mat))))
        }
        split_idx <- split(seq_len(nrow(mat)), rownames(mat))
        out <- lapply(split_idx, function(ix) {
          if (length(ix) == 1L) mat[ix, , drop = FALSE] else {
            med <- matrixStats::colMedians(mat[ix, , drop = FALSE], na.rm = TRUE)
            matrix(med, nrow = 1L, dimnames = list(rownames(mat)[ix][1], colnames(mat)))
          }
        })
        do.call(rbind, out)
      }
      guess_id_column <- function(df) {
        candidates <- c("UID","Gene","Geneid","GENE","symbol","SYMBOL","gene_id","gene_name","feature","id","ID")
        cand <- intersect(candidates, names(df))
        if (length(cand)) return(cand[1])
        ensg_re <- "^ENSG[0-9]+(\\.[0-9]+)?$"
        scores <- vapply(names(df), function(nm) mean(grepl(ensg_re, as.character(df[[nm]])), na.rm = TRUE), numeric(1))
        if (any(scores > 0.5)) return(names(df)[which.max(scores)])
        nn <- names(df)[!vapply(df, is.numeric, logical(1))]
        if (length(nn)) return(nn[1])
        names(df)[1]
      }
      map_ensembl_version_to_symbol <- function(ens_ids, keep_unmapped = TRUE) {
        stripped <- sub("\\..*$", "", ens_ids)
        look <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys   = unique(stripped),
                                      column = "SYMBOL",
                                      keytype= "ENSEMBL",
                                      multiVals = "first")
        sym <- unname(look[stripped])
        if (keep_unmapped) {
          na <- is.na(sym) | sym == ""
          sym[na] <- stripped[na]
        }
        sym
      }
      
      # 1) Try ExpressionSet first
      gse <- tryCatch(getGEO(GEO_number, GSEMatrix = TRUE, AnnotGPL = TRUE), error = function(e) NULL)
      if (!is.null(gse)) {
        esets <- if (is.list(gse)) gse else list(gse)
        dims  <- sapply(esets, function(es) nrow(exprs(es)))
        if (any(dims > 0)) {
          gset <- esets[[ which.max(dims) ]]
          exprs_data <- Biobase::exprs(gset)
          ids <- rownames(exprs_data)
          if (verbose) message("Loaded ExpressionSet: ", nrow(exprs_data), " features × ", ncol(exprs_data), " samples.")
          is_ensg <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids)
          SYMBOL  <- if (any(is_ensg)) map_ensembl_version_to_symbol(ids, keep_unmapped) else ids
          rownames(exprs_data) <- SYMBOL
          gene_mat <- collapse_by_symbol(exprs_data)
          non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
          return(list(
            source            = "GSEMatrix",
            data_with_symbols = NULL,
            gene_mat          = gene_mat,
            non_driver_mat    = non_driver_mat,
            mapping_table     = data.frame(orig_id = ids, SYMBOL = SYMBOL, stringsAsFactors = FALSE)
          ))
        }
      }
      if (!prefer_supp) stop("No non-empty ExpressionSet found and prefer_supp = FALSE.")
      
      # 2) Supplementary files
      base <- file.path(tempdir(), GEO_number)
      dir.create(base, showWarnings = FALSE, recursive = TRUE)
      GEOquery::getGEOSuppFiles(GEO_number, baseDir = base, makeDirectory = TRUE)
      supp_dir <- file.path(base, GEO_number)
      files <- list.files(supp_dir, full.names = TRUE)
      if (!length(files)) stop("No supplementary files found for ", GEO_number)
      
      gz <- grep("\\.gz$", files, value = TRUE)
      for (z in gz) try(R.utils::gunzip(z, overwrite = TRUE), silent = TRUE)
      files <- list.files(supp_dir, full.names = TRUE)
      
      cand <- grep("count|counts|matrix|tpm|fpkm|rpkm|expr|expression|rsem|genes?\\.tsv|\\.txt$|\\.csv$|matrix\\.mtx",
                   files, ignore.case = TRUE, value = TRUE)
      if (!length(cand)) stop("Could not detect an expression-like file in supplements.")
      f <- cand[ which.max(file.info(cand)$size) ]
      if (verbose) message("Reading supplementary file: ", basename(f))
      
      # 10X mtx path
      if (grepl("matrix\\.mtx", basename(f), ignore.case = TRUE)) {
        mm <- Matrix::readMM(f)
        bar <- grep("barcode", files, ignore.case = TRUE, value = TRUE)[1]
        feat<- grep("features|genes", files, ignore.case = TRUE, value = TRUE)[1]
        mat <- as.matrix(mm)
        if (!is.na(bar)) colnames(mat) <- data.table::fread(bar, header = FALSE)[[1]]
        if (!is.na(feat)) rownames(mat) <- data.table::fread(feat, header = FALSE)[[1]]
        ids <- rownames(mat)
        SYMBOL <- map_ensembl_version_to_symbol(ids, keep_unmapped)
        rownames(mat) <- SYMBOL
        gene_mat <- collapse_by_symbol(mat)
        non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
        return(list(
          source            = "supplement_10X",
          data_with_symbols = NULL,
          gene_mat          = gene_mat,
          non_driver_mat    = non_driver_mat,
          mapping_table     = data.frame(orig_id = ids,
                                         ENSEMBL_nover = sub("\\..*$","", ids),
                                         SYMBOL = SYMBOL, stringsAsFactors = FALSE)
        ))
      }
      
      # Delimited table path
      raw_dt <- data.table::fread(f)
      raw_df <- as.data.frame(raw_dt)  # avoid DT column-selection pitfalls
      id_col <- guess_id_column(raw_df)
      ids    <- as.character(raw_df[[id_col]])
      
      # numeric sample columns (coerce if needed)
      num_mask <- vapply(raw_df, is.numeric, logical(1))
      num_mask[ match(id_col, names(raw_df)) ] <- FALSE
      if (!any(num_mask)) {
        non_id <- setdiff(names(raw_df), id_col)
        for (nm in non_id) raw_df[[nm]] <- suppressWarnings(as.numeric(raw_df[[nm]]))
        num_mask <- vapply(raw_df, is.numeric, logical(1))
        num_mask[ match(id_col, names(raw_df)) ] <- FALSE
      }
      if (!any(num_mask)) stop("No numeric sample columns detected in supplementary table.")
      
      num_names <- names(raw_df)[num_mask]
      exprs_all <- as.matrix(raw_df[, num_names, drop = FALSE])
      colnames(exprs_all) <- make.unique(colnames(exprs_all), sep = "_")
      
      SYMBOL <- map_ensembl_version_to_symbol(ids, keep_unmapped)
      rownames(exprs_all) <- SYMBOL
      
      gene_mat <- collapse_by_symbol(exprs_all)
      non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
      
      data_with_symbols <- raw_df
      data_with_symbols$GeneSymbol <- SYMBOL
      ord <- c(id_col, "GeneSymbol", setdiff(names(data_with_symbols), c(id_col,"GeneSymbol")))
      data_with_symbols <- data_with_symbols[, ord]
      
      list(
        source            = "supplement_table",
        data_with_symbols = data_with_symbols,
        gene_mat          = gene_mat,
        non_driver_mat    = non_driver_mat,
        mapping_table     = data.frame(orig_id = ids,
                                       ENSEMBL_nover = sub("\\..*$","", ids),
                                       SYMBOL = SYMBOL, stringsAsFactors = FALSE)
      )
    }
    
    # Prefer using your helper; otherwise do a minimal inline fallback
    if (exists("run_step1_step2")) {
      out <- run_step1_step2(GEO_number, genes_of_interest = genes_of_interest)
      
      probes <- out$probes
      # ---- OUTPUT 2 (exact line preserved) ----
      probes <- probes[!is.na(probes$PROBEID), ]
      exprs_data <- out$gene_mat
      pheno_data <- out$pheno_data
      
    } else if (exists("build_symbol_matrix_from_GEO")) {
      # Minimal inline fallback if run_step1_step2() isn't defined
      res <- build_symbol_matrix_from_GEO(
        GEO_number        = GEO_number,
        genes_of_interest = genes_of_interest,
        prefer_supp       = TRUE,
        keep_unmapped     = TRUE,
        verbose           = TRUE
      )
      
      exprs_data <- res$gene_mat
      pheno_data <- data.frame(sample = colnames(exprs_data), row.names = colnames(exprs_data))
      
      # Build a PROBE-like table from drivers (SYMBOL/ENSEMBL)
      suppressPackageStartupMessages({
        library(org.Hs.eg.db)
        library(dplyr)
      })
      drivers_in <- as.character(genes_of_interest)
      is_ensg    <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", drivers_in)
      ensg_stripped <- ifelse(is_ensg, sub("\\..*$", "", drivers_in), NA_character_)
      sym_map <- rep(NA_character_, length(drivers_in))
      if (any(is_ensg)) {
        look <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys      = unique(ensg_stripped[is_ensg]),
                                      column    = "SYMBOL",
                                      keytype   = "ENSEMBL",
                                      multiVals = "first")
        sym_map[is_ensg] <- unname(look[ensg_stripped[is_ensg]])
      }
      sym_map[!is_ensg] <- drivers_in[!is_ensg]
      
      probes <- tibble::tibble(
        INPUT   = drivers_in,
        ENSEMBL = ifelse(is_ensg, ensg_stripped, NA_character_),
        SYMBOL  = sym_map
      ) %>%
        mutate(
          PROBEID = if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL,
                            if_else(!is.na(ENSEMBL) & ENSEMBL != "", ENSEMBL, NA_character_))
        ) %>%
        distinct(PROBEID, .keep_all = TRUE)
      
      # ---- OUTPUT 2 (exact line preserved) ----
      probes <- probes[!is.na(probes$PROBEID), ]
      
    } else {
      stop("Neither run_step1_step2() nor build_symbol_matrix_from_GEO() is defined for Statement 2.")
    }
  }
  
  # At this point you have:
  # - exprs_data (matrix): microarray exprs (Statement 1) OR SYMBOL-collapsed matrix (Statement 2)
  # - pheno_data (data.frame)
  # - probes (data.frame with PROBEID, SYMBOL, etc., and NA PROBEIDs removed)
  
  
  # -----------------------------
  # Step 2: Prepare expression data for CME simulations
  # -----------------------------
  sentence2 <- "Prepare expression data for CME simulations!" # @@@@@@@@@@@@
  print(sentence2) # @@@@@@@@@@@@
  
  exprs_sub <- exprs_data[rownames(exprs_data) %in% probes$PROBEID, ]
  rownames(exprs_sub) <- probes$SYMBOL[match(rownames(exprs_sub), probes$PROBEID)]
  
  exprs_sub[is.na(exprs_sub)] <- 0.01
  min_val <- min(exprs_sub)
  if(min_val <= 0){
    exprs_sub <- exprs_sub - min_val + 0.01
  }
  
  set.seed(123)
  n_samples <- number_sample
  x0_list <- lapply(rownames(exprs_sub), function(gene) {
    x0 <- exprs_sub[gene, ]
    x0 <- x0[!is.na(x0) & x0 > 0]
    sampled <- sample(x0, size = n_samples, replace = TRUE)
    sampled <- sampled * runif(n_samples, 0.9, 1.1)
    sampled
  })
  names(x0_list) <- rownames(exprs_sub)
  
  # -----------------------------
  # Step 4: Keep sample IDs for each resampled x0 (needed later)
  # -----------------------------
  sentence3 <- "Record sample IDs for each resampled x0!" # @@@@@@@@@@@@
  print(sentence3) # @@@@@@@@@@@@
  
  set.seed(123)
  n_samples <- number_sample
  
  # For each driver gene, sample x0 values *and* remember which sample column they came from
  x0_df <- lapply(rownames(exprs_sub), function(gene) {
    x0 <- exprs_sub[gene, ]
    ok <- !is.na(x0) & x0 > 0
    x0 <- x0[ok]
    idx <- which(ok)
    
    picked <- sample(seq_along(x0), size = n_samples, replace = TRUE)
    data.frame(
      sim       = seq_len(n_samples),
      x0        = as.numeric(x0[picked]) * runif(n_samples, 0.9, 1.1),
      source_id = idx[picked],                  # column index in exprs_sub/exprs_data
      stringsAsFactors = FALSE
    )
  })
  names(x0_df) <- rownames(exprs_sub)
  
  # -----------------------------
  # Step 4: Define CME functions
  # -----------------------------
  sentence4 <- "Define CME functions!" # @@@@@@@@@@@@
  print(sentence4) # @@@@@@@@@@@@
  
  p_xt_analytic <- function(x, t, x0, gamma2, kr, mu){
    a_t <- exp(-gamma2 * t)
    alpha_t <- x0 * a_t
    c <- kr / gamma2
    z_bar <- mu * (1 - a_t)
    x_shift <- x - alpha_t
    x_shift[x_shift <= 0] <- 1e-10
    hyper_vals <- vapply(x_shift, function(z){
      val <- tryCatch(
        hypergeo(c + 1, 2, z * z_bar / a_t),
        error = function(e) 0
      )
      as.numeric(val)
    }, numeric(1))
    term1 <- 1 / x_shift
    term2 <- c * hyper_vals
    p_val <- exp(-mu * x_shift / a_t) / a_t * z_bar * (term1 + term2)
    return(p_val)
  }
  
  simulate_CME <- function(x0, gamma2, kr, mu, t_max, dt = 0.1, method = c("stochastic", "analytic")){
    method <- match.arg(method)
    time <- seq(0, t_max, by = dt)
    x <- numeric(length(time))
    x[1] <- x0
    if(method == "stochastic"){
      for(i in 2:length(time)){
        decay <- gamma2 * x[i-1] * dt
        burst <- rpois(1, lambda = kr * dt) * rexp(1, rate = mu)
        x[i] <- max(x[i-1] - decay + burst, 0)
      }
    } else {
      for(i in 2:length(time)){
        t <- time[i]
        x_vals <- seq(0, max(10*x0, 50), length.out = 200)
        p_vals <- p_xt_analytic(x_vals, t, x0, gamma2, kr, mu)
        p_vals[p_vals < 0] <- 0
        p_vals <- p_vals / sum(p_vals)
        x[i] <- sample(x_vals, size = 1, prob = p_vals)
      }
    }
    data.frame(time = time, x = x)
  }
  
  # -----------------------------
  # Step 5: Estimate CME parameters
  # -----------------------------
  sentence5 <- "Estimate CME parameters!" # @@@@@@@@@@@@
  print(sentence5) # @@@@@@@@@@@@
  
  empirical_stats <- lapply(x0_list, function(vals) {
    c(mean = mean(vals), var = var(vals))
  })
  
  loss_function <- function(params, gene, x0_vec, t_max = 50){
    gamma2 <- params[1]; kr <- params[2]; mu <- params[3]
    sims <- lapply(x0_vec, function(x0) simulate_CME(x0, gamma2, kr, mu, t_max))
    sim_final <- sapply(sims, function(df) tail(df$x, 1))
    sim_stats <- c(mean = mean(sim_final), var = var(sim_final))
    emp_stats <- empirical_stats[[gene]]
    sum((sim_stats - emp_stats)^2)
  }
  
  estimate_parameters <- function(gene, x0_vec){
    init_params <- c(gamma2=0.1, kr=1, mu=1)
    opt <- optim(
      par = init_params, fn = loss_function,
      gene = gene, x0_vec = x0_vec,
      method = "L-BFGS-B",
      lower = c(1e-4,1e-4,1e-4), upper = c(10,50,50)
    )
    return(opt$par)
  }
  
  param_estimates <- lapply(names(x0_list), function(gene){
    estimate_parameters(gene, x0_list[[gene]])
  })
  names(param_estimates) <- names(x0_list)
  
  # -----------------------------
  # Step 6: Re-simulate CME trajectories
  # -----------------------------
  sentence6 <- "Re-simulate CME trajectories!" # @@@@@@@@@@@@
  print(sentence6) # @@@@@@@@@@@@
  
  t_max <- 50
  gene_simulations <- list()
  for(gene in names(x0_list)){
    x0_vec <- x0_list[[gene]]
    est <- param_estimates[[gene]]
    sims <- lapply(x0_vec, function(x0) simulate_CME(x0, est[1], est[2], est[3], t_max))
    sim_df <- do.call(rbind, lapply(seq_along(sims), function(i){
      df <- sims[[i]]; df$sim <- i; df$gene <- gene; df
    }))
    rownames(sim_df) <- NULL
    gene_simulations[[gene]] <- sim_df
  }
  
  # -----------------------------
  # Step 7: Optimal time via tail-variance + mean-separation score
  # -----------------------------
  sentence7 <- "Find the optimal time t_opt per gene!" # @@@@@@@@@@@@
  print(sentence7) # @@@@@@@@@@@@
  
  score_time <- function(df, q_rare = 0.05, w_tail = 0.5, w_mean = 0.5) {
    # df: one gene's simulations (cols: time, x, sim, gene [, source_id])
    by_t <- df %>%
      group_by(time) %>%
      summarize(
        q_low  = quantile(x, q_rare),
        q_high = quantile(x, 1 - q_rare),
        var_tail = {xx <- x[x <= quantile(x, q_rare) | x >= quantile(x, 1 - q_rare)];
        if (length(xx) > 1) var(xx) else 0},
        mean_diff = mean(x[x >= quantile(x, 1 - q_rare)], na.rm = TRUE) -
          mean(x[x <= quantile(x, q_rare)],     na.rm = TRUE),
        .groups = "drop"
      )
    
    # Normalize to [0,1]
    rng_tail <- range(by_t$var_tail, na.rm = TRUE);  den1 <- diff(rng_tail); if (den1 == 0) den1 <- 1
    rng_mean <- range(by_t$mean_diff, na.rm = TRUE); den2 <- diff(rng_mean); if (den2 == 0) den2 <- 1
    
    by_t <- by_t %>%
      mutate(TailVar_n = (var_tail - rng_tail[1]) / den1,
             MeanDiff_n = (mean_diff - rng_mean[1]) / den2,
             Score = w_tail * TailVar_n + w_mean * MeanDiff_n)
    
    t_opt <- by_t$time[ which.max(by_t$Score) ]
    list(t_opt = as.numeric(t_opt), score_table = by_t)
  }
  
  # Compute t_opt for each gene
  t_opt_list <- lapply(names(gene_simulations), function(g) {
    score_time(gene_simulations[[g]])
  })
  names(t_opt_list) <- names(gene_simulations)
  
  # -----------------------------
  # Step 8: Classify subclones at t_opt
  # -----------------------------
  sentence8 <- "Classify low/high/normal subclones at optimal time!" # @@@@@@@@@@@@
  print(sentence8) # @@@@@@@@@@@@
  
  classify_at_topt <- function(df, t_opt) {
    sub <- df %>% dplyr::filter(abs(time - t_opt) < 1e-9)
    mu <- mean(sub$x); sdv <- sd(sub$x)
    sub$type <- ifelse(sub$x <= mu - sdv, "low",
                       ifelse(sub$x >= mu + sdv, "high", "normal"))
    sub
  }
  
  # Ensure we carry source_id from Step 3b
  # Rebuild gene_simulations to include source_id (light patch if missing)
  for (g in names(gene_simulations)) {
    if (!"source_id" %in% names(gene_simulations[[g]])) {
      # Map sim -> source_id from x0_df
      map_tbl <- x0_df[[g]][, c("sim","source_id")]
      gene_simulations[[g]] <- gene_simulations[[g]] %>%
        left_join(map_tbl, by = "sim")
    }
  }
  
  subclone_labels <- lapply(names(gene_simulations), function(g) {
    classify_at_topt(gene_simulations[[g]], t_opt_list[[g]]$t_opt) %>%
      mutate(gene = g)
  })
  names(subclone_labels) <- names(gene_simulations)
  
  # -----------------------------
  # Step 9: Collapse probes -> gene symbols for ALL probes
  # -----------------------------
  suppressPackageStartupMessages({
    library(GEOquery)
    library(AnnotationDbi)
    library(matrixStats)
  })
  
  # --- Decide which statement to use (probe-based ExpressionSet vs ENSEMBL/supplement) ---
  use_stmt1 <- FALSE
  gse_try <- try(GEOquery::getGEO(GEO_number, GSEMatrix = TRUE, AnnotGPL = TRUE), silent = TRUE)
  
  if (!inherits(gse_try, "try-error") && !is.null(gse_try)) {
    esets <- if (is.list(gse_try)) gse_try else list(gse_try)
    dims  <- sapply(esets, function(es) nrow(Biobase::exprs(es)))
    if (any(dims > 0)) {
      gset <- esets[[ which.max(dims) ]]
      ids  <- rownames(Biobase::exprs(gset))
      plat <- annotation(gset)  # e.g. "GPL570"
      
      affy_like <- length(ids) > 0 && any(grepl("_at$", ids, ignore.case = TRUE))
      have_pkg  <- requireNamespace("hgu133plus2.db", quietly = TRUE)
      
      # Use Statement 1 when it looks like an Affy probe matrix and the mapping package is present
      use_stmt1 <- affy_like && have_pkg
      if (use_stmt1) assign(".__gset_stmt1", gset, envir = .GlobalEnv)  # stash for reuse below
    }
  }
  
  if (use_stmt1) {
    # =================================================
    # ================ Statement 1 ====================
    # =================================================
    sentence9 <- "Collapse all probes to gene symbols for association!" # @@@@@@@@@@@@
    print(sentence9)
    
    gset <- get(".__gset_stmt1", envir = .GlobalEnv)  # from detection step
    exprs_data <- Biobase::exprs(gset)
    stopifnot(nrow(exprs_data) > 0)
    platform_id <- annotation(gset)
    message("Loaded platform: ", platform_id, " with ", nrow(exprs_data), " features and ", ncol(exprs_data), " samples.")
    
    # helper
    collapse_by_symbol <- function(mat) {
      if (is.null(mat) || nrow(mat) == 0) {
        return(matrix(numeric(0), nrow = 0, ncol = ncol(mat),
                      dimnames = list(character(0), colnames(mat))))
      }
      split_idx <- split(seq_len(nrow(mat)), rownames(mat))
      out <- lapply(split_idx, function(ix) {
        if (length(ix) == 1L) mat[ix, , drop = FALSE] else {
          med <- matrixStats::colMedians(mat[ix, , drop = FALSE], na.rm = TRUE)
          matrix(med, nrow = 1L, dimnames = list(rownames(mat)[ix][1], colnames(mat)))
        }
      })
      do.call(rbind, out)
    }
    
    ids <- rownames(exprs_data)
    is_probe <- grepl("_at$", ids, ignore.case = TRUE)
    is_ensg  <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids)
    is_sym   <- !is_probe & !is_ensg
    
    # probe -> SYMBOL (Affy GPL570)
    map_probe <- NULL
    if (any(is_probe)) {
      suppressPackageStartupMessages(library(hgu133plus2.db))
      probe_ids <- unique(ids[is_probe])
      map_probe <- AnnotationDbi::select(
        hgu133plus2.db,
        keys    = probe_ids,
        columns = c("PROBEID","SYMBOL"),
        keytype = "PROBEID"
      )
      map_probe <- subset(map_probe, !is.na(PROBEID) & !is.na(SYMBOL))
      if (nrow(map_probe)) {
        map_probe <- map_probe[!duplicated(map_probe$PROBEID), c("PROBEID","SYMBOL")]
        names(map_probe) <- c("orig_id","SYMBOL")
      }
    }
    
    # ENSEMBL(.version) -> SYMBOL
    map_ensg <- NULL
    if (any(is_ensg)) {
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      ensg_full     <- ids[is_ensg]
      ensg_stripped <- sub("\\..*$", "", ensg_full)
      ens_map <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = unique(ensg_stripped),
        columns = c("ENSEMBL","SYMBOL"),
        keytype = "ENSEMBL"
      )
      ens_map <- subset(ens_map, !is.na(ENSEMBL) & !is.na(SYMBOL))
      if (nrow(ens_map)) {
        map_ensg <- merge(
          data.frame(orig_id = ensg_full, ENSEMBL = ensg_stripped, stringsAsFactors = FALSE),
          unique(ens_map[, c("ENSEMBL","SYMBOL")]),
          by = "ENSEMBL",
          all.x = TRUE
        )
        map_ensg <- subset(map_ensg, !is.na(SYMBOL))
        map_ensg <- unique(map_ensg[, c("orig_id","SYMBOL")])
      }
    }
    
    # already SYMBOL-ish
    map_sym <- if (any(is_sym)) data.frame(orig_id = ids[is_sym], SYMBOL = ids[is_sym], stringsAsFactors = FALSE) else NULL
    
    # combine & keep all features (fallback: keep original ID)
    all_map <- do.call(rbind, Filter(NROW, list(map_probe, map_ensg, map_sym)))
    all_map <- all_map[!duplicated(all_map$orig_id), , drop = FALSE]
    unmapped <- setdiff(ids, all_map$orig_id)
    if (length(unmapped)) {
      message("Unmapped feature IDs retained as-is: ", length(unmapped))
      all_map <- rbind(all_map,
                       data.frame(orig_id = unmapped, SYMBOL = unmapped, stringsAsFactors = FALSE))
    }
    all_map <- all_map[match(ids, all_map$orig_id), , drop = FALSE]
    
    exprs_all <- exprs_data[all_map$orig_id, , drop = FALSE]
    rownames(exprs_all) <- all_map$SYMBOL
    
    # ---- OUTPUT 1 ----
    gene_mat       <- collapse_by_symbol(exprs_all)
    non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
    
  } else {
    # =================================================
    # ================ Statement 2 ====================
    # =================================================
    build_symbol_matrix_from_GEO <- function(GEO_number,
                                             genes_of_interest = character(0),
                                             prefer_supp = TRUE,
                                             keep_unmapped = TRUE,
                                             verbose = TRUE) {
      # deps (load quietly)
      suppressPackageStartupMessages({
        library(GEOquery)
        library(AnnotationDbi)
        library(org.Hs.eg.db)
        library(data.table)
        library(matrixStats)
        library(R.utils)
        library(Matrix)
      })
      
      # helpers (scoped inside so you don't clutter your workspace)
      collapse_by_symbol <- function(mat) {
        if (is.null(mat) || nrow(mat) == 0) {
          return(matrix(numeric(0), nrow = 0, ncol = ncol(mat),
                        dimnames = list(character(0), colnames(mat))))
        }
        split_idx <- split(seq_len(nrow(mat)), rownames(mat))
        out <- lapply(split_idx, function(ix) {
          if (length(ix) == 1L) mat[ix, , drop = FALSE] else {
            med <- matrixStats::colMedians(mat[ix, , drop = FALSE], na.rm = TRUE)
            matrix(med, nrow = 1L, dimnames = list(rownames(mat)[ix][1], colnames(mat)))
          }
        })
        do.call(rbind, out)
      }
      guess_id_column <- function(df) {
        candidates <- c("UID","Gene","Geneid","GENE","symbol","SYMBOL","gene_id","gene_name","feature","id","ID")
        cand <- intersect(candidates, names(df))
        if (length(cand)) return(cand[1])
        ensg_re <- "^ENSG[0-9]+(\\.[0-9]+)?$"
        scores <- vapply(names(df), function(nm) mean(grepl(ensg_re, as.character(df[[nm]])), na.rm = TRUE), numeric(1))
        if (any(scores > 0.5)) return(names(df)[which.max(scores)])
        nn <- names(df)[!vapply(df, is.numeric, logical(1))]
        if (length(nn)) return(nn[1])
        names(df)[1]
      }
      map_ensembl_version_to_symbol <- function(ens_ids, keep_unmapped = TRUE) {
        stripped <- sub("\\..*$", "", ens_ids)
        look <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys   = unique(stripped),
                                      column = "SYMBOL",
                                      keytype= "ENSEMBL",
                                      multiVals = "first")
        sym <- unname(look[stripped])
        if (keep_unmapped) {
          na <- is.na(sym) | sym == ""
          sym[na] <- stripped[na]
        }
        sym
      }
      
      # 1) Try ExpressionSet first
      gse <- tryCatch(getGEO(GEO_number, GSEMatrix = TRUE, AnnotGPL = TRUE), error = function(e) NULL)
      if (!is.null(gse)) {
        esets <- if (is.list(gse)) gse else list(gse)
        dims  <- sapply(esets, function(es) nrow(exprs(es)))
        if (any(dims > 0)) {
          gset <- esets[[ which.max(dims) ]]
          exprs_data <- Biobase::exprs(gset)
          ids <- rownames(exprs_data)
          if (verbose) message("Loaded ExpressionSet: ", nrow(exprs_data), " features × ", ncol(exprs_data), " samples.")
          is_ensg <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids)
          SYMBOL  <- if (any(is_ensg)) map_ensembl_version_to_symbol(ids, keep_unmapped) else ids
          rownames(exprs_data) <- SYMBOL
          gene_mat <- collapse_by_symbol(exprs_data)
          non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
          return(list(
            source            = "GSEMatrix",
            data_with_symbols = NULL,
            gene_mat          = gene_mat,
            non_driver_mat    = non_driver_mat,
            mapping_table     = data.frame(orig_id = ids, SYMBOL = SYMBOL, stringsAsFactors = FALSE)
          ))
        }
      }
      if (!prefer_supp) stop("No non-empty ExpressionSet found and prefer_supp = FALSE.")
      
      # 2) Supplementary files
      base <- file.path(tempdir(), GEO_number)
      dir.create(base, showWarnings = FALSE, recursive = TRUE)
      GEOquery::getGEOSuppFiles(GEO_number, baseDir = base, makeDirectory = TRUE)
      supp_dir <- file.path(base, GEO_number)
      files <- list.files(supp_dir, full.names = TRUE)
      if (!length(files)) stop("No supplementary files found for ", GEO_number)
      
      gz <- grep("\\.gz$", files, value = TRUE)
      for (z in gz) try(R.utils::gunzip(z, overwrite = TRUE), silent = TRUE)
      files <- list.files(supp_dir, full.names = TRUE)
      
      cand <- grep("count|counts|matrix|tpm|fpkm|rpkm|expr|expression|rsem|genes?\\.tsv|\\.txt$|\\.csv$|matrix\\.mtx",
                   files, ignore.case = TRUE, value = TRUE)
      if (!length(cand)) stop("Could not detect an expression-like file in supplements.")
      f <- cand[ which.max(file.info(cand)$size) ]
      if (verbose) message("Reading supplementary file: ", basename(f))
      
      # 10X mtx path
      if (grepl("matrix\\.mtx", basename(f), ignore.case = TRUE)) {
        mm <- Matrix::readMM(f)
        bar <- grep("barcode", files, ignore.case = TRUE, value = TRUE)[1]
        feat<- grep("features|genes", files, ignore.case = TRUE, value = TRUE)[1]
        mat <- as.matrix(mm)
        if (!is.na(bar)) colnames(mat) <- data.table::fread(bar, header = FALSE)[[1]]
        if (!is.na(feat)) rownames(mat) <- data.table::fread(feat, header = FALSE)[[1]]
        ids <- rownames(mat)
        SYMBOL <- map_ensembl_version_to_symbol(ids, keep_unmapped)
        rownames(mat) <- SYMBOL
        gene_mat <- collapse_by_symbol(mat)
        non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
        return(list(
          source            = "supplement_10X",
          data_with_symbols = NULL,
          gene_mat          = gene_mat,
          non_driver_mat    = non_driver_mat,
          mapping_table     = data.frame(orig_id = ids,
                                         ENSEMBL_nover = sub("\\..*$","", ids),
                                         SYMBOL = SYMBOL, stringsAsFactors = FALSE)
        ))
      }
      
      # Delimited table path
      raw_dt <- data.table::fread(f)
      raw_df <- as.data.frame(raw_dt)  # avoid DT column-selection pitfalls
      id_col <- guess_id_column(raw_df)
      ids    <- as.character(raw_df[[id_col]])
      
      # numeric sample columns (coerce if needed)
      num_mask <- vapply(raw_df, is.numeric, logical(1))
      num_mask[ match(id_col, names(raw_df)) ] <- FALSE
      if (!any(num_mask)) {
        non_id <- setdiff(names(raw_df), id_col)
        for (nm in non_id) raw_df[[nm]] <- suppressWarnings(as.numeric(raw_df[[nm]]))
        num_mask <- vapply(raw_df, is.numeric, logical(1))
        num_mask[ match(id_col, names(raw_df)) ] <- FALSE
      }
      if (!any(num_mask)) stop("No numeric sample columns detected in supplementary table.")
      
      num_names <- names(raw_df)[num_mask]
      exprs_all <- as.matrix(raw_df[, num_names, drop = FALSE])
      colnames(exprs_all) <- make.unique(colnames(exprs_all), sep = "_")
      
      SYMBOL <- map_ensembl_version_to_symbol(ids, keep_unmapped)
      rownames(exprs_all) <- SYMBOL
      
      gene_mat <- collapse_by_symbol(exprs_all)
      non_driver_mat <- gene_mat[setdiff(rownames(gene_mat), genes_of_interest), , drop = FALSE]
      
      data_with_symbols <- raw_df
      data_with_symbols$GeneSymbol <- SYMBOL
      ord <- c(id_col, "GeneSymbol", setdiff(names(data_with_symbols), c(id_col,"GeneSymbol")))
      data_with_symbols <- data_with_symbols[, ord]
      
      list(
        source            = "supplement_table",
        data_with_symbols = data_with_symbols,
        gene_mat          = gene_mat,
        non_driver_mat    = non_driver_mat,
        mapping_table     = data.frame(orig_id = ids,
                                       ENSEMBL_nover = sub("\\..*$","", ids),
                                       SYMBOL = SYMBOL, stringsAsFactors = FALSE)
      )
    }
    if (!exists("build_symbol_matrix_from_GEO")) {
      stop("build_symbol_matrix_from_GEO() is not defined. Source your Statement 2 loader first.")
    }
    
    sentence9 <- "Collapse all probes to gene symbols for association!" # @@@@@@@@@@@@
    print(sentence9)
    
    res <- build_symbol_matrix_from_GEO(
      GEO_number        = GEO_number,
      genes_of_interest = genes_of_interest,
      prefer_supp       = TRUE,
      keep_unmapped     = TRUE,
      verbose           = TRUE
    )
    
    # ---- OUTPUT 2 ----
    gene_mat       <- res$gene_mat
    non_driver_mat <- res$non_driver_mat
  }
  
  # -----------------------------
  # Step 10: Gene-wise association (two-sample t-test + fold-change)
  # -----------------------------
  sentence10 <- "Run gene-wise association for each driver gene!" # @@@@@@@@@@@@
  print(sentence10) # @@@@@@@@@@@@
  
  assoc_for_driver <- function(driver_gene, labels_df, non_driver_mat) {
    # labels_df: rows for t_opt with columns sim, x, type, source_id, gene
    R_ids <- labels_df$source_id[labels_df$type %in% c("low","high")]
    O_ids <- setdiff(seq_len(ncol(non_driver_mat)), R_ids)
    
    # Means by group
    mean_R <- rowMeans(non_driver_mat[, R_ids, drop = FALSE])
    mean_O <- rowMeans(non_driver_mat[, O_ids, drop = FALSE])
    FC <- mean_R / mean_O
    
    # Two-sample t-test per gene (vectorized over rows)
    t_pvals <- apply(non_driver_mat, 1, function(v) {
      x <- v[R_ids]; y <- v[O_ids]
      if (length(unique(c(x, y))) < 2) return(1.0)
      t.test(x, y, var.equal = FALSE)$p.value
    })
    padj <- p.adjust(t_pvals, method = "BH")
    
    out <- data.frame(
      driver = driver_gene,
      gene   = rownames(non_driver_mat),
      mean_R = mean_R,
      mean_O = mean_O,
      FC     = FC,
      pval   = t_pvals,
      padj   = padj,
      status = ifelse(padj < 0.05 & FC >= 1.2, "enriched",
                      ifelse(padj < 0.05 & FC <= 0.83, "depleted", "ns")),
      stringsAsFactors = FALSE
    )
    out[order(out$padj, -abs(log2(out$FC))), ]
  }
  
  assoc_results <- lapply(names(subclone_labels), function(g) {
    assoc_for_driver(g, subclone_labels[[g]], non_driver_mat)
  })
  names(assoc_results) <- names(subclone_labels)
  
  # Convenience: bind all drivers together
  assoc_all <- dplyr::bind_rows(assoc_results)
  
  # -----------------------------
  # Step 11 (robust): Filter subclonal genes by CME ribbon using new variables
  # -----------------------------
  sentence11 <- "Filter subclonal genes by CME ribbon (robust)!"  # @@@@@@@@@@@@
  print(sentence11)  # @@@@@@@@@@@@
  
  #@suppressPackageStartupMessages({
  #@  library(dplyr)
  #@})
  
  filtered_subclonal_genes_list <- list()   # driver -> character vector of SYMBOLs
  filtered_subclonal_genes_df   <- list()   # driver -> data.frame with stats
  
  for (driver in genes_of_interest) {
    
    # --- Guardrails ---
    if (!driver %in% names(gene_simulations)) {
      warning("No simulations found for driver: ", driver); next
    }
    if (!driver %in% names(t_opt_list)) {
      warning("No t_opt entry found for driver: ", driver); next
    }
    if (!driver %in% names(subclone_labels)) {
      warning("No subclone labels found for driver: ", driver); next
    }
    if (!driver %in% names(assoc_results)) {
      warning("No association results found for driver: ", driver); next
    }
    
    # --- 11.1 CME ribbon at t_opt for this driver (force numeric 'time') ---
    df_sim <- gene_simulations[[driver]]
    band_df <- df_sim %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(
        lower = quantile(x, 0.025, na.rm = TRUE),
        upper = quantile(x, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(time = as.numeric(time))
    
    t_opt_num <- as.numeric(t_opt_list[[driver]]$t_opt)
    idx <- which.min(abs(band_df$time - t_opt_num))
    # base subset to avoid slice()/S4 method dispatch
    lower <- band_df$lower[idx]
    upper <- band_df$upper[idx]
    
    # --- 11.2 Rare sample IDs (source indices) at t_opt from Step 8 ---
    labs <- subclone_labels[[driver]]
    # keep only exact t_opt row if present; otherwise it already contained only t_opt in earlier code
    if ("time" %in% names(labs)) {
      labs <- labs[which.min(abs(as.numeric(labs$time) - t_opt_num)) == 1, , drop = FALSE]
    }
    R_ids <- unique(as.integer(labs$source_id[labs$type %in% c("low", "high")]))
    R_ids <- R_ids[!is.na(R_ids)]
    if (length(R_ids) == 0) {
      warning("No rare sample IDs for driver: ", driver)
      filtered_subclonal_genes_list[[driver]] <- character(0)
      next
    }
    
    # --- 11.3 Candidate subclonal genes from associations (SYMBOL level) ---
    cand_tbl <- assoc_results[[driver]] %>%
      dplyr::mutate(gene = as.character(gene)) %>%
      dplyr::filter(status %in% c("enriched", "depleted"))
    
    if (nrow(cand_tbl) == 0) {
      filtered_subclonal_genes_list[[driver]] <- character(0)
      next
    }
    
    # --- 11.4 Keep genes whose expression in rare samples lies within driver's ribbon ---
    keep_tbl <- cand_tbl %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        within_ribbon = {
          gsym <- as.character(gene)
          if (!gsym %in% rownames(gene_mat)) {
            NA
          } else {
            vals <- as.numeric(gene_mat[gsym, R_ids, drop = TRUE])
            vals <- vals[is.finite(vals)]
            length(vals) > 0 && all(vals >= lower & vals <= upper)
          }
        }
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(isTRUE(within_ribbon))
    
    filtered_subclonal_genes_list[[driver]] <- keep_tbl$gene
    
    # tidy per-driver table (even if empty, keep consistent columns)
    keep_out <- keep_tbl %>%
      dplyr::mutate(
        driver       = driver,
        ribbon_lower = lower,
        ribbon_upper = upper,
        t_opt        = t_opt_num,
        rare_n       = length(R_ids)
      ) %>%
      dplyr::relocate(driver, gene)
    
    filtered_subclonal_genes_df[[driver]] <- keep_out
  }
  
  # Bind per-driver data frames into one tidy table (may be empty)
  filtered_subclonal_genes_all <- dplyr::bind_rows(filtered_subclonal_genes_df)
  
  # Optional: print a compact summary
  message("Kept subclonal genes per driver:")
  #@@@print(lapply(filtered_subclonal_genes_list, length))
  
  filtered_subclonal_genes <- lapply(filtered_subclonal_genes_list, length)
  #@@@print(filtered_subclonal_genes)
  
  # ===============================
  # Step 13 Dynamic calibration for ALL genes + faceted plots
  # ===============================
  sentence13 <- "Dynamic calibration via split simulations (all genes)!" # @@@@@@@@@@@@
  print(sentence13) # @@@@@@@@@@@@
  
  #@  suppressPackageStartupMessages({
  #@   library(dplyr)
  #@   library(tidyr)
  #@   library(ggplot2)
  #@  library(scoringRules)
  #@ })
  
  # ---- Core routine (returns per-time coverage, tails, CRPS + CI band) ----
  dynamic_calibration <- function(theta_vec, x0_vec, times,
                                  S_pred = 200, S_pseudo = 200, dt = diff(times)[1],
                                  seed = 42) {
    set.seed(seed)
    
    # 1) Predictive draws (used to form the band)
    pred_sims <- lapply(seq_len(S_pred), function(i) {
      x0 <- sample(x0_vec, 1, replace = TRUE)
      simulate_CME(x0, theta_vec["gamma2"], theta_vec["kr"], theta_vec["mu"],
                   t_max = max(times), dt = dt, method = "stochastic")
    })
    pred_df <- dplyr::bind_rows(Map(function(df, i) dplyr::mutate(df, sim = i), pred_sims, seq_along(pred_sims)))
    
    # 2) Independent pseudo-observations
    obs_sims <- lapply(seq_len(S_pseudo), function(i) {
      x0 <- sample(x0_vec, 1, replace = TRUE)
      simulate_CME(x0, theta_vec["gamma2"], theta_vec["kr"], theta_vec["mu"],
                   t_max = max(times), dt = dt, method = "stochastic")
    })
    obs_df <- dplyr::bind_rows(Map(function(df, i) dplyr::mutate(df, sim = i), obs_sims, seq_along(obs_sims)))
    
    # 3) Per-time predictive band and draws list
    by_t_pred <- pred_df %>%
      group_by(time) %>%
      summarise(q_low = quantile(x, 0.025),
                q_hi  = quantile(x, 0.975),
                pred_draws = list(x),
                .groups = "drop")
    
    # 4) Pseudo-observations grouped per time
    obs_by_t <- obs_df %>% group_by(time) %>%
      summarise(x_obs = list(x), n_obs = dplyr::n(), .groups = "drop")
    
    # 5) Coverage + CRPS by time
    res <- by_t_pred %>%
      left_join(obs_by_t, by = "time") %>%
      rowwise() %>%
      mutate(
        coverage = mean(unlist(x_obs) >= q_low & unlist(x_obs) <= q_hi),
        tail_low = mean(unlist(x_obs) < q_low),
        tail_hi  = mean(unlist(x_obs) > q_hi),
        crps     = {
          pd <- unlist(pred_draws)
          y  <- unlist(x_obs)
          scoringRules::crps_sample(y = y,
                                    dat = matrix(rep(pd, each = length(y)), nrow = length(y))
          ) %>% mean()
        },
        # Binomial 95% band around nominal 0.95 for visual reference
        se_nom   = sqrt(0.95 * 0.05 / n_obs),
        cov_lo   = 0.95 - 1.96 * se_nom,
        cov_hi   = 0.95 + 1.96 * se_nom
      ) %>%
      ungroup()
    
    res
  }
  
  # ---- Run for all genes and build faceted plots ----
  times_vec <- sort(unique(gene_simulations[[1]]$time))
  genes_run <- intersect(genes_of_interest, names(param_estimates))
  
  # Collect results in one data frame with a 'gene' column
  cal_all <- dplyr::bind_rows(lapply(genes_run, function(g) {
    theta <- setNames(param_estimates[[g]], c("gamma2","kr","mu"))
    out <- dynamic_calibration(theta_vec = theta,
                               x0_vec    = x0_df[[g]]$x0,
                               times     = times_vec,
                               S_pred    = 200,
                               S_pseudo  = 200,
                               dt        = diff(times_vec)[1],
                               seed      = 100 + which(genes_run == g))
    dplyr::mutate(out, gene = g)
  }))
  
  # ---- Faceted coverage plot (with nominal 95% band ribbon) ----
  p_cov <- ggplot(cal_all, aes(time, coverage)) +
    geom_ribbon(aes(ymin = cov_lo, ymax = cov_hi), fill = "#DDAA33", alpha = 0.4) +
    geom_hline(yintercept = 0.95, linetype = 2, color = "#8B0000", linewidth = 1.5) +
    geom_line(color = "#2E8B57") +
    facet_wrap(~ gene, ncol = gene_column_number) +
    labs(title = "Dynamic coverage across genes (95% predictive band)",
         #@@ subtitle = "Shaded ribbon = binomial 95% band around nominal 0.95",
         x = "Time", y = "Coverage") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color="black"),
          panel.border = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold")
    )
  
  #@@ print(p_cov)
  
  # ---- Faceted CRPS plot ----
  p_crps <- ggplot(cal_all, aes(time, crps)) +
    geom_line(color = "#2E8B57") +
    facet_wrap(~ gene, ncol = gene_column_number, scales = "free_y") +
    labs(title = "CRPS over time across genes",
         x = "Time", y = "CRPS (lower is better)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color="black"),
          panel.border = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold")
    )
  
  #@@@print(p_crps)
  
  
  # -----------------------------
  # Step 14: Biological plausibility checks (combined figure)
  # -----------------------------
  sentence14 <- "Report t1/2 and steady-state moments; flag outliers + visualize!"  # @@@@@@@@@@@@
  print(sentence14)  # @@@@@@@@@@@@
  
  #@  suppressPackageStartupMessages({
  #@    library(dplyr)
  #@    library(tidyr)
  #@    library(ggplot2)
  #@    library(scales)
  #@    library(patchwork)   # for combining plots
  #@  })
  
  # ---- sizing controls ----
  pt_size   <- 5.8   # point size for genes (was 2)
  line_w    <- 4.2   # general line width (was 0.6–0.7)
  seg_w     <- 3.3   # lollipop stem width
  hline_w   <- 3.1   # dashed threshold lines
  text_size <- 4.2   # label size for ggrepel (if used)
  
  # ---- 14.1 Build plausibility table from fitted parameters -------------------
  half_life_bounds <- c(0.05, 200)
  
  bio_df <- dplyr::bind_rows(lapply(names(param_estimates), function(g) {
    th <- setNames(param_estimates[[g]], c("gamma2","kr","mu"))
    data.frame(
      gene    = g,
      t_half  = log(2) / th["gamma2"],
      mean_ss = (th["kr"] * th["mu"]) / th["gamma2"],
      var_ss  = (th["kr"] * (th["mu"]^2)) / th["gamma2"],  # Exp bursts
      stringsAsFactors = FALSE
    )
  }))
  
  bio_df <- bio_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      t_half  = stats::median(t_half,  na.rm = TRUE),
      mean_ss = stats::median(mean_ss, na.rm = TRUE),
      var_ss  = stats::median(var_ss,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      flag_half = t_half < half_life_bounds[1] | t_half > half_life_bounds[2],
      flag_mean = !is.finite(mean_ss) | mean_ss <= 0,
      flag_any  = flag_half | flag_mean
    )
  
  bio_summary <- bio_df %>%
    dplyr::summarise(
      n_genes     = dplyr::n(),
      n_flagged   = sum(flag_any),
      med_t_half  = stats::median(t_half),
      iqr_t_half  = IQR(t_half),
      med_mean_ss = stats::median(mean_ss),
      med_var_ss  = stats::median(var_ss)
    )
  #@@@print(bio_summary)
  #@@@print(bio_df)
  
  # ---- 14.2 Helpers -----------------------------------------------------------
  si_lab <- function(...) scales::label_number(scale_cut = scales::cut_si(""), ...)
  
  phi_hat <- stats::median(pmax(bio_df$var_ss, 1e-12) / pmax(bio_df$mean_ss, 1e-12)^2, na.rm = TRUE)
  m_rng   <- range(pmax(bio_df$mean_ss, 1e-12), na.rm = TRUE)
  x_line  <- exp(seq(log(m_rng[1]), log(m_rng[2]), length.out = 200))
  ref_lines <- tibble::tibble(
    mean_ss = x_line,
    poisson = x_line,
    nb_like = x_line + phi_hat * x_line^2
  ) %>%
    tidyr::pivot_longer(c(poisson, nb_like), names_to = "model", values_to = "var_ref")
  
  # Common theme tweak
  base_theme <- theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color="black", linewidth = 0.4),
          panel.border = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  # ---- 14.3 (A) Half-life vs steady-state mean (log–log) ----------------------
  p14_scatter <- ggplot(
    bio_df %>% dplyr::mutate(mean_ss_p = pmax(mean_ss, 1e-12), t_half_p = pmax(t_half, 1e-12)),
    aes(mean_ss_p, t_half_p, color = flag_any)
  ) +
    geom_point(size = pt_size, alpha = 0.9) +
    geom_hline(yintercept = half_life_bounds, linetype = 2, color = "#A52A2A", linewidth = hline_w) +
    scale_x_log10(labels = si_lab()) +
    scale_y_log10(labels = si_lab()) +
    scale_color_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                       labels = c(`FALSE` = "OK", `TRUE` = "Flagged")) +
    labs(title = "Half-life vs steady-state mean",
         x = "Steady-state mean (m∞)", y = "Half-life t½", color = "Plausibility") +
    base_theme
  
  if (any(bio_df$flag_any) && requireNamespace("ggrepel", quietly = TRUE)) {
    p14_scatter <- p14_scatter +
      ggrepel::geom_text_repel(
        data = dplyr::filter(bio_df, flag_any),
        aes(label = gene), size = text_size, min.segment.length = 0,
        max.overlaps = Inf, show.legend = FALSE
      )
  }
  
  # ---- 14.4 (B) Variance–mean with reference curves ---------------------------
  p14_varmean <- ggplot() +
    geom_point(data = bio_df %>%
                 dplyr::mutate(mean_ss_p = pmax(mean_ss, 1e-12),
                               var_ss_p  = pmax(var_ss,  1e-12)),
               aes(mean_ss_p, var_ss_p, color = flag_any),
               size = pt_size, alpha = 0.9) +
    geom_line(data = ref_lines, aes(mean_ss, var_ref, linetype = model),
              linewidth = line_w, color = "#A52A2A") +
    scale_x_log10(labels = si_lab()) +
    scale_y_log10(labels = si_lab()) +
    scale_color_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                       labels = c(`FALSE` = "OK", `TRUE` = "Flagged")) +
    scale_linetype_manual(values = c(poisson = "solid", nb_like = "dashed"),
                          labels = c(poisson = "Poisson: Var = Mean",
                                     nb_like  = paste0("NB-like: Var = Mean + ",
                                                       signif(phi_hat, 2), "·Mean²"))) +
    labs(title = "Variance–mean with reference curves",
         x = "Steady-state mean (m∞)", y = "Steady-state variance (v∞)",
         color = "Plausibility", linetype = "Reference") +
    base_theme
  
  # ---- 14.5 (C) Half-life lollipop by gene -----------------------------------
  bio_df_ord <- bio_df %>%
    dplyr::arrange(dplyr::desc(t_half)) %>%
    dplyr::mutate(gene_f = factor(gene, levels = rev(gene)))
  
  p14_lollipop <- ggplot(
    bio_df_ord %>% dplyr::mutate(t_half_p = pmax(t_half, 1e-12)),
    aes(gene_f, t_half_p, color = flag_any)
  ) +
    geom_segment(aes(xend = gene_f, y = min(t_half_p)/2, yend = t_half_p),
                 linewidth = seg_w, alpha = 0.7) +
    geom_point(size = pt_size) +
    geom_hline(yintercept = half_life_bounds, linetype = 2, color = "#A52A2A", linewidth = hline_w) +
    scale_y_log10(labels = si_lab()) +
    scale_color_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                       labels = c(`FALSE` = "OK", `TRUE` = "Flagged")) +
    coord_flip() +
    labs(title = "Transcript half-life (log scale)",
         x = NULL, y = "t½ (log scale)", color = "Plausibility") +
    base_theme
  
  # ---- 14.6 Combine into one figure ------------------------------------------
  # Layout: [A | B] over [C]
  #@@combined_p14 <- (p14_scatter | p14_varmean) / p14_lollipop +
  #@@  plot_layout(guides = "collect", heights = c(1, 1)) &
  #@@  theme(legend.position = "right")
  
  #@@print(combined_p14)
  
  # Optional: save combined figure
  # ggsave("step14_biological_plausibility_combined.png", combined_p14, width = 12, height = 9, dpi = 300)
  
  
  # ggsave("step14_half_life_vs_mean.pdf", p14_scatter, width = 7, height = 5)
  # ggsave("step14_var_vs_mean.pdf",     p14_varmean, width = 7, height = 5)
  # ggsave("step14_lollipop_half_life.pdf", p14_lollipop, width = 8, height = 6)
  
  # ===============================
  # Faceted CME overlay visualization (ggplot)
  # ===============================
  #@ library(dplyr)
  #@ library(ggplot2)
  
  plot_cme_overlays_faceted <- function(gene_sims,
                                        subclone_labs,
                                        t_opt_list,
                                        genes,
                                        alpha = 0.05,
                                        ncol = 3) {
    # keep genes we actually have
    genes <- intersect(genes, intersect(names(gene_sims), names(subclone_labs)))
    stopifnot(length(genes) > 0)
    
    # 95% band + mean per time, per gene
    bands <- dplyr::bind_rows(lapply(genes, function(g) {
      gene_sims[[g]] %>%
        group_by(time) %>%
        summarise(
          mean = mean(x),
          low  = quantile(x, alpha/2),
          hi   = quantile(x, 1 - alpha/2),
          .groups = "drop"
        ) %>%
        mutate(gene = g)
    }))
    
    # points at t_opt (already labeled low/high/normal in subclone_labs)
    pts <- dplyr::bind_rows(lapply(genes, function(g) {
      df <- subclone_labs[[g]]
      # ensure 'gene' column exists
      if (!"gene" %in% names(df)) df$gene <- g
      df
    }))
    
    # vertical line at t_opt per facet
    vlines <- tibble::tibble(
      gene  = genes,
      t_opt = vapply(genes, function(g) t_opt_list[[g]]$t_opt, numeric(1))
    )
    
    ggplot() +
      geom_ribbon(data = bands, aes(x = time, ymin = low, ymax = hi),
                  alpha = 0.25, fill = "#DDAA33") +
      geom_line(data = bands, aes(x = time, y = mean),
                color = "#444444", linewidth = 0.6) +
      geom_vline(data = vlines, aes(xintercept = t_opt),
                 linetype = 2, color = "grey40") +
      geom_point(data = pts, aes(x = time, y = x, color = type),
                 alpha = 0.7, size = 2.4) +
      scale_color_manual(values = c(low = "blue", normal = "grey60", high = "red")) +
      facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
      labs(
        title = "CME envelopes with subclone labels at t_opt",
        x = "Time", y = "Expression", color = "Subclone"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color="black"),
            panel.border = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(face = "bold")
      )
  }
  
  # ---------- Usage ----------
  # Make sure genes_of_interest exists (as in your script)
  cme_overlays_faceted <- plot_cme_overlays_faceted(
    gene_sims      = gene_simulations,
    subclone_labs  = subclone_labels,
    t_opt_list     = t_opt_list,
    genes          = genes_of_interest,  # or specify a subset, e.g., c("TP53","BRCA1",...)
    alpha          = 0.05,
    ncol           = gene_column_number
  )
  # save to file
  ggsave("cme_envelopes_subclones_faceted.pdf", cme_overlays_faceted, width = 10, height = 7)
  
  #@@@print(cme_overlays_faceted)
  #@@@@@@@@@@@@@@@
  
  plot_cme <- function(gene_sims,
                       subclone_labs,
                       t_opt_list,
                       genes,
                       alpha = 0.05,
                       ncol = 3) {
    # keep genes we actually have
    genes <- intersect(genes, intersect(names(gene_sims), names(subclone_labs)))
    stopifnot(length(genes) > 0)
    
    # 95% band + mean per time, per gene
    bands <- dplyr::bind_rows(lapply(genes, function(g) {
      gene_sims[[g]] %>%
        group_by(time) %>%
        summarise(
          mean = mean(x),
          low  = quantile(x, alpha/2),
          hi   = quantile(x, 1 - alpha/2),
          .groups = "drop"
        ) %>%
        mutate(gene = g)
    }))
    
    # points at t_opt (already labeled low/high/normal in subclone_labs)
    pts <- dplyr::bind_rows(lapply(genes, function(g) {
      df <- subclone_labs[[g]]
      # ensure 'gene' column exists
      if (!"gene" %in% names(df)) df$gene <- g
      df
    }))
    
    # vertical line at t_opt per facet
    vlines <- tibble::tibble(
      gene  = genes,
      t_opt = vapply(genes, function(g) t_opt_list[[g]]$t_opt, numeric(1))
    )
    
    ggplot() +
      geom_ribbon(data = bands, aes(x = time, ymin = low, ymax = hi),
                  alpha = 0.25, fill = "#DDAA33") +
      #@@geom_line(data = bands, aes(x = time, y = mean),
      #@@  color = "#444444", linewidth = 0.6) +
      #@@ geom_vline(data = vlines, aes(xintercept = t_opt),
      #@@           linetype = 2, color = "grey40") +
      #@@geom_point(data = pts, aes(x = time, y = x, color = type),
      #@@       alpha = 0.7, size = 2.4) +
      #@@scale_color_manual(values = c(low = "blue", normal = "grey60", high = "red")) +
      facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
      labs(
        title = "CME envelopes with subclone labels at t_opt",
        x = "Time", y = "Expression", color = "Subclone"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color="black"),
            panel.border = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(face = "bold")
      )
  }
  
  # Make sure genes_of_interest exists (as in your script)
  cme_dynamics <- plot_cme(
    gene_sims      = gene_simulations,
    subclone_labs  = subclone_labels,
    t_opt_list     = t_opt_list,
    genes          = genes_of_interest,  # or specify a subset, e.g., c("TP53","BRCA1",...)
    alpha          = 0.05,
    ncol           = gene_column_number
  )
  #@@print(cme_dynamics)
  
  # ===============================
  # Filtered subclonal genes all
  # ===============================
  sentence15 <- "Filtered subclonal genes all!"  # @@@@@@@@@@@@
  print(sentence15)  # @
  
  #@ suppressPackageStartupMessages({
  #@  library(dplyr)
  #@  library(tidyr)
  #@  library(ggplot2)
  #@  library(stringr)
  #@})
  
  plot_cme_overlays_with_filter <- function(gene_sims,
                                            subclone_labs,
                                            t_opt_list,
                                            genes,
                                            filtered_df,            # filtered_subclonal_genes_all from Step 11
                                            alpha = 0.05,
                                            ncol = 3,
                                            max_labels = 6,
                                            rect_width = NULL,
                                            mean_linewidth = 1.0,
                                            vline_width = 0.9,
                                            point_size = 3.2,
                                            rect_alpha = 0.15) {
    
    # keep genes we actually have in all inputs
    genes <- intersect(genes, Reduce(intersect, list(names(gene_sims), names(subclone_labs), names(t_opt_list))))
    stopifnot(length(genes) > 0)
    
    # --------- 1) Standardize/repair filtered_df -------------------------------
    if (is.null(filtered_df) || !is.data.frame(filtered_df) || nrow(filtered_df) == 0) {
      filtered_df <- data.frame(driver = character(), gene = character(), status = character(), stringsAsFactors = FALSE)
    } else {
      # normalize driver column name
      if (!"driver" %in% names(filtered_df)) {
        d_alt <- intersect(c("Driver","DRIVER","driver_gene"), names(filtered_df))
        if (length(d_alt)) names(filtered_df)[names(filtered_df) == d_alt[1]] <- "driver"
      }
      # normalize gene column name
      if (!"gene" %in% names(filtered_df)) {
        g_alt <- intersect(c("GENE","SYMBOL","symbol","Gene"), names(filtered_df))
        if (length(g_alt)) names(filtered_df)[names(filtered_df) == g_alt[1]] <- "gene"
      }
      # create status if absent
      if (!"status" %in% names(filtered_df)) filtered_df$status <- NA_character_
      # cast
      filtered_df$driver <- as.character(filtered_df$driver)
      filtered_df$gene   <- as.character(filtered_df$gene)
      filtered_df$status <- as.character(filtered_df$status)
    }
    
    # --------- 2) 95% band + mean per time, per gene ---------------------------
    bands <- dplyr::bind_rows(lapply(genes, function(g) {
      df <- gene_sims[[g]]
      df$time <- as.numeric(df$time)
      df %>%
        dplyr::group_by(time) %>%
        dplyr::summarise(
          mean = mean(x),
          low  = stats::quantile(x, alpha/2, na.rm = TRUE),
          hi   = stats::quantile(x, 1 - alpha/2, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(gene = g)
    }))
    
    # --------- 3) subclone points (from Step 8 labels) -------------------------
    pts <- dplyr::bind_rows(lapply(genes, function(g) {
      df <- subclone_labs[[g]]
      if (!"gene" %in% names(df)) df$gene <- g
      if (!"type" %in% names(df)) df$type <- "normal"
      df$time <- as.numeric(df$time)
      df
    }))
    
    # --------- 4) t_opt + rectangle around t_opt (auto width) ------------------
    if (is.null(rect_width)) {
      step_med <- dplyr::bind_rows(lapply(genes, function(g) {
        tt <- sort(unique(as.numeric(gene_sims[[g]]$time)))
        data.frame(dt = diff(tt), gene = g)
      })) %>%
        dplyr::summarise(h = stats::median(dt, na.rm = TRUE)) %>%
        dplyr::pull(h)
      if (!is.finite(step_med) || step_med <= 0) step_med <- 0.02 * diff(range(bands$time))
      rect_width <- 2 * step_med
    }
    
    ribbon_topt <- dplyr::bind_rows(lapply(genes, function(g) {
      t_opt <- as.numeric(t_opt_list[[g]]$t_opt)
      bd <- dplyr::filter(bands, gene == g)
      idx <- which.min(abs(bd$time - t_opt))
      data.frame(
        gene = g,
        t_opt = t_opt,
        lower = bd$low[idx],
        upper = bd$hi[idx],
        xmin = t_opt - rect_width/2,
        xmax = t_opt + rect_width/2,
        ymin = bd$low[idx],
        ymax = bd$hi[idx]
      )
    }))
    
    # --------- 5) annotations: kept subclonal genes & counts -------------------
    ann <- filtered_df %>%
      dplyr::filter(driver %in% genes) %>%
      dplyr::group_by(driver) %>%
      dplyr::summarise(
        kept_total = dplyr::n_distinct(gene),
        kept_enr   = sum(status == "enriched", na.rm = TRUE),
        kept_dep   = sum(status == "depleted", na.rm = TRUE),
        kept_list  = {
          u <- unique(gene)
          if (length(u) == 0) "(none)" else paste(head(u, max_labels), collapse = ", ")
        },
        .groups = "drop"
      ) %>%
      dplyr::right_join(data.frame(driver = genes), by = "driver") %>%
      dplyr::mutate(
        kept_total = tidyr::replace_na(kept_total, 0),
        kept_enr   = tidyr::replace_na(kept_enr,   0),
        kept_dep   = tidyr::replace_na(kept_dep,   0),
        kept_list  = ifelse(is.na(kept_list), "(none)", kept_list),
        gene = driver
      ) %>%
      dplyr::select(gene, kept_total, kept_enr, kept_dep, kept_list)
    
    # position annotation slightly above the band
    ytops <- bands %>% dplyr::group_by(gene) %>% dplyr::summarise(y_top = max(hi, na.rm = TRUE), .groups = "drop")
    xtops <- bands %>% dplyr::group_by(gene) %>% dplyr::summarise(x_min = min(time), .groups = "drop")
    ann <- ann %>%
      dplyr::left_join(ytops, by = "gene") %>%
      dplyr::left_join(xtops, by = "gene") %>%
      dplyr::mutate(
        y_lab = y_top * 1.02,
        x_lab = x_min,
        label = paste0("Kept: ", kept_total, "  (E=", kept_enr, ", D=", kept_dep, ")\n",
                       "Genes: ", kept_list)
      )
    
    vlines <- ribbon_topt %>% dplyr::select(gene, t_opt)
    
    # --------- 6) plot ----------------------------------------------------------
    p <- ggplot() +
      # band over all time
      geom_ribbon(data = bands, aes(x = time, ymin = low, ymax = hi),
                  alpha = 0.25, fill = "#DDAA33") +
      geom_line(data = bands, aes(x = time, y = mean),
                color = "#444444", linewidth = mean_linewidth) +
      # shaded rectangle around t_opt
      geom_rect(data = ribbon_topt,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "#DDAA33", alpha = rect_alpha) +
      geom_vline(data = vlines, aes(xintercept = t_opt),
                 linetype = 2, color = "grey40", linewidth = vline_width) +
      # points from subclone labels
      geom_point(data = pts, aes(x = time, y = x, color = type),
                 alpha = 0.8, size = point_size) +
      scale_color_manual(values = c(low = "blue", normal = "grey60", high = "red")) +
      facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
      labs(
        title = "CME envelopes with subclone labels at t_opt",
        x = "Time", y = "Expression", color = "Subclone"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color="black"),
            panel.border = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"),
            plot.margin = margin(8, 40, 8, 8)) +
      coord_cartesian(clip = "off") +
      # driver-level annotation (kept counts and gene list)
      geom_text(data = ann, aes(x = x_lab, y = y_lab, label = label),
                hjust = 0, vjust = 0, size = 3, inherit.aes = FALSE)
    
    p
  }
  
  filteredSubclonalGenes <- plot_cme_overlays_with_filter(
    gene_sims     = gene_simulations,
    subclone_labs = subclone_labels,
    t_opt_list    = t_opt_list,
    genes         = genes_of_interest,
    filtered_df   = filtered_subclonal_genes_all,  # can lack 'status' now
    alpha         = 0.05,
    ncol          = gene_column_number,
    max_labels    = 6
  )
  
  #@@print(filteredSubclonalGenes)
  
  # -----------------------------
  # Step 13 (updated): Compute heterogeneity metrics per driver (new variables)
  # -----------------------------
  sentence16 <- "Compute heterogeneity metrics per driver!"  # @@@@@@@@@@@@
  print(sentence16)  # @@@@@@@@@@@@
  
  # helper: Shannon entropy per row (expects nonnegative values)
  row_entropy <- function(mat, eps = 1e-12) {
    rs <- rowSums(mat, na.rm = TRUE)
    P  <- mat / pmax(rs, eps)
    P[P <= 0] <- 1                       # makes p*log(p)=0 where p==0
    -rowSums(P * log(P), na.rm = TRUE)
  }
  
  # Ensure filtered_subclonal_genes_all exists and has expected columns
  if (!exists("filtered_subclonal_genes_all") || !is.data.frame(filtered_subclonal_genes_all)) {
    filtered_subclonal_genes_all <- tibble::tibble(driver=character(), gene=character(), status=character())
  } else {
    # normalize common alt names
    if (!"driver" %in% names(filtered_subclonal_genes_all)) {
      d_alt <- intersect(c("Driver","DRIVER","driver_gene"), names(filtered_subclonal_genes_all))
      if (length(d_alt)) names(filtered_subclonal_genes_all)[names(filtered_subclonal_genes_all)==d_alt[1]] <- "driver"
    }
    if (!"gene" %in% names(filtered_subclonal_genes_all)) {
      g_alt <- intersect(c("GENE","SYMBOL","symbol","Gene"), names(filtered_subclonal_genes_all))
      if (length(g_alt)) names(filtered_subclonal_genes_all)[names(filtered_subclonal_genes_all)==g_alt[1]] <- "gene"
    }
    if (!"status" %in% names(filtered_subclonal_genes_all)) filtered_subclonal_genes_all$status <- NA_character_
    filtered_subclonal_genes_all <- filtered_subclonal_genes_all %>%
      mutate(driver = as.character(driver), gene = as.character(gene), status = as.character(status))
  }
  
  # -----------------------------
  # Robust accessors for subclone_labels
  # -----------------------------
  # return a data.frame of labels for a given driver, or NULL
  .get_labels <- function(driver, subclone_labels) {
    # list-of-dfs (preferred)
    if (is.list(subclone_labels) && !is.data.frame(subclone_labels)) {
      if (!driver %in% names(subclone_labels)) return(NULL)
      df <- subclone_labels[[driver]]
      if (!is.data.frame(df) || nrow(df) == 0) return(NULL)
      return(df)
    }
    # single data.frame with a 'driver' column
    if (is.data.frame(subclone_labels)) {
      if (!"driver" %in% names(subclone_labels)) return(NULL)
      df <- subclone_labels[subclone_labels$driver == driver, , drop = FALSE]
      if (!nrow(df)) return(NULL)
      return(df)
    }
    NULL
  }
  
  # Which drivers are actually available in subclone_labels?
  .available_drivers <- function(genes_of_interest, subclone_labels) {
    if (is.list(subclone_labels) && !is.data.frame(subclone_labels)) {
      return(intersect(genes_of_interest, names(subclone_labels)))
    } else if (is.data.frame(subclone_labels) && "driver" %in% names(subclone_labels)) {
      return(intersect(genes_of_interest, unique(subclone_labels$driver)))
    } else {
      return(character(0))
    }
  }
  
  drivers <- .available_drivers(genes_of_interest, subclone_labels)
  if (!length(drivers)) {
    message("No overlapping drivers between genes_of_interest and subclone_labels. Returning empty heterogeneity table.")
    tumor_heterogeneity_df <- tibble::tibble(
      driver      = character(), GENE = character(),
      n_rare      = integer(),   n_other = integer(),
      mean_R      = double(),    var_R   = double(),
      sd_R        = double(),    cv_R    = double(),
      entropy_R   = double(),
      mean_O      = double(),    var_O   = double(),
      sd_O        = double(),    cv_O    = double(),
      FC_R_over_O = double(),
      status      = character()
    )
  } else {
    # -----------------------------
    # Build tumor_heterogeneity (robust; no [[g]] usage)
    # -----------------------------
    tumor_heterogeneity <- lapply(drivers, function(driver) {
      # --- Labels for this driver ---
      labs <- .get_labels(driver, subclone_labels)
      if (is.null(labs)) {
        message("[", driver, "] no labels found in subclone_labels; skipping.")
        return(NULL)
      }
      # need source_id + type
      req_cols <- c("source_id","type")
      if (!all(req_cols %in% names(labs))) {
        message("[", driver, "] labels missing required columns {source_id, type}; skipping.")
        return(NULL)
      }
      
      # Rare sample IDs at t_opt (low/high)
      R_ids <- labs$source_id[labs$type %in% c("low","high")]
      R_ids <- suppressWarnings(as.integer(as.character(R_ids)))
      R_ids <- unique(R_ids[!is.na(R_ids)])
      if (!length(R_ids)) {
        message("[", driver, "] no rare sample IDs (low/high); skipping.")
        return(NULL)
      }
      
      # Complement set
      if (!exists("gene_mat")) {
        stop("gene_mat is missing; cannot compute heterogeneity.")
      }
      O_ids <- setdiff(seq_len(ncol(gene_mat)), R_ids)
      
      # Kept non-driver genes for this driver after Step 11
      kept_genes <- filtered_subclonal_genes_all %>%
        dplyr::filter(.data$driver == .env$driver) %>%
        dplyr::pull(gene) %>%
        unique()
      
      kept_genes <- intersect(kept_genes, rownames(gene_mat))
      if (!length(kept_genes)) {
        message("[", driver, "] no kept genes present in gene_mat; skipping.")
        return(NULL)
      }
      
      # Expression matrices (SYMBOL rows)
      expr_R <- gene_mat[kept_genes, R_ids, drop = FALSE]
      expr_O <- if (length(O_ids)) gene_mat[kept_genes, O_ids, drop = FALSE] else NULL
      
      # --- Heterogeneity metrics within rare set ---
      mean_R <- rowMeans(expr_R, na.rm = TRUE)
      var_R  <- matrixStats::rowVars(expr_R, na.rm = TRUE)
      sd_R   <- sqrt(var_R)
      cv_R   <- sd_R / pmax(mean_R, 1e-12)
      ent_R  <- row_entropy(pmax(expr_R, 0))
      
      # Optional contrasts vs. others
      if (!is.null(expr_O) && ncol(expr_O) > 0) {
        mean_O <- rowMeans(expr_O, na.rm = TRUE)
        var_O  <- matrixStats::rowVars(expr_O, na.rm = TRUE)
        sd_O   <- sqrt(var_O)
        cv_O   <- sd_O / pmax(mean_O, 1e-12)
        fc_RO  <- mean_R / pmax(mean_O, 1e-12)
      } else {
        mean_O <- var_O <- sd_O <- cv_O <- fc_RO <- rep(NA_real_, length(kept_genes))
      }
      
      # Bring along status (enriched/depleted) from Step 10/11 (if present)
      status_map <- filtered_subclonal_genes_all %>%
        dplyr::filter(.data$driver == .env$driver) %>%
        dplyr::select(gene, status) %>%
        dplyr::distinct()
      
      out <- data.frame(
        driver      = driver,
        GENE        = kept_genes,
        n_rare      = length(R_ids),
        n_other     = length(O_ids),
        mean_R      = as.numeric(mean_R),
        var_R       = as.numeric(var_R),
        sd_R        = as.numeric(sd_R),
        cv_R        = as.numeric(cv_R),
        entropy_R   = as.numeric(ent_R),
        mean_O      = as.numeric(mean_O),
        var_O       = as.numeric(var_O),
        sd_O        = as.numeric(sd_O),
        cv_O        = as.numeric(cv_O),
        FC_R_over_O = as.numeric(fc_RO),
        stringsAsFactors = FALSE,
        row.names = NULL
      ) %>%
        dplyr::left_join(status_map, by = c("GENE" = "gene"))
      
      out
    })
    
    tumor_heterogeneity_df <- dplyr::bind_rows(tumor_heterogeneity)
    
    if (!nrow(tumor_heterogeneity_df)) {
      message("No heterogeneity rows constructed for available drivers. Returning empty table.")
      tumor_heterogeneity_df <- tibble::tibble(
        driver      = character(), GENE = character(),
        n_rare      = integer(),   n_other = integer(),
        mean_R      = double(),    var_R   = double(),
        sd_R        = double(),    cv_R    = double(),
        entropy_R   = double(),
        mean_O      = double(),    var_O   = double(),
        sd_O        = double(),    cv_O    = double(),
        FC_R_over_O = double(),
        status      = character()
      )
    }
  }
  
  # (optional) inspect
  # str(tumor_heterogeneity_df)
  
  # ===============================
  # 12.x Tumor heterogeneity visualization (robust with fallbacks)
  # ===============================
  #@ suppressPackageStartupMessages({
  #@  library(dplyr)
  #@  library(tidyr)
  #@  library(ggplot2)
  #@  library(pheatmap)
  #@  library(RColorBrewer)
  #@})
  
  # ------------------------------
  # ALWAYS-ON heterogeneity heatmap
  # ------------------------------
  # ================================
  # STEP: Tumor heterogeneity heatmap
  #  - Compact legends (no overlap)
  #  - Driver legend named 'matrix_25'
  #  - Robust to list/data.frame subclone_labels
  #  - Safe fallbacks if filtered/assoc sets are empty
  # ================================
  suppressPackageStartupMessages({
    library(dplyr)
    library(pheatmap)
    library(RColorBrewer)
    library(matrixStats)
  })
  
  always_heatmap_tumor_heterogeneity <- function(gene_mat,
                                                 genes_of_interest,
                                                 subclone_labels,                 # list-of-dfs or one data.frame with 'driver'
                                                 filtered_subclonal_genes_all=NULL,
                                                 assoc_results=NULL,              # optional per-driver table (gene, FC, padj, status)
                                                 topK_assoc=10,
                                                 n_topvar_driver=15,
                                                 n_global_topvar=50,
                                                 driver_legend_name="matrix_25",
                                                 fontsize=11,
                                                 fontsize_row=7) {
    stopifnot(is.matrix(gene_mat) || is.data.frame(gene_mat))
    gene_mat <- as.matrix(gene_mat)
    if (is.null(rownames(gene_mat))) rownames(gene_mat) <- paste0("gene_", seq_len(nrow(gene_mat)))
    if (is.null(colnames(gene_mat))) colnames(gene_mat) <- paste0("sample_", seq_len(ncol(gene_mat)))
    
    # --- helpers for subclone_labels shapes ---
    .drivers_available <- function(goi, lbls) {
      if (is.list(lbls) && !is.data.frame(lbls)) {
        intersect(goi, names(lbls))
      } else if (is.data.frame(lbls) && "driver" %in% names(lbls)) {
        intersect(goi, unique(lbls$driver))
      } else character(0)
    }
    .labels_for <- function(g, lbls) {
      if (is.list(lbls) && !is.data.frame(lbls)) {
        if (!g %in% names(lbls)) return(NULL)
        lbls[[g]]
      } else if (is.data.frame(lbls) && "driver" %in% names(lbls)) {
        lbls[lbls$driver == g, , drop=FALSE]
      } else NULL
    }
    
    drivers <- .drivers_available(genes_of_interest, subclone_labels)
    
    # --- 1) Rare sample indices by driver; fallback to ALL samples
    rare_by_driver <- setNames(vector("list", length(drivers)), drivers)
    for (g in drivers) {
      labs <- .labels_for(g, subclone_labels)
      if (is.data.frame(labs) && all(c("source_id","type") %in% names(labs))) {
        ids <- labs$source_id[labs$type %in% c("low","high")]
        ids <- suppressWarnings(as.integer(as.character(ids)))
        rare_by_driver[[g]] <- unique(ids[!is.na(ids)])
      } else {
        rare_by_driver[[g]] <- integer(0)
      }
    }
    rare_samples_all <- sort(unique(unlist(rare_by_driver)))
    if (!length(rare_samples_all)) {
      message("No rare samples found; using ALL samples for columns.")
      rare_samples_all <- seq_len(ncol(gene_mat))
    }
    
    # --- 2) Build kept_by_driver with fallbacks (filtered -> assoc -> topvar per-driver)
    fsg <- NULL
    if (!is.null(filtered_subclonal_genes_all) && nrow(filtered_subclonal_genes_all) > 0) {
      fsg <- filtered_subclonal_genes_all
      if (!"driver" %in% names(fsg)) {
        d_alt <- intersect(c("Driver","DRIVER","driver_gene"), names(fsg))
        if (length(d_alt)) names(fsg)[names(fsg)==d_alt[1]] <- "driver"
      }
      if (!"gene" %in% names(fsg)) {
        g_alt <- intersect(c("GENE","SYMBOL","symbol","Gene"), names(fsg))
        if (length(g_alt)) names(fsg)[names(fsg)==g_alt[1]] <- "gene"
      }
    }
    
    kept_by_driver <- setNames(vector("list", length(drivers)), drivers)
    kept_source    <- list()
    
    for (g in drivers) {
      g_filtered <- character(0)
      if (!is.null(fsg)) {
        g_filtered <- fsg %>%
          dplyr::filter(.data$driver == .env$g) %>%
          dplyr::pull(gene) %>% unique() %>%
          intersect(rownames(gene_mat))
      }
      
      g_assoc <- character(0)
      if (!length(g_filtered) && !is.null(assoc_results) && !is.null(assoc_results[[g]])) {
        g_assoc <- assoc_results[[g]] %>%
          dplyr::filter(status %in% c("enriched","depleted")) %>%
          dplyr::arrange(padj, dplyr::desc(abs(log2(FC)))) %>%
          dplyr::slice_head(n = topK_assoc) %>%
          dplyr::pull(gene) %>% unique() %>%
          intersect(rownames(gene_mat))
      }
      
      g_topv <- character(0)
      if (!length(g_filtered) && !length(g_assoc)) {
        cols <- if (length(rare_by_driver[[g]])) rare_by_driver[[g]] else seq_len(ncol(gene_mat))
        mat  <- gene_mat[, cols, drop = FALSE]
        if (g %in% rownames(mat)) mat <- mat[setdiff(rownames(mat), g), , drop = FALSE]
        if (nrow(mat) > 0) {
          vrs <- matrixStats::rowVars(mat, na.rm = TRUE)
          ord <- order(vrs, decreasing = TRUE)
          idx <- ord[seq_len(min(n_topvar_driver, length(ord)))]
          g_topv <- rownames(mat)[idx]
        }
      }
      
      genes_g <- unique(c(g_filtered, g_assoc, g_topv))
      kept_by_driver[[g]] <- genes_g
      
      if (length(genes_g)) {
        src <- rep("filtered", length(genes_g))
        if (length(g_assoc)) src[genes_g %in% g_assoc] <- "assoc"
        if (length(g_topv))  src[genes_g %in% g_topv]  <- "topvar_driver"
        kept_source[[g]] <- data.frame(gene = genes_g, source = src, driver = g, stringsAsFactors = FALSE)
      } else {
        kept_source[[g]] <- data.frame(gene = character(0), source = character(0), driver = character(0))
      }
    }
    
    # union across drivers; if still empty, use global top-variance genes
    all_kept_genes <- sort(unique(unlist(kept_by_driver)))
    if (!length(all_kept_genes)) {
      message("No kept/assoc/topvar-per-driver genes; using GLOBAL top-variance genes.")
      vrs <- matrixStats::rowVars(gene_mat[, rare_samples_all, drop = FALSE], na.rm = TRUE)
      ord <- order(vrs, decreasing = TRUE)
      idx <- ord[seq_len(min(n_global_topvar, length(ord)))]
      all_kept_genes <- rownames(gene_mat)[idx]
      kept_source[["__global__"]] <- data.frame(
        gene = all_kept_genes, source = "topvar_global", driver = "(global)", stringsAsFactors = FALSE
      )
    }
    
    # --- 3) Build matrix + COMPACT annotations (so legends don’t overlap) ----
    expr_sub_all <- gene_mat[all_kept_genes, rare_samples_all, drop = FALSE]
    expr_scaled  <- t(scale(t(expr_sub_all))); expr_scaled[is.na(expr_scaled)] <- 0
    
    # Row annotation: collapse composite origins -> single class
    source_df <- dplyr::bind_rows(kept_source)
    if (!nrow(source_df)) source_df <- data.frame(gene = all_kept_genes, source="(none)", driver="(none)")
    
    gene_to_source <- source_df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        origin_raw = paste(sort(unique(source[source != "(none)"])), collapse = ", "),
        n_drivers  = dplyr::n_distinct(driver[driver != "(global)"]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Origin = dplyr::case_when(
          grepl("\\bfiltered\\b", origin_raw)      ~ "filtered",
          grepl("\\bassoc\\b",    origin_raw)      ~ "assoc",
          grepl("topvar_driver",  origin_raw)      ~ "topvar_driver",
          grepl("topvar_global",  origin_raw)      ~ "topvar_global",
          origin_raw == ""                          ~ "(none)",
          TRUE                                      ~ "mixed"
        )
      )
    
    annotation_row <- gene_to_source %>%
      dplyr::select(gene, Origin, n_drivers) %>%
      dplyr::mutate(row.names = gene) %>% as.data.frame()
    rownames(annotation_row) <- annotation_row$row.names; annotation_row$row.names <- NULL
    
    # Column annotation: legend named EXACTLY 'matrix_25' (driver names + None/Multiple)
    sample_to_drivers <- lapply(seq_along(rare_samples_all), function(j) {
      s <- rare_samples_all[j]
      names(rare_by_driver)[vapply(rare_by_driver, function(v) s %in% v, logical(1))]
    })
    driver_label <- vapply(sample_to_drivers, function(v) {
      if (length(v) == 0) "None" else if (length(v) == 1) v else "Multiple"
    }, character(1))
    annotation_col <- data.frame(
      matrix_25 = factor(driver_label, levels = c("None", drivers, "Multiple")),
      row.names = colnames(expr_sub_all),
      check.names = FALSE
    )
    
    # ---- Colors that MATCH PRESENT LEVELS EXACTLY ----
    # Origin colors (subset to levels that appear)
    annotation_row$Origin <- droplevels(factor(annotation_row$Origin))
    origin_master <- c(filtered = "#1b9e77",
                       assoc = "#d95f02",
                       topvar_driver = "#7570b3",
                       topvar_global = "#636363",
                       mixed = "#a6cee3",
                       `(none)` = "#bbbbbb")
    origin_present <- levels(annotation_row$Origin)
    origin_colors  <- origin_master[origin_present]
    annotation_row$Origin <- factor(annotation_row$Origin, levels = names(origin_colors))
    
    # Driver legend ('matrix_25') palette (subset to present levels)
    drv_base <- RColorBrewer::brewer.pal(8, "Set2")
    drv_pal  <- if (length(drivers) <= length(drv_base)) drv_base[seq_along(drivers)]
    else colorRampPalette(drv_base)(length(drivers))
    names(drv_pal) <- drivers
    m25_levels <- levels(annotation_col$matrix_25)
    m25_colors <- c(`None` = "#e0e0e0", drv_pal, `Multiple` = "#444444")[m25_levels]
    
    ann_colors <- list(
      Origin    = origin_colors,    # factor
      matrix_25 = m25_colors        # factor
      # numeric 'n_drivers' is left without a palette -> pheatmap uses a continuous legend
    )
    
    heat_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101)
    
    # Use readable names for the numeric driver count
    annotation_row_use <- annotation_row[, c("Origin","n_drivers"), drop = FALSE]
    colnames(annotation_row_use)[2] <- "Drivers (#)"
    
    # --- 4) Draw heatmap (compact, non-overlapping legends) ---
    p <- pheatmap::pheatmap(
      expr_scaled,
      cluster_rows = TRUE, cluster_cols = TRUE,
      show_rownames = TRUE, show_colnames = FALSE,
      annotation_row = annotation_row_use,
      annotation_col = annotation_col,
      annotation_colors = list(
        Origin       = ann_colors$Origin,
        matrix_25    = ann_colors$matrix_25
        # 'Drivers (#)' is numeric; no discrete colors supplied
      ),
      color = heat_pal,
      main = "Tumor heterogeneity: subclonal gene expression (row-scaled)",
      fontsize = fontsize,
      fontsize_row = fontsize_row,
      legend = TRUE,
      annotation_legend = TRUE,
      border_color = NA
    )
    
    return(p)
  }
  
  # ---------- Run the step ----------
  Heatmap_tumorHeterogeneity <- always_heatmap_tumor_heterogeneity(
    gene_mat                  = gene_mat,
    genes_of_interest         = genes_of_interest,
    subclone_labels           = subclone_labels,
    filtered_subclonal_genes_all = filtered_subclonal_genes_all,  # can be NULL/empty
    assoc_results             = assoc_results,                     # can be NULL/empty
    topK_assoc                = 10,
    n_topvar_driver           = 15,
    n_global_topvar           = 50,
    driver_legend_name        = "matrix_25",
    fontsize                  = 11,
    fontsize_row              = 7
  )
  
  # View
  #@@@print(Heatmap_tumorHeterogeneity)
  
  tumor_heterogeneity_df <- dplyr::bind_rows(tumor_heterogeneity)
  rownames(tumor_heterogeneity_df) <- NULL
  
  # Optional quick peek:
  # print(dplyr::arrange(tumor_heterogeneity_df, driver, dplyr::desc(cv_R)))
  
  
  # ---------- Barplot: fraction of low vs high among rare per driver ----------
  heterogeneity_summary_df <- dplyr::bind_rows(lapply(drivers, function(g) {
    labs <- subclone_labels[[g]]
    if (is.null(labs) || !nrow(labs)) {
      return(data.frame(driver = g, low = 0, high = 0, n_rare = 0))
    }
    n_low  <- sum(labs$type == "low",  na.rm = TRUE)
    n_high <- sum(labs$type == "high", na.rm = TRUE)
    n_rare <- n_low + n_high
    if (n_rare == 0) data.frame(driver = g, low = 0, high = 0, n_rare = 0) else
      data.frame(driver = g, low = n_low / n_rare, high = n_high / n_rare, n_rare = n_rare)
  }))
  
  heterogeneity_summary_long <- heterogeneity_summary_df %>%
    tidyr::pivot_longer(cols = c("low","high"), names_to = "subclone_type", values_to = "fraction")
  
  Barplot_tumorHeterogeneity <- ggplot(heterogeneity_summary_long, aes(x = driver, y = fraction, fill = subclone_type)) +
    geom_col(position = "stack") +
    geom_text(aes(label = scales::percent(fraction, accuracy = 1)),
              position = position_stack(vjust = 0.5), size = 3, color = "white") +
    scale_fill_manual(values = c(low = "blue", high = "red"),
                      name = "Subclone Type", labels = c("High", "Low")) +
    labs(title = "Fraction of low/high subclones at t\u2092\u209a per driver",
         subtitle = "Proportions are within the (low+high) rare set for each driver",
         x = "Driver gene", y = "Fraction") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color="black"),
          panel.border = element_blank())
  
  #@@@ print(Barplot_tumorHeterogeneity)
  
  # =============================================================================
  # Step 15: Compute subclone likelihood & detect high-risk (robust to missing cols)
  # =============================================================================
  sentence17 <- "Compute subclone likelihood & detect high-risk!"  # @@@@@@@@@@@@
  print(sentence17)  # @@@@@@@@@@@@
  
  #@suppressPackageStartupMessages({
  #@  library(dplyr)
  #@  library(tidyr)
  #@ library(ggplot2)
  #@ library(matrixStats)
  #@})
  
  .top_assoc <- function(g, topK = 10) {
    if (!exists("assoc_results") || is.null(assoc_results[[g]])) return(character(0))
    assoc_results[[g]] %>%
      dplyr::filter(status %in% c("enriched","depleted")) %>%
      dplyr::arrange(padj, dplyr::desc(abs(log2(FC)))) %>%
      dplyr::slice_head(n = topK) %>%
      dplyr::pull(gene) %>%
      unique()
  }
  
  ##########@@@@@@@@@@@@@@@@@@@
  #@@suppressPackageStartupMessages({
  #@@  library(dplyr)
  #@@  library(matrixStats)
  #@@})
  
  # ---------- CONFIG (you can tune) ----------
  fallback_top_var_genes <- 20   # if no kept/assoc genes, pick top-variance genes in R_ids
  driver_tail_q          <- 0.05 # for fallback rare IDs from driver gene expression (5%/95%)
  
  # ---------- Utilities ----------
  .safe_unique_int <- function(x) {
    x <- suppressWarnings(as.integer(as.character(x)))
    unique(x[!is.na(x)])
  }
  
  .get_rare_ids <- function(driver, subclone_labels, gene_mat, driver_tail_q = 0.05) {
    # 1) subclone_labels low/high at t_opt
    if (!is.null(subclone_labels[[driver]]) &&
        all(c("source_id","type") %in% names(subclone_labels[[driver]]))) {
      labs <- subclone_labels[[driver]]
      R <- .safe_unique_int(labs$source_id[labs$type %in% c("low","high")])
      if (length(R)) return(list(ids = R, how = "subclone_labels"))
    }
    # 2) driver_rare_samples if present
    if (exists("driver_rare_samples") && !is.null(driver_rare_samples[[driver]])) {
      R <- .safe_unique_int(driver_rare_samples[[driver]])
      if (length(R)) return(list(ids = R, how = "driver_rare_samples"))
    }
    # 3) Data-driven tails from driver’s own expression (if present in gene_mat)
    if (!is.null(gene_mat) && driver %in% rownames(gene_mat)) {
      v <- as.numeric(gene_mat[driver, ])
      lo <- as.numeric(quantile(v, probs = driver_tail_q, na.rm = TRUE))
      hi <- as.numeric(quantile(v, probs = 1 - driver_tail_q, na.rm = TRUE))
      R <- which(v <= lo | v >= hi)
      R <- .safe_unique_int(R)
      if (length(R)) return(list(ids = R, how = "driver_gene_tails"))
    }
    list(ids = integer(0), how = "none")
  }
  
  .get_kept_genes <- function(driver, filtered_df, assoc_results, gene_mat, R_ids, topK_var = 20) {
    # 1) kept after ribbon filter
    kept <- character(0)
    if (!is.null(filtered_df) && is.data.frame(filtered_df)) {
      # standardize column names
      df <- filtered_df
      if (!"driver" %in% names(df)) {
        d_alt <- intersect(c("Driver","DRIVER","driver_gene"), names(df))
        if (length(d_alt)) names(df)[names(df) == d_alt[1]] <- "driver"
      }
      if (!"gene" %in% names(df)) {
        g_alt <- intersect(c("GENE","SYMBOL","symbol","Gene"), names(df))
        if (length(g_alt)) names(df)[names(df) == g_alt[1]] <- "gene"
      }
      if (all(c("driver","gene") %in% names(df))) {
        kept <- df %>% dplyr::filter(.data$driver == .env$driver) %>% dplyr::pull(gene) %>% unique()
      }
    }
    if (length(kept)) return(list(genes = kept, how = "filtered"))
    
    # 2) top assoc_results fallback
    if (exists("assoc_results") && !is.null(assoc_results[[driver]])) {
      kept <- .top_assoc(driver, topK = max(10, topK_var))
      if (length(kept)) return(list(genes = kept, how = "assoc_results"))
    }
    
    # 3) top-variance genes within rare set (exclude the driver row if present)
    if (!is.null(gene_mat) && length(R_ids)) {
      mat <- gene_mat[, R_ids, drop = FALSE]
      if (driver %in% rownames(mat)) mat <- mat[setdiff(rownames(mat), driver), , drop = FALSE]
      if (nrow(mat) > 0) {
        vrs <- matrixStats::rowVars(as.matrix(mat), na.rm = TRUE)
        ord <- order(vrs, decreasing = TRUE)
        top_idx <- ord[seq_len(min(topK_var, length(ord)))]
        kept <- rownames(mat)[top_idx]
        if (length(kept)) return(list(genes = kept, how = "top_variance"))
      }
    }
    list(genes = character(0), how = "none")
  }
  
  # ---------- 15.1 Build/standardize tumor_heterogeneity_df (robust) ----------
  needed_cols <- c("driver","GENE","mean_R","mean_O","cv_R","cv_O","entropy_R","FC_R_over_O","n_rare","n_other")
  
  build_tumor_heterogeneity <- function() {
    drivers <- intersect(genes_of_interest, names(subclone_labels))
    if (length(drivers) == 0) {
      message("No drivers available to build tumor_heterogeneity_df."); 
      # return empty tibble with required cols
      return(tibble::tibble(
        driver = character(), GENE = character(), n_rare = integer(), n_other = integer(),
        mean_R = double(), mean_O = double(), cv_R = double(), cv_O = double(),
        entropy_R = double(), FC_R_over_O = double()
      ))
    }
    
    diag_msgs <- list()
    hetero_list <- lapply(drivers, function(driver) {
      # Rare IDs with fallbacks
      rid <- .get_rare_ids(driver, subclone_labels, gene_mat, driver_tail_q = driver_tail_q)
      R_ids <- rid$ids
      if (!length(R_ids)) {
        diag_msgs[[length(diag_msgs) + 1]] <<- paste0("[", driver, "] no rare samples (how=", rid$how, "). Skipping.")
        return(NULL)
      }
      O_ids <- setdiff(seq_len(ncol(gene_mat)), R_ids)
      
      # Kept genes with fallbacks
      gkp <- .get_kept_genes(driver, filtered_subclonal_genes_all, assoc_results, gene_mat,
                             R_ids, topK_var = fallback_top_var_genes)
      kept_genes <- intersect(unique(gkp$genes), rownames(gene_mat))
      if (!length(kept_genes)) {
        diag_msgs[[length(diag_msgs) + 1]] <<- paste0("[", driver, "] no kept genes (how=", gkp$how, "). Skipping.")
        return(NULL)
      }
      
      # Expression matrices (SYMBOL rows)
      expr_R <- gene_mat[kept_genes, R_ids, drop = FALSE]
      expr_O <- if (length(O_ids)) gene_mat[kept_genes, O_ids, drop = FALSE] else NULL
      
      # Metrics within rare set
      mean_R <- rowMeans(expr_R, na.rm = TRUE)
      var_R  <- matrixStats::rowVars(expr_R, na.rm = TRUE)
      sd_R   <- sqrt(var_R)
      cv_R   <- sd_R / pmax(mean_R, 1e-12)
      ent_R  <- row_entropy(pmax(expr_R, 0))
      
      # Contrasts vs. others
      if (!is.null(expr_O) && ncol(expr_O) > 0) {
        mean_O <- rowMeans(expr_O, na.rm = TRUE)
        var_O  <- matrixStats::rowVars(expr_O, na.rm = TRUE)
        sd_O   <- sqrt(var_O)
        cv_O   <- sd_O / pmax(mean_O, 1e-12)
        fc_RO  <- mean_R / pmax(mean_O, 1e-12)
      } else {
        mean_O <- var_O <- sd_O <- cv_O <- fc_RO <- rep(NA_real_, length(kept_genes))
      }
      
      tibble::tibble(
        driver      = driver,
        GENE        = kept_genes,
        n_rare      = length(R_ids),
        n_other     = length(O_ids),
        mean_R      = as.numeric(mean_R),
        mean_O      = as.numeric(mean_O),
        cv_R        = as.numeric(cv_R),
        cv_O        = as.numeric(cv_O),
        entropy_R   = as.numeric(ent_R),
        FC_R_over_O = as.numeric(fc_RO)
      )
    })
    
    out <- dplyr::bind_rows(hetero_list)
    
    # Print concise diagnostics
    if (length(diag_msgs)) message(paste(diag_msgs, collapse = "\n"))
    
    # Return empty tibble with expected columns instead of stopping
    if (is.null(out) || !nrow(out)) {
      message("No heterogeneity rows constructed. Returning empty table with required columns.")
      out <- tibble::tibble(
        driver = character(), GENE = character(), n_rare = integer(), n_other = integer(),
        mean_R = double(), mean_O = double(), cv_R = double(), cv_O = double(),
        entropy_R = double(), FC_R_over_O = double()
      )
    }
    out
  }
  
  # Use existing if valid; otherwise rebuild
  if (exists("tumor_heterogeneity_df") && is.data.frame(tumor_heterogeneity_df)) {
    th_df <- tumor_heterogeneity_df
    if (!"driver" %in% names(th_df)) {
      d_alt <- intersect(c("Driver","DRIVER","driver_gene","gene_driver"), names(th_df))
      if (length(d_alt)) names(th_df)[names(th_df) == d_alt[1]] <- "driver"
    }
    if (!"GENE" %in% names(th_df)) {
      g_alt <- intersect(c("gene","Gene","SYMBOL","symbol"), names(th_df))
      if (length(g_alt)) names(th_df)[names(th_df) == g_alt[1]] <- "GENE"
    }
    if (!all(needed_cols %in% names(th_df))) {
      message("tumor_heterogeneity_df lacks required columns; rebuilding from current objects.")
      th_df <- build_tumor_heterogeneity()
    }
  } else {
    message("tumor_heterogeneity_df missing; building from current objects.")
    th_df <- build_tumor_heterogeneity()
  }
  
  # Final standardized table
  tumor_heterogeneity_df <- th_df %>%
    mutate(driver = as.character(driver), GENE = as.character(GENE))
  
  ##########@@@@@@@@@@@@@@@@@@@
  
  # Ensure types
  th_df <- th_df %>%
    dplyr::mutate(driver = as.character(driver), GENE = as.character(GENE))
  
  # Ensure filtered_subclonal_genes_all exists and is standardized
  if (!exists("filtered_subclonal_genes_all") || !is.data.frame(filtered_subclonal_genes_all)) {
    filtered_subclonal_genes_all <- tibble::tibble(driver = character(), gene = character(), status = character())
  } else {
    fsg <- filtered_subclonal_genes_all
    if (!"driver" %in% names(fsg)) {
      d_alt <- intersect(c("Driver","DRIVER","driver_gene"), names(fsg))
      if (length(d_alt)) names(fsg)[names(fsg) == d_alt[1]] <- "driver"
    }
    if (!"gene" %in% names(fsg)) {
      g_alt <- intersect(c("GENE","SYMBOL","symbol","Gene"), names(fsg))
      if (length(g_alt)) names(fsg)[names(fsg) == g_alt[1]] <- "gene"
    }
    if (!"status" %in% names(fsg)) fsg$status <- NA_character_
    filtered_subclonal_genes_all <- fsg %>%
      dplyr::mutate(driver = as.character(driver), gene = as.character(gene), status = as.character(status)) %>%
      dplyr::distinct(driver, gene, status)
  }
  
  # Join status and engineer features
  advantage_df <- th_df %>%
    dplyr::left_join(
      filtered_subclonal_genes_all %>% dplyr::select(driver, gene, status),
      by = c("driver" = "driver", "GENE" = "gene")
    ) %>%
    dplyr::mutate(
      fc          = pmax(FC_R_over_O, 1e-12),
      log2_fc     = log2(fc),
      cv_ratio    = cv_R / pmax(cv_O, 1e-12),
      log2_cv     = log2(cv_ratio),
      ent_norm    = entropy_R / pmax(log(n_rare + 1e-12), 1e-12),
      init_expr   = mean_O,
      final_expr  = mean_R
    )
  
  # ---------- 15.2 Likelihood score ----------
  compute_likelihood_score <- function(adv_df, w_fc = 0.6, w_cv = 0.3, w_ent = 0.1) {
    if (is.null(adv_df) || !nrow(adv_df)) return(adv_df)
    mm <- function(x) {
      x <- as.numeric(x); ok <- is.finite(x)
      if (!any(ok)) return(rep(0, length(x)))
      rng <- range(x[ok]); if (diff(rng) == 0) return(rep(0.5, length(x)))
      (x - rng[1]) / diff(rng)
    }
    adv_df %>%
      dplyr::mutate(
        s_fc   = mm(log2_fc),
        s_cv   = mm(log2_cv),
        s_ent  = mm(ent_norm),
        raw_score = w_fc * s_fc + w_cv * s_cv + w_ent * s_ent,
        likelihood_score = mm(raw_score)
      )
  }
  
  advantage_scored <- compute_likelihood_score(advantage_df)
  
  # ---------- 15.3 Detect high-risk subclones ----------
  detect_high_risk_subclones <- function(adv_df, likelihood_threshold = 0.7,
                                         per_driver_percentile = FALSE, percentile = 0.7) {
    if (is.null(adv_df) || !nrow(adv_df)) return(list(flagged_subclones = adv_df, gene_summary = NULL))
    if (per_driver_percentile) {
      thr_tbl <- adv_df %>% dplyr::group_by(driver) %>%
        dplyr::summarize(thr = stats::quantile(likelihood_score, probs = percentile, na.rm = TRUE), .groups = "drop")
      adv_flag <- adv_df %>% dplyr::left_join(thr_tbl, by = "driver") %>% dplyr::mutate(high_risk = likelihood_score >= thr)
    } else {
      adv_flag <- adv_df %>% dplyr::mutate(high_risk = likelihood_score >= likelihood_threshold)
    }
    
    gene_summary <- adv_flag %>% dplyr::group_by(driver) %>%
      dplyr::summarize(
        n_subclones = dplyr::n(),
        n_high_risk = sum(high_risk, na.rm = TRUE),
        prop_high_risk = ifelse(n_subclones > 0, n_high_risk / n_subclones, 0),
        median_score = stats::median(likelihood_score, na.rm = TRUE),
        .groups = "drop"
      ) %>% dplyr::arrange(dplyr::desc(prop_high_risk), dplyr::desc(median_score))
    
    list(flagged_subclones = adv_flag, gene_summary = gene_summary)
  }
  
  risk_results <- detect_high_risk_subclones(
    adv_df = advantage_scored,
    likelihood_threshold = 0.7,
    per_driver_percentile = FALSE,
    percentile = 0.7
  )
  flagged_subclones_df <- risk_results$flagged_subclones
  gene_summary_df       <- risk_results$gene_summary
  #@@@print(gene_summary_df)
  
  # ---------- 15.4 Visualization: high-risk at t_opt over CME envelopes ----------
  visualize_high_risk_subclones <- function(gene_simulations,
                                            flagged_subclones_df,
                                            t_opt_list,
                                            alpha_band = 0.05,
                                            max_labels = 5,
                                            ncol = 3,
                                            label_size = 5.2,          # <- bigger text for subclone labels
                                            label_fontface = "bold",   # <- bold labels
                                            use_label_boxes = TRUE) {  # <- draw boxes behind text
    drivers <- intersect(names(gene_simulations), names(t_opt_list))
    if (is.null(flagged_subclones_df) || !nrow(flagged_subclones_df)) {
      message("No subclones to visualize."); return(ggplot() + theme_void())
    }
    
    bands <- dplyr::bind_rows(lapply(drivers, function(g) {
      df <- gene_simulations[[g]]; df$time <- as.numeric(df$time)
      df %>% dplyr::group_by(time) %>%
        dplyr::summarize(
          mean = mean(x),
          low  = stats::quantile(x, alpha_band/2, na.rm = TRUE),
          hi   = stats::quantile(x, 1 - alpha_band/2, na.rm = TRUE),
          .groups = "drop"
        ) %>% dplyr::mutate(driver = g)
    }))
    
    vlines <- tibble::tibble(
      driver = drivers,
      t_opt  = vapply(drivers, function(g) as.numeric(t_opt_list[[g]]$t_opt), numeric(1))
    )
    
    pts <- flagged_subclones_df %>%
      dplyr::mutate(driver = as.character(driver)) %>%
      dplyr::filter(driver %in% drivers) %>%
      dplyr::left_join(vlines, by = "driver") %>%
      dplyr::mutate(
        Risk   = ifelse(high_risk, "High-risk", "Other"),
        status = ifelse(is.na(status), "enriched", status),
        sizeA  = pmin(8, pmax(2.5, 3 + abs(log2(pmax(FC_R_over_O, 1e-12))))) # slightly larger points
      )
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      lab_tbl <- pts %>% dplyr::group_by(driver) %>%
        dplyr::arrange(dplyr::desc(likelihood_score)) %>%
        dplyr::slice_head(n = max_labels) %>%
        dplyr::ungroup()
    } else lab_tbl <- pts[0, ]
    
    p <- ggplot() +
      geom_ribbon(data = bands, aes(x = time, ymin = low, ymax = hi),
                  alpha = 0.22, fill = "#DDAA33") +
      geom_line(data = bands, aes(x = time, y = mean),
                color = "#444444", linewidth = 0.9) +
      geom_vline(data = vlines, aes(xintercept = t_opt),
                 linetype = 2, color = "black", linewidth = 1.1) +
      geom_point(data = pts, aes(x = t_opt, y = mean_R,
                                 color = Risk, shape = status, size = sizeA),
                 alpha = 0.9, stroke = 0.3) +
      scale_color_manual(values = c("High-risk" = "#D62728", "Other" = "#1F77B4")) +
      scale_shape_manual(values = c(enriched = 16, depleted = 17)) +
      scale_size_identity() +
      facet_wrap(~ driver, scales = "free_y", ncol = ncol) +
      labs(title = "High-risk subclones at t\u2092\u209a overlaid on CME envelopes",
           subtitle = "Point size ∝ |log2 FC_R/O|; shape = enriched/depleted; color = risk flag",
           x = "Time", y = "Expression (rare-cohort mean for kept genes)",
           color = "Risk", shape = "Status") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.line  = element_line(color = "black"),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"))
    
    # Bigger, clearer labels
    if (nrow(lab_tbl)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        if (isTRUE(use_label_boxes)) {
          p <- p + ggrepel::geom_label_repel(
            data = lab_tbl,
            aes(x = t_opt, y = mean_R, label = GENE, color = Risk),
            size = label_size, fontface = label_fontface,
            min.segment.length = 0, box.padding = 0.15, point.padding = 0.25,
            label.size = 0.2, show.legend = FALSE
          )
        } else {
          p <- p + ggrepel::geom_text_repel(
            data = lab_tbl,
            aes(x = t_opt, y = mean_R, label = GENE, color = Risk),
            size = label_size, fontface = label_fontface,
            min.segment.length = 0, point.padding = 0.25,
            show.legend = FALSE
          )
        }
      } else {
        # fallback without ggrepel
        p <- p + geom_text(
          data = lab_tbl,
          aes(x = t_opt, y = mean_R, label = GENE, color = Risk),
          size = label_size, fontface = label_fontface, show.legend = FALSE, vjust = -0.6
        )
      }
    }
    
    p
  }
  
  # Example call (only change is optional label_size if you want even larger text)
  ncol_facets <- if (exists("gene_column_number")) gene_column_number else 3
  viz_subclone <- visualize_high_risk_subclones(
    gene_simulations     = gene_simulations,
    flagged_subclones_df = flagged_subclones_df,
    t_opt_list           = t_opt_list,
    alpha_band           = 0.05,
    max_labels           = 5,
    ncol                 = ncol_facets,
    label_size           = 5.5,     # tweak here
    label_fontface       = "bold",
    use_label_boxes      = TRUE
  )
  
  #@@@print(viz_subclone)
  # vector of unique high-risk genes (handles possible NA in high_risk)
  ## keep original order (duplicates kept)
  high_risk_genes_in_order <- flagged_subclones_df$GENE[
    flagged_subclones_df$high_risk %in% TRUE
  ]
  
  ## keep original order but drop duplicates (first occurrence only)
  high_risk_genes_in_order_unique <- high_risk_genes_in_order[!duplicated(high_risk_genes_in_order)]
  
  #@@@high_risk_genes_in_order
  #@@@high_risk_genes_in_order_unique
  
  
  return(list(cme_dynamics = cme_dynamics,
              #combined_p14 = combined_p14,
              gene_simulations = gene_simulations,
              subclone_labels = subclone_labels,
              t_opt_list = t_opt_list,
              cme_overlays_faceted = cme_overlays_faceted,
              p_cov = p_cov,
              p_crps = p_crps,
              p14_scatter = p14_scatter,
              p14_varmean = p14_varmean,
              p14_lollipop = p14_lollipop,
              filtered_subclonal_genes = filtered_subclonal_genes,
              filteredSubclonalGenes = filteredSubclonalGenes,
              Heatmap_tumorHeterogeneity = Heatmap_tumorHeterogeneity,
              Barplot_tumorHeterogeneity = Barplot_tumorHeterogeneity,
              viz_subclone = viz_subclone,
              high_risk_genes_in_order_unique = high_risk_genes_in_order_unique))
}
