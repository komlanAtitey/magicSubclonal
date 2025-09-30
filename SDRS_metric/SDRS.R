library(stats)

# === OA-only clinical scorer with survival OR proxy outcomes ===
compute_SDRS_OA <- function(
    X,                         # genes x samples matrix (rownames = gene symbols)
    drivers,                   # driver gene symbols
    S,                         # tool's predicted genes (ranked; we'll use top-k after filtering)
    outcome,                   # list describing the outcome; see below
    k                 = NULL,  # evaluate top-k genes (default = all available after filtering)
    pathways          = NULL,  # if NULL, will call build_pathways_for_drivers(drivers)
    filter_to_driver_pathways = TRUE,
    directional       = FALSE, # if TRUE: penalize protective (β<0) / reward adverse (β>0)
    cap_decibels      = 10,    # cap for -log10(p) mapping to [0,1]
    covariates_df     = NULL,  # optional (nrow = ncol(X)) clinical covariates
    seed              = 1
){
  set.seed(seed)
  
  # --- helpers ---
  in_universe <- function(g) g[!is.na(g) & nzchar(g)]
  safe_intersect <- function(a, b) unique(intersect(a, b))
  row_zscore <- function(M){
    M <- as.matrix(M)
    mu <- rowMeans(M, na.rm = TRUE)
    sdv <- apply(M, 1, sd, na.rm = TRUE)
    sdv[sdv == 0 | is.na(sdv)] <- 1
    sweep(sweep(M, 1, mu, "-"), 1, sdv, "/")
  }
  scale_p_to_unit <- function(p, cap = 10){
    if (is.na(p)) return(NA_real_)
    val <- -log10(max(p, .Machine$double.xmin))
    pmin(1, val / cap)
  }
  
  # --- inputs ---
  X <- as.matrix(X)
  if (is.null(rownames(X))) stop("X must have rownames (gene symbols).")
  G <- rownames(X); N <- ncol(X)
  
  drivers <- unique(in_universe(drivers))
  S       <- unique(in_universe(S))
  S_in_X  <- safe_intersect(S, G)
  
  # --- pathways from drivers (used to prioritize genes) ---
  if (is.null(pathways)) {
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      pathways <- tryCatch(
        build_pathways_for_drivers(drivers),
        error = function(e) { warning("build_pathways_for_drivers() failed: ", conditionMessage(e)); list() }
      )
    } else {
      warning("msigdbr not installed; proceeding without pathway filtering.")
      pathways <- list()
    }
  }
  T_set <- character(0)
  if (length(pathways) > 0) {
    T_set <- unique(unlist(pathways, use.names = FALSE))
    T_set <- safe_intersect(T_set, G)
  }
  
  # --- choose the set to score (prefer driver-pathway overlap if available) ---
  if (isTRUE(filter_to_driver_pathways) && length(T_set) > 0) {
    S_pref <- safe_intersect(S_in_X, T_set)
    S_use  <- if (length(S_pref) > 0) S_pref else S_in_X
  } else {
    S_use <- S_in_X
  }
  if (length(S_use) == 0) {
    warning("No genes from S found in X (or after pathway filtering).")
    return(list(
      SDRS = NA_real_,
      subscores = list(OA = list(score = NA_real_, p_value = NA_real_, effect = NA_real_, effect_name = NA_character_)),
      details = list(k_effective = 0, n_genes_in_X = length(S_in_X), in_driver_pathways = 0)
    ))
  }
  if (is.null(k)) k <- length(S_use)
  S_use <- head(S_use, k)
  
  # --- per-sample gene-set score ---
  Zx    <- row_zscore(X)
  score <- colMeans(Zx[S_use, , drop = FALSE], na.rm = TRUE)
  
  # --- assemble modeling frame ---
  df <- data.frame(score = as.numeric(score))
  if (!is.null(covariates_df)) {
    if (nrow(covariates_df) != N) stop("covariates_df must have nrow = ncol(X).")
    df <- cbind(df, covariates_df)
  }
  
  # --- outcome dispatcher ---
  type <- outcome$type
  s_OA <- NA_real_; pval <- NA_real_; effect <- NA_real_; effect_name <- NA_character_
  mapping_used <- NA_character_
  
  if (identical(type, "survival")) {
    if (!requireNamespace("survival", quietly = TRUE))
      stop('Outcome type "survival" requires the survival package.')
    time   <- outcome$time
    status <- outcome$status
    if (length(time) != N || length(status) != N) stop("time/status length must match ncol(X).")
    df$time <- time; df$status <- status
    
    fmla <- as.formula(paste0("survival::Surv(time, status) ~ score",
                              if (!is.null(covariates_df)) " + ." else ""))
    fit <- tryCatch(survival::coxph(fmla, data = df, ties = "efron"), error = function(e) NULL)
    if (!is.null(fit)) {
      s   <- summary(fit)
      beta <- unname(coef(fit)["score"])
      pval <- suppressWarnings(s$coefficients["score","Pr(>|z|)"])
      effect <- unname(exp(beta)); effect_name <- "HR"
      s_OA <- scale_p_to_unit(pval, cap = cap_decibels)
      mapping_used <- "-log10(p)→[0,1]"
      if (isTRUE(outcome$directional) || isTRUE(directional)) {
        s_OA <- (1 + sign(beta) * s_OA) / 2
      }
    }
    
  } else if (identical(type, "binary")) {
    y <- outcome$y
    if (length(y) != N) stop("Binary outcome y must have length ncol(X).")
    # coerce to 0/1
    if (is.factor(y)) y <- as.numeric(y) - 1
    y <- as.integer(y)
    if (!all(y %in% c(0L,1L))) stop("Binary y must be 0/1 or a 2-level factor.")
    df$y <- y
    
    fmla <- as.formula(paste0("y ~ score", if (!is.null(covariates_df)) " + ." else ""))
    fit <- tryCatch(glm(fmla, data = df, family = binomial()), error = function(e) NULL)
    if (!is.null(fit)) {
      co <- summary(fit)$coefficients
      beta <- unname(co["score","Estimate"])
      se   <- unname(co["score","Std. Error"])
      zval <- beta / se
      pval <- 2 * (1 - pnorm(abs(zval)))
      effect <- unname(exp(beta)); effect_name <- "OR"
      s_OA <- scale_p_to_unit(pval, cap = cap_decibels)  # or map AUC if you prefer
      mapping_used <- "-log10(p)→[0,1]"
      if (isTRUE(outcome$directional) || isTRUE(directional)) {
        s_OA <- (1 + sign(beta) * s_OA) / 2
      }
    }
    
  } else if (identical(type, "continuous")) {
    y <- outcome$y
    if (length(y) != N) stop("Continuous outcome y must have length ncol(X).")
    df$y <- as.numeric(y)
    
    fmla <- as.formula(paste0("y ~ score", if (!is.null(covariates_df)) " + ." else ""))
    fit <- tryCatch(lm(fmla, data = df), error = function(e) NULL)
    if (!is.null(fit)) {
      co <- summary(fit)$coefficients
      beta <- unname(co["score","Estimate"])
      pval <- unname(co["score","Pr(>|t|)"])
      effect <- beta; effect_name <- "β (per 1-SD score)"
      s_OA <- scale_p_to_unit(pval, cap = cap_decibels)
      mapping_used <- "-log10(p)→[0,1]"
      if (isTRUE(outcome$directional) || isTRUE(directional)) {
        s_OA <- (1 + sign(beta) * s_OA) / 2
      }
    }
    
  } else {
    stop('Provide outcome$type as one of: "binary", "continuous", or "survival".')
  }
  
  path_hits <- character(0)
  if (length(pathways) > 0) {
    path_hits <- names(pathways)[vapply(pathways, function(gs) any(gs %in% S_use), logical(1))]
  }
  
  list(
    SDRS = s_OA,   # clinical SDRS (OA-only)
    subscores = list(
      OA = list(score = s_OA, p_value = pval, effect = effect, effect_name = effect_name)
    ),
    details = list(
      k_effective = length(S_use),
      genes_used = S_use,
      n_genes_in_X = length(S_in_X),
      n_genes_in_driver_pathways = length(safe_intersect(S_in_X, T_set)),
      driver_pathway_genes_total = length(T_set),
      pathways_hit_by_S = path_hits,
      mapping = mapping_used,
      directional = isTRUE(outcome$directional) || isTRUE(directional),
      outcome_type = type
    )
  )
}

build_pathways_for_drivers <- function(
    drivers,
    species  = "Homo sapiens",
    keep_cat = c("H","C2","C5"),
    keep_sub = c("", "CP:KEGG", "CP:REACTOME")
){
  msig <- msigdbr::msigdbr(species = species)
  nm <- names(msig)
  
  col_cat  <- if ("gs_cat" %in% nm) "gs_cat" else if ("gs_collection" %in% nm) "gs_collection" else stop("No category column")
  col_sub  <- if ("gs_subcat" %in% nm) "gs_subcat" else if ("gs_subcollection" %in% nm) "gs_subcollection" else NA_character_
  col_name <- if ("gs_name" %in% nm) "gs_name" else stop("No pathway-name column")
  col_gene <- if ("gene_symbol" %in% nm) "gene_symbol" else stop("No gene-symbol column")
  
  if (is.na(col_sub)) { msig$..subcat <- NA_character_; col_sub <- "..subcat" }
  
  keep <- (msig[[col_cat]] %in% keep_cat) & (is.na(msig[[col_sub]]) | msig[[col_sub]] %in% keep_sub)
  msig2 <- msig[keep, , drop = FALSE]
  
  pathways <- lapply(split(msig2[[col_gene]], msig2[[col_name]]), unique)
  pathways_driver <- pathways[vapply(pathways, function(gs) any(gs %in% drivers), logical(1))]
  
  # If none found, widen once: include all subcollections within the kept categories
  if (length(pathways_driver) == 0) {
    keep <- (msig[[col_cat]] %in% keep_cat)
    msig2 <- msig[keep, , drop = FALSE]
    pathways <- lapply(split(msig2[[col_gene]], msig2[[col_name]]), unique)
    pathways_driver <- pathways[vapply(pathways, function(gs) any(gs %in% drivers), logical(1))]
  }
  pathways_driver
}
