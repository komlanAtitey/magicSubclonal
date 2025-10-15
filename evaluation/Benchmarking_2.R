setwd("/Users/atiteyk2/Documents/magicSubclonal") 
getwd()


#################################################################################################
#################################################################################################
#################################################################################################
## =========================================================
## 0) Libraries
## =========================================================
suppressWarnings({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(mclust)    # EM for Gaussians
  library(nnls)      # nonnegative least squares
  library(precrec)   # ROC/PR metrics
})

set.seed(1)

## =========================================================
## 1) Simulate bulk RNA-seq with K subclones + markers
## =========================================================
simulate_bulk <- function(G = 400, M_per_clone = 12, K = 3, N = 240,
                          baseline_mu = 5, marker_lfc = 1.2, noise_sd = 0.20) {
  # True subclone signatures (genes x K) on linear scale
  S0 <- matrix(rexp(G, rate = 1 / baseline_mu), nrow = G, ncol = K)
  markers <- vector("list", K)
  for (k in 1:K) {
    idx <- ((k - 1) * M_per_clone + 1):(k * M_per_clone)
    markers[[k]] <- idx
    S0[idx, k] <- S0[idx, k] * exp(marker_lfc)  # up-weight clone-specific markers
  }
  rownames(S0) <- paste0("g", seq_len(G))
  colnames(S0) <- paste0("clone", seq_len(K))
  marker_map <- tibble(
    gene  = rownames(S0)[unique(unlist(markers))],
    clone = rep(paste0("clone", seq_len(K)), each = M_per_clone)
  )
  
  # Sample mixture weights W (K x N)
  rdirichlet <- function(n, alpha) {
    out <- matrix(NA_real_, n, length(alpha))
    for (i in 1:n) {
      g <- rgamma(length(alpha), shape = alpha, rate = 1)
      out[i, ] <- g / sum(g)
    }
    out
  }
  W <- t(rdirichlet(N, alpha = rep(2, K)))
  rownames(W) <- paste0("clone", seq_len(K))
  colnames(W) <- paste0("s", seq_len(N))
  
  # Expression with multiplicative log-normal noise
  X_mean <- S0 %*% W                    # (G x N)
  Eps    <- matrix(rnorm(G * N, 0, noise_sd), G, N)
  X      <- X_mean * exp(Eps)
  rownames(X) <- rownames(S0)
  colnames(X) <- colnames(W)
  
  # Ground truth labels for (sample, gene): 1 if gene's clone equals dominant clone in sample
  dominant_idx <- apply(W, 2, which.max)               # length N
  names(dominant_idx) <- colnames(W)
  # gene -> clone index (1..K)
  markers_lookup <- setNames(as.integer(sub("clone", "", marker_map$clone)), marker_map$gene)
  
  # Build (sample, gene) frame only for marker genes
  mgenes <- marker_map$gene
  grid   <- expand.grid(sample_id = colnames(X), gene = mgenes, stringsAsFactors = FALSE)
  grid$label <- as.integer(markers_lookup[grid$gene] == dominant_idx[grid$sample_id])
  
  list(X = X, W_true = W, S_true = S0,
       markers = marker_map, gt = as_tibble(grid), K = K,
       markers_lookup = markers_lookup)
}

sim <- simulate_bulk(G = 400, M_per_clone = 12, K = 3, N = 240)
X <- sim$X; K <- sim$K; markers <- sim$markers; gt <- sim$gt
markers_lookup <- sim$markers_lookup
X_log <- log1p(X)

## =========================================================
## 2) Folds
## =========================================================
make_folds <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  idx <- sample(seq_len(n))
  split(idx, rep(1:K, length.out = n))
}
Kfold <- 5
folds <- make_folds(ncol(X_log), K = Kfold)

## =========================================================
## 3) Utilities
## =========================================================
# Isotonic regression calibrator: returns f(x) clamped to [0,1]
iso_fit <- function(scores, labels) {
  df <- data.frame(x = as.numeric(scores), y = as.numeric(labels))
  df <- df[is.finite(df$x) & !is.na(df$y), , drop = FALSE]
  if (nrow(df) < 2) return(function(x) rep(0.5, length(x)))
  o   <- order(df$x); x <- df$x[o]; y <- df$y[o]
  agg <- aggregate(y, by = list(x), FUN = mean)
  xs  <- as.numeric(agg$Group.1)
  ys  <- as.numeric(agg$x)
  iso <- stats::isoreg(xs, ys)
  f   <- stats::approxfun(iso$x, iso$yf, method = "linear", rule = 2)
  function(x) pmin(pmax(f(as.numeric(x)), 0), 1)
}

# Safe extractor aligning a matrix to a (sample_id, gene) tibble
.extract_scores <- function(pred_mat, pair_df) {
  stopifnot(!is.null(rownames(pred_mat)), !is.null(colnames(pred_mat)))
  ridx <- match(pair_df$sample_id, rownames(pred_mat))
  cidx <- match(pair_df$gene,      colnames(pred_mat))
  as.numeric(pred_mat[cbind(ridx, cidx)])
}

# Normalize binary labels to integer {0,1}
.norm_labels01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  x <- suppressWarnings(as.integer(x))
  if (!all(x %in% c(0L, 1L), na.rm = TRUE)) stop("Labels must be {0,1}.")
  x
}

# Scale vector to [0,1] (unused here, left for convenience)
.scale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

## =========================================================
## 4) Predictors
## =========================================================
# 4a) Improved "magicSubclonal": per-gene EM (2-Gaussian) + mean-shrinkage
magic_em_predict <- function(X_train_log, X_test_log,
                             shrink_lambda = 0.10,
                             models = c("V","E"),
                             eps_sd = 1e-3) {
  G   <- nrow(X_train_log)
  out <- matrix(NA_real_, nrow = ncol(X_test_log), ncol = G,
                dimnames = list(colnames(X_test_log), rownames(X_train_log)))
  for (g in seq_len(G)) {
    xtr <- as.numeric(X_train_log[g, ])
    fit <- try(Mclust(xtr, G = 2, modelNames = models, verbose = FALSE), silent = TRUE)
    if (!inherits(fit, "try-error") && !is.null(fit$parameters$mean)) {
      mu <- as.numeric(fit$parameters$mean)
      sg <- sqrt(as.numeric(fit$parameters$variance$sigmasq))
      pi <- as.numeric(fit$parameters$pro)
    } else {
      km <- kmeans(xtr, centers = 2, nstart = 5)
      mu <- tapply(xtr, km$cluster, mean)
      sg <- tapply(xtr, km$cluster, sd)
      pi <- as.numeric(table(km$cluster)) / length(xtr)
    }
    # shrink means toward gene mean
    m0 <- mean(xtr); mu <- (1 - shrink_lambda) * mu + shrink_lambda * m0
    sg <- pmax(sg, eps_sd)
    ord <- order(mu); mu <- mu[ord]; sg <- sg[ord]; pi <- pi[ord]
    
    xte <- as.numeric(X_test_log[g, ])
    ll1 <- dnorm(xte, mean = mu[1], sd = sg[1], log = TRUE) + log(pi[1])
    ll2 <- dnorm(xte, mean = mu[2], sd = sg[2], log = TRUE) + log(pi[2])
    out[, g] <- plogis(ll2 - ll1)  # posterior of "high" component
  }
  out
}

# 4b) ss-Deconv predictor: per-sample NNLS to get clone weights, map to genes via markers
#     Returns (N_test x G) scores: for gene g, score = weight of its clone (fallback: mean weight)
ss_deconv_predict <- function(S_hat, X_test, markers_lookup) {
  stopifnot(nrow(S_hat) == nrow(X_test))
  K <- ncol(S_hat)
  N <- ncol(X_test)
  G <- nrow(X_test)
  out <- matrix(NA_real_, nrow = N, ncol = G,
                dimnames = list(colnames(X_test), rownames(X_test)))
  for (j in seq_len(N)) {
    fit <- nnls::nnls(S_hat, X_test[, j])
    w   <- as.numeric(coef(fit)); w <- pmax(w, 0); s <- sum(w); if (s > 0) w <- w / s
    k_idx <- markers_lookup[colnames(out)]   # gene -> clone index
    sc <- rep(mean(w), G)                    # default fallback = mean weight
    sel <- !is.na(k_idx)
    sc[sel] <- w[k_idx[sel]]
    out[j, ] <- sc
  }
  out
}

# 4c) Mixture model on samples (train PCA -> GMM with K comps); gene prob = resp of gene's clone
mm_fit_predict <- function(Xtr_log, Xte_log, K) {
  pr  <- prcomp(t(Xtr_log), center = TRUE, scale. = FALSE)
  q   <- max(1, min(10, ncol(pr$x)))
  Ztr <- pr$x[, 1:q, drop = FALSE]
  fit <- Mclust(Ztr, G = K, verbose = FALSE)
  Zte <- scale(t(Xte_log), center = pr$center, scale = FALSE) %*% pr$rotation[, 1:q, drop = FALSE]
  pr_tr <- predict(fit, newdata = Ztr)$z
  pr_te <- predict(fit, newdata = Zte)$z
  list(model = fit, resp_train = pr_tr, resp_test = pr_te)
}

# Map responsibilities to gene-level probabilities using gene->clone lookup
mm_marker_prob <- function(resp_mat, model_means_unused, keep_unused,
                           gene_id, all_genes, markers_lookup = NULL) {
  if (is.null(markers_lookup)) return(rep(0.5, nrow(resp_mat)))
  k0 <- markers_lookup[gene_id]
  if (is.na(k0)) return(rep(0.5, nrow(resp_mat)))
  as.numeric(resp_mat[, k0])
}

## =========================================================
## 5) 5-fold CV: build dictionary, fit, predict, calibrate
## =========================================================
pool <- list(ss = list(), magic = list(), mm = list())

for (f in seq_along(folds)) {
  message(sprintf("Fold %d / %d", f, length(folds)))
  te_idx  <- folds[[f]]
  tr_idx  <- setdiff(seq_len(ncol(X_log)), te_idx)
  
  Xtr     <- X[, tr_idx, drop = FALSE]
  Xte     <- X[, te_idx, drop = FALSE]
  Xtr_log <- X_log[, tr_idx, drop = FALSE]
  Xte_log <- X_log[, te_idx, drop = FALSE]
  
  # ----- Build a crude dictionary S_hat from training -----
  # SVD proxy to get nonnegative "Htr" (Ntr x K), then S_hat = (Xtr %*% Htr) / colSums(Htr)
  sv  <- svd(scale(t(Xtr), center = TRUE, scale = FALSE))   # on samples
  Kp  <- K
  Htr <- sv$u[, 1:Kp, drop = FALSE]                         # Ntr x K
  Htr <- pmax(Htr, 0)                                       # nonnegative
  Htr <- sweep(Htr, 2, pmax(colSums(Htr), 1e-8), "/")       # normalize columns
  SH  <- Xtr %*% Htr                                        # (G x K)
  denom <- pmax(colSums(Htr), 1e-8)                         # length K
  S_hat <- sweep(SH, 2, denom, "/")                         # <-- FIX: column-wise divide
  S_hat[S_hat < 0] <- 0
  colnames(S_hat) <- paste0("clone", seq_len(K))
  
  # ----- Predictions (TEST) -----
  pred_ss <- ss_deconv_predict(S_hat, Xte, markers_lookup)  # (Nte x G)
  pred_ms <- magic_em_predict(Xtr_log, Xte_log)             # (Nte x G)
  
  mm <- mm_fit_predict(Xtr_log, Xte_log, K = K)
  all_genes <- rownames(X)
  pred_mm <- matrix(0, nrow = ncol(Xte_log), ncol = nrow(Xte_log),
                    dimnames = list(colnames(Xte_log), rownames(Xte_log)))
  for (g in seq_along(all_genes)) {
    pred_mm[, g] <- mm_marker_prob(mm$resp_test, NULL, NULL,
                                   gene_id = all_genes[g],
                                   all_genes = all_genes,
                                   markers_lookup = markers_lookup)
  }
  
  # ----- Training labels/pairs for calibration -----
  mgenes <- intersect(markers$gene, rownames(X))
  gt_tr <- gt %>% filter(sample_id %in% colnames(X)[tr_idx], gene %in% mgenes) %>%
    arrange(sample_id, gene)
  gt_te <- gt %>% filter(sample_id %in% colnames(X)[te_idx], gene %in% mgenes) %>%
    arrange(sample_id, gene)
  if (nrow(gt_tr) == 0 || nrow(gt_te) == 0) next
  
  # ----- TRAIN predictions for calibration (aligned to gt_tr) -----
  pred_ss_tr <- ss_deconv_predict(S_hat, Xtr, markers_lookup)
  pred_ms_tr <- magic_em_predict(Xtr_log, Xtr_log)
  
  pr_tr <- mm$resp_train
  pred_mm_tr <- matrix(0, nrow = ncol(Xtr_log), ncol = nrow(Xtr_log),
                       dimnames = list(colnames(Xtr_log), rownames(Xtr_log)))
  for (g in seq_along(all_genes)) {
    pred_mm_tr[, g] <- mm_marker_prob(pr_tr, NULL, NULL,
                                      gene_id = all_genes[g],
                                      all_genes = all_genes,
                                      markers_lookup = markers_lookup)
  }
  
  # ----- Fit isotonic calibrators per method -----
  iso_ss <- iso_fit(.extract_scores(pred_ss_tr[, mgenes, drop = FALSE], gt_tr), gt_tr$label)
  iso_ms <- iso_fit(.extract_scores(pred_ms_tr[, mgenes, drop = FALSE], gt_tr), gt_tr$label)
  iso_mm <- iso_fit(.extract_scores(pred_mm_tr[, mgenes, drop = FALSE], gt_tr), gt_tr$label)
  
  # ----- Apply to TEST (aligned to gt_te) -----
  score_ss_te <- .extract_scores(pred_ss[, mgenes, drop = FALSE], gt_te)
  score_ms_te <- .extract_scores(pred_ms[, mgenes, drop = FALSE], gt_te)
  score_mm_te <- .extract_scores(pred_mm[, mgenes, drop = FALSE], gt_te)
  
  cal_ss_te <- iso_ss(score_ss_te)
  cal_ms_te <- iso_ms(score_ms_te)
  cal_mm_te <- iso_mm(score_mm_te)
  
  # ----- Assemble OOF rows for this fold -----
  te_df <- gt_te %>% mutate(
    score_ss = as.numeric(cal_ss_te),
    score_ms = as.numeric(cal_ms_te),
    score_mm = as.numeric(cal_mm_te)
  )
  
  # Long format per method
  pool$ss[[f]] <- te_df %>%
    transmute(sample_id, gene, label, score = score_ss, method = "ss-Deconv")
  
  pool$magic[[f]] <- te_df %>%
    transmute(sample_id, gene, label, score = score_ms, method = "magicSubclonal(EM)")
  
  pool$mm[[f]] <- te_df %>%
    transmute(sample_id, gene, label, score = score_mm, method = "MixtureModel")
}

## =========================================================
## 6) Combine predictions (robustly) and metrics
## =========================================================
preds_raw <- bind_rows(c(pool$ss, pool$magic, pool$mm))
if ("gene" %in% names(preds_raw)) {
  preds <- preds_raw %>%
    dplyr::rename(marker = gene) %>%
    dplyr::mutate(label = .norm_labels01(label))
} else if ("marker" %in% names(preds_raw)) {
  preds <- preds_raw %>%
    dplyr::mutate(label = .norm_labels01(label))
} else {
  stop("Neither 'gene' nor 'marker' column found in preds_raw.")
}


# ---- Metrics per method ----

# ---- Define helper function ----
eval_by_method <- function(df) {
  df <- df %>%
    dplyr::filter(is.finite(score), !is.na(label)) %>%
    dplyr::mutate(label = .norm_labels01(label))
  
  if (nrow(df) < 2 || length(unique(df$label)) < 2) {
    return(tibble(
      AUC_ROC = NA_real_,
      AUC_PRC = NA_real_,
      Brier = NA_real_,
      Calib_Slope = NA_real_,
      Calib_Intercept = NA_real_
    ))
  }
  
  # Compute ROC/PRC AUCs
  mm <- precrec::mmdata(
    scores = list(df$score),
    labels = list(df$label),
    modnames = unique(df$method),
    posclass = 1L
  )
  
  ev <- precrec::evalmod(mm)
  au <- precrec::auc(ev) %>%
    dplyr::select(modnames, curvetypes, aucs) %>%
    tidyr::pivot_wider(
      names_from = curvetypes,
      values_from = aucs,
      names_prefix = "AUC_"
    )
  
  # Compute Brier score
  brier <- mean((df$score - df$label)^2)
  
  # Compute calibration intercept and slope
  s  <- pmin(pmax(df$score, 1e-6), 1 - 1e-6)
  lo <- qlogis(s)
  
  cal <- tryCatch({
    fit <- glm(df$label ~ lo, family = binomial())
    co  <- coef(summary(fit))
    tibble(
      Calib_Intercept = unname(co[1, 1]),
      Calib_Slope     = unname(co[2, 1])
    )
  }, error = function(e) {
    tibble(Calib_Intercept = NA_real_, Calib_Slope = NA_real_)
  })
  
  dplyr::bind_cols(
    au %>% dplyr::select(dplyr::starts_with("AUC_")),
    tibble(Brier = brier),
    cal
  )
}

# ---- Compute leaderboard ----
leaderboard <- preds %>%
  dplyr::group_split(method) %>%
  purrr::map_dfr(~ eval_by_method(.x) %>%
                   dplyr::mutate(method = unique(.x$method))) %>%
  dplyr::select(method, AUC_ROC, AUC_PRC, Brier, Calib_Slope, Calib_Intercept) %>%
  dplyr::arrange(dplyr::desc(AUC_ROC))

# ---- Inspect output ----
print(leaderboard)


## =========================================================
## 7) ROC / PR plots across methods
## =========================================================
## =========================================================
## ROC / PR with consistent labels across methods (precrec)
## =========================================================

# Helper to enforce 0/1 integer labels
.norm_labels01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  x <- suppressWarnings(as.integer(x))
  if (!all(x %in% c(0L, 1L), na.rm = TRUE)) stop("Labels must be {0,1}.")
  x
}

# 1) Clean and align: one row per (sample_id, marker) with all methods as columns
methods_present <- sort(unique(preds$method))

wide <- preds %>%
  dplyr::filter(is.finite(score)) %>%
  dplyr::mutate(label = .norm_labels01(label),
                method = as.character(method)) %>%
  dplyr::select(sample_id, marker, label, method, score) %>%
  tidyr::pivot_wider(names_from = method, values_from = score) %>%
  # keep only rows where every method has a score
  tidyr::drop_na(dplyr::all_of(methods_present)) %>%
  dplyr::arrange(sample_id, marker)

if (nrow(wide) == 0L) stop("No common (sample, marker) rows found across all methods after alignment.")

# 2) Build lists of aligned score vectors (one numeric vector per method)
label_vec    <- .norm_labels01(wide$label)
scores_list  <- lapply(methods_present, function(m) as.numeric(wide[[m]]))
names(scores_list) <- methods_present

# 3) Feed to precrec — single dataset (dsid = 1), multiple models
mm_all <- precrec::mmdata(
  scores   = scores_list,
  labels   = label_vec,
  modnames = methods_present,
  posclass = 1L
)
ev_all <- precrec::evalmod(mm_all)


######################################

# Extract data from precrec object
roc_df <- fortify(ev_all, curvetypes = "ROC")
prc_df <- fortify(ev_all, curvetypes = "PRC")

# ROC plot with thicker lines
p_roc <- ggplot(roc_df, aes(x = 1 - x, y = y, color = modname)) +
  geom_line(linewidth = 1.8) +  # Increase line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray70") +
  theme_minimal(base_size = 13) +
  labs(
    title = "ROC Curves by Method (Aligned Labels)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Method"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# PRC plot with thicker lines
p_pr <- ggplot(prc_df, aes(x = x, y = y, color = modname)) +
  geom_line(linewidth = 1.8) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Precision–Recall Curves by Method (Aligned Labels)",
    x = "Recall",
    y = "Precision",
    color = "Method"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# Display
p_roc
p_pr

library(precrec)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- 1. Compute AUC table ----
auc_tbl <- precrec::auc(ev_all) %>%
  dplyr::select(modnames, curvetypes, aucs) %>%
  tidyr::pivot_wider(
    names_from = curvetypes,
    values_from = aucs,
    names_prefix = "AUC_"
  ) %>%
  dplyr::rename(method = modnames) %>%
  dplyr::arrange(dplyr::desc(AUC_ROC))

print(auc_tbl)

# ---- 2. Prepare data for plotting ----
auc_long <- auc_tbl %>%
  tidyr::pivot_longer(cols = starts_with("AUC_"),
                      names_to = "Metric",
                      values_to = "AUC")

# ---- 3. Plot ----
p_auc <- ggplot(auc_long, aes(x = reorder(method, AUC), y = AUC, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(AUC, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "AUC (ROC & PRC) by Method",
    x = "Method",
    y = "AUC Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  )

# ---- 4. Show the plot ----
p_auc


