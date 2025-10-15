setwd("/Users/atiteyk2/Documents/Master_Equation")
getwd()

# ============================================================
# magicSubclonal — Microenvironment-aware checks (CLEAN SCRIPT)
# 1) Stability (raw vs residualized)
# 2) Microenvironment surrogacy checks
# 3) Purity-stratified clinical consistency
# Also: MCP-counter fractions + ESTIMATE purity (via immunedeconv)
# ============================================================

# ---------- Packages ----------
need_cran <- c(
  "dplyr","tibble","stringr","MCPcounter","immunedeconv",
  "ppcor","survival","broom","Hmisc"
)
to_install <- setdiff(need_cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(MCPcounter)
  library(immunedeconv)
  library(ppcor)
  library(survival)
  library(broom)
  library(Hmisc)   # cut2()
})

# ---------- Load your data ----------
# Expect: data_GSE9891_Research_6 with a GeneSymbol column and expression columns

load("data/data_GSE9891_Research_6.rdata")
stopifnot(exists("data_GSE9891_Research_6"))
expr_df <- as.data.frame(data_GSE9891_Research_6, stringsAsFactors = FALSE)

load("data/data_GSE171415_Research_2.rdata")
stopifnot(exists("data_GSE171415_Research_2"))
expr_df <- as.data.frame(data_GSE171415_Research_2, stringsAsFactors = FALSE)

load("data/data_GSE21422_Research_7.rdata")
stopifnot(exists("data_GSE21422_Research_7"))
expr_df <- as.data.frame(data_GSE21422_Research_7, stringsAsFactors = FALSE)

# ---------- Build genes × samples matrix (HGNC symbols) ----------
# Ensure GeneSymbol is a column
if (!"GeneSymbol" %in% names(expr_df)) {
  expr_df <- tibble::rownames_to_column(expr_df, "GeneSymbol")
}

# Keep gene + numeric columns
num_cols <- names(expr_df)[vapply(expr_df, is.numeric, logical(1))]
stopifnot(length(num_cols) > 0)
expr_df <- expr_df[, c("GeneSymbol", num_cols), drop = FALSE]

# Clean symbols (uppercase, strip multi-maps like "A1BG /// A1BG-AS1" -> "A1BG")
clean_symbol_cell <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  x <- toupper(trimws(x))
  # split by common delimiters: /// ; , |
  parts <- unlist(strsplit(x, "\\s*(///|;|,|\\|)\\s*"))
  parts <- parts[parts != ""]
  if (length(parts) == 0) return(NA_character_)
  parts[1]
}
expr_df$GeneSymbol <- vapply(expr_df$GeneSymbol, clean_symbol_cell, character(1))
expr_df <- expr_df %>% filter(!is.na(GeneSymbol), GeneSymbol != "")

# Collapse duplicates by highest variance across numeric columns
expr_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  mutate(.var = apply(as.matrix(pick(where(is.numeric))), 1, var, na.rm = TRUE)) %>%
  slice_max(order_by = .var, n = 1, with_ties = FALSE) %>%
  ungroup()
expr_df$.var <- NULL

# Matrix
expr_mat <- as.matrix(expr_df[, num_cols, drop = FALSE])
rownames(expr_mat) <- expr_df$GeneSymbol
mode(expr_mat) <- "numeric"
expr_mat[is.na(expr_mat)] <- 0

samples <- colnames(expr_mat)
genes   <- rownames(expr_mat)
message(sprintf("Expression matrix: %d genes x %d samples", nrow(expr_mat), ncol(expr_mat)))

# ---------- MCP-counter: Immune / Stromal / Endothelial fractions ----------
mcp_scores <- MCPcounter.estimate(expr_mat, featuresType = "HUGO_symbols")

immune_pops <- c("T cells","CD8 T cells","Cytotoxic lymphocytes","NK cells",
                 "B lineage","Monocytic lineage","Myeloid dendritic cells","Neutrophils")
req_pops <- c(immune_pops, "Endothelial cells", "Fibroblasts")

# Resolve function conflicts up-front
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::setdiff)
  conflicted::conflicts_prefer(base::union)
  conflicted::conflicts_prefer(base::intersect)
}


missing_pops <- setdiff(req_pops, rownames(mcp_scores))
if (length(missing_pops)) stop("Missing MCP-counter populations: ", paste(missing_pops, collapse=", "))

immune_signal      <- colSums(mcp_scores[immune_pops, , drop = FALSE], na.rm = TRUE)
stromal_signal     <- mcp_scores["Fibroblasts", ]
endothelial_signal <- mcp_scores["Endothelial cells", ]

to_pos <- function(x) x - min(x, na.rm = TRUE)
num_df <- data.frame(
  Immune      = to_pos(immune_signal),
  Stromal     = to_pos(stromal_signal),
  Endothelial = to_pos(endothelial_signal),
  row.names   = samples
)
den <- rowSums(num_df); den[den == 0] <- NA
fractions_df <- as.data.frame(num_df / den)
fractions_df[] <- lapply(fractions_df, function(v) pmin(pmax(v, 0), 1))
stopifnot(identical(rownames(fractions_df), samples))

# ---------- ESTIMATE (via immunedeconv) → Purity ----------
expr_df_for_id <- data.frame(GeneSymbol = rownames(expr_mat), expr_mat, check.names = FALSE)
est_tbl <- immunedeconv::deconvolute(expr_df_for_id, method = "estimate", tumor = TRUE)

# Robust fetch of ESTIMATEScore regardless of row labels
get_ESTIMATEScore <- function(mat) {
  if (is.data.frame(mat) && "cell_type" %in% names(mat)) {
    rn <- mat$cell_type
    M  <- as.matrix(mat[, base::setdiff(names(mat), "cell_type"), drop = FALSE])
    rownames(M) <- rn
  } else {
    M <- as.matrix(mat)
  }
  norm_key <- function(x) tolower(gsub("[^a-z]", "", x))
  rk <- norm_key(rownames(M))
  is_est <- rk %in% c("estimatescore","estimate","estimatescoretcga")
  is_imm <- rk %in% c("immunescore","immune")
  is_str <- rk %in% c("stromalscore","stromal","stroma")
  if (any(is_est)) {
    es <- M[which(is_est)[1], , drop = TRUE]
  } else if (any(is_imm) && any(is_str)) {
    es <- M[which(is_imm)[1], , drop = TRUE] + M[which(is_str)[1], , drop = TRUE]
  } else {
    stop("Could not find ESTIMATE/Immune/Stromal rows. Available rows: ",
         paste(rownames(M), collapse = ", "))
  }
  es <- as.numeric(es); names(es) <- colnames(M); es
}
ESTIMATEScore <- get_ESTIMATEScore(est_tbl)

# Yoshihara cosine transform → purity in [0,1]
purity <- cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
purity <- pmin(pmax(purity, 0), 1)
purity <- purity[samples]; names(purity) <- samples

message("Purity summary:"); print(summary(purity))

# ---------- Residualization helper (fractions + purity) ----------
residualize_expression <- function(expr, fractions_df, purity) {
  stopifnot(identical(colnames(expr), rownames(fractions_df)))
  stopifnot(identical(colnames(expr), names(purity)))
  covar <- cbind(fractions_df, purity = purity[rownames(fractions_df)])
  design <- model.matrix(~ ., data = as.data.frame(covar))
  resids <- matrix(NA_real_, nrow = nrow(expr), ncol = ncol(expr),
                   dimnames = dimnames(expr))
  for (i in seq_len(nrow(expr))) {
    y <- as.numeric(expr[i, ])
    fit <- tryCatch(lm.fit(x = design, y = y), error = function(e) NULL)
    if (!is.null(fit)) {
      yhat <- design %*% coef(fit)
      resids[i, ] <- y - as.numeric(yhat)
    }
  }
  resids
}

expr_resid <- residualize_expression(expr_mat, fractions_df, purity)

# ---------- Driver-timed labels (placeholder) ----------
# Choose top-variance genes as proxy drivers (adjust 'k' as desired)
k <- 4
gene_var <- apply(expr_mat, 1, var, na.rm = TRUE)
# ----------------------------
# 0) Check driver coverage now
# ----------------------------
#drivers <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")
drivers <- c("EGFR", "TP53", "CHEK2", "UBA1")
drivers <- c("LTBP4", "BMP4", "GREM1") 
present_drivers <- intersect(drivers, rownames(expr_mat))
missing_drivers <- setdiff(drivers, rownames(expr_mat))
cat("Present drivers: ", paste(present_drivers, collapse=", "), "\n")
if (length(missing_drivers)) {
  cat("Missing drivers (not in rownames(expr_mat)): ",
      paste(missing_drivers, collapse=", "), "\n")
}

# If none present, stop early to avoid all-NA downstream
stopifnot(length(present_drivers) > 0)

# ------------------------------------------------
# 1) Robust labeler: rank-based top-k (NA-tolerant)
#    - Guarantees at least k positives when data exist
#    - Replaces non-finite values with median
# ------------------------------------------------
compute_labels_topk <- function(expr_mat, driver_genes, top_prop = 0.10, min_pos = 1) {
  stopifnot(all(driver_genes %in% rownames(expr_mat)))
  n <- ncol(expr_mat)
  k <- max(min_pos, floor(top_prop * n))  # at least 'min_pos' positives
  
  labs <- sapply(driver_genes, function(g) {
    x <- as.numeric(expr_mat[g, ])
    # fix non-finite values
    if (!all(is.finite(x))) {
      med <- median(x[is.finite(x)], na.rm = TRUE)
      x[!is.finite(x)] <- med
    }
    # if still all the same value, create an all-zero label (can't rank)
    if (length(unique(x)) <= 1) return(rep(0L, n))
    
    # rank largest -> smallest; take top-k
    r <- rank(-x, ties.method = "min", na.last = "keep")
    lab <- as.integer(r <= k)
    # if everything NA, return NA vector to be caught by Jaccard
    if (all(is.na(lab))) lab[] <- NA_integer_
    lab
  })
  labs <- as.data.frame(labs, check.names = FALSE)
  rownames(labs) <- colnames(expr_mat)
  labs
}

# Recompute labels with the robust labeler (use same drivers for raw & resid)
labels_raw <- compute_labels_topk(expr_mat,   present_drivers, top_prop = 0.10, min_pos = 1)
labels_res <- compute_labels_topk(expr_resid, present_drivers, top_prop = 0.10, min_pos = 1)

# -----------------------------------------
# 2) NA-safe Jaccard (as you already added)
# -----------------------------------------
# --- 0) Quick sanity on covariates (fractions + purity) ---
cat("\nNA counts (fractions):\n"); print(colSums(is.na(fractions_df)))
cat("\nNA count (purity):\n"); print(sum(is.na(purity)))

# If any NA, median-impute and re-normalize fractions to sum to 1
if (anyNA(fractions_df) || anyNA(purity)) {
  frac_fix <- as.data.frame(fractions_df)
  for (j in names(frac_fix)) {
    m <- median(frac_fix[[j]], na.rm = TRUE)
    frac_fix[[j]][!is.finite(frac_fix[[j]]) | is.na(frac_fix[[j]])] <- m
  }
  rs <- rowSums(frac_fix)
  rs[rs == 0 | !is.finite(rs)] <- 1
  frac_fix <- frac_fix / rs
  fractions_df <- frac_fix
  
  pur_fix <- purity
  pm <- median(pur_fix, na.rm = TRUE)
  pur_fix[!is.finite(pur_fix) | is.na(pur_fix)] <- pm
  purity <- pur_fix
}

stopifnot(identical(rownames(fractions_df), colnames(expr_mat)))
stopifnot(identical(names(purity),          colnames(expr_mat)))

# --- 1) Robust residualization: ridge projection (handles collinearity, no NA propagation) ---
residualize_ridge <- function(expr, fractions_df, purity, lambda = 1e-6) {
  stopifnot(identical(colnames(expr), rownames(fractions_df)))
  stopifnot(identical(colnames(expr), names(purity)))
  
  Z <- cbind(fractions_df, purity = purity[rownames(fractions_df)])
  X <- model.matrix(~ 0 + Immune + Stromal + Endothelial + purity, data = as.data.frame(Z))
  # center columns to improve conditioning
  Xc <- scale(X, center = TRUE, scale = FALSE)
  XtX <- crossprod(Xc)
  p <- ncol(Xc)
  # ridge to ensure invertible (very small; doesn’t change fits materially)
  beta_hat <- function(y) {
    XtY <- crossprod(Xc, y)
    solve(XtX + diag(lambda, p), XtY, tol = 1e-12)
  }
  
  resids <- matrix(NA_real_, nrow = nrow(expr), ncol = ncol(expr), dimnames = dimnames(expr))
  for (i in seq_len(nrow(expr))) {
    y <- as.numeric(expr[i, ])
    y[!is.finite(y)] <- median(y[is.finite(y)], na.rm = TRUE)
    b  <- beta_hat(y)
    yh <- as.numeric(Xc %*% b)
    resids[i, ] <- y - yh
  }
  resids
}

expr_resid <- residualize_ridge(expr_mat, fractions_df, purity)

# Diagnostics: residual variance & finiteness
cat("\nResidual NA fraction per driver:\n")
for (d in c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")) {
  x <- as.numeric(expr_resid[d, ])
  cat(sprintf("%s  finite=%d/%d  sd=%.6g\n", d, sum(is.finite(x)), length(x), sd(x, na.rm = TRUE)))
}

# --- 2) Deterministic top-K labeler that ALWAYS returns K positives (handles ties) ---
compute_labels_force_topk <- function(expr_mat, driver_genes, k = 10, tie_eps = 1e-12) {
  stopifnot(all(driver_genes %in% rownames(expr_mat)))
  labs <- sapply(driver_genes, function(g) {
    x <- as.numeric(expr_mat[g, ])
    # replace non-finite by median
    if (any(!is.finite(x))) {
      med <- median(x[is.finite(x)], na.rm = TRUE)
      x[!is.finite(x)] <- med
    }
    # if perfectly flat, add tiny deterministic jitter to break ties
    if (isTRUE(all.equal(var(x), 0))) {
      idx <- seq_along(x)
      x <- x + (idx - mean(idx)) * tie_eps
    }
    ord <- order(x, decreasing = TRUE)
    lab <- integer(length(x))
    take <- ord[seq_len(min(k, length(ord)))]
    lab[take] <- 1L
    lab
  })
  labs <- as.data.frame(labs, check.names = FALSE)
  rownames(labs) <- colnames(expr_mat)
  labs
}

#drivers <- c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN")
drivers <- c("EGFR", "TP53", "CHEK2", "UBA1")
drivers <- c("LTBP4", "BMP4", "GREM1") 
K <- 10

labels_raw <- compute_labels_force_topk(expr_mat,   drivers, k = K)
labels_res <- compute_labels_force_topk(expr_resid, drivers, k = K)

# --- 3) Verify positives now exist in residuals ---
pos_counts <- tibble::tibble(
  driver  = drivers,
  raw_pos = vapply(drivers, \(d) sum(labels_raw[[d]] == 1, na.rm = TRUE), numeric(1)),
  res_pos = vapply(drivers, \(d) sum(labels_res[[d]] == 1, na.rm = TRUE), numeric(1))
)
print(pos_counts)

# --- 4) Jaccard overlap (order-agnostic, NA-safe) ---
top_ids <- function(v) names(v)[!is.na(v) & v == 1]
jaccard_sets <- function(a_ids, b_ids) {
  u <- base::union(a_ids, b_ids); if (length(u) == 0) return(NA_real_)
  length(base::intersect(a_ids, b_ids)) / length(u)
}

label_stability <- tibble::tibble(
  driver  = drivers,
  jaccard = vapply(
    drivers,
    \(d) jaccard_sets(top_ids(labels_raw[[d]]), top_ids(labels_res[[d]])),
    numeric(1)
  )
) |>
  dplyr::mutate(stable_0p85 = !is.na(jaccard) & jaccard >= 0.85)

label_stability_summary <- label_stability |>
  dplyr::summarise(
    median_jaccard      = median(jaccard, na.rm = TRUE),
    IQR_low             = quantile(jaccard, 0.25, na.rm = TRUE),
    IQR_high            = quantile(jaccard, 0.75, na.rm = TRUE),
    pct_drivers_ge_0p85 = mean(stable_0p85, na.rm = TRUE) * 100
  )

cat("\n--- Jaccard stability (raw vs residual) ---\n"); print(label_stability)
cat("\n--- Summary ---\n"); print(label_stability_summary)





#@@@@@@@@@@
jaccard_binary <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (!any(ok)) return(NA_real_)
  a <- as.integer(a[ok] == 1)
  b <- as.integer(b[ok] == 1)
  num <- sum(a == 1 & b == 1)
  den <- sum((a == 1) | (b == 1))
  if (den == 0) return(NA_real_)
  num / den
}

label_stability <- tibble(
  driver  = colnames(labels_raw),
  jaccard = vapply(colnames(labels_raw),
                   function(d) jaccard_binary(labels_raw[[d]], labels_res[[d]]),
                   numeric(1))
) %>%
  mutate(stable_0p85 = !is.na(jaccard) & jaccard >= 0.85)

label_stability_summary <- label_stability %>%
  summarise(
    median_jaccard      = suppressWarnings(median(jaccard, na.rm = TRUE)),
    IQR_low             = suppressWarnings(quantile(jaccard, 0.25, na.rm = TRUE)),
    IQR_high            = suppressWarnings(quantile(jaccard, 0.75, na.rm = TRUE)),
    pct_drivers_ge_0p85 = mean(stable_0p85, na.rm = TRUE) * 100
  )

print(label_stability)
print(label_stability_summary)

# -----------------------------------------
# 3) Optional: see which drivers are NA/empty
# -----------------------------------------
problem_drivers <- label_stability %>% filter(is.na(jaccard)) %>% pull(driver)
if (length(problem_drivers)) {
  cat("Drivers with undefined/NA Jaccard (likely missing or all-constant after residualization): ",
      paste(problem_drivers, collapse=", "), "\n")
  # quick diagnostics:
  for (d in problem_drivers) {
    if (d %in% rownames(expr_mat)) {
      xr <- expr_mat[d, ]; xs <- expr_resid[d, ]
      cat(sprintf("\n[%s] raw: range=%.3f..%.3f  resid: range=%.3f..%.3f  NA_raw=%d NA_resid=%d\n",
                  d, min(xr, na.rm=TRUE), max(xr, na.rm=TRUE),
                  min(xs, na.rm=TRUE), max(xs, na.rm=TRUE),
                  sum(!is.finite(xr)), sum(!is.finite(xs))))
    } else {
      cat(sprintf("\n[%s] not present in expr_mat rownames.\n", d))
    }
  }
}


# ---------- (2) Microenvironment surrogacy: partial correlations ----------
partial_cor_results <- {
  res <- list()
  for (drv in colnames(labels_raw)) {
    y <- labels_raw[[drv]]
    for (frac in colnames(fractions_df)) {
      pc <- ppcor::pcor.test(
        x = as.numeric(y),
        y = as.numeric(fractions_df[[frac]]),
        z = data.frame(purity = purity)
      )
      res[[length(res) + 1]] <- tibble(
        driver = drv, fraction = frac,
        partial_r = pc$estimate, p_value = pc$p.value
      )
    }
  }
  bind_rows(res) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"),
           pass = (abs(partial_r) < 0.2) & (q_value >= 0.05))
}

partial_cor_summary <- partial_cor_results %>%
  group_by(driver) %>%
  summarise(
    pct_fractions_passing = mean(pass, na.rm = TRUE) * 100,
    max_abs_r             = max(abs(partial_r), na.rm = TRUE),
    min_q                 = min(q_value, na.rm = TRUE),
    .groups               = "drop"
  )

cat("\n=== (2) Surrogacy (partial correlations) ===\n")
print(partial_cor_results)
print(partial_cor_summary)

# ---------- (3) Purity-stratified clinical consistency ----------
# Build toy outcomes if your data doesn't already include them.
# (Replace with your real survival/binary outcomes if available.)
set.seed(42)
# Simulate modest association with labels & purity for demonstration
base_haz <- 0.02
any_label <- as.integer(rowMeans(labels_raw) > 0)
haz <- base_haz * ifelse(any_label == 1, 1.6, 1.0) * (1 + 0.5*(1 - purity))
time <- rexp(length(samples), rate = haz)
status <- as.integer(runif(length(samples)) > 0.15)
surv_df <- data.frame(time = time, status = status, row.names = samples)

binary_outcome <- rbinom(length(samples), 1,
                         plogis(-0.5 + 0.6*any_label - 0.3*purity))
names(binary_outcome) <- samples

purity_tertile <- cut2(purity, g = 3)

fit_cox_safe <- function(df) {
  if (!any(df$status == 1) || length(unique(df$label)) < 2) {
    return(tibble(HR = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  fit <- coxph(Surv(time, status) ~ label, data = df)
  out <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |> filter(term == "label")
  tibble(HR = out$estimate, conf.low = out$conf.low, conf.high = out$conf.high, p.value = out$p.value)
}

fit_logit_safe <- function(df) {
  if (length(unique(df$y)) < 2 || length(unique(df$label)) < 2) {
    return(tibble(OR = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  fit <- glm(y ~ label, data = df, family = binomial())
  out <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |> filter(term == "label")
  tibble(OR = out$estimate, conf.low = out$conf.low, conf.high = out$conf.high, p.value = out$p.value)
}

# Use one representative driver for risk models (or loop over drivers)
drv_use <- colnames(labels_raw)[1]

# Cox by purity tertile + interaction
cox_results <- {
  lab <- labels_raw[[drv_use]]
  df <- tibble(sample = samples, time = surv_df$time, status = surv_df$status,
               label = lab, purity = purity, tertile = factor(purity_tertile))
  by_tertile <- df %>% group_by(tertile) %>% group_modify(~ fit_cox_safe(.x)) %>% ungroup()
  fit_int  <- coxph(Surv(time, status) ~ label * purity, data = df)
  int_p    <- tidy(fit_int) %>% filter(term == "label:purity") %>% pull(p.value)
  list(by_tertile = by_tertile, interaction_p = ifelse(length(int_p), int_p, NA_real_))
}
cat("\n=== (3a) Cox (by purity tertile) for driver:", drv_use, "===\n")
print(cox_results$by_tertile)
cat("Interaction p (label × purity):", cox_results$interaction_p, "\n")

# Logistic by purity tertile + interaction
logit_results <- {
  lab <- labels_raw[[drv_use]]
  df <- tibble(sample = samples, y = binary_outcome[sample], label = lab,
               purity = purity, tertile = factor(purity_tertile))
  by_tertile <- df %>% group_by(tertile) %>% group_modify(~ fit_logit_safe(.x)) %>% ungroup()
  fit_int <- glm(y ~ label * purity, data = df, family = binomial())
  int_p   <- tidy(fit_int) %>% filter(term == "label:purity") %>% pull(p.value)
  list(by_tertile = by_tertile, interaction_p = ifelse(length(int_p), int_p, NA_real_))
}
cat("\n=== (3b) Logistic (by purity tertile) for driver:", drv_use, "===\n")
print(logit_results$by_tertile)
cat("Interaction p (label × purity):", logit_results$interaction_p, "\n")

# ---------- Quick summaries to quote ----------
cat("\n--- Quick summaries ---\n")
cat(sprintf("Label stability: median Jaccard = %.2f; %%drivers ≥0.85 = %.1f%%\n",
            label_stability_summary$median_jaccard,
            label_stability_summary$pct_drivers_ge_0p85))
cat(sprintf("Microenv surrogacy: overall %% driver–fraction pairs passing = %.1f%%\n",
            mean(partial_cor_results$pass) * 100))
cat(sprintf("Purity interaction (Cox) p = %.3g | (Logit) p = %.3g\n",
            cox_results$interaction_p, logit_results$interaction_p))

