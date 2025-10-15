setwd("/Users/atiteyk2/Documents/Master_Equation")
getwd()

# Install if needed
pkgs <- c("tidyverse", "pROC", "precrec", "cowplot", "scales")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(tidyverse)
library(pROC)     # ROC AUC + DeLong CI (optional)
library(precrec)  # ROC & PR curves/AUCs from pooled predictions
library(cowplot)  # plot grids
library(scales)   # pretty breaks

# ====== Benchmarking magicSubclonal vs baselines (ROC, PR, Calibration, Brier) ======
# Inputs expected (long format, one row per sample×method×marker):
#   sample_id    : unique sample identifier
#   fold         : CV fold ID; each sample appears in exactly ONE test fold
#   method       : "magicSubclonal", "NMF", "sciClone", etc.
#   score        : predicted probability or score (higher = more "subclonal")
#   marker_name  : gene name for the known subclonal marker (e.g., "TP53")
#   marker_label : 0/1 label indicating whether THIS SAMPLE truly has that marker-positive subclone
#
# You can read this from a CSV produced by your pipeline (see READ-IN section).
# If you already have a per-sample pooled label (any-marker), see the "ALTERNATIVE" note below.

# ---------------------------
# Known subclonal markers (from manuscript Table 1)
# ---------------------------
known_markers <- c(
  "TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN",
  "EGFR","CHEK2","UBA1",
  "LTBP4","BMP4","GREM1",
  "GATA3","GPNMB","BHLHE41"
)

# ---------------------------
# Packages
# ---------------------------
pkgs <- c("tidyverse", "pROC", "precrec", "cowplot", "scales")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------
# READ-IN (edit path as needed)
# ---------------------------
# Expected columns in CSV: sample_id, fold, method, score, marker_name, marker_label
# Scores must be OUT-OF-FOLD (i.e., from held-out fold); if not, split beforehand.
# Example:
# preds_raw <- readr::read_csv("predictions_long.csv")

# For demonstration, here's a small synthetic generator. Replace with your real data import.
set.seed(42)
n_samp <- 300
k_folds <- 5
methods <- c("magicSubclonal", "NMF", "sciClone")
samples <- tibble(sample_id = sprintf("S%03d", 1:n_samp),
                  fold = sample(1:k_folds, n_samp, replace = TRUE))
# Simulate per-marker truth table (sparse positives for known markers)
truth_long <- samples %>%
  tidyr::crossing(marker_name = known_markers) %>%
  mutate(marker_label = rbinom(n(), 1, prob = 0.10))  # ~10% positives per marker

# Simulate method scores with signal aligned to truth (magic strongest)
gen_scores <- function(base, signal) pmin(pmax(base + signal, 0), 1)
preds_raw <- truth_long %>%
  tidyr::crossing(method = methods) %>%
  left_join(samples, by = c("sample_id","fold")) %>%
  mutate(
    base = runif(n()),
    signal = case_when(
      method == "magicSubclonal" ~ 0.35*marker_label,
      method == "NMF"            ~ 0.20*marker_label,
      TRUE                       ~ 0.15*marker_label
    ),
    score = gen_scores(base, signal)
  ) %>% 
  select(sample_id, fold, method, score, marker_name, marker_label)

# ---------------------------
# (Optional) Isotonic calibration (leave-one-fold-out, LOFO)
# ---------------------------
use_isotonic <- TRUE  # set FALSE to skip calibration

eps <- 1e-6
clamp01 <- function(p) pmin(pmax(p, eps), 1 - eps)

iso_cv_calibrate <- function(df) {
  # df requires: sample_id, fold, method, score, marker_label
  df %>%
    group_by(method) %>%
    group_modify(function(d, key) {
      folds <- sort(unique(d$fold))
      pred_iso <- numeric(nrow(d))
      for (f in folds) {
        train_idx <- which(d$fold != f)
        test_idx  <- which(d$fold == f)
        # Fit isotonic on training partitions of CURRENT method and CURRENT marker_name pooling
        fit <- isoreg(x = d$score[train_idx], y = d$marker_label[train_idx])
        iso_fun <- approxfun(x = fit$x, y = fit$y, rule = 2, ties = mean)
        pred_iso[test_idx] <- iso_fun(d$score[test_idx])
      }
      d$score_iso <- clamp01(pred_iso)
      d
    }) %>%
    ungroup()
}

# Apply isotonic per method, pooled across markers (no leakage because of folds)
preds_cal <- if (use_isotonic) iso_cv_calibrate(preds_raw) else preds_raw %>% mutate(score_iso = score)

# For downstream use the calibrated score
preds_long <- preds_cal %>%
  mutate(score = ifelse(is.na(score_iso), score, score_iso)) %>%
  select(sample_id, fold, method, score, marker_name, marker_label)

# ---------------------------
# Build two evaluation targets:
#  (A) Per-marker evaluation (each marker_name separately)
#  (B) Pooled "any-marker positive" per sample (label 1 if ANY known marker positive)
# ---------------------------
# (A) Per-marker stays as is (preds_long).
# (B) Collapse truth to "any positive marker" per sample; combine multiple rows for the same sample×method.


any_marker_truth <- preds_long %>%
  group_by(sample_id) %>%
  summarize(any_marker_label = as.integer(any(marker_label == 1)),
            .groups = "drop")

preds_any <- preds_long %>%
  group_by(sample_id, fold, method) %>%
  summarize(score = mean(score), .groups = "drop") %>%  # average across markers for pooled decision; you can change to max()
  left_join(any_marker_truth, by = "sample_id") %>%
  rename(marker = any_marker_label)

# ---------------------------
# Helper functions for metrics and curves
# ---------------------------
calib_stats <- function(y, p) {
  p <- clamp01(p); logit_p <- qlogis(p)
  m <- suppressWarnings(glm(y ~ logit_p, family = binomial()))
  co <- coef(summary(m))
  tibble(
    calib_intercept = unname(co[1, "Estimate"]),
    calib_intercept_se = unname(co[1, "Std. Error"]),
    calib_slope = unname(co[2, "Estimate"]),
    calib_slope_se = unname(co[2, "Std. Error"])
  )
}

reliability_curve <- function(y, p, bins = 10) {
  p <- clamp01(p)
  q <- cut_number(p, n = bins)
  tibble(y = y, p = p, bin = q) %>%
    group_by(bin) %>%
    summarize(pred_mean = mean(p), obs_rate = mean(y), n_bin = n(), .groups = "drop") %>%
    arrange(pred_mean)
}

brier_score <- function(y, p) mean((clamp01(p) - y)^2)

compute_fold_metrics <- function(df) {
  # df: one METHOD across folds with columns: fold, score, marker
  df %>%
    group_by(fold) %>%
    group_modify(function(d, key) {
      # ROC AUC
      roc_obj <- pROC::roc(d$marker, d$score, quiet = TRUE, direction = "<")
      auc_roc <- as.numeric(pROC::auc(roc_obj))
      # PR AUC
      mm <- precrec::mmdata(scores = d$score, labels = d$marker, modnames = "m", dsids = 1)
      ev <- precrec::evalmod(mm)
      aucs <- precrec::auc(ev)
      pr_auc <- aucs %>% filter(curvetypes == "PRC") %>% pull(aucs)
      # Brier
      bs <- brier_score(d$marker, d$score)
      # Calibration
      cs <- calib_stats(d$marker, d$score)
      tibble(
        auc_roc = auc_roc,
        auc_pr  = pr_auc,
        brier   = bs,
        calib_intercept = cs$calib_intercept,
        calib_slope     = cs$calib_slope
      )
    }) %>%
    ungroup()
}

pooled_curves_and_aucs <- function(df) {
  # df: columns sample_id, method, score, marker
  df %>%
    group_by(method) %>%
    group_modify(function(d, key) {
      mm <- precrec::mmdata(scores = d$score, labels = d$marker,
                            modnames = unique(d$method), dsids = 1)
      ev <- precrec::evalmod(mm)
      roc_df <- as_tibble(precrec::fortify(ev, "ROC"))
      pr_df  <- as_tibble(precrec::fortify(ev, "PRC"))
      aucs <- precrec::auc(ev) %>%
        select(curvetypes, aucs) %>%
        pivot_wider(names_from = curvetypes, values_from = aucs, names_prefix = "AUC_")
      list(roc = roc_df, pr = pr_df, aucs = aucs)
    }) %>% 
    ungroup()
}

# ---------------------------
# (A) Per-marker evaluation table
# ---------------------------
# ---------------------------
# (A) Per-marker evaluation table — FIXED
# ---------------------------
per_marker_results <- preds_long %>%
  filter(marker_name %in% known_markers) %>%
  group_by(marker_name, method) %>%
  group_modify(function(d, key) {
    # pooled curves
    mm <- precrec::mmdata(scores = d$score, labels = d$marker_label, modnames = "m", dsids = 1)
    ev <- precrec::evalmod(mm)
    aucs <- precrec::auc(ev)
    auc_roc <- aucs %>% dplyr::filter(curvetypes == "ROC") %>% dplyr::pull(aucs)
    auc_pr  <- aucs %>% dplyr::filter(curvetypes == "PRC") %>% dplyr::pull(aucs)
    
    # fold metrics
    ftab <- d %>%
      dplyr::rename(marker = marker_label) %>%
      dplyr::select(fold, score, marker) %>%
      compute_fold_metrics()
    
    # return a one-row tibble (not summarize())
    tibble::tibble(
      AUC_ROC_pooled = auc_roc,
      AUC_PR_pooled  = auc_pr,
      AUC_ROC_mean = mean(ftab$auc_roc), AUC_ROC_sd = sd(ftab$auc_roc),
      AUC_PR_mean  = mean(ftab$auc_pr),  AUC_PR_sd  = sd(ftab$auc_pr),
      Brier_mean   = mean(ftab$brier),   Brier_sd   = sd(ftab$brier),
      CalibSlope_mean = mean(ftab$calib_slope), CalibSlope_sd = sd(ftab$calib_slope),
      CalibIntercept_mean = mean(ftab$calib_intercept), CalibIntercept_sd = sd(ftab$calib_intercept)
    )
  }) %>%
  ungroup() %>%
  arrange(marker_name, desc(AUC_ROC_mean))


# ---------------------------
# (B) Pooled "any-marker" benchmarking
# ---------------------------
# Rebuild preds_any so method is preserved
preds_any <- preds_long %>%
  dplyr::group_by(sample_id, fold, method) %>%   # <-- keep method here
  dplyr::summarize(score = mean(score), .groups = "drop") %>%
  dplyr::left_join(
    preds_long %>%
      dplyr::group_by(sample_id) %>%
      dplyr::summarize(any_marker_label = as.integer(any(marker_label == 1)), .groups = "drop"),
    by = "sample_id"
  ) %>%
  dplyr::rename(marker = any_marker_label)

# Quick sanity check
if (!all(c("sample_id","method","score","marker") %in% names(preds_any))) {
  stop("preds_any is missing required columns. It has: ",
       paste(names(preds_any), collapse = ", "))
}

pooled_curves_and_aucs <- function(df) {
  # df must have: sample_id, method, score, marker
  req <- c("sample_id","method","score","marker")
  if (!all(req %in% names(df))) {
    stop("Input df must have columns: ", paste(req, collapse=", "),
         ". Found: ", paste(names(df), collapse=", "))
  }
  # split by method to avoid dplyr's tidy-eval warnings
  sp <- split(df, df$method)
  out <- lapply(names(sp), function(m) {
    d <- sp[[m]]
    # guard against empty method slice
    if (nrow(d) == 0) return(NULL)
    mm <- precrec::mmdata(scores = d$score, labels = d$marker,
                          modnames = m, dsids = 1)
    ev <- precrec::evalmod(mm)
    roc_df <- tibble::as_tibble(ggplot2::fortify(ev, "ROC"))
    pr_df  <- tibble::as_tibble(ggplot2::fortify(ev, "PRC"))
    aucs <- precrec::auc(ev) %>%
      dplyr::select(curvetypes, aucs) %>%
      tidyr::pivot_wider(names_from = curvetypes, values_from = aucs, names_prefix = "AUC_")
    tibble::tibble(method = m, roc = list(roc_df), pr = list(pr_df), aucs = list(aucs))
  })
  dplyr::bind_rows(out)
}


# This select will now succeed because method exists
pool_any <- pooled_curves_and_aucs(
  preds_any %>% dplyr::select(sample_id, method, score, marker)
)

roc_any <- pool_any %>%
  dplyr::select(method, roc) %>%
  tidyr::unnest(roc)

pr_any <- pool_any %>%
  dplyr::select(method, pr) %>%
  tidyr::unnest(pr)

auc_any <- pool_any %>%
  dplyr::select(method, aucs) %>%
  tidyr::unnest(aucs) %>%
  dplyr::rename(AUC_ROC = AUC_ROC, AUC_PR = AUC_PRC)

print(auc_any)


# Fold-wise summaries
fold_summ_any <- preds_any %>%
  group_by(method) %>%
  group_modify(~ compute_fold_metrics(.x)) %>%
  ungroup()

fold_summary_table <- fold_summ_any %>%
  group_by(method) %>%
  summarize(
    AUC_ROC_mean = mean(auc_roc), AUC_ROC_sd = sd(auc_roc),
    AUC_PR_mean  = mean(auc_pr),  AUC_PR_sd  = sd(auc_pr),
    Brier_mean   = mean(brier),   Brier_sd   = sd(brier),
    CalibSlope_mean = mean(calib_slope), CalibSlope_sd = sd(calib_slope),
    CalibIntercept_mean = mean(calib_intercept), CalibIntercept_sd = sd(calib_intercept),
    .groups = "drop"
  )

bench_any <- auc_any %>%
  left_join(fold_summary_table, by = "method") %>%
  arrange(desc(AUC_ROC_mean))

print(bench_any)

# ---------------------------
# Calibration (pooled) for any-marker
# ---------------------------
calib_curves_any <- preds_any %>%
  group_by(method) %>%
  group_modify(function(d, key) {
    rel <- reliability_curve(d$marker, d$score, bins = 10)
    stats <- calib_stats(d$marker, d$score)
    rel$calib_intercept <- stats$calib_intercept
    rel$calib_slope     <- stats$calib_slope
    rel
  }) %>%
  ungroup()

# ---------------------------
# PLOTS (save to disk)
# ---------------------------
dir.create("benchmark_out", showWarnings = FALSE)

p_roc <- ggplot(roc_any, aes(x = x, y = y, color = method)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "ROC Curves (pooled across test folds, any-marker)",
       color = "Method") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  )
p_roc

p_pr <- ggplot(pr_any, aes(x = x, y = y, color = method)) +
  geom_line(size = 1) +
  labs(x = "Recall", y = "Precision",
       title = "Precision–Recall Curves (pooled, any-marker)",
       color = "Method") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  )
p_pr

p_calib <- calib_curves_any %>%
  mutate(lbl = sprintf("%s  (α=%.2f, β=%.2f)", method, calib_intercept, calib_slope)) %>%
  ggplot(aes(x = pred_mean, y = obs_rate, color = method)) +
  geom_point(aes(size = n_bin), alpha = 0.9) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_size_continuous(name = "Bin size") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Mean predicted probability", y = "Observed event rate",
       title = "Calibration (Reliability) Curves (any-marker)",
       color = "Method") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  )
p_calib

