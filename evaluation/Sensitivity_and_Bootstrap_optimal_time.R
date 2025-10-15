setwd("/Users/atiteyk2/Documents/Master_Equation")
getwd()

# -----------------------------
# Sensitivity & Bootstrap of t_opt under CME
# -----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stats)
})

#genes_of_interest <- c("LTBP4", "BMP4", "GREM1") 
#genes_of_interest <- c("GATA3", "GPNMB", "BHLHE41") 

# ---- User inputs (replace df_params with your estimates) ----
set.seed(1)
df_params <- tibble::tibble(
  #@@gene  = c("TP53","PIK3CA","GATA3","GPNMB"),
  #@@kr    = c(0.30, 0.25, 0.15, 0.12),   # burst-initiation rate
  #@@mu    = c(2.0,  1.6,  1.2,  1.1),    # burst size (mean transcripts per burst)
  #@@gamma = c(0.20, 0.22, 0.18, 0.16),   # decay rate
  #@@x0    = c(1.0,  0.8,  0.7,  0.6)     # initial mean expression at t=0 (can be sample mean)
  
  #@@gene = c("TP53","BRCA1","BRCA2","ARID1A","PIK3CA","PTEN"),
  #@@kr    = c(0.30, 0.25, 0.15, 0.12, 0.10, 0.8),   # burst-initiation rate
  #@@mu    = c(2.0,  1.6,  1.2,  1.1, 0.8, 0.6),    # burst size (mean transcripts per burst)
  #@@gamma = c(0.20, 0.22, 0.18, 0.16, 0.14, 0.12),   # decay rate
  #@@x0    = c(1.0,  0.8,  0.7,  0.6, 0.5, 0.4)     # initial mean expression at t=0 (can be sample mean)
  
  #@@gene =  c("LTBP4", "BMP4", "GREM1") ,
  #@@kr    = c(0.30, 0.25, 0.15),   # burst-initiation rate
  #@@mu    = c(2.0,  1.6,  1.2),    # burst size (mean transcripts per burst)
  #@@gamma = c(0.20, 0.22, 0.18),   # decay rate
  #@@x0    = c(1.0,  0.8,  0.7)     # initial mean expression at t=0 (can be sample mean)
  
  gene =  c("GATA3", "GPNMB", "BHLHE41"),
  kr    = c(0.32, 0.24, 0.14),   # burst-initiation rate
  mu    = c(2.2,  1.4,  1.3),    # burst size (mean transcripts per burst)
  gamma = c(0.22, 0.21, 0.17),   # decay rate
  x0    = c(1.0,  0.85,  0.76)     # initial mean expression at t=0 (can be sample mean)
  
)

# ---- Core CME-derived helpers ----

# Time-dependent mean under CME: E[X](t) = E_inf + (x0 - E_inf) exp(-gamma t)
cme_mean <- function(t, kr, mu, gamma, x0) {
  Einf <- kr * mu / gamma
  Einf + (x0 - Einf) * exp(-gamma * t)
}

# Approximate variance at time t:
# We anchor to the stationary variance structure Var = m * (1 + mu*gamma/(kr+gamma))
# and let it track the time-dependent mean m(t). This preserves burst-driven overdispersion.
cme_var <- function(m_t, kr, mu, gamma) {
  m_t * (1 + (mu * gamma) / (kr + gamma))
}

# Draw ensemble at time t via Negative Binomial with matched mean/variance
# NB parameterization: var = m + m^2/size  => size = m^2/(var - m) when var > m
rnb_m_var <- function(n, m, v) {
  m <- pmax(m, 1e-9)
  v <- pmax(v, m + 1e-9)
  size <- (m^2) / (v - m)
  prob <- size / (size + m)
  rnbinom(n, size = size, prob = prob)
}

# Compute TailVar and MeanDiff at time t from an ensemble vector x
tail_metrics <- function(x, q_low = 0.1, q_high = 0.9) {
  lo <- quantile(x, q_low, names = FALSE, type = 7)
  hi <- quantile(x, q_high, names = FALSE, type = 7)
  upper <- x[x >= hi]
  lower <- x[x <= lo]
  tailvar <- var(c(upper, lower))
  meandiff <- mean(upper) - mean(lower)
  c(TailVar = tailvar, MeanDiff = meandiff)
}

# Normalize a vector to [0,1]
norm01 <- function(z) {
  zmin <- min(z); zmax <- max(z)
  if (abs(zmax - zmin) < 1e-12) return(rep(0, length(z)))
  (z - zmin) / (zmax - zmin)
}

# Compute t_opt for one gene given parameters
compute_topt <- function(kr, mu, gamma, x0,
                         n_sims = 500,
                         w_tail = 1, w_mean = 1,
                         q_low = 0.1, q_high = 0.9,
                         t_grid = NULL) {
  # Time grid: by default, cover 0 to 3 half-lives
  if (is.null(t_grid)) {
    tau_half <- log(2) / gamma
    t_grid <- seq(0, 3 * tau_half, length.out = 100)
  }
  # Simulate ensembles across time
  mets <- matrix(NA_real_, nrow = length(t_grid), ncol = 2,
                 dimnames = list(NULL, c("TailVar","MeanDiff")))
  for (i in seq_along(t_grid)) {
    t <- t_grid[i]
    m_t <- cme_mean(t, kr, mu, gamma, x0)
    v_t <- cme_var(m_t, kr, mu, gamma)
    x   <- rnb_m_var(n_sims, m_t, v_t)
    mm  <- tail_metrics(x, q_low = q_low, q_high = q_high)
    mets[i, ] <- mm
  }
  TVn <- norm01(mets[, "TailVar"])
  MDn <- norm01(mets[, "MeanDiff"])
  score <- w_tail * TVn + w_mean * MDn
  i_max <- which.max(score)
  list(
    t_opt = t_grid[i_max],
    score = score,
    t_grid = t_grid,
    TailVar = mets[, "TailVar"],
    MeanDiff = mets[, "MeanDiff"]
  )
}

# ---- Sensitivity analysis: ±10% perturbations in kr, mu, gamma ----
perturb_params <- function(kr, mu, gamma, pct = 0.10) {
  list(
    list(name = "kr_minus10",  kr = kr*(1 - pct), mu = mu,         gamma = gamma),
    list(name = "kr_plus10",   kr = kr*(1 + pct), mu = mu,         gamma = gamma),
    list(name = "mu_minus10",  kr = kr,          mu = mu*(1-pct),  gamma = gamma),
    list(name = "mu_plus10",   kr = kr,          mu = mu*(1+pct),  gamma = gamma),
    list(name = "gamma_minus10", kr = kr,        mu = mu,          gamma = gamma*(1-pct)),
    list(name = "gamma_plus10",  kr = kr,        mu = mu,          gamma = gamma*(1+pct))
  )
}

analyze_gene_sensitivity <- function(row,
                                     n_sims = 500,
                                     w_tail = 1, w_mean = 1,
                                     q_low = 0.1, q_high = 0.9) {
  kr <- row$kr; mu <- row$mu; gamma <- row$gamma; x0 <- row$x0
  base <- compute_topt(kr, mu, gamma, x0, n_sims = n_sims,
                       w_tail = w_tail, w_mean = w_mean,
                       q_low = q_low, q_high = q_high)
  base_t <- base$t_opt
  perts <- perturb_params(kr, mu, gamma, pct = 0.10)
  out <- purrr::map_dfr(perts, function(p) {
    tp <- compute_topt(p$kr, p$mu, p$gamma, x0, n_sims = n_sims,
                       w_tail = w_tail, w_mean = w_mean,
                       q_low = q_low, q_high = q_high)$t_opt
    tibble::tibble(
      gene = row$gene, perturb = p$name,
      t_opt_base = base_t, t_opt_pert = tp,
      rel_dev = abs(tp - base_t) / pmax(base_t, 1e-12)
    )
  })
  out
}

# ---- Bootstrap analysis of t_opt (re-simulate ensembles B times) ----
bootstrap_topt <- function(row, B = 500,
                           n_sims = 500,
                           w_tail = 1, w_mean = 1,
                           q_low = 0.1, q_high = 0.9) {
  kr <- row$kr; mu <- row$mu; gamma <- row$gamma; x0 <- row$x0
  vals <- numeric(B)
  for (b in seq_len(B)) {
    vals[b] <- compute_topt(kr, mu, gamma, x0, n_sims = n_sims,
                            w_tail = w_tail, w_mean = w_mean,
                            q_low = q_low, q_high = q_high)$t_opt
  }
  tibble::tibble(
    gene = row$gene,
    mean_topt = mean(vals),
    sd_topt = sd(vals),
    cv_topt = ifelse(mean(vals) > 0, sd(vals)/mean(vals), NA_real_)
  )
}

# ==============================
# Run analyses
# ==============================

# 1) Sensitivity to ±10% parameter perturbations
sens_results <- df_params %>%
  split(.$gene) %>%
  map_dfr(~ analyze_gene_sensitivity(.x[1, ])) %>%
  arrange(gene, perturb)

# Per-gene mean relative deviation
sens_summary_by_gene <- sens_results %>%
  group_by(gene) %>%
  summarise(mean_rel_dev = mean(rel_dev), .groups = "drop")

# Overall mean relative deviation across genes
overall_mean_rel_dev <- mean(sens_summary_by_gene$mean_rel_dev)

# 2) Bootstrap stability of t_opt (B = 500)
boot_results <- df_params %>%
  split(.$gene) %>%
  map_dfr(~ bootstrap_topt(.x[1, ], B = 500))

median_cv <- median(boot_results$cv_topt, na.rm = TRUE)

# ==============================
# Report the two statements
# ==============================
cat("\n--- Sensitivity summary ---\n")
print(sens_summary_by_gene)
cat(sprintf("\nAcross all genes, the mean relative deviation in t_opt under ±10%% parameter perturbations was %.3f.\n",
            overall_mean_rel_dev))

cat("\n--- Bootstrap summary (B = 500) ---\n")
print(boot_results)
cat(sprintf("\nThe median coefficient of variation (CV) of t_opt across genes was %.3f.\n",
            median_cv))

# ==============================
# Visualization of t_opt sensitivity & bootstrap stability
# ==============================
pkgs <- c("ggplot2","cowplot","scales")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(scales)
})

dir.create("topt_viz_out", showWarnings = FALSE)

# ---- 1) Heatmap of relative deviation under ±10% perturbations ----
sens_hmap_df <- sens_results %>%
  dplyr::mutate(
    perturb = factor(
      perturb,
      levels = c("kr_minus10","kr_plus10","mu_minus10","mu_plus10","gamma_minus10","gamma_plus10"),
      labels = c("k[r] -10%","k[r] +10%","mu -10%","mu +10%","gamma -10%","gamma +10%")
    )
  )

p_heat <- ggplot(sens_hmap_df, aes(x = gene, y = perturb, fill = rel_dev)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Relative\ndeviation",
                      labels = label_percent(accuracy = 1)) +
  labs(title = expression(paste("Sensitivity of ", t[opt], " to ", "\u00B110% parameter perturbations")),
       x = "Gene", y = "Perturbation") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  ) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(margin = margin(r = 8)))
p_heat
#ggsave("topt_viz_out/fig_heatmap_topt_relative_deviation.png", p_heat, width = 8.5, height = 4.8, dpi = 300)

# ---- 2) Faceted bar charts of relative deviation (per gene) ----
p_bars <- ggplot(sens_hmap_df, aes(x = perturb, y = rel_dev)) +
  geom_col(width = 0.65, fill = "grey30") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = expression(paste("Relative deviation of ", t[opt], " under \u00B110% perturbations")),
       x = NULL, y = "Relative deviation") +
  coord_flip() +
  facet_wrap(~ gene, ncol = 2, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  ) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))
p_bars
#ggsave("topt_viz_out/fig_bars_topt_relative_deviation.png", p_bars, width = 8.5, height = 6.5, dpi = 300)

# ---- 3) Violin + boxplot of bootstrap t_opt distributions (B = 500) ----
# Build full bootstrap draws for nice distributions (if you only computed summaries)
bootstrap_topt_vals <- function(row, B = 500,
                                n_sims = 500,
                                w_tail = 1, w_mean = 1,
                                q_low = 0.1, q_high = 0.9) {
  kr <- row$kr; mu <- row$mu; gamma <- row$gamma; x0 <- row$x0
  vals <- numeric(B)
  for (b in seq_len(B)) {
    vals[b] <- compute_topt(kr, mu, gamma, x0, n_sims = n_sims,
                            w_tail = w_tail, w_mean = w_mean,
                            q_low = q_low, q_high = q_high)$t_opt
  }
  tibble::tibble(gene = row$gene, t_opt = vals)
}

boot_vals <- df_params %>%
  split(.$gene) %>%
  purrr::map_dfr(~ bootstrap_topt_vals(.x[1, ], B = 500))

# Merge in CVs for informative subtitles/labels
boot_cv <- boot_results %>% dplyr::select(gene, cv_topt)
boot_vals_plot <- dplyr::left_join(boot_vals, boot_cv, by = "gene")

p_violin <- ggplot(boot_vals_plot, aes(x = gene, y = t_opt)) +
  geom_violin(fill = "grey80", color = "grey40", linewidth = 0.3, width = 0.9) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, fill = "white") +
  labs(title = expression(paste("Bootstrap distribution of ", t[opt], " (B = 500)")),
       subtitle = "Box shows IQR; violin shows density. CV per gene reported in caption/table.",
       x = "Gene", y = expression(t[opt])) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  )
p_violin
#@@ggsave("topt_viz_out/fig_violin_bootstrap_topt.png", p_violin, width = 8.5, height = 5.0, dpi = 300)

# ---- 4) Overlay of Score(t) curves for a selected gene (baseline vs. perturbations) ----
# Helper to compute Score(t) series (re-using your TailVar/MeanDiff definition via compute_topt)
score_series <- function(label, kr, mu, gamma, x0,
                         n_sims = 500, q_low = 0.1, q_high = 0.9,
                         w_tail = 1, w_mean = 1, t_grid = NULL) {
  o <- compute_topt(kr, mu, gamma, x0, n_sims = n_sims,
                    w_tail = w_tail, w_mean = w_mean,
                    q_low = q_low, q_high = q_high, t_grid = t_grid)
  tibble::tibble(
    t = o$t_grid,
    Score = o$score,
    label = label,
    t_opt = o$t_opt
  )
}

# Choose a gene to visualize (first by default)
gene_to_plot <- df_params$gene[1]
row <- df_params %>% dplyr::filter(gene == gene_to_plot) %>% dplyr::slice(1)

# Baseline + ±10% per-parameter overlays
plist <- perturb_params(row$kr, row$mu, row$gamma, pct = 0.10)
score_df <- dplyr::bind_rows(
  score_series("baseline", row$kr, row$mu, row$gamma, row$x0),
  purrr::map_dfr(plist, ~ score_series(.x$name, .x$kr, .x$mu, .x$gamma, row$x0))
)

# Extract unique t_opt per label for vertical markers
# Make sure 'vdf' is robust if any NA slipped into t_opt
vdf <- score_df %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(t_opt = dplyr::last(stats::na.omit(t_opt)), .groups = "drop")

p_score <- ggplot(score_df, aes(x = t, y = Score, color = label)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_vline(data = vdf, aes(xintercept = t_opt, color = label),
             linetype = 2, linewidth = 0.6, alpha = 0.8, show.legend = FALSE) +
  labs(title = paste0("Score curves for ", gene_to_plot,
                      " (baseline vs. \u00B110% perturbations)"),
       x = "Time (same units as CME grid)", y = "Composite Score(t)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(color = "black"),
        legend.position = "top",
        legend.key.height = unit(12, "pt"),
        strip.text = element_text(face = "bold")
  )+
  theme(legend.position = "bottom", legend.title = element_blank())
p_score
#@@ggsave(sprintf("topt_viz_out/fig_score_overlay_%s.png", gene_to_plot), p_score,
#@@       width = 8.5, height = 5.0, dpi = 300)

# ---- Optional: Assemble a quick dashboard image ----
dash <- cowplot::plot_grid(
  p_heat, p_bars, p_violin, p_score,
  labels = c("A","B","C","D"), ncol = 2, rel_heights = c(1,1)
)
dash
#@@ggsave("topt_viz_out/fig_dashboard_topt_sensitivity_bootstrap.png", dash,
#@@       width = 12, height = 9, dpi = 300)

# ------------------------------
# Printing where files are saved
# ------------------------------
cat("\nSaved plots to 'topt_viz_out/':\n",
    " - fig_heatmap_topt_relative_deviation.png\n",
    " - fig_bars_topt_relative_deviation.png\n",
    " - fig_violin_bootstrap_topt.png\n",
    sprintf(" - fig_score_overlay_%s.png\n", gene_to_plot),
    " - fig_dashboard_topt_sensitivity_bootstrap.png\n", sep = "")
