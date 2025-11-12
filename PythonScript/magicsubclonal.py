"""
magicsubclonal.py
CME-based subclonal analysis & visualization with full data export.

Figures saved:
  - CME_Dynamics_AllDrivers.png
  - HighRiskSubcloneGenes_noOverlap.png
  - TumorHeterogeneity_annotated.png
  - Fraction_LowHigh_Subclones.png
  - RiskClassHeatmap.png

Tables saved (CSV):
  - exprs_data.csv
  - x0_list__<DRIVER>.csv
  - param_estimates.csv
  - gene_simulations__<DRIVER>.csv
  - t_opt.csv
  - subclone_labels__<DRIVER>.csv
  - filtered_subclonal_genes_all.csv
  - heterogeneity_metrics.csv
  - advantage_scores.csv
  - risk_state_matrix_rowZ.csv        (matrix used in RiskClassHeatmap)
  - heterogeneity_matrix_rowZ.csv     (matrix used in heterogeneity heatmap)

Also:
  - manifest.json  (paths to all outputs)
"""

from __future__ import annotations
import os, json
import numpy as np
import pandas as pd
from scipy import optimize
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from matplotlib.patches import Patch

# optional: non-overlapping text labels
try:
    from adjustText import adjust_text
    _HAS_ADJUSTTEXT = True
except Exception:
    _HAS_ADJUSTTEXT = False


# =========
# IO utils
# =========
def _ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def _save_csv(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, index=True)

def _np_to_py(obj):
    """Make dicts/arrays JSON-serializable."""
    if isinstance(obj, dict):
        return {str(k): _np_to_py(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_np_to_py(x) for x in obj]
    if isinstance(obj, (np.integer, )):
        return int(obj)
    if isinstance(obj, (np.floating, )):
        return float(obj)
    if isinstance(obj, (np.ndarray, )):
        return obj.tolist()
    return obj


# =========================
# CME simulation primitives
# =========================
def _simulate_CME(x0: float, gamma2: float, kr: float, mu: float,
                  t_max: float, dt: float = 0.1) -> pd.DataFrame:
    time = np.arange(0, t_max + 1e-9, dt)
    x = np.zeros_like(time, dtype=float)
    x[0] = float(x0)
    for i in range(1, len(time)):
        decay = gamma2 * x[i - 1] * dt
        n_bursts = np.random.poisson(kr * dt)
        burst = np.random.exponential(1.0 / mu, n_bursts).sum() if n_bursts > 0 else 0.0
        x[i] = max(x[i - 1] - decay + burst, 0.0)
    return pd.DataFrame({"time": time, "x": x})


def _loss_function(params: np.ndarray, x0_vec: np.ndarray) -> float:
    gamma2, kr, mu = params
    finals = []
    for x0 in x0_vec:
        df = _simulate_CME(float(x0), gamma2, kr, mu, t_max=50.0)
        finals.append(df["x"].iloc[-1])
    finals = np.asarray(finals, float)
    if np.any(~np.isfinite(finals)):
        return 1e9
    return float(np.var(finals))


def _score_time(df: pd.DataFrame,
                q_rare: float = 0.05,
                w_tail: float = 0.5,
                w_mean: float = 0.5) -> float:
    rows = []
    for t, grp in df.groupby("time"):
        x = grp["x"].to_numpy(dtype=float)
        lo, hi = np.quantile(x, [q_rare, 1 - q_rare])
        mask = (x <= lo) | (x >= hi)
        var_tail = np.var(x[mask], ddof=1) if mask.any() else 0.0
        mean_diff = (x[x >= hi].mean() - x[x <= lo].mean()) if mask.any() else 0.0
        rows.append((float(t), var_tail, mean_diff))
    stab = pd.DataFrame(rows, columns=["time", "var_tail", "mean_diff"])
    if stab.empty:
        return 0.0
    mm = lambda v: (v - v.min()) / (v.max() - v.min() + 1e-12)
    stab["score"] = w_tail * mm(stab["var_tail"]) + w_mean * mm(stab["mean_diff"])
    return float(stab.loc[stab["score"].idxmax(), "time"])


# =========================
# Figure 0: CME dynamics
# =========================
def _plot_cme_dynamics_all(gene_sims: Dict[str, pd.DataFrame],
                           t_opt: Dict[str, float]) -> plt.Figure:
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = sns.color_palette("tab10", len(gene_sims))
    for (g, df), col in zip(gene_sims.items(), colors):
        bt = df.groupby("time")["x"].agg(
            mean="mean",
            low=lambda z: np.quantile(z, 0.025),
            hi=lambda z: np.quantile(z, 0.975)
        ).reset_index()
        ax.plot(bt["time"], bt["mean"], color=col, lw=1.6, label=g)
        ax.fill_between(bt["time"], bt["low"], bt["hi"], color=col, alpha=0.25)
        ax.axvline(t_opt.get(g, np.nan), ls="--", color=col, lw=1.2, alpha=0.5)
    ax.set_title("CME Dynamics across driver genes", fontsize=15, fontweight="bold")
    ax.set_xlabel("Time"); ax.set_ylabel("Expression")
    ax.legend(title="Driver")
    fig.tight_layout()
    return fig


# =======================================
# Figure 1: High-risk subclone gene labels
# =======================================
def _plot_high_risk_subclone_genes_panel(
    gene_sims: Dict[str, pd.DataFrame],
    t_opt_map: Dict[str, float],
    filtered_df: pd.DataFrame,
    exprs_mat: pd.DataFrame,
    subclone_labels: Dict[str, pd.DataFrame],
    drivers: List[str],
    ncol: int = 3,
    max_labels_per_driver: int = 10,
) -> plt.Figure:

    sns.set_style("whitegrid")

    def _rare0_for(driver: str) -> np.ndarray:
        labs = subclone_labels.get(driver)
        if isinstance(labs, pd.DataFrame) and not labs.empty and {"type", "source_id"}.issubset(labs.columns):
            r1 = pd.to_numeric(labs.loc[labs["type"].isin(["low", "high"]), "source_id"],
                               errors="coerce").dropna().astype(int).unique()
            return np.array([i - 1 for i in r1 if 1 <= i <= exprs_mat.shape[1]], dtype=int)
        return np.array([], dtype=int)

    def _band_for(driver: str) -> pd.DataFrame:
        df = gene_sims[driver]
        return df.groupby("time")["x"].agg(
            low=lambda z: np.quantile(z, 0.025),
            hi=lambda z: np.quantile(z, 0.975),
            mean="mean"
        ).reset_index()

    facet = [d for d in drivers if d in gene_sims and d in t_opt_map]
    n = len(facet)
    ncol = max(1, int(ncol))
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.3 * ncol, 4.5 * nrow), squeeze=False)

    for i, drv in enumerate(facet):
        r, c = divmod(i, ncol)
        ax = axes[r, c]
        bt = _band_for(drv)
        ax.fill_between(bt["time"], bt["low"], bt["hi"], color="#EBD2A6", alpha=0.45)
        ax.plot(bt["time"], bt["mean"], color="#3C3C3C", lw=1.4)
        tstar = float(t_opt_map[drv])
        ax.axvline(tstar, ls="--", color="black", lw=1.3)

        sub = pd.DataFrame(columns=["gene", "status"])
        if isinstance(filtered_df, pd.DataFrame) and not filtered_df.empty:
            sub = filtered_df[filtered_df["driver"] == drv][["gene", "status"]].drop_duplicates()
        sub = pd.concat([sub[sub["status"] == "enriched"],
                         sub[sub["status"] == "depleted"]], ignore_index=True).head(max_labels_per_driver)

        R0 = _rare0_for(drv)
        texts = []
        for _, row in sub.iterrows():
            gene = str(row["gene"])
            status = str(row["status"])
            color = "#D62728" if status == "enriched" else "#1f77b4"
            if gene in exprs_mat.index and R0.size:
                vals = exprs_mat.loc[gene].to_numpy(dtype=float)[R0]
                vals = vals[np.isfinite(vals)]
                y = float(np.mean(vals)) if vals.size else np.nan
            else:
                y = np.nan
            if not np.isfinite(y):
                y = float(bt.iloc[(bt["time"] - tstar).abs().argmin()]["mean"])
            ax.scatter([tstar], [y], s=30, color=color, zorder=3)
            txt = ax.text(tstar + 0.6, y, gene, fontsize=10, color=color, fontweight="bold",
                          bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=color, lw=1.2, alpha=0.95),
                          va="center")
            texts.append(txt)

        if texts:
            if _HAS_ADJUSTTEXT:
                adjust_text(texts, ax=ax,
                            expand_points=(1.2, 1.2), expand_text=(1.1, 1.4),
                            arrowprops=dict(arrowstyle="-", color="gray", lw=0.6))
            else:
                ys = [t.get_position()[1] for t in texts]
                if len(ys) > 1:
                    jitter = np.linspace(-0.6, 0.6, len(ys))
                    for t, dy in zip(texts, jitter):
                        x, y0 = t.get_position()
                        t.set_position((x, y0 + dy))

        ax.set_title(drv, fontsize=13, fontweight="bold")
        ax.set_xlabel("Time"); ax.set_ylabel("")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    for j in range(n, nrow * ncol):
        rr, cc = divmod(j, ncol); axes[rr, cc].axis("off")

    fig.suptitle(r"High risk subclones at optimal time $t_{opt}$ on CME envelopes",
                 fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout()
    return fig


# ==========================================
# Heterogeneity heatmap (robust)
# ==========================================
def _plot_tumor_heterogeneity_heatmap_upgraded(
    exprs_mat: pd.DataFrame,
    subclone_labels: Dict[str, pd.DataFrame],
    filtered_df: pd.DataFrame,
    drivers: List[str],
    n_top_mad: int = 60,
    z_clip: float = 2.5
) -> Tuple[plt.Figure, pd.DataFrame]:
    rare_ids = []
    for d, labs in subclone_labels.items():
        if labs is None or labs.empty:
            continue
        if {"type", "source_id"}.issubset(labs.columns):
            r = pd.to_numeric(labs.loc[labs["type"].isin(["low", "high"]), "source_id"],
                              errors="coerce").dropna().astype(int).tolist()
            rare_ids.extend(r)
    rare_ids = sorted(set([i for i in rare_ids if 1 <= i <= exprs_mat.shape[1]]))
    if rare_ids:
        cols0 = [i - 1 for i in rare_ids]
        Xall = exprs_mat.iloc[:, cols0].copy()
    else:
        Xall = exprs_mat.copy()

    keep = set()
    if isinstance(filtered_df, pd.DataFrame) and not filtered_df.empty:
        keep.update(filtered_df["gene"].astype(str).tolist())

    med = Xall.median(axis=1)
    mad = (Xall.sub(med, axis=0)).abs().median(axis=1)
    top_mad = mad.sort_values(ascending=False).head(min(n_top_mad, mad.size)).index.tolist()
    keep.update(top_mad)

    X = Xall.loc[Xall.index.intersection(sorted(keep))].copy()
    if X.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off"); ax.text(0.5, 0.5, "No heatmap data", ha="center", va="center")
        return fig, pd.DataFrame(index=[], columns=[])

    med = X.median(axis=1)
    mad = (X.sub(med, axis=0)).abs().median(axis=1).replace(0, np.nan)
    Xz = X.sub(med, axis=0).div(mad, axis=0)  # robust row z
    Xz = Xz.replace([np.inf, -np.inf], np.nan).fillna(0.0).clip(-z_clip, z_clip)

    def col_label(j0):
        j1 = j0 + 1
        owners = []
        for d, labs in subclone_labels.items():
            if labs is None or labs.empty:
                continue
            if {"type", "source_id"}.issubset(labs.columns):
                r = pd.to_numeric(labs.loc[labs["type"].isin(["low", "high"]), "source_id"],
                                  errors="coerce").dropna().astype(int).tolist()
                if j1 in r:
                    owners.append(d)
        if len(owners) == 0: return "None"
        if len(owners) == 1: return owners[0]
        return "Multiple"

    col_lab = pd.Series([col_label(j) for j in range(Xz.shape[1])],
                        index=Xz.columns, name="Driver label")

    sns.set(font_scale=0.85)
    base = sns.color_palette("Set2", len(drivers))
    drv_colors = {"None": "#cfcfcf", "Multiple": "#444444"}
    drv_colors.update({d: base[i % len(base)] for i, d in enumerate(drivers)})
    cmap = sns.color_palette("vlag", as_cmap=True)

    g = sns.clustermap(
        Xz, cmap=cmap, method="average", metric="euclidean",
        linewidths=0, xticklabels=False, yticklabels=True,
        col_colors=col_lab.map(drv_colors),
        figsize=(14, 9)
    )
    g.ax_heatmap.set_title("Tumor heterogeneity: sub-clonal gene expression (row-scaled, robust)",
                           fontsize=14, fontweight="bold", pad=20)
    leg = [Patch(facecolor=c, label=l) for l, c in drv_colors.items()]
    g.ax_heatmap.legend(handles=leg, title="Driver label",
                        loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return g.fig, Xz


# =======================================================
# Figure 3: Fraction of Low/High subclones per driver gene
# =======================================================
def _plot_fraction_low_high(subclone_labels: Dict[str, pd.DataFrame],
                            drivers: List[str]) -> plt.Figure:
    rows = []
    for d in drivers:
        labs = subclone_labels.get(d)
        if labs is None or labs.empty:
            continue
        n_low = int((labs["type"] == "low").sum())
        n_high = int((labs["type"] == "high").sum())
        total = max(n_low + n_high, 1)
        rows.append({"driver": d, "Low": n_low / total, "High": n_high / total})
    df = pd.DataFrame(rows)
    if df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off"); ax.text(0.5, 0.5, "No subclone labels to plot", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(df["driver"], df["Low"], label="Low", color="blue")
    ax.bar(df["driver"], df["High"], bottom=df["Low"], label="High", color="red")
    for i, row in df.iterrows():
        ax.text(i, row["Low"]/2, f"{int(row['Low']*100)}%", color="white",
                ha="center", va="center", fontsize=10)
        ax.text(i, row["Low"] + row["High"]/2, f"{int(row['High']*100)}%",
                color="white", ha="center", va="center", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction")
    ax.set_title("Fraction of Low/High subclone per driver", fontsize=14, fontweight="bold")
    ax.legend(title="Subclone Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    return fig


# =======================================================
# Figure 4: Risk-State heatmap (styled like heterogeneity)
# =======================================================
def _plot_risk_state_heatmap(
    exprs_mat: pd.DataFrame,
    subclone_labels: Dict[str, pd.DataFrame],
    drivers: List[str],
    base_gene_set: List[str] | None = None,
    n_top_mad_fallback: int = 60,
    z_clip: float = 2.5
) -> Tuple[plt.Figure, pd.DataFrame]:
    def class_indices(labs: pd.DataFrame, state: str) -> np.ndarray:
        sel = labs.loc[labs["type"] == state, "source_id"]
        ids = pd.to_numeric(sel, errors="coerce").dropna().astype(int).tolist()
        return np.array([i - 1 for i in ids if 1 <= i <= exprs_mat.shape[1]], dtype=int)

    blocks = []
    meta = []
    for drv in drivers:
        labs = subclone_labels.get(drv)
        if labs is None or labs.empty or "type" not in labs.columns or "source_id" not in labs.columns:
            continue
        for state in ["Low", "Normal", "High"]:
            idx = class_indices(labs, state.lower())
            col = exprs_mat.iloc[:, idx].mean(axis=1) if idx.size else pd.Series(np.nan, index=exprs_mat.index)
            blocks.append(col.rename((drv, state)))
            meta.append((drv, state))

    if not blocks:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off"); ax.text(0.5, 0.5, "No risk-state columns to plot", ha="center", va="center")
        return fig, pd.DataFrame(index=[], columns=[])

    M = pd.concat(blocks, axis=1)
    M.columns = pd.MultiIndex.from_tuples(M.columns, names=["Driver", "State"])

    # choose rows
    if base_gene_set:
        rows = sorted(set(base_gene_set))
        X = M.loc[M.index.intersection(rows)].copy()
    else:
        med = M.median(axis=1)
        mad = (M.sub(med, axis=0)).abs().median(axis=1)
        rows = mad.sort_values(ascending=False).head(min(n_top_mad_fallback, mad.size)).index.tolist()
        X = M.loc[rows].copy()

    # robust row z
    med = X.median(axis=1)
    mad = (X.sub(med, axis=0)).abs().median(axis=1).replace(0, np.nan)
    Xz = X.sub(med, axis=0).div(mad, axis=0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    Xz = Xz.clip(-z_clip, z_clip)

    flat_cols = [f"{drv} | {st}" for drv, st in Xz.columns.to_list()]
    Xz.columns = flat_cols
    drv_series = pd.Series([c.split(" | ")[0] for c in flat_cols], index=Xz.columns, name="Driver")
    st_series  = pd.Series([c.split(" | ")[1] for c in flat_cols], index=Xz.columns, name="State")

    base = sns.color_palette("Set2", len(drivers))
    drv_colors = {d: base[i % len(base)] for i, d in enumerate(drivers)}
    state_colors = {"Low": "#1f77b4", "Normal": "#7f7f7f", "High": "#d62728"}

    col_colors = pd.DataFrame({
        "Driver": drv_series.map(drv_colors),
        "State": st_series.map(state_colors)
    }, index=Xz.columns)

    sns.set(font_scale=0.85)
    g = sns.clustermap(
        Xz,
        cmap=sns.color_palette("vlag", as_cmap=True),
        method="average", metric="euclidean",
        linewidths=0, xticklabels=True, yticklabels=True,
        col_colors=col_colors,
        figsize=(14, 9)
    )
    g.ax_heatmap.set_title("Risk-state heatmap: sub-clonal gene expression (row-scaled, robust)",
                           fontsize=14, fontweight="bold", pad=20)

    drv_legend = [Patch(facecolor=c, label=l) for l, c in drv_colors.items()]
    st_legend  = [Patch(facecolor=c, label=l) for l, c in state_colors.items()]
    g.ax_heatmap.legend(handles=drv_legend, title="Driver label",
                        loc="upper left", bbox_to_anchor=(1.02, 1.0))
    g.ax_heatmap.legend(handles=st_legend, title="State",
                        loc="upper left", bbox_to_anchor=(1.02, 0.72))
    return g.fig, Xz


# ==================
# Main API function
# ==================
def magicsubclonal(input_csv: str,
                   genes_of_interest: List[str],
                   number_sample: int,
                   gene_column_number: int,
                   output_dir: str = "magicsubclonal_outputs") -> Dict:
    # prepare output folders
    figs_dir = _ensure_dir(os.path.join(output_dir, "figures"))
    tbls_dir = _ensure_dir(os.path.join(output_dir, "tables"))

    # 1) Load data
    df = pd.read_csv(input_csv)
    if "GeneSymbol" not in df.columns:
        raise ValueError("The CSV must have a 'GeneSymbol' column.")
    df = df.dropna(subset=["GeneSymbol"]).set_index("GeneSymbol")
    exprs_data = df.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.01)

    # save input expression matrix
    _save_csv(exprs_data, os.path.join(tbls_dir, "exprs_data.csv"))

    # 2) Build x0 (and save per-driver)
    x0_list: Dict[str, np.ndarray] = {}
    for g in genes_of_interest:
        if g not in exprs_data.index:
            print(f"Warning: {g} not found; skipping.")
            continue
        vals = exprs_data.loc[g].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size == 0:
            x0_list[g] = np.array([])
            continue
        idx = np.random.choice(vals.size, size=number_sample, replace=True)
        sampled = vals[idx] * np.random.uniform(0.9, 1.1, number_sample)
        x0_list[g] = sampled
        pd.DataFrame({"sim": np.arange(1, sampled.size + 1), "x0": sampled}) \
          .to_csv(os.path.join(tbls_dir, f"x0_list__{g}.csv"), index=False)

    # 3) Fit CME parameters
    param_estimates: Dict[str, np.ndarray] = {}
    for g, x0v in x0_list.items():
        if x0v.size == 0:
            continue
        init = np.array([0.1, 1.0, 1.0], float)
        res = optimize.minimize(_loss_function, init, args=(x0v,),
                                bounds=[(1e-4, 10), (1e-4, 50), (1e-4, 50)],
                                method="L-BFGS-B")
        param_estimates[g] = res.x
    if param_estimates:
        pd.DataFrame([
            {"driver": g, "gamma2": float(v[0]), "kr": float(v[1]), "mu": float(v[2])}
            for g, v in param_estimates.items()
        ]).to_csv(os.path.join(tbls_dir, "param_estimates.csv"), index=False)

    # 4) Simulate trajectories (and save per-driver)
    gene_simulations: Dict[str, pd.DataFrame] = {}
    for g, x0v in x0_list.items():
        if g not in param_estimates or x0v.size == 0:
            continue
        gamma2, kr, mu = map(float, param_estimates[g])
        sims = []
        for i, x0 in enumerate(x0v, start=1):
            sim = _simulate_CME(float(x0), gamma2, kr, mu, t_max=50.0)
            sim["sim"] = i
            sim["gene"] = g
            sims.append(sim)
        df_sim = pd.concat(sims, ignore_index=True)
        gene_simulations[g] = df_sim
        df_sim.to_csv(os.path.join(tbls_dir, f"gene_simulations__{g}.csv"), index=False)

    # 5) t_opt + labels at t_opt (save t_opt.csv and labels per driver)
    t_opt: Dict[str, float] = {}
    subclone_labels: Dict[str, pd.DataFrame] = {}
    for g, df_sim in gene_simulations.items():
        tstar = _score_time(df_sim)
        t_opt[g] = tstar
        at_t = df_sim[np.isclose(df_sim["time"], tstar)].copy()
        mu, sdv = at_t["x"].mean(), at_t["x"].std()
        at_t["type"] = np.where(at_t["x"] <= mu - sdv, "low",
                         np.where(at_t["x"] >= mu + sdv, "high", "normal"))
        at_t["source_id"] = at_t["sim"]  # sim -> sample proxy
        subclone_labels[g] = at_t
        at_t.to_csv(os.path.join(tbls_dir, f"subclone_labels__{g}.csv"), index=False)

    if t_opt:
        pd.DataFrame([{"driver": g, "t_opt": float(t)} for g, t in t_opt.items()]) \
          .to_csv(os.path.join(tbls_dir, "t_opt.csv"), index=False)

    # 6) Demo filtered subclonal genes (replace with your selection if available)
    rng = np.random.default_rng(42)
    filtered_rows = []
    for drv in [d for d in genes_of_interest if d in gene_simulations]:
        picks = rng.choice(exprs_data.index, size=min(6, exprs_data.shape[0]), replace=False)
        for fg in picks:
            filtered_rows.append({
                "driver": drv,
                "gene": str(fg),
                "status": rng.choice(["enriched", "depleted"], p=[0.6, 0.4])
            })
    filtered_subclonal_genes_all = pd.DataFrame(filtered_rows)
    filtered_subclonal_genes_all.to_csv(os.path.join(tbls_dir, "filtered_subclonal_genes_all.csv"), index=False)

    # =========
    # FIGURES
    # =========
    fig_dynamics = _plot_cme_dynamics_all(gene_simulations, t_opt)
    fig_dynamics.savefig(os.path.join(figs_dir, "CME_Dynamics_AllDrivers.png"),
                         dpi=300, bbox_inches="tight")

    fig_highrisk = _plot_high_risk_subclone_genes_panel(
        gene_sims=gene_simulations,
        t_opt_map=t_opt,
        filtered_df=filtered_subclonal_genes_all,
        exprs_mat=exprs_data,
        subclone_labels=subclone_labels,
        drivers=[d for d in genes_of_interest if d in gene_simulations],
        ncol=gene_column_number,
        max_labels_per_driver=12
    )
    fig_highrisk.savefig(os.path.join(figs_dir, "HighRiskSubcloneGenes_noOverlap.png"),
                         dpi=300, bbox_inches="tight")

    fig_heatmap, hetero_mat = _plot_tumor_heterogeneity_heatmap_upgraded(
        exprs_mat=exprs_data,
        subclone_labels=subclone_labels,
        filtered_df=filtered_subclonal_genes_all,
        drivers=[d for d in genes_of_interest if d in gene_simulations],
        n_top_mad=60,
        z_clip=2.5
    )
    fig_heatmap.savefig(os.path.join(figs_dir, "TumorHeterogeneity_annotated.png"),
                        dpi=300, bbox_inches="tight")
    if isinstance(hetero_mat, pd.DataFrame) and not hetero_mat.empty:
        _save_csv(hetero_mat, os.path.join(tbls_dir, "heterogeneity_matrix_rowZ.csv"))

    fig_frac = _plot_fraction_low_high(subclone_labels,
                                       [d for d in genes_of_interest if d in gene_simulations])
    fig_frac.savefig(os.path.join(figs_dir, "Fraction_LowHigh_Subclones.png"),
                     dpi=300, bbox_inches="tight")

    # NEW: Risk-state heatmap (styled like heterogeneity) + export its matrix
    base_set = (filtered_subclonal_genes_all["gene"].tolist()
                if not filtered_subclonal_genes_all.empty else None)
    fig_riskstate, risk_mat = _plot_risk_state_heatmap(
        exprs_mat=exprs_data,
        subclone_labels=subclone_labels,
        drivers=[d for d in genes_of_interest if d in gene_simulations],
        base_gene_set=base_set,
        n_top_mad_fallback=60,
        z_clip=2.5
    )
    fig_riskstate.savefig(os.path.join(figs_dir, "RiskClassHeatmap.png"),
                          dpi=300, bbox_inches="tight")
    if isinstance(risk_mat, pd.DataFrame) and not risk_mat.empty:
        _save_csv(risk_mat, os.path.join(tbls_dir, "risk_state_matrix_rowZ.csv"))

    # 8) Heterogeneity + simple risk scoring (export both)
    heterogeneity = []
    for gene in exprs_data.index:
        v = exprs_data.loc[gene].values
        mean_v = float(np.mean(v))
        var_v = float(np.var(v))
        sd_v = float(np.sqrt(var_v))
        cv_v = float(sd_v / max(mean_v, 1e-12))
        P = v / max(np.sum(v), 1e-12)
        P = np.where(P > 0, P, 1.0)
        entropy_v = float(-np.sum(P * np.log(P)))
        heterogeneity.append([gene, mean_v, var_v, sd_v, cv_v, entropy_v])
    heterogeneity_df = pd.DataFrame(
        heterogeneity, columns=["GENE", "mean", "var", "sd", "cv", "entropy"]
    )
    _save_csv(heterogeneity_df, os.path.join(tbls_dir, "heterogeneity_metrics.csv"))

    adv = heterogeneity_df.copy()
    adv["log2_cv"] = np.log2(np.maximum(adv["cv"], 1e-12))
    adv["entropy_n"] = adv["entropy"] / max(np.log(exprs_data.shape[1]), 1e-12)
    def mm(v):
        v = np.asarray(v, dtype=float)
        if not np.any(np.isfinite(v)): return np.zeros_like(v)
        lo, hi = np.nanmin(v), np.nanmax(v)
        return np.zeros_like(v) if hi - lo == 0 else (v - lo) / (hi - lo)
    adv["score"] = 0.6 * mm(adv["log2_cv"]) + 0.4 * mm(adv["entropy_n"])
    adv["high_risk"] = adv["score"] >= 0.7
    adv.to_csv(os.path.join(tbls_dir, "advantage_scores.csv"), index=False)

    # manifest.json (paths to everything)
    manifest = {
        "figures": {
            "CME_Dynamics_AllDrivers": os.path.join(figs_dir, "CME_Dynamics_AllDrivers.png"),
            "HighRiskSubcloneGenes_noOverlap": os.path.join(figs_dir, "HighRiskSubcloneGenes_noOverlap.png"),
            "TumorHeterogeneity_annotated": os.path.join(figs_dir, "TumorHeterogeneity_annotated.png"),
            "Fraction_LowHigh_Subclones": os.path.join(figs_dir, "Fraction_LowHigh_Subclones.png"),
            "RiskClassHeatmap": os.path.join(figs_dir, "RiskClassHeatmap.png"),
        },
        "tables": {
            "exprs_data": os.path.join(tbls_dir, "exprs_data.csv"),
            "param_estimates": os.path.join(tbls_dir, "param_estimates.csv") if param_estimates else None,
            "t_opt": os.path.join(tbls_dir, "t_opt.csv") if t_opt else None,
            "filtered_subclonal_genes_all": os.path.join(tbls_dir, "filtered_subclonal_genes_all.csv"),
            "heterogeneity_metrics": os.path.join(tbls_dir, "heterogeneity_metrics.csv"),
            "advantage_scores": os.path.join(tbls_dir, "advantage_scores.csv"),
            "heterogeneity_matrix_rowZ": os.path.join(tbls_dir, "heterogeneity_matrix_rowZ.csv") if not (hetero_mat.empty if isinstance(hetero_mat, pd.DataFrame) else True) else None,
            "risk_state_matrix_rowZ": os.path.join(tbls_dir, "risk_state_matrix_rowZ.csv") if not (risk_mat.empty if isinstance(risk_mat, pd.DataFrame) else True) else None,
            "x0_list_per_driver": {g: os.path.join(tbls_dir, f"x0_list__{g}.csv") for g in x0_list.keys()},
            "gene_simulations_per_driver": {g: os.path.join(tbls_dir, f"gene_simulations__{g}.csv") for g in gene_simulations.keys()},
            "subclone_labels_per_driver": {g: os.path.join(tbls_dir, f"subclone_labels__{g}.csv") for g in subclone_labels.keys()},
        }
    }
    with open(os.path.join(output_dir, "manifest.json"), "w") as f:
        json.dump(_np_to_py(manifest), f, indent=2)

    # return everything
    figures = {
        "cme_dynamics_all": fig_dynamics,
        "high_risk_panel": fig_highrisk,
        "heterogeneity_heatmap": fig_heatmap,
        "fraction_low_high": fig_frac,
        "risk_class_heatmap": fig_riskstate,
    }
    return {
        "exprs_data": exprs_data,
        "x0_list": x0_list,
        "param_estimates": param_estimates,
        "gene_simulations": gene_simulations,
        "t_opt": t_opt,
        "subclone_labels": subclone_labels,
        "filtered_subclonal_genes_all": filtered_subclonal_genes_all,
        "heterogeneity": heterogeneity_df,
        "advantage_scores": adv,
        "figures": figures,
        "output_dir": output_dir,
        "tables_dir": tbls_dir,
        "figures_dir": figs_dir,
        "manifest": manifest
    }


if __name__ == "__main__":
    print("This module is intended to be called from run_magicsubclonal.py")

