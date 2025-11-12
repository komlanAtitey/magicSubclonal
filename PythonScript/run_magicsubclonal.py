#!/usr/bin/env python3
# run_magicsubclonal.py

import argparse
import os
import re
import sys
import json
import time
import matplotlib
import matplotlib.pyplot as plt

from magicsubclonal import magicsubclonal

def parse_genes_arg(s: str) -> list[str]:
    if not s:
        return []
    import re
    s = re.sub(r"[;, \n\t]+", ",", s.strip())
    return [g for g in (x.strip() for x in s.split(",")) if g]

def load_genes(args) -> list[str]:
    genes = []
    if args.genes:
        genes = parse_genes_arg(args.genes)
    if args.genes_file:
        with open(args.genes_file, "r") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    genes.append(line)
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            out.append(g); seen.add(g)
    return out

def main():
    parser = argparse.ArgumentParser(description="Run magicsubclonal pipeline.")
    parser.add_argument("--csv", required=True, help="Path to ovarian_data.csv (first column must be GeneSymbol)")
    parser.add_argument("--genes", help='Inline list, e.g. --genes "TP53,BRCA1,BRCA2"', default=None)
    parser.add_argument("--genes-file", help="File with one gene per line", default=None)
    parser.add_argument("--samples", type=int, default=200, help="number_sample for resampling (default 200)")
    parser.add_argument("--ncol", type=int, default=3, help="facet columns for figures")
    parser.add_argument("--outdir", default="magicsubclonal_figures", help="Directory to save figures/outputs")
    parser.add_argument("--no-show", action="store_true", help="Do not display figures interactively")
    parser.add_argument("--save-csvs", action="store_true", help="Save key tables as CSV")
    parser.add_argument("--save-json", action="store_true", help="Save a compact JSON summary")
    parser.add_argument("--topk", type=int, default=10, help="Top-K genes to print per driver")
    args = parser.parse_args()

    genes = load_genes(args)
    if not genes:
        print("ERROR: No genes provided. Use --genes or --genes-file.", file=sys.stderr)
        sys.exit(2)

    # --- FIX: resolve CSV path BEFORE changing directories ---
    csv_abs = os.path.abspath(args.csv)

    # Make and enter output directory (artifacts will land here)
    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    print("Running magicsubclonal ...")
    # Pass the absolute CSV path so it still works after chdir
    results = magicsubclonal(csv_abs, genes, number_sample=args.samples, gene_column_number=args.ncol)

    print("\n✅ Analysis complete!")
    print("\nTop estimated CME parameters:")
    for g, p in results["param_estimates"].items():
        print(f"  {g}: gamma2={p[0]:.3f}, kr={p[1]:.3f}, mu={p[2]:.3f}")

    fsg_all = results.get("filtered_subclonal_genes_all")
    if isinstance(fsg_all, dict):
        for d, lst in fsg_all.items():
            print(f"\nDriver {d} — kept (ribbon): {', '.join(lst[:args.topk]) if lst else '(none)'}")
    else:
        print("\nKept subclonal genes (ribbon) per driver:")
        if fsg_all is not None and not fsg_all.empty:
            for d in genes:
                sub = fsg_all[fsg_all["driver"] == d]["gene"].unique().tolist()
                print(f"  {d}: {', '.join(sub[:args.topk]) if sub else '(none)'}")
        else:
            print("  (none)")

    if args.save_csvs:
        results["exprs_data"].to_csv("exprs_data.csv")
        if results.get("filtered_subclonal_genes_all") is not None:
            results["filtered_subclonal_genes_all"].to_csv("filtered_subclonal_genes_all.csv", index=False)
        if results.get("heterogeneity") is not None:
            results["heterogeneity"].to_csv("heterogeneity.csv", index=False)
        print("💾 Saved CSVs")

    if args.save_json:
        summary = {
            "genes": genes,
            "param_estimates": {g: [float(x) for x in p] for g, p in results["param_estimates"].items()},
            "t_opt": {g: float(t) for g, t in results["t_opt"].items()},
        }
        with open("summary.json", "w") as fh:
            import json
            json.dump(summary, fh, indent=2)
        print("💾 Saved summary.json")

    figs = results.get("figures", {})
    if not args.no_show and len(figs):
        try:
            matplotlib.get_backend()
        except Exception:
            os.environ["MPLBACKEND"] = "TkAgg"
        for name, fig in figs.items():
            print(f"🖼️ Displaying: {name} ... (close window to continue)")
            fig.show()
            plt.show(block=True)
            time.sleep(0.2)

    print("\n🎉 All done!")
    print(f"Artifacts saved in: {os.path.abspath(os.getcwd())}")

if __name__ == "__main__":
    main()

