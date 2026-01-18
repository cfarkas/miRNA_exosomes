#!/usr/bin/env python3
# ===============================================================
# tf_immune_heatmaps_ranked_table.py  (additional/legacy utility)
#
# Post-process g:Profiler outputs (GO_Terms_*.csv)
#   - TF heatmap (>= N miRNA presence)
#   - Immune terms heatmap (>= N miRNA presence) using a STRICT synonym list
#   - Sankey: TF -> Immune terms (simple co-occurrence across miRNAs)
#   - Ranked immune pathways table (GO:BP + REAC) with wide regex filter
#
# Usage:
#   python tf_immune_heatmaps_ranked_table.py --input <GO_Terms_dir> --outdir <outdir>
#
# Example:
#   python tf_immune_heatmaps_ranked_table.py \
#     --input ./figures/without_sample2/enrichment/gprofiler/GO_Terms \
#     --outdir ./figures/without_sample2/enrichment/summary
#
# Note:
#   This script is included because you asked to "join this script too".
#   For a *more legible* Sankey (top-N nodes, fixed left/right layout, counts in labels),
#   prefer: tf_immune_analysis.py (also included in this package).
# ===============================================================

import os
import glob
import re
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except Exception:
    HAS_PLOTLY = False


# -----------------------------
# TF DICTIONARY OF SYNONYMS
# -----------------------------
TRANSCRIPTION_FACTORS = {
    "FOXP3":  ["FOXP3", "Foxp3", "FKH10", "scurfin"],
    "GATA3":  ["GATA3", "GATA-3"],
    "TBX21":  ["T-bet", "Tbet", "TBX21"],
    "RORC":   ["RORC", "RORγt", "RORgt", "ROR-gamma", "NR1F3"],
    "BCL6":   ["BCL6", "B-cell lymphoma 6"],
    "MAF":    ["c-Maf", "MAF"],
    "BACH2":  ["BACH2"],
    "PRDM1":  ["BLIMP1", "Blimp-1", "PRDM1"],
    "AHR":    ["AHR", "aryl hydrocarbon receptor", "AhR"],
    "EGR2":   ["EGR2", "EGR-2"],
    "EGR3":   ["EGR3", "EGR-3"],
    "IKZF2":  ["Helios", "IKZF2"],
    "IKZF1":  ["Ikaros", "IKZF1"],
    "BCL11B": ["BCL-11B", "BCL11B"],
    "NFIL3":  ["NFIL3", "E4BP4"],
    "FOXO1":  ["FOXO1", "Fkhr", "FKHRL1", "FOXO1A"],
    "IRF4":   ["IRF4"],
    "STAT5":  ["STAT5", "Stat5"],
    "NFAT":   ["NFAT", "Nuclear factor of activated T-cells"],
    "AP-1":   ["AP-1", "AP1", "FOS", "JUN"],
    "NF-κB":  ["NF-kB", "NF-kappaB", "RelA", "RelB", "c-Rel", "NFκB", "NFKB"],
}

TF_PATTERNS = {}
for tf, syns in TRANSCRIPTION_FACTORS.items():
    escaped = [re.escape(s) for s in syns]
    TF_PATTERNS[tf] = re.compile("|".join(escaped), re.IGNORECASE)


def find_matched_tfs(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    matched = []
    for tf, pat in TF_PATTERNS.items():
        if pat.search(text):
            matched.append(tf)
    return matched


def make_tf_label(matched_list, original_term_id):
    if not matched_list:
        return original_term_id
    return "TF: " + ", ".join(matched_list)


# -----------------------------
# IMMUNE DICTIONARY (strict list)
# -----------------------------
IMMUNE_SYNONYMS_STRICT = [
    "positive regulation of innate immune response",
    "somatic diversification of immune receptors",
    "lymphocyte activation",
    "lymphocyte differentiation",
    "regulation of lymphocyte differentiation",
    "regulation of lymphocyte activation",
    "positive regulation of lymphocyte activation",
    "positive regulation of lymphocyte differentiation",
    "appendix; lymphoid tissue[≥Low]",
    "appendix; lymphoid tissue[≥Medium]",
    "appendix; lymphoid tissue[High]",
    "skin 1; lymphocytes[≥Low]",
    "skin 1; lymphocytes[≥Medium]",
    "skin 1; lymphocytes[≥High]",
    "skin 2; lymphocytes[≥Low]",
    "skin 2; lymphocytes[≥Medium]",
    "skin 2; lymphocytes[≥High]",
    "rectum; mucosal lymphoid cells[High]",
    "miR targeted genes in lymphocytes",
]
IMMUNE_PATTERN_STRICT = re.compile("|".join(re.escape(s) for s in IMMUNE_SYNONYMS_STRICT), re.IGNORECASE)


def row_matches_immune_strict(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    return IMMUNE_PATTERN_STRICT.search(text) is not None


# -----------------------------
# IMMUNE REGEX (wide) for ranked table
# -----------------------------
IMMUNE_REGEX_WIDE = (
    r"(Th17|Th1|Th2|T[\s-]?cell|B[\s-]?cell|Treg|regulatory T|"
    r"lymphocyte|leukocyte|immune|immun|inflamm|"
    r"cytokine|chemokine|interleukin|interferon|"
    r"NF[-\s]?kappaB|NF[-\s]?κB|NFκB|NFKB|"
    r"toll[-\s]?like|\bTLR\b|JAK[-\s]?STAT|\bTNF\b|"
    r"antigen|MHC|"
    r"IL[-\s]?\d+)"
)
IMMUNE_PATTERN_WIDE = re.compile(IMMUNE_REGEX_WIDE, re.IGNORECASE)


def row_matches_immune_wide(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    return IMMUNE_PATTERN_WIDE.search(text) is not None


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Directory containing GO_Terms_*.csv")
    ap.add_argument("--outdir", required=True, help="Output directory for plots/tables")
    ap.add_argument("--min_mirnas", type=int, default=2, help="Min miRNA presence threshold (default: 2)")
    ap.add_argument("--p_cutoff", type=float, default=0.05, help="p_value cutoff for ranked immune table (default: 0.05)")
    ap.add_argument("--valid_sources", default="GO:BP,GO:MF,TF,REAC",
                    help="Comma-separated list of sources to keep for heatmaps (default: GO:BP,GO:MF,TF,REAC)")
    args = ap.parse_args()

    in_dir = args.input
    outdir = args.outdir
    min_mirnas = args.min_mirnas
    p_cutoff = args.p_cutoff
    valid_sources = [x.strip() for x in args.valid_sources.split(",") if x.strip()]

    ensure_dir(outdir)
    out_heatmaps = ensure_dir(os.path.join(outdir, "heatmaps"))
    out_sankey = ensure_dir(os.path.join(outdir, "sankey"))
    out_tables = ensure_dir(os.path.join(outdir, "tables"))
    out_combined = ensure_dir(os.path.join(outdir, "combined"))

    pattern = os.path.join(in_dir, "GO_Terms_*.csv")
    file_list = sorted(glob.glob(pattern))
    if not file_list:
        raise SystemExit(f"No files found matching: {pattern}")

    all_dfs = []
    for fname in file_list:
        base = os.path.splitext(os.path.basename(fname))[0]
        mirna = base.replace("GO_Terms_", "")
        df = pd.read_csv(fname)
        df["miRNA"] = mirna
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_path = os.path.join(out_combined, "combined_gprofiler_terms.csv")
    combined_df.to_csv(combined_path, index=False)

    required_cols = {"source", "term_id", "term_name", "p_value", "miRNA"}
    missing = required_cols - set(combined_df.columns)
    if missing:
        raise ValueError(f"Missing columns in combined_df: {missing}.")

    filtered_df = combined_df[combined_df["source"].isin(valid_sources)].copy()

    # ---- Heatmap 1: TF
    filtered_df["MatchedTFs"] = filtered_df.apply(
        lambda r: find_matched_tfs(r.get("term_name"), r.get("term_id")),
        axis=1
    )
    focused_df = filtered_df[filtered_df["MatchedTFs"].apply(lambda x: len(x) > 0)].copy()
    focused_df["TFlabel"] = focused_df.apply(
        lambda r: make_tf_label(r["MatchedTFs"], r.get("term_id")),
        axis=1
    )
    focused_df["presence"] = 1
    focused_df.drop_duplicates(subset=["TFlabel", "miRNA"], inplace=True)

    pivot_tf = focused_df.pivot_table(index="TFlabel", columns="miRNA",
                                      values="presence", fill_value=0, aggfunc="max")
    pivot_tf["row_sum"] = pivot_tf.sum(axis=1)
    subset_tf = pivot_tf[pivot_tf["row_sum"] >= min_mirnas].drop(columns=["row_sum"])

    tf_heatmap_path = os.path.join(out_heatmaps, "heatmap_with_TF_names.pdf")
    if subset_tf.empty:
        print(f"No TF terms found in >= {min_mirnas} miRNAs. Skipping TF heatmap.")
    else:
        sns.set(style="white", font_scale=1.3)
        g = sns.clustermap(subset_tf, linewidths=0.5, linecolor="gray",
                           figsize=(16, 8), cbar_kws={"label": "Presence (1) / Absence (0)"})
        g.fig.suptitle("Detected TFs (presence/absence)", y=1.02)
        g.savefig(tf_heatmap_path, format="pdf", bbox_inches="tight")
        plt.close()

    # ---- Heatmap 2: immune strict
    filtered_df["ImmuneMatchStrict"] = filtered_df.apply(
        lambda r: row_matches_immune_strict(r.get("term_name"), r.get("term_id")),
        axis=1
    )
    immune_df = filtered_df[filtered_df["ImmuneMatchStrict"]].copy()
    immune_df["presence"] = 1
    immune_df.drop_duplicates(subset=["term_id", "miRNA"], inplace=True)

    pivot_immune = immune_df.pivot_table(index="term_id", columns="miRNA",
                                         values="presence", fill_value=0, aggfunc="max")
    pivot_immune["row_sum"] = pivot_immune.sum(axis=1)
    subset_immune = pivot_immune[pivot_immune["row_sum"] >= min_mirnas].drop(columns=["row_sum"])

    term2name = {}
    for _, row in immune_df.iterrows():
        tid = row.get("term_id")
        tname = row.get("term_name")
        term2name[tid] = f"{tid} {tname}" if pd.notnull(tid) and pd.notnull(tname) else tid
    subset_immune.index = [term2name.get(t, t) for t in subset_immune.index]

    immune_heatmap_path = os.path.join(out_heatmaps, "heatmap_immune_terms.pdf")
    if subset_immune.empty:
        print(f"No strict immune terms found in >= {min_mirnas} miRNAs. Skipping immune heatmap.")
    else:
        sns.set(style="white", font_scale=1.3)
        cmap_discrete = ListedColormap(["lightgrey", "darkred"])
        norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap_discrete.N)
        g2 = sns.clustermap(subset_immune, cmap=cmap_discrete, norm=norm,
                            linewidths=0.5, linecolor="gray",
                            figsize=(16, 8),
                            cbar_kws={"label": "Presence (1) / Absence (0)"})
        g2.fig.suptitle("Detected Immune Terms (strict list)", y=1.02)
        g2.savefig(immune_heatmap_path, format="pdf", bbox_inches="tight")
        plt.close()

    # ---- Sankey (simple)
    sankey_path = os.path.join(out_sankey, "tf_immune_sankey.html")
    if not HAS_PLOTLY:
        print("Plotly not available; skipping Sankey.")
    elif subset_tf.empty or subset_immune.empty:
        print("Cannot build Sankey; no TF or immune data found.")
    else:
        tf_labels = list(subset_tf.index)
        immune_labels = list(subset_immune.index)
        pair_counts = defaultdict(int)
        common_cols = set(subset_tf.columns).intersection(set(subset_immune.columns))
        for col in common_cols:
            tf_present = subset_tf[subset_tf[col] == 1].index
            immune_present = subset_immune[subset_immune[col] == 1].index
            for tf in tf_present:
                for im in immune_present:
                    pair_counts[(tf, im)] += 1

        all_nodes = tf_labels + immune_labels
        node_index = {label: i for i, label in enumerate(all_nodes)}
        links = [(node_index[tf], node_index[im], v) for (tf, im), v in pair_counts.items() if v > 0]

        if links:
            fig = go.Figure(data=[go.Sankey(
                node=dict(pad=15, thickness=20, line=dict(color="black", width=0.5), label=all_nodes),
                link=dict(source=[s for s, _, _ in links],
                          target=[t for _, t, _ in links],
                          value=[v for _, _, v in links])
            )])
            fig.update_layout(title_text="Sankey: TF (left) -> Immune Terms (right)", font_size=14)
            fig.write_html(sankey_path)
            print(f"Sankey saved: {sankey_path}")
        else:
            print("No TF -> Immune co-occurrences found. Skipping Sankey.")

    # ---- Ranked immune pathways table
    enrich_sources_for_table = {"GO:BP", "REAC"}
    df_sig = combined_df[
        (combined_df["source"].isin(enrich_sources_for_table)) &
        (pd.to_numeric(combined_df["p_value"], errors="coerce") <= p_cutoff)
    ].copy()

    df_sig["ImmuneMatchWide"] = df_sig.apply(
        lambda r: row_matches_immune_wide(r.get("term_name"), r.get("term_id")),
        axis=1
    )
    immune_hits = df_sig[df_sig["ImmuneMatchWide"]].copy()
    ranked_path = os.path.join(out_tables, "immune_pathways_ranked_GO_REAC.csv")

    if not immune_hits.empty:
        ranked = (immune_hits.groupby(["source", "term_id", "term_name"], as_index=False)
                  .agg(best_p_value=("p_value", "min"),
                       n_miRNAs=("miRNA", "nunique"),
                       miRNAs=("miRNA", lambda x: "; ".join(sorted(set(map(str, x)))))))
        best = pd.to_numeric(ranked["best_p_value"], errors="coerce").fillna(1.0)
        best = np.maximum(best, np.finfo(float).tiny)
        ranked["minus_log10_p"] = -np.log10(best)
        ranked = ranked.sort_values(["best_p_value", "n_miRNAs"], ascending=[True, False]).reset_index(drop=True)
        ranked.insert(0, "rank", ranked.index + 1)
        ranked.to_csv(ranked_path, index=False)
        print(f"Ranked immune pathway table saved to: {ranked_path}")
    else:
        print(f"No immune pathways found with p_value <= {p_cutoff} in GO:BP/REAC.")

    print("DONE.")


if __name__ == "__main__":
    main()
