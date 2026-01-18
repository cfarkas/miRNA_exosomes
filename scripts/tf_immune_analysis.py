#!/usr/bin/env python3
# =============================================================================
# tf_immune_analysis.py
#
# Reads g:Profiler per-miRNA enrichment CSVs (GO_Terms_*.csv) and produces:
#   - TF presence/absence heatmap (blank separators, shrunken dendrogram & color key)
#   - Immune-term presence/absence heatmap (same aesthetics)
#   - Sankey TF -> Immune co-occurrence (HTML)
#   - Ranked, enriched immune pathways table (GO:BP + Reactome) in a CSV
#
# Usage:
#   python tf_immune_analysis.py --input <GO_Terms_dir> --outdir <outdir>
#
# Example:
#   python tf_immune_analysis.py \
#     --input results/without_sample2/06_enrichment/GO_Terms \
#     --outdir results/without_sample2/07_tf_immune
#
# =============================================================================

import argparse
import glob
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

try:
    import plotly.graph_objects as go
    PLOTLY_OK = True
except Exception:
    PLOTLY_OK = False

# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser(description="TF + Immune term analysis from g:Profiler GO_Terms_*.csv files")
    p.add_argument("--input", required=True, help="Directory containing GO_Terms_*.csv files")
    p.add_argument("--outdir", required=True, help="Output directory for plots/tables")
    p.add_argument("--min_presence", type=int, default=2, help="Min # miRNAs a term must appear in to keep (default: 2)")
    p.add_argument("--pval_cutoff", type=float, default=0.05, help="p_value cutoff for ranked immune table (default: 0.05)")

    # Sankey controls (for legibility)
    p.add_argument("--sankey_top_tfs", type=int, default=10,
                   help="Show top N TF nodes in Sankey by miRNA-count (default: 10). Use 0 for all.")
    p.add_argument("--sankey_top_terms", type=int, default=8,
                   help="Show top N Immune-term nodes in Sankey by miRNA-count (default: 8). Use 0 for all.")
    p.add_argument("--sankey_min_link", type=int, default=1,
                   help="Minimum co-occurrence count to draw a link (default: 1). Increase to declutter.")
    p.add_argument("--sankey_scale", type=float, default=1.0,
                   help="Multiply link widths by this constant (default: 1.0) for visibility.")
    p.add_argument("--sankey_palette", default="Set3",
                   help="Plotly qualitative palette for RIGHT nodes (e.g. Set3, Plotly, D3, Dark2, Pastel).")
    p.add_argument("--sankey_left_color", default="#d978c2",
                   help="Hex color for LEFT TF nodes (default: pink).")
    p.add_argument("--sankey_alpha", type=float, default=0.55,
                   help="Link transparency 0-1 (default: 0.55).")
    p.add_argument("--sankey_width", type=int, default=1500, help="Sankey figure width in px (default: 1500).")
    p.add_argument("--sankey_height", type=int, default=900, help="Sankey figure height in px (default: 900).")

    return p.parse_args()

# -----------------------------
# TF dictionary
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
    "NF-κB":  ["NF-kB", "NF-kappaB", "RelA", "RelB", "c-Rel", "NFκB", "NFKB"]
}

TF_PATTERNS = {tf: re.compile("|".join(re.escape(s) for s in syns), re.IGNORECASE)
               for tf, syns in TRANSCRIPTION_FACTORS.items()}

def find_matched_tfs(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    return [tf for tf, pat in TF_PATTERNS.items() if pat.search(text)]

def make_tf_label(matched_list):
    if not matched_list:
        return None
    return "TF: " + ", ".join(sorted(set(matched_list)))

# -----------------------------
# Immune terms dictionary (for presence/absence heatmap)
# -----------------------------
IMMUNE_SYNONYMS = [
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
    "miR targeted genes in lymphocytes"
]
IMMUNE_PATTERN = re.compile("|".join(re.escape(s) for s in IMMUNE_SYNONYMS), re.IGNORECASE)

def row_matches_immune(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    return IMMUNE_PATTERN.search(text) is not None

# Wide immune regex for ranked table
IMMUNE_REGEX_WIDE = re.compile(
    r"(Th17|Th1|Th2|T[\s-]?cell|B[\s-]?cell|Treg|regulatory T|"
    r"lymphocyte|leukocyte|immune|immun|inflamm|"
    r"cytokine|chemokine|interleukin|interferon|"
    r"NF[-\s]?kappaB|NF[-\s]?κB|NFκB|NFKB|"
    r"toll[-\s]?like|\bTLR\b|JAK[-\s]?STAT|\bTNF\b|"
    r"antigen|MHC|"
    r"IL[-\s]?\d+)",
    re.IGNORECASE
)

def row_matches_immune_wide(term_name, term_id):
    text = f"{term_name or ''} {term_id or ''}"
    return IMMUNE_REGEX_WIDE.search(text) is not None

# -----------------------------
# Plot helpers
# -----------------------------
def save_clustermap(df, out_pdf, title, cmap, norm=None, discrete=False):
    sns.set(style="white", font_scale=1.4)

    # Clustering requires at least 2 observations on an axis; otherwise scipy linkage can fail.
    row_cluster = df.shape[0] >= 2
    col_cluster = df.shape[1] >= 2

    # blank separators: thick white grid lines
    g = sns.clustermap(
        df,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        cmap=cmap,
        norm=norm,
        linewidths=1.6,
        linecolor="white",
        figsize=(16, 8),
        dendrogram_ratio=(0.08, 0.12),       # shrink dendrogram
        cbar_pos=(0.02, 0.83, 0.03, 0.12),   # shrink color key
        cbar_kws={"label": "Presence (1) / Absence (0)"}
    )
    g.fig.suptitle(title, y=1.02)
    g.savefig(out_pdf, bbox_inches="tight")
    plt.close()

def main():
    args = parse_args()
    in_dir = args.input
    outdir = args.outdir

    # Output subfolders
    out_heatmaps = os.path.join(outdir, "heatmaps")
    out_tables = os.path.join(outdir, "tables")
    out_sankey = os.path.join(outdir, "sankey")
    for d in (out_heatmaps, out_tables, out_sankey):
        os.makedirs(d, exist_ok=True)

    files = sorted(glob.glob(os.path.join(in_dir, "GO_Terms_*.csv")))
    if not files:
        raise FileNotFoundError(f"No GO_Terms_*.csv files found in: {in_dir}")

    all_dfs = []
    for f in files:
        base = os.path.splitext(os.path.basename(f))[0]
        mirna = base.replace("GO_Terms_", "")
        df = pd.read_csv(f)
        df["miRNA"] = mirna
        all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)

    # Some g:Profiler exports can contain NA term_id/term_name
    combined["term_id"] = combined["term_id"].astype(str)
    combined["term_name"] = combined["term_name"].astype(str)

    VALID_SOURCES = {"GO:BP", "GO:MF", "TF", "REAC"}
    df = combined[combined["source"].isin(VALID_SOURCES)].copy()

    # ------------------------------------------------------------
    # TF heatmap (presence/absence)
    # ------------------------------------------------------------
    df["MatchedTFs"] = df.apply(lambda r: find_matched_tfs(r.get("term_name"), r.get("term_id")), axis=1)
    tf_rows = df[df["MatchedTFs"].apply(lambda x: len(x) > 0)].copy()
    tf_rows["TFlabel"] = tf_rows["MatchedTFs"].apply(make_tf_label)
    tf_rows["presence"] = 1
    tf_rows = tf_rows.drop_duplicates(subset=["TFlabel", "miRNA"])

    tf_pivot = tf_rows.pivot_table(index="TFlabel", columns="miRNA", values="presence", fill_value=0, aggfunc="max")
    tf_pivot["row_sum"] = tf_pivot.sum(axis=1)
    tf_subset = tf_pivot[tf_pivot["row_sum"] >= args.min_presence].drop(columns=["row_sum"])

    if tf_subset.empty:
        print("No TF rows found (>=min_presence). Skipping TF heatmap.")
    else:
        cmap_tf = sns.light_palette("navy", as_cmap=True, reverse=True)
        save_clustermap(
            tf_subset,
            out_pdf=os.path.join(out_heatmaps, "Detected_TFs_heatmap.pdf"),
            title="Detected Transcription Factors",
            cmap=cmap_tf
        )

    # ------------------------------------------------------------
    # Immune heatmap (presence/absence)
    # ------------------------------------------------------------
    df["ImmuneMatch"] = df.apply(lambda r: row_matches_immune(r.get("term_name"), r.get("term_id")), axis=1)
    immune_rows = df[df["ImmuneMatch"]].copy()
    immune_rows["presence"] = 1
    immune_rows = immune_rows.drop_duplicates(subset=["term_id", "miRNA"])

    immune_pivot = immune_rows.pivot_table(index="term_id", columns="miRNA", values="presence", fill_value=0, aggfunc="max")
    immune_pivot["row_sum"] = immune_pivot.sum(axis=1)
    immune_subset = immune_pivot[immune_pivot["row_sum"] >= args.min_presence].drop(columns=["row_sum"])

    # Map term_id -> "term_id term_name"
    term2name = {}
    for _, r in immune_rows.iterrows():
        tid = r.get("term_id")
        tname = r.get("term_name")
        if pd.notnull(tid) and pd.notnull(tname):
            term2name[tid] = f"{tid} {tname}"
        else:
            term2name[tid] = tid
    immune_subset.index = [term2name.get(i, i) for i in immune_subset.index]

    if immune_subset.empty:
        print("No immune rows found (>=min_presence). Skipping immune heatmap.")
    else:
        cmap_discrete = ListedColormap(["lightgrey", "darkred"])
        norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap_discrete.N)
        save_clustermap(
            immune_subset,
            out_pdf=os.path.join(out_heatmaps, "Detected_ImmuneTerms_heatmap.pdf"),
            title="Detected Immune Terms",
            cmap=cmap_discrete,
            norm=norm
        )

    # ------------------------------------------------------------
    # Sankey TF -> Immune co-occurrence (optional)
    # ------------------------------------------------------------
    if PLOTLY_OK and (not tf_subset.empty) and (not immune_subset.empty):
        tf_labels = list(tf_subset.index)
        immune_labels = list(immune_subset.index)

        pair_counts = defaultdict(int)
        common_cols = set(tf_subset.columns).intersection(set(immune_subset.columns))

        for col in common_cols:
            tf_present = tf_subset[tf_subset[col] == 1].index
            immune_present = immune_subset[immune_subset[col] == 1].index
            for tf_lab in tf_present:
                for im_lab in immune_present:
                    pair_counts[(tf_lab, im_lab)] += 1

        nodes = tf_labels + immune_labels
        node_index = {lab: i for i, lab in enumerate(nodes)}
        links = [(node_index[a], node_index[b], v) for (a, b), v in pair_counts.items() if v > 0]

        if links:
                        # -----------------------------
            # Improved Sankey (more legible):
            #   - top-N TFs and top-N Immune terms
            #   - labels include counts on a new line (like the reference figure)
            #   - fixed left/right positioning, larger canvas
            #   - link colors match the RIGHT node colors (with alpha)
            # -----------------------------
            import plotly.express as px
            import plotly.io as pio

            def _evenly_spaced(n, start=0.05, end=0.95):
                if n <= 1:
                    return [0.5]
                step = (end - start) / (n - 1)
                return [start + i * step for i in range(n)]

            def _hex_to_rgba(hex_color, alpha=0.55):
                h = hex_color.lstrip("#")
                if len(h) != 6:
                    return f"rgba(120,120,120,{alpha})"
                r = int(h[0:2], 16)
                g = int(h[2:4], 16)
                b = int(h[4:6], 16)
                return f"rgba({r},{g},{b},{alpha})"

            # Node counts (miRNA presence counts)
            tf_counts = tf_subset.sum(axis=1).sort_values(ascending=False)
            immune_counts = immune_subset.sum(axis=1).sort_values(ascending=False)

            tf_keep = list(tf_counts.index)
            immune_keep = list(immune_counts.index)

            if args.sankey_top_tfs and args.sankey_top_tfs > 0:
                tf_keep = tf_keep[:args.sankey_top_tfs]
            if args.sankey_top_terms and args.sankey_top_terms > 0:
                immune_keep = immune_keep[:args.sankey_top_terms]

            # Subset for Sankey
            sank_tf = tf_subset.loc[tf_keep, tf_subset.columns]
            sank_im = immune_subset.loc[immune_keep, immune_subset.columns]

            # Recompute link counts (TF, immune) by miRNA co-occurrence
            pair_counts2 = defaultdict(int)
            common_cols2 = set(sank_tf.columns).intersection(set(sank_im.columns))
            for col in common_cols2:
                tfs_present = sank_tf[sank_tf[col] == 1].index
                ims_present = sank_im[sank_im[col] == 1].index
                for tf in tfs_present:
                    for im in ims_present:
                        pair_counts2[(tf, im)] += 1

            # Filter links
            links2 = [(a, b, v) for (a, b), v in pair_counts2.items() if v >= args.sankey_min_link]
            if not links2:
                print("No TF->Immune co-occurrences after filtering. Skipping Sankey.")
            else:
                # Order nodes by counts (desc) for readability
                tf_keep = sorted(tf_keep, key=lambda x: tf_counts.get(x, 0), reverse=True)
                immune_keep = sorted(immune_keep, key=lambda x: immune_counts.get(x, 0), reverse=True)

                # Build labels: include count on a new line (like reference)
                left_labels = [f"{tf}\n{int(tf_counts.get(tf, 0))}" for tf in tf_keep]
                right_labels = [f"{im}\n{int(immune_counts.get(im, 0))}" for im in immune_keep]

                nodes = left_labels + right_labels
                n_tf = len(left_labels)
                n_im = len(right_labels)

                node_index = {label: i for i, label in enumerate(nodes)}
                # Map original names -> labeled names
                tf_label_map = {tf: left_labels[i] for i, tf in enumerate(tf_keep)}
                im_label_map = {im: right_labels[i] for i, im in enumerate(immune_keep)}

                # Palette for RIGHT nodes
                pal_name = str(args.sankey_palette).strip().lower()
                palette_dict = {
                    "set3": px.colors.qualitative.Set3,
                    "plotly": px.colors.qualitative.Plotly,
                    "d3": px.colors.qualitative.D3,
                    "dark2": px.colors.qualitative.Dark2,
                    "pastel": px.colors.qualitative.Pastel,
                    "pastel1": px.colors.qualitative.Pastel1,
                    "bold": px.colors.qualitative.Bold,
                    "prism": px.colors.qualitative.Prism,
                    "safe": px.colors.qualitative.Safe,
                    "vivid": px.colors.qualitative.Vivid
                }
                right_base = palette_dict.get(pal_name, px.colors.qualitative.Set3)
                right_colors = [right_base[i % len(right_base)] for i in range(n_im)]
                left_color = args.sankey_left_color

                node_colors = [left_color] * n_tf + right_colors

                # Build link arrays with colors matching the RIGHT node (target)
                src_idx, tgt_idx, val_arr, link_colors = [], [], [], []
                for tf, im, v in links2:
                    if tf not in tf_keep or im not in immune_keep:
                        continue
                    s_lab = tf_label_map[tf]
                    t_lab = im_label_map[im]
                    s = node_index[s_lab]
                    t = node_index[t_lab]
                    src_idx.append(s)
                    tgt_idx.append(t)
                    val_arr.append(float(v) * float(args.sankey_scale))

                    # Target node color
                    im_pos = immune_keep.index(im)
                    link_colors.append(_hex_to_rgba(right_colors[im_pos], alpha=args.sankey_alpha))

                # Fixed positions (left / right)
                x = [0.01] * n_tf + [0.99] * n_im
                y = _evenly_spaced(n_tf, 0.05, 0.95) + _evenly_spaced(n_im, 0.05, 0.95)

                fig = go.Figure(data=[go.Sankey(
                    arrangement="fixed",
                    node=dict(
                        pad=18,
                        thickness=22,
                        line=dict(color="black", width=0.6),
                        label=nodes,
                        color=node_colors,
                        x=x,
                        y=y
                    ),
                    link=dict(
                        source=src_idx,
                        target=tgt_idx,
                        value=val_arr,
                        color=link_colors
                    )
                )])

                fig.update_layout(
                    title_text="Sankey: TF (left) → Immune GO terms (right)",
                    font_size=18,
                    width=args.sankey_width,
                    height=args.sankey_height,
                    margin=dict(l=40, r=340, t=80, b=40)
                )

                out_html = os.path.join(out_sankey, "tf_immune_sankey.html")
                fig.write_html(out_html)
                print(f"Sankey saved to: {out_html}")

                # Optional static export (requires kaleido)
                try:
                    out_png = os.path.join(out_sankey, "tf_immune_sankey.png")
                    fig.write_image(out_png, scale=2)
                    print(f"Sankey PNG saved to: {out_png}")
                except Exception:
                    print("Static PNG export skipped (install kaleido: pip install -U kaleido).")
        else:
            print("No TF->Immune co-occurrences found. Skipping Sankey.")
    else:
        if not PLOTLY_OK:
            print("Plotly not available; skipping Sankey.")
        else:
            print("TF or Immune subset empty; skipping Sankey.")

    # ------------------------------------------------------------
    # Ranked immune pathways table (GO:BP + REAC)
    # ------------------------------------------------------------
    required_cols = {"source", "term_id", "term_name", "p_value", "miRNA"}
    missing = required_cols - set(combined.columns)
    if missing:
        raise ValueError(f"Missing columns in combined data: {missing}")

    combined["p_value"] = pd.to_numeric(combined["p_value"], errors="coerce")

    ranked_df = combined[
        combined["source"].isin(["GO:BP", "REAC"]) &
        (combined["p_value"] <= args.pval_cutoff)
    ].copy()

    ranked_df["ImmuneMatchWide"] = ranked_df.apply(lambda r: row_matches_immune_wide(r.get("term_name"), r.get("term_id")), axis=1)
    immune_hits = ranked_df[ranked_df["ImmuneMatchWide"]].copy()

    if immune_hits.empty:
        print("No immune pathways found for ranked table with current filters.")
    else:
        ranked = (
            immune_hits
            .groupby(["source", "term_id", "term_name"], as_index=False)
            .agg(
                best_p=("p_value", "min"),
                n_miRNAs=("miRNA", "nunique"),
                miRNAs=("miRNA", lambda x: "; ".join(sorted(set(x))))
            )
        )
        ranked["minus_log10_p"] = -np.log10(ranked["best_p"].astype(float))
        ranked = ranked.sort_values(["best_p", "n_miRNAs"], ascending=[True, False]).reset_index(drop=True)
        ranked.insert(0, "rank", ranked.index + 1)

        ranked.to_csv(os.path.join(out_tables, "immune_pathways_ranked_GO_REAC.csv"), index=False)
        print("Saved ranked immune pathways table:",
              os.path.join(out_tables, "immune_pathways_ranked_GO_REAC.csv"))

    print("DONE. Outputs written under:", outdir)

if __name__ == "__main__":
    main()
