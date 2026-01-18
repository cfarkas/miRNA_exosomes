# miRNA differential expression pipeline (RUVSeq + edgeR)

End-to-end pipeline to analyze miRNA count matrices (e.g., small RNA-seq counts) comparing **Control (Ctrl)** vs **Oligomycin (Olig)**.

## What this repository provides

- **Normalization** with **RUVSeq** (`RUVr`) to estimate/remove unwanted variation.
- **Differential expression** with **edgeR**:
  - non-paired exact test (always)
  - paired GLM (optional; automatically enabled when pairs can be inferred from sample name suffixes)
- **Single publication-ready heatmap** (default **viridis**) with unclipped labels.
  - A *paired-order* heatmap exists but is **OFF by default**.
- **Volcano plots** using **EnhancedVolcano** (labels restricted to significant miRNAs + selected top hits).
- **Target retrieval** with **multiMiR**.
- **Functional enrichment** per miRNA using **g:Profiler** (`gprofiler2`) for **GO:BP / GO:MF / TF / REAC**.
- **Seed (m8) conservation plot** vs mouse miRNAs (plots **only positives**, auto-scaled height).
- **Python post-processing**: TF/immune presence–absence heatmaps + improved Sankey + ranked immune pathways table.

> Note: multiMiR and g:Profiler require an internet connection.

## Repository structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── env/
│   ├── RUVSeq_env.yml
│   └── requirements.txt
├── scripts/
│   ├── run_mirna_pipeline.R
│   ├── tf_immune_analysis.py
│   ├── tf_immune_heatmaps_ranked_table.py
│   └── cross_species_mouse_seed_enrichment.R
├── resources/
│   └── miR_Family_Info.txt.zip
├── data/
│   └── miRNA-counts_example.csv
└── docs/
    ├── Commands_2026.txt
    └── miRNA_without_sample2_report.docx
```

## 1) Installation

### 1.1 Create the conda environment (R + Bioconductor)

```bash
git clone https://github.com/cfarkas/miRNA_exosomes.git
cd miRNA_exosomes

conda env create -f env/RUVSeq_env.yml
conda activate RUVSeq_env
```

### 1.2 Python dependencies (optional, for TF/immune post-processing)

```bash
pip install -r env/requirements.txt
```

## 2) Input format

`--input` must be a CSV where:
- rows = miRNAs
- columns = samples
- entries = raw counts
- the first column is used as rownames (miRNA IDs)

### Sample naming conventions

The pipeline infers groups from sample names:
- sample names must contain **Ctrl** or **Olig** (case-insensitive)

Paired analysis is attempted when sample names end with a numeric suffix:
- `Ctrl_1`, `Olig_1`, `Ctrl_2`, `Olig_2`, ...

## 3) Run the analysis

### 3.1 Recommended run for this project: **without sample2**

This reproduces the analysis excluding `Ctrl_2` and `Olig_2`. Inside ```miRNA_exosomes``` repository

```bash
mkdir -p results

Rscript scripts/run_mirna_pipeline.R \
  --input data/miRNA-counts_example.csv \
  --outdir results/without_sample2 \
  --drop_samples Ctrl_2,Olig_2
```

Key outputs:
- `results/without_sample2/03_heatmaps/DEG_heatmap_grouped.pdf`  (default **viridis**, unclipped labels)
- `results/without_sample2/04_volcano/Volcano_EnhancedVolcano_nonpaired.pdf`
- `results/without_sample2/04_volcano/Volcano_EnhancedVolcano_paired.pdf` (only if pairs can be inferred)
- `results/without_sample2/07_seed_conservation/Seed_m8_conservation_up_positive.pdf`

### 3.2 All samples. 

Inside ```miRNA_exosomes``` repository: 

```bash
Rscript scripts/run_mirna_pipeline.R \
  --input data/miRNA-counts_example.csv \
  --outdir results/all_samples
```

## 4) Python post-processing (TF/immune heatmaps + Sankey)

Inside ```miRNA_exosomes``` repository:

```bash
python3 scripts/tf_immune_analysis.py \
  --input results/without_sample2/06_enrichment/GO_Terms \
  --outdir results/without_sample2/07_tf_immune \
  --sankey_top_tfs 10 \
  --sankey_top_terms 8 \
  --sankey_min_link 2 \
  --sankey_palette Set3
```

Outputs:
- `results/without_sample2/07_tf_immune/heatmaps/*.pdf`
- `results/without_sample2/07_tf_immune/sankey/tf_immune_sankey.html`
- `results/without_sample2/07_tf_immune/tables/immune_pathways_ranked_GO_REAC.csv`

## 5) Optional: cross-species mouse validation (seed conservation + mouse enrichment)

This script supports cross-species validation by mapping human DE miRNAs to a TargetScan-style miRNA family table and then running a seed-scan in mouse 3'UTRs. Inside ```miRNA_exosomes``` repository

```bash
Rscript scripts/cross_species_mouse_seed_enrichment.R \
  --results_dir results/without_sample2 \
  --design non_paired \
  --direction up \
  --fdr 0.05 \
  --logfc 1 \
  --outdir results/without_sample2/08_cross_species_mouse_seedscan \
  --target_mode seedscan_utr \
  --mir_family_table resources/miR_Family_Info.txt.zip \
  --mouse_genome mm10 \
  --site_types 8mer,7mer-m8 \
  --min_sites 1
```

## 6) Useful options

### Heatmap options
- `--color` : heatmap palette (default: `viridis`)
- `--heatmap_cex_row` : row label size (default tuned to avoid clipping)
- `--heatmap_cex_col` : column label size
- `--paired_heatmap true` : (optional) generate a paired-order heatmap

Examples:
```bash
# viridis family
--color viridis
--color magma
--color inferno

# custom ramp
--color "navy,white,firebrick3"
```

### DE thresholds
- `--pCutoff` : FDR cutoff for volcano categorization (default: 0.05)
- `--FCcutoff` : |log2FC| cutoff (default: 1)
- `--fdr_heatmap` : FDR cutoff for heatmap inclusion (default: 0.05)

## License

MIT (see `LICENSE`).
