#!/usr/bin/env Rscript

# =============================================================================
# run_mirna_pipeline.R
#
# End-to-end miRNA pipeline:
#   - Load miRNA count matrix
#   - Optional sample dropping (e.g., Ctrl_2,Olig_2)
#   - RUVSeq (RUVr) normalization
#   - edgeR differential expression (non-paired + paired if possible)
#   - DEG heatmaps (blank group separations + shrunken dendrogram + shrunken key)
#   - Volcano plot via EnhancedVolcano (Up/Down/No diff legend) with custom labels
#   - multiMiR target retrieval
#   - g:Profiler enrichment per miRNA (GO:BP, GO:MF, TF, REAC) with robust CSV export
#
# Usage examples:
#   # All samples
#   Rscript run_mirna_pipeline.R --input miRNA-counts.csv --outdir results/all_samples
#
#   # Without sample2 (drops Ctrl_2 and Olig_2)
#   Rscript run_mirna_pipeline.R --input miRNA-counts.csv --outdir results/without_sample2 --drop_samples Ctrl_2,Olig_2
#
# Notes:
#   - Requires packages: RUVSeq, edgeR, EnhancedVolcano, multiMiR, gprofiler2, gplots, viridis
#   - Recommended: run inside your conda env (RUVSeq_env)
# =============================================================================

suppressPackageStartupMessages({
  library(RUVSeq)
  library(edgeR)
  library(gplots)
  library(viridis)
  library(dplyr)
  library(ggplot2)
  library(EnhancedVolcano)
  library(ggrepel)
  library(multiMiR)
  library(gprofiler2)
})

# -----------------------------
# Minimal CLI parser (no optparse dependency)
# -----------------------------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

input_file <- get_arg("--input", default = NA_character_)
outdir     <- get_arg("--outdir", default = "results")
drop_str   <- get_arg("--drop_samples", default = "")
fdr_heatmap <- as.numeric(get_arg("--fdr_heatmap", default = "0.05"))
pCutoff     <- as.numeric(get_arg("--pCutoff", default = "0.05"))
FCcutoff    <- as.numeric(get_arg("--FCcutoff", default = "1.0"))
k_ruv       <- as.integer(get_arg("--k_ruv", default = "1"))
label_str   <- get_arg("--label_mirnas", default =
                         "miR-3681-5p,miR-365a-5p,miR-6081,miR-193a-3p,miR-4684-5p,miR-622,miR-335-5p,miR-6883-3p,miR-3064-5p,miR-4478")
# Volcano aesthetics
volcano_width  <- as.numeric(get_arg("--volcano_width", default = "9.2"))
volcano_height <- as.numeric(get_arg("--volcano_height", default = "8.6"))
label_top_n_each_side <- as.integer(get_arg("--label_top_n", default = "3"))
max_labels_total <- as.integer(get_arg("--max_labels_total", default = "14"))
# Heatmap aesthetics (defaults tuned to match the reference viridis heatmap)
# - viridis palette by default
# - modest text sizes so labels don't get clipped
heatmap_cex_row  <- as.numeric(get_arg("--heatmap_cex_row", default = "0.60"))
heatmap_cex_col  <- as.numeric(get_arg("--heatmap_cex_col", default = "1.40"))
# Backward-compatible CLI arg (separators are no longer drawn by default)
heatmap_sepwidth <- as.numeric(get_arg("--heatmap_sepwidth", default = "0"))
heatmap_color    <- get_arg("--color", default = "viridis")


# Optional paired-order heatmap (OFF by default; enable with --paired_heatmap true)
use_paired_heatmap <- tolower(get_arg("--paired_heatmap", default = "false")) %in% c("true","t","1","yes","y")

if (is.na(input_file) || input_file == "") {
  stop("ERROR: --input is required (path to miRNA-counts.csv).")
}

# -----------------------------
# Output folder structure
# -----------------------------
dir_qc       <- file.path(outdir, "00_qc")
dir_ruv      <- file.path(outdir, "01_ruvseq")
dir_edger    <- file.path(outdir, "02_edger")
dir_heatmap  <- file.path(outdir, "03_heatmaps")
dir_volcano  <- file.path(outdir, "04_volcano")
dir_targets  <- file.path(outdir, "05_targets")
dir_enrich   <- file.path(outdir, "06_enrichment")
dir_go_terms <- file.path(dir_enrich, "GO_Terms")
dir_seed     <- file.path(outdir, "07_seed_conservation")

for (d in c(dir_qc, dir_ruv, dir_edger, dir_heatmap, dir_volcano, dir_targets, dir_enrich, dir_go_terms, dir_seed)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# -----------------------------
# Helpers
# -----------------------------
clean_mirna_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("^hsa-", "", x, ignore.case = TRUE)
  x <- gsub("^mmu-", "", x, ignore.case = TRUE)
  x <- gsub("^rno-", "", x, ignore.case = TRUE)
  x
}

extract_pair_id <- function(sample_names) {
  # Returns numeric vector if names end with _NUMBER, else NA
  nums <- suppressWarnings(as.numeric(sub(".*_(\\d+)$", "\\1", sample_names)))
  if (all(is.na(nums))) return(rep(NA_real_, length(sample_names)))
  nums
}

order_by_suffix_number <- function(x) {
  nums <- suppressWarnings(as.numeric(sub(".*_(\\d+)$", "\\1", x)))
  ord <- order(is.na(nums), nums, x)
  x[ord]
}

make_group_factor <- function(sample_names) {
  grp <- ifelse(grepl("Ctrl", sample_names, ignore.case = TRUE), "Ctrl",
                ifelse(grepl("Olig", sample_names, ignore.case = TRUE), "Olig", NA))
  if (any(is.na(grp))) {
    stop("Could not infer group for some samples. Ensure column names contain 'Ctrl' or 'Olig'.")
  }
  factor(grp, levels = c("Ctrl","Olig"))
}

plot_qq <- function(df_counts, out_pdf) {
  pdf(out_pdf, width = 25, height = 10)
  n <- ncol(df_counts)
  par(mfrow = c(2, max(2, ceiling(n/2))))
  for (i in 1:ncol(df_counts)) {
    qqnorm(df_counts[[i]], main = paste("Q-Q plot for", names(df_counts)[i]))
    qqline(df_counts[[i]], col = "steelblue")
  }
  dev.off()
}

plot_rle_pca <- function(seqset, group_factor, out_pdf, title_prefix="") {
  colors <- RColorBrewer::brewer.pal(6, "Set2")
  pdf(out_pdf, width = 6, height = 6)
  plotRLE(seqset, outline = FALSE, ylim = c(-4, 4), col = colors[group_factor],
          main = paste0(title_prefix, "RLE"))
  plotPCA(seqset, col = colors[group_factor], cex = 1.2,
          main = paste0(title_prefix, "PCA"))
  dev.off()
}

get_heatmap_colors <- function(palette_name, n = 256) {
  pal_raw <- ifelse(is.null(palette_name) || is.na(palette_name) || palette_name == "", "greenred", palette_name)
  pal <- tolower(pal_raw)

  # ------------------------------------------------------------
  # 0) User-defined ramps: comma-separated colors
  #    Example: --color "navy,white,firebrick3"
  #             --color "#440154,#21908C,#FDE725"
  # ------------------------------------------------------------
  if (grepl(",", pal_raw)) {
    cols <- trimws(unlist(strsplit(pal_raw, ",")))
    cols <- cols[cols != ""]
    if (length(cols) >= 2) {
      return(colorRampPalette(cols)(n))
    }
  }

  # Allow suffixes to reverse palettes: _r, -r, _rev, -rev
  reverse_flag <- FALSE
  if (grepl("(_r|-r|_rev|-rev)$", pal)) {
    reverse_flag <- TRUE
    pal <- gsub("(_r|-r|_rev|-rev)$", "", pal)
  }

  # ------------------------------------------------------------
  # 1) Common diverging ramps (high-contrast, publication-friendly)
  # ------------------------------------------------------------
  if (pal %in% c("greenred", "gbr", "green-black-red", "green_black_red")) {
    cols <- colorRampPalette(c("green", "black", "red"))(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("redgreen", "rbg", "red-black-green", "red_black_green")) {
    cols <- colorRampPalette(c("red", "black", "green"))(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("bluewhitered", "bwr", "blue-white-red")) {
    cols <- colorRampPalette(c("blue", "white", "red"))(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("redwhiteblue", "rwb", "red-white-blue")) {
    cols <- colorRampPalette(c("red", "white", "blue"))(n)
    return(if (reverse_flag) rev(cols) else cols)
  }

  # ------------------------------------------------------------
  # 2) Base R palettes
  # ------------------------------------------------------------
  if (pal %in% c("heat", "heatcolors", "heat.colors")) {
    cols <- heat.colors(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("terrain", "terraincolors", "terrain.colors")) {
    cols <- terrain.colors(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("topo", "topocolors", "topo.colors")) {
    cols <- topo.colors(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("cm", "cmcolors", "cm.colors")) {
    cols <- cm.colors(n)
    return(if (reverse_flag) rev(cols) else cols)
  }

  # ------------------------------------------------------------
  # 3) Viridis palettes (+ turbo via viridisLite)
  #    Options: viridis, magma, inferno, plasma, cividis, turbo
  # ------------------------------------------------------------
  if (pal %in% c("viridis", "magma", "inferno", "plasma", "cividis")) {
    cols <- get(pal, asNamespace("viridis"))(n)
    return(if (reverse_flag) rev(cols) else cols)
  }
  if (pal %in% c("turbo")) {
    # turbo is in viridisLite
    cols <- viridisLite::turbo(n)
    return(if (reverse_flag) rev(cols) else cols)
  }

  # ------------------------------------------------------------
  # 4) RColorBrewer palettes (any name in brewer.pal.info)
  #    Example: --color RdBu, Spectral, BrBG, PiYG, YlGnBu, ...
  #    Reverse with --color RdBu_r
  # ------------------------------------------------------------
  pal2 <- pal
  if (startsWith(pal2, "brewer:")) pal2 <- sub("^brewer:", "", pal2)

  suppressWarnings({
    brewer_info <- try(RColorBrewer::brewer.pal.info, silent = TRUE)
  })
  if (!inherits(brewer_info, "try-error")) {
    brewer_names <- rownames(brewer_info)
    brewer_match <- brewer_names[tolower(brewer_names) == pal2]
    if (length(brewer_match) == 1) {
      maxc <- brewer_info[brewer_match, "maxcolors"]
      base_cols <- RColorBrewer::brewer.pal(maxc, brewer_match)
      cols <- colorRampPalette(base_cols)(n)
      return(if (reverse_flag) rev(cols) else cols)
    }
  }

  warning("Unknown --color palette: ", palette_name,
          ". Using green-black-red (greenred).",
          "\nTip: you can pass comma-separated colors, e.g. --color 'navy,white,firebrick3'.")
  cols <- colorRampPalette(c("green", "black", "red"))(n)
  return(if (reverse_flag) rev(cols) else cols)
}

plot_deg_heatmap <- function(mat_log2, out_pdf, main_title, colsep_positions,
                             cex_row = 1.10, cex_col = 1.60,
                             sepwidth_val = 0.05,
                             color_mode = "greenred") {
  stopifnot(is.matrix(mat_log2))
  if (nrow(mat_log2) < 2 || ncol(mat_log2) < 2) {
    stop("Heatmap matrix must be at least 2x2. Current dim: ",
         paste(dim(mat_log2), collapse="x"))
  }

  # Layout to shrink dendrogram/key and keep key on top-right
  lmat <- rbind(c(0, 3, 4),
                c(2, 1, 0))

  # Make dendrogram & key thinner to reduce blank space
  lwid <- c(0.70, 6.6, 0.65)  # row dendro | heatmap | key
  lhei <- c(0.50, 6.6)        # top row (col dendro+key) | heatmap

  hm_cols <- get_heatmap_colors(color_mode, n = 256)

  pdf(out_pdf, width = 17, height = 16)

  # Build heatmap.2 args so we can optionally disable group separators when sepwidth_val == 0
  hm_args <- list(
    mat_log2,
    main = main_title,
    density.info = "none",
    trace = "none",
    margins = c(8, 10),
    cexRow = cex_row,
    cexCol = cex_col,
    col = hm_cols,

    # --- Shrink key ---
    key = TRUE,
    keysize = 0.60,
    key.title = "Color key",
    key.xlab = "Value",
    key.par = list(mar = c(1.5, 1.5, 0.8, 0.8)),

    # --- Shrink dendrogram/key via layout ---
    lmat = lmat, lwid = lwid, lhei = lhei
  )

  # Only draw separators if explicitly requested (sepwidth_val > 0)
  if (!is.null(colsep_positions) && length(colsep_positions) > 0 && is.finite(sepwidth_val) && sepwidth_val > 0) {
    hm_args$colsep <- colsep_positions
    hm_args$sepcolor <- "white"
    hm_args$sepwidth <- c(sepwidth_val, sepwidth_val)
  }

  do.call(heatmap.2, hm_args)
  dev.off()
}

plot_volcano_enhancedvolcano <- function(res_table,
                                         out_pdf,
                                         plot_title,
                                         pCutoff,
                                         FCcutoff,
                                         always_label,
                                         label_top_n_each_side = 3,
                                         max_labels_total = 14,
                                         pointSize = 2.1,
                                         labSize = 8.0,
                                         xlim_left = -5,
                                         xlim_right = NA_real_,
                                         width = 9.2,
                                         height = 8.6) {

  stopifnot(is.data.frame(res_table))
  stopifnot(all(c("logFC","FDR") %in% colnames(res_table)))

  res <- res_table
  res$miRNA <- rownames(res)

  # Standardize miRNA labels for plotting and matching
  res$miRNA <- clean_mirna_id(res$miRNA)
  res$FDR <- as.numeric(res$FDR)

  # Collapse duplicates (prevents repeated labels)
  res <- res %>%
    arrange(FDR, desc(abs(logFC))) %>%
    distinct(miRNA, .keep_all = TRUE)

  res <- res %>%
    mutate(
      FDR_safe = pmax(FDR, .Machine$double.xmin),
      negLog10FDR = -log10(FDR_safe),
      Regulation = case_when(
        FDR < pCutoff & logFC >  FCcutoff ~ "Up",
        FDR < pCutoff & logFC < -FCcutoff ~ "Down",
        TRUE                              ~ "No diff."
      )
    )
  res$Regulation <- factor(res$Regulation, levels = c("Up","No diff.","Down"))

  # sanitize always_label list
  always_label <- clean_mirna_id(always_label)
  always_label <- always_label[!is.na(always_label) & always_label != ""]
  always_label <- unique(always_label)

  top_up <- res %>%
    filter(Regulation == "Up") %>%
    arrange(FDR, desc(logFC)) %>%
    head(label_top_n_each_side) %>%
    pull(miRNA)

  top_down <- res %>%
    filter(Regulation == "Down") %>%
    arrange(FDR, logFC) %>%
    head(label_top_n_each_side) %>%
    pull(miRNA)

  always_present <- res %>%
    filter(miRNA %in% always_label, Regulation != "No diff.") %>%
    pull(miRNA)

  # Build label set and cap total labels to avoid clutter
  select_labels <- unique(c(always_present, top_up, top_down))
  if (length(select_labels) > max_labels_total) {
    remaining <- setdiff(select_labels, always_present)
    rem_ranked <- res %>%
      filter(miRNA %in% remaining) %>%
      arrange(FDR, desc(abs(logFC))) %>%
      pull(miRNA)
    select_labels <- unique(c(always_present, head(rem_ranked, max(0, max_labels_total - length(always_present)))))
  }

  label_df <- res %>% filter(miRNA %in% select_labels)

  # X/Y limits for nicer aesthetics
  max_x <- max(res$logFC, na.rm = TRUE)
  if (is.na(xlim_right)) {
    # extra room for right-side labels
    xlim_right <- min(12, max(6, ceiling(max_x + 2.0)))
  }
  max_y <- max(res$negLog10FDR, na.rm = TRUE)
  ylim <- c(0, max(3, max_y + 0.6))

  # Base EnhancedVolcano WITHOUT its own labels (we add ggrepel labels ourselves)
  p <- EnhancedVolcano(
    res,
    lab = rep("", nrow(res)),   # <- suppress EV labels
    x = "logFC",
    y = "FDR",
    title = plot_title,
    subtitle = "",
    xlab = "Log2(fold change)",
    ylab = "-Log10(FDR)",
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    cutoffLineType = "dashed",
    cutoffLineCol = "black",
    xlim = c(xlim_left, xlim_right),
    ylim = ylim,
    pointSize = pointSize,
    labSize = labSize,
    col = c("grey70","grey70","grey70","grey70"),
    legendPosition = "none"
  )

  # Overlay colored points + legend Up/No diff/Down
  p <- p +
    geom_point(
      data = res,
      aes(x = logFC, y = negLog10FDR, color = Regulation),
      inherit.aes = FALSE,
      size = pointSize,
      alpha = 0.95
    ) +
    scale_color_manual(
      name = "Regulation",
      values = c("Up" = "red2", "No diff." = "grey70", "Down" = "royalblue")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 15),
      plot.title   = element_text(hjust = 0.5, size = 20),
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 16),
      plot.margin  = margin(10, 45, 10, 10)  # extra right margin for labels
    ) +
    coord_cartesian(clip = "off")

  # Add repelled labels that do NOT cover red dots:
  # - move labels horizontally away from the volcano (left labels to the left, right labels to the right)
  # - constrain movement mostly along y to avoid drifting back onto the points
  if (nrow(label_df) > 0) {
    label_pos <- label_df %>% filter(logFC >= 0)
    label_neg <- label_df %>% filter(logFC <  0)

    if (nrow(label_pos) > 0) {
      p <- p + ggrepel::geom_label_repel(
        data = label_pos,
        aes(x = logFC, y = negLog10FDR, label = miRNA),
        inherit.aes = FALSE,
        nudge_x = 0.9,
        direction = "y",
        box.padding = 0.45,
        point.padding = 0.40,
        segment.size = 0.3,
        segment.alpha = 0.8,
        min.segment.length = 0,
        seed = 123,
        size = labSize/3,         # labSize was in points; convert to ggplot "size" scale
        label.size = 0.25,
        color = "black",
        fill = "white"
      )
    }

    if (nrow(label_neg) > 0) {
      p <- p + ggrepel::geom_label_repel(
        data = label_neg,
        aes(x = logFC, y = negLog10FDR, label = miRNA),
        inherit.aes = FALSE,
        nudge_x = -0.9,
        direction = "y",
        box.padding = 0.45,
        point.padding = 0.40,
        segment.size = 0.3,
        segment.alpha = 0.8,
        min.segment.length = 0,
        seed = 123,
        size = labSize/3,
        label.size = 0.25,
        color = "black",
        fill = "white"
      )
    }
  }

  ggsave(out_pdf, plot = p, width = width, height = height, device = "pdf")
}

# -----------------------------
# Seed (m8) conservation: human DE miRNAs vs mouse miRNAs
#   - Uses miRBaseConverter built-in miRBase sequences (if installed)
#   - Plots ONLY the positive/conserved ones (n_mouse_seed_matches > 0)
#     to keep the figure compact and readable.
# -----------------------------
seed_m8_from_seq <- function(seq) {
  if (is.null(seq) || is.na(seq)) return(NA_character_)
  seq <- toupper(as.character(seq))
  seq <- gsub("T", "U", seq)  # just in case
  if (nchar(seq) < 8) return(NA_character_)
  substr(seq, 2, 8)  # nt2-nt8 (7mer-m8)
}

plot_seed_m8_conservation <- function(res_table,
                                     out_pdf,
                                     selection = c("up", "down"),
                                     fdr_cut = 0.05,
                                     fc_cut = 1.0,
                                     targetVersion = "v22",
                                     positive_only = TRUE,
                                     width = 7.0,
                                     height = NA_real_) {

  selection <- match.arg(selection)

  if (!requireNamespace("miRBaseConverter", quietly = TRUE)) {
    message("Seed conservation plot skipped: miRBaseConverter not installed.")
    return(invisible(NULL))
  }

  stopifnot(is.data.frame(res_table))
  stopifnot(all(c("logFC", "FDR") %in% colnames(res_table)))

  df <- res_table
  df$miRNA <- clean_mirna_id(rownames(df))
  df$FDR <- as.numeric(df$FDR)
  df$logFC <- as.numeric(df$logFC)

  # Select DE miRNAs
  if (selection == "up") {
    sel <- df %>%
      filter(!is.na(FDR)) %>%
      filter(FDR < fdr_cut, logFC > fc_cut)
  } else {
    sel <- df %>%
      filter(!is.na(FDR)) %>%
      filter(FDR < fdr_cut, logFC < -fc_cut)
  }

  if (nrow(sel) == 0) {
    message("Seed conservation plot: no DE miRNAs for selection=", selection,
            " (FDR<", fdr_cut, ", |logFC|>", fc_cut, ")")
    return(invisible(NULL))
  }

  # Retrieve mature miRNAs for human and mouse (built-in miRBase)
  hsa_all <- miRBaseConverter::getAllMiRNAs(version = targetVersion, type = "mature", species = "hsa")
  mmu_all <- miRBaseConverter::getAllMiRNAs(version = targetVersion, type = "mature", species = "mmu")

  # Compute seeds (m8) for mouse, then count by seed
  mmu_all$seed_m8 <- vapply(mmu_all$Sequence, seed_m8_from_seq, FUN.VALUE = character(1))
  mmu_seed_counts <- mmu_all %>%
    filter(!is.na(seed_m8) & seed_m8 != "") %>%
    count(seed_m8, name = "n_mouse")

  # Map selected human miRNAs to their mature sequences to get the seed
  # NOTE: our pipeline stores miRNAs WITHOUT the 'hsa-' prefix.
  sel <- sel %>%
    mutate(hsa_full = paste0("hsa-", miRNA),
           hsa_full_lower = tolower(hsa_full))

  hsa_map <- hsa_all %>%
    mutate(Name_lower = tolower(Name),
           seed_m8 = vapply(Sequence, seed_m8_from_seq, FUN.VALUE = character(1))) %>%
    select(Name_lower, seed_m8)

  sel <- sel %>%
    left_join(hsa_map, by = c("hsa_full_lower" = "Name_lower")) %>%
    left_join(mmu_seed_counts, by = "seed_m8") %>%
    mutate(
      n_mouse = ifelse(is.na(n_mouse), 0L, as.integer(n_mouse)),
      seed_conserved = n_mouse > 0
    )

  # Keep only the positive/conserved ones if requested
  if (positive_only) {
    sel <- sel %>% filter(n_mouse > 0)
  }

  if (nrow(sel) == 0) {
    message("Seed conservation plot: no conserved seeds for selection=", selection)
    return(invisible(NULL))
  }

  # Order by significance (FDR) like the rest of the pipeline
  sel <- sel %>% arrange(FDR, desc(abs(logFC)))
  sel$miRNA <- factor(sel$miRNA, levels = rev(sel$miRNA))

  # Dynamic figure height so labels don't collide (smaller when fewer miRNAs)
  if (is.na(height)) {
    height <- max(2.8, min(9.5, 1.8 + 0.15 * nrow(sel)))
  }

  # Compact bar plot (only positives)
  p <- ggplot(sel, aes(x = miRNA, y = n_mouse)) +
    geom_col(fill = "forestgreen", color = "black", width = 0.85) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0, max(sel$n_mouse), by = 1), limits = c(0, max(sel$n_mouse))) +
    labs(
      title = "Seed (m8) conservation: human DE miRNAs vs mouse miRNAs",
      subtitle = paste0("Selection: ", selection, " (FDR<", fdr_cut, ", |logFC|>", fc_cut, ")"),
      x = NULL,
      y = "# mouse miRNAs with identical seed"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 7),
      legend.position = "none"
    )

  # Save PDF + PNG + CSV table
  ggsave(out_pdf, plot = p, width = width, height = height, device = "pdf")
  ggsave(sub("\\.pdf$", ".png", out_pdf), plot = p, width = width, height = height, dpi = 300, bg = "white")
  write.csv(sel, file = sub("\\.pdf$", ".csv", out_pdf), row.names = FALSE)

  invisible(sel)
}

# -----------------------------
# Load data
# -----------------------------
counts_raw <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Keep only Ctrl/Olig columns
counts <- counts_raw[, grepl("Ctrl|Olig", colnames(counts_raw), ignore.case = TRUE), drop = FALSE]

# Drop samples if requested
drop_samples <- character(0)
if (!is.null(drop_str) && drop_str != "") {
  drop_samples <- trimws(unlist(strsplit(drop_str, ",")))
  keep <- !colnames(counts) %in% drop_samples
  counts <- counts[, keep, drop = FALSE]
}

if (ncol(counts) < 4) {
  stop("Not enough samples after filtering/dropping. Found ", ncol(counts), " columns.")
}

# QC: QQ plots (raw)
plot_qq(counts, file.path(dir_qc, "QQ_plots_raw.pdf"))

# Infer group
group <- make_group_factor(colnames(counts))

# -----------------------------
# RUVSeq: build SeqExpressionSet
# -----------------------------
set <- newSeqExpressionSet(as.matrix(counts),
                           phenoData = data.frame(group = group, row.names = colnames(counts)))

plot_rle_pca(set, group, file.path(dir_qc, "RLE_PCA_no_normalization.pdf"), title_prefix = "")

# Estimate unwanted factors using residuals (RUVr)
design <- model.matrix(~group, data = pData(set))
y <- DGEList(counts = counts(set), group = group)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
residuals_fit <- residuals(fit, type="deviance")

genes <- rownames(counts)[grep("", rownames(counts))]
set_ruv <- RUVr(set, genes, k = k_ruv, residuals_fit)

plot_rle_pca(set_ruv, group, file.path(dir_qc, "RLE_PCA_RUVr.pdf"), title_prefix = "RUVr ")

# Save normalized counts
norm_file <- file.path(dir_ruv, "Ctrl_vs_Olig_empirical_norm.tsv")
write.table(normCounts(set_ruv), file = norm_file, quote = FALSE, sep = "\t")

# Load normalized counts
norm_counts <- read.table(norm_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# QC: QQ plots (normalized)
plot_qq(as.data.frame(norm_counts), file.path(dir_qc, "QQ_plots_RUVr.pdf"))

# -----------------------------
# edgeR DE: non-paired exactTest
# -----------------------------
dge <- DGEList(counts = norm_counts, group = group)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

et <- exactTest(dge)
res_nonpaired <- topTags(et, n = nrow(dge))$table
write.csv(res_nonpaired, file.path(dir_edger, "edgeR_results_non_paired.csv"), row.names = TRUE)

# -----------------------------
# edgeR DE: paired GLM (if pair IDs can be inferred)
# -----------------------------
pair_ids <- extract_pair_id(colnames(norm_counts))
paired_ok <- !all(is.na(pair_ids)) && length(unique(pair_ids[!is.na(pair_ids)])) >= 2

if (paired_ok) {
  pair_factor <- factor(pair_ids)
  group_glm <- group

  y_glm <- DGEList(counts = norm_counts)
  y_glm <- calcNormFactors(y_glm)
  design_glm <- model.matrix(~ pair_factor + group_glm)
  y_glm <- estimateDisp(y_glm, design_glm)
  fit_glm <- glmFit(y_glm, design_glm)

  coef_name <- "group_glmOlig"
  if (!(coef_name %in% colnames(design_glm))) {
    # fallback: pick last column
    coef_idx <- ncol(design_glm)
  } else {
    coef_idx <- which(colnames(design_glm) == coef_name)
  }

  lrt <- glmLRT(fit_glm, coef = coef_idx)
  res_paired <- topTags(lrt, n = Inf)$table
  write.csv(res_paired, file.path(dir_edger, "edgeR_results_paired.csv"), row.names = TRUE)
} else {
  message("Paired analysis skipped: could not infer pair IDs from sample names.")
}

# -----------------------------
# Heatmap of DEGs (grouped order Ctrl then Olig) with blank group separation
# -----------------------------
res_for_heat <- read.csv(file.path(dir_edger, "edgeR_results_non_paired.csv"), row.names = 1, check.names = FALSE)
res_for_heat$miRNA <- clean_mirna_id(rownames(res_for_heat))
res_for_heat$FDR <- as.numeric(res_for_heat$FDR)

# Filter by FDR; fallback if <2
sig_mirs <- res_for_heat %>%
  filter(!is.na(FDR)) %>%
  filter(FDR < fdr_heatmap) %>%
  arrange(FDR) %>%
  pull(miRNA)

if (length(sig_mirs) < 2) {
  message("WARNING: <2 DE miRNAs with FDR < ", fdr_heatmap, ". Using TOP 50 by FDR for heatmap.")
  sig_mirs <- res_for_heat %>%
    filter(!is.na(FDR)) %>%
    arrange(FDR) %>%
    head(50) %>%
    pull(miRNA)
}

rownames(norm_counts) <- clean_mirna_id(rownames(norm_counts))

# Column order grouped by condition
samples <- colnames(norm_counts)
ctrl_samples <- order_by_suffix_number(samples[grepl("Ctrl", samples, ignore.case = TRUE)])
olig_samples <- order_by_suffix_number(samples[grepl("Olig", samples, ignore.case = TRUE)])
ordered_samples <- c(ctrl_samples, olig_samples)

n_ctrl <- length(ctrl_samples)
colsep_positions <- if (n_ctrl > 0 && n_ctrl < length(ordered_samples)) n_ctrl else NULL

mat <- norm_counts[rownames(norm_counts) %in% sig_mirs, ordered_samples, drop = FALSE]
mat_log <- log2(as.matrix(mat) + 1)

plot_deg_heatmap(
  mat_log2 = mat_log,
  out_pdf = file.path(dir_heatmap, "DEG_heatmap_grouped.pdf"),
  main_title = "Control vs Olig (Log Counts)",
  colsep_positions = colsep_positions,
  cex_row = heatmap_cex_row,
  cex_col = heatmap_cex_col,
  sepwidth_val = heatmap_sepwidth,
  color_mode = heatmap_color
)

# Optional paired-order heatmap: Ctrl_i next to Olig_i, with separators between pairs
if (use_paired_heatmap && paired_ok) {
  pair_ids2 <- extract_pair_id(colnames(norm_counts))
  pair_levels <- sort(unique(pair_ids2[!is.na(pair_ids2)]))
  paired_order <- c()
  for (pid in pair_levels) {
    c_s <- colnames(norm_counts)[grepl("Ctrl", colnames(norm_counts), ignore.case = TRUE) & pair_ids2 == pid]
    o_s <- colnames(norm_counts)[grepl("Olig", colnames(norm_counts), ignore.case = TRUE) & pair_ids2 == pid]
    paired_order <- c(paired_order, c_s, o_s)
  }
  paired_order <- paired_order[paired_order %in% colnames(norm_counts)]
  # separators after each pair block (every 2 columns)
  if (length(paired_order) >= 4) {
    colsep_pairs <- seq(2, length(paired_order) - 2, by = 2)
  } else {
    colsep_pairs <- NULL
  }

  mat2 <- norm_counts[rownames(norm_counts) %in% sig_mirs, paired_order, drop = FALSE]
  mat2_log <- log2(as.matrix(mat2) + 1)

  plot_deg_heatmap(
    mat_log2 = mat2_log,
    out_pdf = file.path(dir_heatmap, "DEG_heatmap_paired_order.pdf"),
    main_title = "Control vs Olig (Log Counts) - paired order",
    colsep_positions = colsep_pairs,
    cex_row = heatmap_cex_row,
    cex_col = heatmap_cex_col,
    sepwidth_val = heatmap_sepwidth,
    color_mode = heatmap_color
  )
}

# -----------------------------
# Volcano plots (EnhancedVolcano): non-paired AND paired (if available)
# -----------------------------
labels_to_show <- unique(trimws(unlist(strsplit(label_str, ","))))
labels_to_show <- labels_to_show[labels_to_show != ""]

# 1) Non-paired volcano
res_volcano_nonpaired <- read.csv(
  file.path(dir_edger, "edgeR_results_non_paired.csv"),
  row.names = 1,
  check.names = FALSE
)

plot_volcano_enhancedvolcano(
  res_table = res_volcano_nonpaired,
  out_pdf = file.path(dir_volcano, "Volcano_EnhancedVolcano_nonpaired.pdf"),
  plot_title = paste0("Control vs Oligomycin (", basename(outdir), ") - non-paired"),
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  always_label = labels_to_show,
  label_top_n_each_side = label_top_n_each_side,
  max_labels_total = max_labels_total,
  width = volcano_width,
  height = volcano_height,
  xlim_left = -5
)

# 2) Paired volcano (only if the paired results were produced)
paired_csv <- file.path(dir_edger, "edgeR_results_paired.csv")
if (file.exists(paired_csv)) {

  res_volcano_paired <- read.csv(
    paired_csv,
    row.names = 1,
    check.names = FALSE
  )

  plot_volcano_enhancedvolcano(
    res_table = res_volcano_paired,
    out_pdf = file.path(dir_volcano, "Volcano_EnhancedVolcano_paired.pdf"),
    plot_title = paste0("Control vs Oligomycin (", basename(outdir), ") - paired GLM"),
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    always_label = labels_to_show,
    label_top_n_each_side = label_top_n_each_side,
    max_labels_total = max_labels_total,
    width = volcano_width,
    height = volcano_height,
    xlim_left = -5
  )

} else {
  message("Paired volcano skipped: paired results file not found at: ", paired_csv)
}



# -----------------------------
# Seed (m8) conservation plots (human DE miRNAs vs mouse miRNAs)
#   - Plots only the positive / conserved ones (n_mouse > 0)
#   - Uses smaller, auto-scaled plot height
# -----------------------------
tryCatch({
  plot_seed_m8_conservation(
    res_table  = res_volcano_nonpaired,
    out_pdf    = file.path(dir_seed, "Seed_m8_conservation_up_positive.pdf"),
    selection  = "up",
    fdr_cut    = pCutoff,
    fc_cut     = FCcutoff,
    positive_only = TRUE,
    width      = 7.0
  )

  plot_seed_m8_conservation(
    res_table  = res_volcano_nonpaired,
    out_pdf    = file.path(dir_seed, "Seed_m8_conservation_down_positive.pdf"),
    selection  = "down",
    fdr_cut    = pCutoff,
    fc_cut     = FCcutoff,
    positive_only = TRUE,
    width      = 7.0
  )
}, error = function(e) {
  message("Seed conservation plotting failed: ", conditionMessage(e))
})
# -----------------------------
# Target retrieval (multiMiR) + g:Profiler enrichment (GO+REAC)
# -----------------------------
# NOTE: Update mirna list as desired. Here we keep your original list (requires 'hsa-' prefix).
mirnas <- paste0("hsa-", c(
  "miR-4434", "miR-7112-3p", "miR-4684-5p", "miR-3064-5p",
  "miR-4478", "miR-4534", "miR-8069", "miR-365a-5p",
  "miR-3180", "miR-6887-5p", "miR-4535", "miR-4257",
  "miR-6736-5p", "miR-758-5p", "miR-4470", "miR-5682",
  "miR-4302", "miR-2467-5p", "miR-6071", "miR-3166",
  "miR-298", "miR-3681-5p", "miR-4267", "miR-135b-3p",
  "miR-1178-5p", "miR-4303", "miR-6842-3p", "miR-6838-3p",
  "miR-367-5p", "miR-2114-5p", "miR-646", "miR-933",
  "miR-623", "miR-4697-5p", "miR-208a-3p", "miR-6762-3p",
  "miR-939-5p", "miR-218-1-3p", "miR-4296", "miR-6720-5p",
  "miR-193a-3p", "miR-6750-5p", "miR-33b-5p", "miR-3714",
  "miR-3677-5p", "miR-4519", "miR-6080", "miR-4746-5p",
  "miR-497-3p", "miR-648", "miR-4266", "miR-5008-5p",
  "miR-6848-3p", "miR-6770-3p", "miR-619-3p", "miR-572",
  "miR-6715b-3p", "miR-3169", "miR-6774-3p", "miR-662",
  "miR-7703", "miR-6816-3p", "miR-127-5p", "miR-4690-3p",
  "miR-1199-3p", "miR-6893-3p", "miR-6510-3p", "miR-4673",
  "miR-6081", "miR-6821-3p", "miR-3124-5p", "miR-6791-3p",
  "miR-3917", "miR-935", "miR-4707-5p", "miR-4515"
))

message("Querying multiMiR targets (this can take time)...")
mm <- get_multimir(mirna = mirnas, summary = TRUE)

miRNA_targets <- mm@data %>%
  dplyr::select(mature_mirna_id, target_symbol, database) %>%
  distinct()

write.csv(miRNA_targets, file.path(dir_targets, "miRNA_Target_Genes.csv"), row.names = FALSE)

message("Running g:Profiler enrichment per miRNA (GO:BP/GO:MF/TF/REAC)...")
unique_mirnas <- unique(miRNA_targets$mature_mirna_id)

for (mirna in unique_mirnas) {

  target_genes_vec <- miRNA_targets %>%
    filter(mature_mirna_id == mirna) %>%
    pull(target_symbol) %>%
    unique() %>%
    na.omit()

  if (length(target_genes_vec) == 0) next

  mirna_safe <- gsub("[^A-Za-z0-9_.-]", "_", mirna)

  go_results <- tryCatch(
    gost(
      query = target_genes_vec,
      organism = "hsapiens",
      correction_method = "fdr",
      sources = c("GO:BP", "GO:MF", "TF", "REAC")
    ),
    error = function(e) {
      message("gost() failed for ", mirna, " : ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(go_results) || is.null(go_results$result) || nrow(go_results$result) == 0) next

  go_results_df <- as.data.frame(go_results$result)

  # Convert list-columns to CSV-safe strings
  is_list_col <- vapply(go_results_df, is.list, logical(1))
  if (any(is_list_col)) {
    go_results_df[is_list_col] <- lapply(go_results_df[is_list_col], function(col) {
      vapply(col, function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_character_)
        paste(unlist(x), collapse = ",")
      }, FUN.VALUE = NA_character_)
    })
  }

  # Rename intersection -> intersection_genes if present
  if ("intersection" %in% colnames(go_results_df)) {
    colnames(go_results_df)[colnames(go_results_df) == "intersection"] <- "intersection_genes"
  }

  write.csv(go_results_df,
            file = file.path(dir_go_terms, paste0("GO_Terms_", mirna_safe, ".csv")),
            row.names = FALSE)

  write.csv(data.frame(miRNA = mirna, target_genes = paste(target_genes_vec, collapse=";")),
            file = file.path(dir_go_terms, paste0("Target_Genes_", mirna_safe, ".csv")),
            row.names = FALSE)
}

message("DONE. Outputs written under: ", outdir)
