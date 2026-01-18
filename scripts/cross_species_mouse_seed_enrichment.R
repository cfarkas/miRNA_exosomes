#!/usr/bin/env Rscript
# ===============================================================
# Cross-species validation: human EV miRNAs -> mouse targeting
#   (Seed-scan mode using TargetScan-style miR family table)
#
# FIX v2:
#   - Robust parser for miR_Family_Info.txt (plain/zip) WITHOUT assuming exact column names.
#   - Handles BOM, extra header/comment lines, and name variations ("miR family" vs "miR.family" etc.)
#
# What this script does (for ONE DE contrast/run):
#   1) Reads an edgeR DE table from an existing pipeline /results folder
#   2) Selects a miRNA subset (default = significantly UP miRNAs)
#   3) Fetches HUMAN mature sequences + Seed+m8 from a miR family table (e.g. TargetScan miR_Family_Info)
#      - very relaxed name matching (ignores hsa- prefixes, case, underscores, dot-versions, etc.)
#   4) Computes seed conservation vs MOUSE miRNAs using the SAME miR family table (species_id=10090)
#   5) Option A (recommended): scans mouse 3'UTRs for canonical seed sites (8mer, 7mer-m8, optional)
#   6) Runs functional enrichment in mouse (GO:BP + Reactome) via g:Profiler (organism = mmusculus)
#   7) Writes tables + PDF + summary
#
# Requirements (R packages):
#   - dplyr, readr, stringr, ggplot2, tibble
#   - gprofiler2
#   - AnnotationDbi, org.Mm.eg.db
#   - Biostrings, GenomicFeatures
#   - Either TxDb+BSgenome for mm10/mm39 OR --mouse_utr_fasta + --mouse_tx2gene
# ===============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(gprofiler2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(Biostrings)
  library(GenomicFeatures)
})

# ---------------------------
# Arg parsing (no optparse)
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}
has_flag <- function(flag) any(args == flag)

if (has_flag("--help") || length(args) == 0) {
  cat("
cross_species_mouse_seed_enrichment.R

Required:
  --results_dir          Pipeline run folder containing 02_edger/edgeR_results_*.csv
  --mir_family_table     miR family table (TSV or ZIP containing TSV) with columns like:
                         miR family, Seed+m8, Species ID, MiRBase ID, Mature sequence, ...

Optional:
  --design               paired | non_paired (default: non_paired)
  --de_table             Path to edgeR table (overrides --results_dir/--design)
  --direction            up | down | all_sig | all (default: up)
  --fdr                  FDR cutoff (default: 0.05)
  --logfc                logFC cutoff (default: 1)
  --outdir               output directory (default: <results_dir>/08_cross_species_mouse)
  --immune_regex         custom regex for immune term filtering (optional)

Targeting:
  --target_mode          seedscan_utr (default) | none
  --mouse_genome         mm10 | mm39 (default: mm10)
  --mouse_utr_fasta      Precomputed mouse 3'UTR FASTA (optional alternative to TxDb/BSgenome)
  --mouse_tx2gene        Transcript->Entrez mapping table (TSV/CSV) if using --mouse_utr_fasta
  --site_types           comma-separated: 8mer,7mer-m8,7mer-A1,6mer (default: 8mer,7mer-m8)
  --min_sites            minimum number of sites per transcript (default: 1)

Example:
  Rscript cross_species_mouse_seed_enrichment.R \\
    --results_dir results/without_sample2 \\
    --design non_paired \\
    --direction up --fdr 0.05 --logfc 1 \\
    --mir_family_table mirna_pipeline_package_v6/miR_Family_Info.txt.zip \\
    --target_mode seedscan_utr --mouse_genome mm10 \\
    --outdir results/without_sample2/08_cross_species_mouse_seedscan
\n")
  quit(status = 0)
}

results_dir <- get_arg("--results_dir", default = NULL)
design      <- get_arg("--design", default = "non_paired")
de_table    <- get_arg("--de_table", default = NULL)

direction   <- get_arg("--direction", default = "up")
fdr_cut     <- as.numeric(get_arg("--fdr", default = "0.05"))
logfc_cut   <- as.numeric(get_arg("--logfc", default = "1"))
outdir      <- get_arg("--outdir", default = NULL)

mir_family_table <- get_arg("--mir_family_table", default = NULL)

target_mode <- get_arg("--target_mode", default = "seedscan_utr")
mouse_genome <- get_arg("--mouse_genome", default = "mm10")
mouse_utr_fasta <- get_arg("--mouse_utr_fasta", default = NULL)
mouse_tx2gene   <- get_arg("--mouse_tx2gene", default = NULL)
site_types <- get_arg("--site_types", default = "8mer,7mer-m8")
min_sites  <- as.integer(get_arg("--min_sites", default = "1"))

species_human <- as.integer(get_arg("--species_human", default = "9606"))
species_mouse <- as.integer(get_arg("--species_mouse", default = "10090"))

immune_regex <- get_arg("--immune_regex", default =
  paste0(
    "(Th17|Th1|Th2|T[ -]?cell|B[ -]?cell|Treg|regulatory T|",
    "lymphocyte|leukocyte|immune|immun|inflamm|",
    "cytokine|chemokine|interleukin|interferon|",
    "NF[- ]?kappaB|NF[- ]?κB|NFκB|NFKB|",
    "toll[- ]?like|\\bTLR\\b|JAK[- ]?STAT|\\bTNF\\b|",
    "antigen|MHC|IL[- ]?\\d+)"
  )
)

if (is.null(outdir)) {
  if (is.null(results_dir)) stop("ERROR: Provide --results_dir or --outdir + --de_table")
  outdir <- file.path(results_dir, "08_cross_species_mouse")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(outdir, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
tab_dir <- file.path(outdir, "tables");  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Helpers
# ---------------------------
fix_empty_colnames <- function(df) {
  cn <- colnames(df)
  bad <- is.na(cn) | cn == ""
  if (any(bad)) {
    cn[bad] <- paste0("V", which(bad))
    colnames(df) <- cn
  }
  df
}

normalize_mir_key <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "")
  x <- str_replace_all(x, "_", "-")
  x <- tolower(x)
  x <- sub("\\s.*$", "", x)
  # Strip 3-letter miRBase species prefixes (e.g. "hsa-", "mmu-")
  # while avoiding the common bug of stripping the core "mir-" prefix
  # from names like "miR-365a-5p". We only strip a 3-letter prefix
  # when the string does NOT already start with "mir-" or "let-".
  x <- sub("^(?!mir-|let-)([a-z]{3})-", "", x, perl = TRUE)
  # Be forgiving about missing hyphen after 'mir'/'let' (e.g. "mir21")
  x <- sub("^(mir|let)(?=[0-9])", "\\1-", x, perl = TRUE)
  # If an upstream preprocessing step already removed the leading "mir-"
  # (e.g. "365a-5p"), re-add it so it can still match the family table.
  x <- ifelse(grepl("^[0-9]", x), paste0("mir-", x), x)
  x
}

strip_dot_version <- function(x) sub("\\.[0-9]+$", "", x)
has_arm <- function(x) grepl("-(5p|3p)$", x, ignore.case = TRUE)

seed_to_site_patterns <- function(seed_m8_rna) {
  seed <- toupper(seed_m8_rna)
  seed_dna <- chartr("U", "T", seed)
  rc7 <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seed_dna)))
  seed6 <- substr(seed_dna, 1, 6)
  rc6 <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seed6)))
  list(
    `7mer-m8` = rc7,
    `8mer`    = paste0(rc7, "A"),
    `7mer-A1` = paste0(rc6, "A"),
    `6mer`    = rc6
  )
}

flatten_list_cols <- function(df) {
  if (!is.data.frame(df)) return(df)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(x) {
        if (is.null(x) || length(x) == 0) return("")
        paste(as.character(x), collapse = ";")
      }, character(1))
    }
  }
  df
}

# ---------------------------
# Locate DE table
# ---------------------------
if (is.null(de_table)) {
  if (is.null(results_dir)) stop("ERROR: Provide --de_table or --results_dir")
  edger_dir <- file.path(results_dir, "02_edger")
  if (!dir.exists(edger_dir)) stop("ERROR: expected folder not found: ", edger_dir)
  de_table <- if (design == "paired") {
    file.path(edger_dir, "edgeR_results_paired.csv")
  } else {
    file.path(edger_dir, "edgeR_results_non_paired.csv")
  }
}
if (!file.exists(de_table)) stop("ERROR: DE table not found: ", de_table)

message("Reading DE table: ", de_table)
de <- read.csv(de_table, stringsAsFactors = FALSE, check.names = FALSE)
de <- fix_empty_colnames(de)

first_col <- colnames(de)[1]
de$miRNA_raw <- de[[first_col]]
if (all(grepl("^[0-9]+$", as.character(de$miRNA_raw))) && ("X" %in% colnames(de))) {
  de$miRNA_raw <- de$X
}
stopifnot(all(c("logFC", "FDR") %in% colnames(de)))

de <- de %>%
  mutate(
    miRNA_raw   = as.character(miRNA_raw),
    miRNA_label = {
      x <- str_trim(miRNA_raw)
      x <- str_replace_all(x, "\\s+", "")
      x
    },
    miRNA_key = normalize_mir_key(miRNA_label),
    miRNA_key_nodot = strip_dot_version(miRNA_key)
  )

# ---------------------------
# Subset selection
# ---------------------------
sel <- de
if (direction == "up") {
  sel <- de %>% filter(FDR <= fdr_cut, logFC >= logfc_cut)
} else if (direction == "down") {
  sel <- de %>% filter(FDR <= fdr_cut, logFC <= -logfc_cut)
} else if (direction == "all_sig") {
  sel <- de %>% filter(FDR <= fdr_cut, abs(logFC) >= logfc_cut)
} else if (direction == "all") {
  sel <- de
} else {
  stop("Unknown --direction: ", direction, " (use up|down|all_sig|all)")
}
message("Selected miRNAs (", direction, "): ", nrow(sel))

# ---------------------------
# Read miR family table (TSV or ZIP), ROBUSTLY
# ---------------------------
if (is.null(mir_family_table)) stop("ERROR: Provide --mir_family_table (TSV or ZIP containing TSV).")
if (!file.exists(mir_family_table)) stop("ERROR: miR family table not found: ", mir_family_table)

read_mir_family_table <- function(path) {
  real_path <- path
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    tmpdir <- tempfile("mirfam_unzip_")
    dir.create(tmpdir)
    utils::unzip(path, exdir = tmpdir)
    cand <- list.files(tmpdir, recursive = TRUE, full.names = TRUE)
    cand <- cand[grepl("\\.(txt|tsv|tab|csv)$", cand, ignore.case = TRUE)]
    if (length(cand) == 0) stop("ERROR: ZIP contains no .txt/.tsv/.tab/.csv file: ", path)
    real_path <- cand[1]
  }

  # Detect header line index (some files may have preamble/comments)
  lines <- readLines(real_path, warn = FALSE)
  hdr_idx <- which(
    grepl("species\\s*id", lines, ignore.case = TRUE) &
      grepl("mirbase", lines, ignore.case = TRUE) &
      grepl("seed", lines, ignore.case = TRUE)
  )[1]
  skip <- if (!is.na(hdr_idx) && hdr_idx > 1) hdr_idx - 1 else 0

  # First try: tab-delimited
  fam <- tryCatch({
    read.delim(real_path,
               sep = "\t", header = TRUE, skip = skip,
               stringsAsFactors = FALSE,
               quote = "", comment.char = "",
               check.names = FALSE, fill = TRUE)
  }, error = function(e) NULL)

  # Fallback: whitespace-delimited if tab parse failed
  if (is.null(fam) || ncol(fam) < 5) {
    fam <- tryCatch({
      read.table(real_path,
                 sep = "", header = TRUE, skip = skip,
                 stringsAsFactors = FALSE,
                 quote = "", comment.char = "",
                 check.names = FALSE, fill = TRUE)
    }, error = function(e) NULL)
  }

  if (is.null(fam) || nrow(fam) == 0) stop("ERROR: Failed to parse miR family table: ", real_path)

  # Clean colnames: remove BOM, normalize spaces
  cn <- colnames(fam)
  cn <- sub("^\ufeff", "", cn)           # remove BOM at start
  cn <- gsub("\ufeff", "", cn, fixed = TRUE)
  cn <- gsub("\\s+", " ", cn)
  cn <- trimws(cn)
  colnames(fam) <- cn

  cn_norm <- tolower(cn)
  cn_norm <- gsub("\\s+", " ", cn_norm)
  cn_norm <- trimws(cn_norm)

  find_col <- function(patterns) {
    idx <- integer(0)
    for (p in patterns) {
      hit <- which(grepl(p, cn_norm, ignore.case = TRUE))
      idx <- c(idx, hit)
    }
    idx <- unique(idx)
    if (length(idx) == 0) return(NA_integer_)
    idx[1]
  }

  idx_family <- find_col(c("^mir\\s*family$", "^mir_family$", "^mir\\.family$"))
  idx_seed   <- find_col(c("^seed\\+m8$", "^seed\\s*\\+\\s*m8$", "^seed\\.?\\+?m8$", "seed\\+m8"))
  idx_species<- find_col(c("^species\\s*id$", "^species_id$", "^taxon\\s*id$"))
  idx_mirbase<- find_col(c("^mirbase\\s*id$", "^mirbase_id$", "^mirbase\\.id$"))
  idx_mature <- find_col(c("^mature\\s*sequence$", "^mature_sequence$", "^mature\\.sequence$"))

  idx_fc     <- find_col(c("^family\\s*conservation\\??$", "family\\s*conservation"))
  idx_acc    <- find_col(c("^mirbase\\s*accession$", "^mirbase_accession$"))

  need_idx <- c(idx_family, idx_seed, idx_species, idx_mirbase, idx_mature)
  if (any(is.na(need_idx))) {
    msg <- paste0(
      "ERROR: Could not identify required columns in miR family table.\n",
      "Parsed columns were:\n  - ", paste(cn, collapse = "\n  - "), "\n\n",
      "Required (approx): miR family | Seed+m8 | Species ID | MiRBase ID | Mature sequence\n",
      "Tip: ensure you are using TargetScan miR_Family_Info.txt (tab-delimited) or adjust parsing.\n"
    )
    stop(msg)
  }

  out <- tibble::tibble(
    mir_family        = as.character(fam[[idx_family]]),
    seed_m8           = as.character(fam[[idx_seed]]),
    species_id        = suppressWarnings(as.integer(fam[[idx_species]])),
    mirbase_id        = as.character(fam[[idx_mirbase]]),
    mature_sequence   = as.character(fam[[idx_mature]]),
    family_conservation = if (!is.na(idx_fc)) as.character(fam[[idx_fc]]) else NA_character_,
    mirbase_accession   = if (!is.na(idx_acc)) as.character(fam[[idx_acc]]) else NA_character_
  ) %>%
    mutate(
      seed_m8 = str_trim(seed_m8),
      mirbase_id = str_trim(mirbase_id),
      mature_sequence = str_trim(mature_sequence),
      mirbase_id_key = normalize_mir_key(mirbase_id),
      mirbase_id_key_nodot = strip_dot_version(mirbase_id_key)
    )

  out
}

message("Reading miR family table: ", mir_family_table)
fam <- read_mir_family_table(mir_family_table)

fam_hsa <- fam %>% filter(species_id == species_human)
fam_mmu <- fam %>% filter(species_id == species_mouse)

message("miR family records: human=", nrow(fam_hsa), " mouse=", nrow(fam_mmu))

# Build lookup maps for human
hsa_lookup <- fam_hsa %>%
  arrange(mirbase_id) %>%
  distinct(mirbase_id_key, .keep_all = TRUE)

hsa_lookup_nodot <- fam_hsa %>%
  arrange(mirbase_id) %>%
  distinct(mirbase_id_key_nodot, .keep_all = TRUE)

# Mouse seed -> list of mouse miRNAs
mmu_seed2ids <- fam_mmu %>%
  filter(!is.na(seed_m8), seed_m8 != "") %>%
  group_by(seed_m8) %>%
  summarise(mouse_mirbase_ids = paste(unique(mirbase_id), collapse = ";"),
            n_mouse = n_distinct(mirbase_id),
            .groups = "drop") %>%
  distinct(seed_m8, .keep_all = TRUE)

mmu_seed_map <- setNames(mmu_seed2ids$mouse_mirbase_ids, mmu_seed2ids$seed_m8)
mmu_seed_n   <- setNames(mmu_seed2ids$n_mouse, mmu_seed2ids$seed_m8)

# Map selected DE miRNAs to human sequences/seeds (very relaxed)
lookup_human_mir <- function(mir_key, mir_key_nodot) {
  cand <- c(mir_key, mir_key_nodot)
  if (!has_arm(mir_key) && !has_arm(mir_key_nodot)) {
    cand <- c(cand,
              paste0(mir_key, "-5p"),
              paste0(mir_key, "-3p"),
              paste0(mir_key_nodot, "-5p"),
              paste0(mir_key_nodot, "-3p"))
  }
  cand <- unique(cand[!is.na(cand) & cand != ""])
  for (k in cand) {
    hit <- hsa_lookup %>% filter(mirbase_id_key == k)
    if (nrow(hit) == 0) hit <- hsa_lookup_nodot %>% filter(mirbase_id_key_nodot == k)
    if (nrow(hit) > 0) {
      hit <- hit[1, , drop = FALSE]
      return(list(
        matched = TRUE,
        mirbase_id = as.character(hit$mirbase_id),
        seed_m8 = as.character(hit$seed_m8),
        mature_sequence = as.character(hit$mature_sequence),
        mir_family = as.character(hit$mir_family),
        mirbase_accession = as.character(hit$mirbase_accession),
        family_conservation = as.character(hit$family_conservation)
      ))
    }
  }
  list(matched = FALSE,
       mirbase_id = NA_character_,
       seed_m8 = NA_character_,
       mature_sequence = NA_character_,
       mir_family = NA_character_,
       mirbase_accession = NA_character_,
       family_conservation = NA_character_)
}

mapped <- lapply(seq_len(nrow(sel)), function(i) {
  lookup_human_mir(sel$miRNA_key[i], sel$miRNA_key_nodot[i])
})


get_chr1 <- function(x, field) {
  v <- x[[field]]
  if (is.null(v) || length(v) == 0) return(NA_character_)
  v <- as.character(v[1])
  if (is.na(v) || v == "") return(NA_character_)
  v
}

sel$mirbase_id_matched <- vapply(mapped, get_chr1, character(1), field = "mirbase_id")
sel$seed_m8 <- vapply(mapped, get_chr1, character(1), field = "seed_m8")
sel$mature_sequence <- vapply(mapped, get_chr1, character(1), field = "mature_sequence")
sel$mir_family <- vapply(mapped, get_chr1, character(1), field = "mir_family")
sel$mirbase_accession <- vapply(mapped, get_chr1, character(1), field = "mirbase_accession")
sel$family_conservation <- vapply(mapped, get_chr1, character(1), field = "family_conservation")
sel$mirbase_match_ok <- !is.na(sel$mirbase_id_matched) & sel$mirbase_id_matched != ""

# Seed conservation (human seed -> mouse miRNAs with same seed)
sel$mouse_miRNAs_same_seed <- vapply(sel$seed_m8, function(seed) {
  if (!is.na(seed) && seed != "" && seed %in% names(mmu_seed_map)) mmu_seed_map[[seed]] else ""
}, character(1))

sel$n_mouse_miRNAs_same_seed <- vapply(sel$seed_m8, function(seed) {
  if (!is.na(seed) && seed != "" && seed %in% names(mmu_seed_n)) as.integer(mmu_seed_n[[seed]]) else 0L
}, integer(1))

sel$has_mouse_seed_match <- sel$n_mouse_miRNAs_same_seed > 0

# Export mapping diagnostics
map_out <- file.path(tab_dir, "human_miRNA_familytable_mapping.tsv")
write_tsv(
  sel %>%
    dplyr::select(miRNA_label, miRNA_key, miRNA_key_nodot, mirbase_id_matched, seed_m8, mature_sequence, mir_family, mirbase_accession, family_conservation),
  map_out
)
message("Wrote: ", map_out)

not_found_out <- file.path(tab_dir, "human_miRNAs_not_found_in_familytable.txt")
writeLines(sel$miRNA_label[!sel$mirbase_match_ok], con = not_found_out)
message("Wrote: ", not_found_out)

# ---------------------------
# Mouse targeting: seedscan in 3'UTRs
# ---------------------------
load_mouse_utrs <- function(mouse_genome, mouse_utr_fasta = NULL, mouse_tx2gene = NULL) {
  if (!is.null(mouse_utr_fasta)) {
    if (!file.exists(mouse_utr_fasta)) stop("ERROR: --mouse_utr_fasta not found: ", mouse_utr_fasta)
    message("Loading mouse 3'UTRs from FASTA: ", mouse_utr_fasta)
    utr_seqs <- Biostrings::readDNAStringSet(mouse_utr_fasta)
    tx_ids <- names(utr_seqs)

    tx2gene <- NULL
    if (!is.null(mouse_tx2gene)) {
      if (!file.exists(mouse_tx2gene)) stop("ERROR: --mouse_tx2gene not found: ", mouse_tx2gene)
      message("Loading transcript->gene mapping: ", mouse_tx2gene)
      tx2gene_df <- tryCatch({
        readr::read_tsv(mouse_tx2gene, show_col_types = FALSE, progress = FALSE)
      }, error = function(e) {
        readr::read_csv(mouse_tx2gene, show_col_types = FALSE, progress = FALSE)
      })
      cn <- colnames(tx2gene_df)
      tx_col <- cn[which(tolower(cn) %in% c("tx","tx_name","transcript","transcript_id","transcriptid","txid","name"))[1]]
      gene_col <- cn[which(tolower(cn) %in% c("gene","gene_id","geneid","entrez","entrezid","entrez_id","symbol"))[1]]
      if (is.na(tx_col) || is.na(gene_col)) {
        stop("ERROR: --mouse_tx2gene must contain transcript and gene columns. Found: ", paste(cn, collapse = ", "))
      }
      tx2gene <- tx2gene_df %>%
        transmute(tx_name = as.character(.data[[tx_col]]),
                  gene_id = as.character(.data[[gene_col]])) %>%
        filter(!is.na(tx_name), tx_name != "") %>%
        distinct(tx_name, gene_id)
    } else {
      warning("No --mouse_tx2gene provided; treating FASTA names as gene IDs (NOT recommended).")
      tx2gene <- tibble(tx_name = tx_ids, gene_id = tx_ids)
    }
    return(list(utr_seqs = utr_seqs, tx2gene = tx2gene))
  }

  txdb_pkg <- if (mouse_genome == "mm10") "TxDb.Mmusculus.UCSC.mm10.knownGene" else if (mouse_genome == "mm39") "TxDb.Mmusculus.UCSC.mm39.knownGene" else NA
  if (is.na(txdb_pkg)) stop("ERROR: Unsupported --mouse_genome: ", mouse_genome, " (use mm10 or mm39)")

  bs_pkg <- if (mouse_genome == "mm10") "BSgenome.Mmusculus.UCSC.mm10" else "BSgenome.Mmusculus.UCSC.mm39"

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    stop("ERROR: Missing package ", txdb_pkg, ". Install it OR provide --mouse_utr_fasta + --mouse_tx2gene.")
  }
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("ERROR: Missing package BSgenome. Install BSgenome + ", bs_pkg, " OR provide --mouse_utr_fasta.")
  }
  if (!requireNamespace(bs_pkg, quietly = TRUE)) {
    stop("ERROR: Missing package ", bs_pkg, ". Install it OR provide --mouse_utr_fasta.")
  }

  suppressPackageStartupMessages(library(txdb_pkg, character.only = TRUE))
  suppressPackageStartupMessages(library(bs_pkg, character.only = TRUE))

  txdb <- getExportedValue(txdb_pkg, txdb_pkg)
  genome <- get(bs_pkg)

  message("Extracting mouse 3'UTRs from ", txdb_pkg, " + ", bs_pkg, " ...")
  utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  utr_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, utr3)

  tx_gr <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  tx2gene <- tibble(tx_name = as.character(tx_gr$tx_name),
                    gene_id = as.character(tx_gr$gene_id)) %>%
    filter(!is.na(tx_name), tx_name != "", !is.na(gene_id), gene_id != "") %>%
    distinct(tx_name, gene_id)

  keep_tx <- intersect(names(utr_seqs), tx2gene$tx_name)
  utr_seqs <- utr_seqs[keep_tx]
  tx2gene <- tx2gene %>% filter(tx_name %in% keep_tx)

  message("Mouse UTRs loaded: transcripts=", length(utr_seqs), " mapped_tx2gene=", nrow(tx2gene))
  list(utr_seqs = utr_seqs, tx2gene = tx2gene)
}

seedscan_targets_for_seed <- function(seed_m8, utr_seqs, site_types_vec, min_sites) {
  if (is.na(seed_m8) || seed_m8 == "") return(character(0))
  pats <- seed_to_site_patterns(seed_m8)

  counts <- list()
  if ("7mer-m8" %in% site_types_vec || "8mer" %in% site_types_vec) {
    counts[["7mer-m8"]] <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["7mer-m8"]]), utr_seqs, fixed = TRUE)
  }
  if ("8mer" %in% site_types_vec) {
    counts[["8mer"]] <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["8mer"]]), utr_seqs, fixed = TRUE)
  }
  if ("7mer-A1" %in% site_types_vec) {
    counts[["7mer-A1"]] <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["7mer-A1"]]), utr_seqs, fixed = TRUE)
  }
  if ("6mer" %in% site_types_vec) {
    counts[["6mer"]] <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["6mer"]]), utr_seqs, fixed = TRUE)
  }

  score <- integer(length(utr_seqs))
  if ("8mer" %in% site_types_vec) score <- score + counts[["8mer"]]
  if ("7mer-m8" %in% site_types_vec) {
    if ("8mer" %in% site_types_vec) {
      score <- score + pmax(counts[["7mer-m8"]] - counts[["8mer"]], 0L)
    } else {
      score <- score + counts[["7mer-m8"]]
    }
  }
  if ("7mer-A1" %in% site_types_vec) score <- score + counts[["7mer-A1"]]
  if ("6mer" %in% site_types_vec) score <- score + counts[["6mer"]]

  names(utr_seqs)[which(score >= min_sites)]
}

site_types_vec <- str_split(site_types, ",", simplify = TRUE) %>% as.character()
site_types_vec <- str_trim(site_types_vec)
site_types_vec <- site_types_vec[site_types_vec != ""]
site_types_vec <- unique(site_types_vec)

sel$mouse_targets_entrez <- ""
sel$n_mouse_targets <- 0L
all_entrez <- character(0)

if (target_mode == "seedscan_utr") {
  usable <- sel %>% filter(!is.na(seed_m8), seed_m8 != "")
  if (nrow(usable) == 0) {
    warning("No selected miRNAs had Seed+m8 available from the miR family table. Skipping targeting/enrichment.")
  } else {
    utr_obj <- load_mouse_utrs(mouse_genome, mouse_utr_fasta, mouse_tx2gene)
    utr_seqs <- utr_obj$utr_seqs
    tx2gene <- utr_obj$tx2gene

    seeds <- unique(usable$seed_m8)
    message("Seed-scan targeting: unique seeds=", length(seeds),
            " site_types=", paste(site_types_vec, collapse = ","),
            " min_sites=", min_sites)

    seed2genes <- list()
    dbg_out <- file.path(tab_dir, "seedscan_site_counts_by_transcript_example.tsv")

    for (si in seq_along(seeds)) {
      seed <- seeds[si]
      tx_hits <- seedscan_targets_for_seed(seed, utr_seqs, site_types_vec, min_sites)
      if (length(tx_hits) == 0) {
        seed2genes[[seed]] <- character(0)
      } else {
        gene_ids <- tx2gene %>% filter(tx_name %in% tx_hits) %>% pull(gene_id) %>% unique()
        seed2genes[[seed]] <- as.character(gene_ids)
      }

      if (si == 1) {
        pats <- seed_to_site_patterns(seed)
        c7 <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["7mer-m8"]]), utr_seqs, fixed = TRUE)
        c8 <- Biostrings::vcountPattern(Biostrings::DNAString(pats[["8mer"]]), utr_seqs, fixed = TRUE)
        dbg <- tibble(
          transcript = names(utr_seqs),
          count_7mer_m8 = as.integer(c7),
          count_8mer = as.integer(c8)
        ) %>%
          arrange(desc(count_8mer), desc(count_7mer_m8)) %>%
          slice_head(n = 200)
        write_tsv(dbg, dbg_out)
        message("Wrote: ", dbg_out)
      }
    }

    sel$mouse_targets_entrez <- vapply(sel$seed_m8, function(seed) {
      if (is.na(seed) || seed == "" || !(seed %in% names(seed2genes))) return("")
      genes <- seed2genes[[seed]]
      if (length(genes) == 0) return("")
      paste(unique(genes), collapse = ";")
    }, character(1))

    sel$n_mouse_targets <- vapply(sel$mouse_targets_entrez, function(x) {
      if (is.na(x) || x == "") return(0L)
      length(unique(strsplit(x, ";", fixed = TRUE)[[1]]))
    }, integer(1))

    all_entrez <- unique(unlist(strsplit(sel$mouse_targets_entrez[sel$mouse_targets_entrez != ""], ";", fixed = TRUE)))
    all_entrez <- unique(all_entrez[!is.na(all_entrez) & all_entrez != ""])
  }

} else if (target_mode == "none") {
  message("Target mode = none. Skipping targeting/enrichment.")
} else {
  stop("ERROR: Unknown --target_mode: ", target_mode, " (use seedscan_utr|none)")
}

# Convert Entrez -> Symbol
sym_map <- tibble(ENTREZID = character(0), SYMBOL = character(0))
mouse_symbols <- character(0)

if (length(all_entrez) > 0) {
  sym_map <- AnnotationDbi::select(org.Mm.eg.db, keys = all_entrez,
                                  columns = c("SYMBOL"), keytype = "ENTREZID") %>%
    as_tibble() %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENTREZID, SYMBOL)
  mouse_symbols <- unique(sym_map$SYMBOL)
}

# ---------------------------
# Write summary tables
# ---------------------------
seed_out <- file.path(tab_dir, "seed_conservation_and_mouse_targets.csv")
sel_out <- sel %>%
  dplyr::select(
    miRNA_label, logFC, FDR,
    mirbase_id_matched, mir_family, mirbase_accession, family_conservation,
    seed_m8, mature_sequence,
    has_mouse_seed_match, n_mouse_miRNAs_same_seed, mouse_miRNAs_same_seed,
    n_mouse_targets, mouse_targets_entrez
  ) %>%
  arrange(FDR, desc(logFC))

write.csv(sel_out, seed_out, row.names = FALSE)
message("Wrote: ", seed_out)

targets_out <- file.path(tab_dir, "mouse_target_genes_used_for_enrichment.csv")
write.csv(sym_map %>% arrange(SYMBOL), targets_out, row.names = FALSE)
message("Wrote: ", targets_out)

# ---------------------------
# Enrichment
# ---------------------------
enrich_full_out <- file.path(tab_dir, "mouse_enrichment_GO_BP_REAC_full.csv")
enrich_immune_out <- file.path(tab_dir, "mouse_enrichment_GO_BP_REAC_immune_only.csv")

enrich_df <- data.frame()
enrich_immune <- data.frame()

if (target_mode != "none") {
  if (length(mouse_symbols) < 10) {
    warning("Too few mouse targets for enrichment (n=", length(mouse_symbols), "). Skipping gost().")
  } else {
    gp <- gost(
      query = mouse_symbols,
      organism = "mmusculus",
      correction_method = "fdr",
      sources = c("GO:BP", "REAC")
    )

    if (!is.null(gp) && "result" %in% names(gp) && !is.null(gp$result)) {
      enrich_df <- as.data.frame(gp$result)
      enrich_df <- flatten_list_cols(enrich_df)
      write.csv(enrich_df, enrich_full_out, row.names = FALSE)
      message("Wrote: ", enrich_full_out)

      enrich_immune <- enrich_df %>%
        filter(str_detect(term_name, regex(immune_regex, ignore_case = TRUE))) %>%
        arrange(p_value, desc(intersection_size))
      write.csv(enrich_immune, enrich_immune_out, row.names = FALSE)
      message("Wrote: ", enrich_immune_out)
    } else {
      warning("gost() returned no results.")
    }
  }
}

# ---------------------------
# Figure PDF
# ---------------------------
pdf_out <- file.path(fig_dir, "cross_species_seed_conservation_and_mouse_enrichment.pdf")

##
## Seed (m8) conservation plot
##   - User requested: show ONLY positive/conserved seeds (no empty/zero rows)
##   - Render all bars in green, no TRUE/FALSE legend
##
p1_df <- sel %>%
  filter(!is.na(n_mouse_miRNAs_same_seed)) %>%
  filter(n_mouse_miRNAs_same_seed > 0)

if (nrow(p1_df) == 0) {
  p1 <- ggplot() +
    theme_void() +
    labs(
      title = "Seed (m8) conservation: human DE miRNAs vs mouse miRNAs",
      subtitle = paste0("Selection: ", direction, " (FDR≤", fdr_cut, ", |logFC|≥", logfc_cut, ")\nNo conserved seeds found (n_mouse = 0 for all selected miRNAs).")
    )
} else {
  # Order by number of mouse miRNAs sharing the seed (then by significance)
  p1_df <- p1_df %>%
    arrange(desc(n_mouse_miRNAs_same_seed), FDR, desc(abs(logFC)))
  p1_df$miRNA_label <- factor(p1_df$miRNA_label, levels = rev(p1_df$miRNA_label))

  p1 <- ggplot(p1_df, aes(x = miRNA_label, y = n_mouse_miRNAs_same_seed)) +
    geom_col(fill = "#2E7D32", color = "black", width = 0.85, linewidth = 0.2) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0, max(p1_df$n_mouse_miRNAs_same_seed, na.rm = TRUE), by = 1)) +
    labs(
      title = "Seed (m8) conservation: human DE miRNAs vs mouse miRNAs",
      subtitle = paste0("Selection: ", direction, " (FDR≤", fdr_cut, ", |logFC|≥", logfc_cut, ")"),
      x = NULL,
      y = "# mouse miRNAs with identical seed"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7)
    )
}

p2_data <- NULL
if (is.data.frame(enrich_immune) && nrow(enrich_immune) > 0) {
  p2_data <- enrich_immune
} else if (is.data.frame(enrich_df) && nrow(enrich_df) > 0) {
  p2_data <- enrich_df
}

if (!is.null(p2_data)) {
  top_n <- min(20, nrow(p2_data))
  p2_top <- p2_data %>%
    slice_head(n = top_n) %>%
    mutate(minus_log10_fdr = -log10(p_value),
           term_label = str_wrap(term_name, width = 45),
           term_label = factor(term_label, levels = rev(unique(term_label))))

  p2 <- ggplot(p2_top, aes(x = minus_log10_fdr, y = term_label, color = source, size = intersection_size)) +
    geom_point(alpha = 0.9) +
    labs(
      title = "Mouse functional enrichment (seed-scan targets in mouse 3'UTRs)",
      subtitle = "Top terms (immune-only if available); g:Profiler organism=mmusculus; sources=GO:BP,REAC",
      x = "-log10(FDR)",
      y = NULL,
      color = "Source",
      size = "Intersection\nsize"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
} else {
  p2 <- ggplot() + theme_void() + labs(title = "No enrichment results to plot.")
}

pdf(pdf_out, width = 11, height = 8.5, onefile = TRUE)
print(p1)
print(p2)
dev.off()
message("Wrote: ", pdf_out)

# ---------------------------
# Summary TXT
# ---------------------------
summary_out <- file.path(outdir, "selection_and_cross_species_summary.txt")

total_selected <- nrow(sel)
mapped_n <- sum(sel$mirbase_match_ok, na.rm = TRUE)
seed_known <- sum(!is.na(sel$seed_m8) & sel$seed_m8 != "", na.rm = TRUE)
seed_conserved <- sum(sel$has_mouse_seed_match, na.rm = TRUE)
usable_targets <- sum(sel$n_mouse_targets > 0, na.rm = TRUE)
target_gene_n <- length(mouse_symbols)

cat(
  "Cross-species validation summary\n",
  "================================\n\n",
  "DE table: ", de_table, "\n",
  "Selection: direction=", direction, ", FDR≤", fdr_cut, ", |logFC|≥", logfc_cut, "\n\n",
  "miR family table: ", mir_family_table, "\n",
  "Human species_id=", species_human, " Mouse species_id=", species_mouse, "\n\n",
  "Targeting mode: ", target_mode, "\n",
  "Mouse genome: ", mouse_genome, "\n",
  "Site types: ", paste(site_types_vec, collapse = ","), " ; min_sites=", min_sites, "\n\n",
  "Counts:\n",
  "  Total selected miRNAs: ", total_selected, "\n",
  "  Selected miRNAs mapped to family table (human): ", mapped_n, "\n",
  "  Selected miRNAs with Seed+m8 available: ", seed_known, "\n",
  "  Selected miRNAs with ≥1 mouse miRNA sharing the same seed: ", seed_conserved, "\n",
  "  Selected miRNAs with ≥1 mouse target gene (seedscan): ", usable_targets, "\n",
  "  Union of mouse target genes (SYMBOL) used for enrichment: ", target_gene_n, "\n\n",
  "Outputs:\n",
  "  - Mapping TSV: ", map_out, "\n",
  "  - Missing miRNAs TXT: ", not_found_out, "\n",
  "  - CSV seed/targets: ", seed_out, "\n",
  "  - CSV target genes: ", targets_out, "\n",
  "  - CSV enrichment full: ", enrich_full_out, "\n",
  "  - CSV enrichment immune-only: ", enrich_immune_out, "\n",
  "  - PDF figure: ", pdf_out, "\n",
  sep = "",
  file = summary_out
)

message("Wrote: ", summary_out)
message("DONE.")
