# Merges the fold-change objects of mwtab_fc on opt_global_inchi_key, with one row
# per key and one column per pair of study variables of each object, saved as
# `merged_fc` in mwtab_fc/merged_fold_change.RData. Rows with an NA key are dropped,
# as are rows with another kind of value in the column, such as the KEGG IDs of
# ST000044 and ST000061, and objects left without rows contribute no columns.
#
# Rows of one object that share a key are collapsed to the geometric median of
# their values, which keeps fold-changes reciprocal when a pair is reversed. A key
# with one row in an object keeps its value as is. The n_values assay counts the
# values behind each cell, so `n_values > 1` marks the collapsed ones, and
# metadata(merged_fc)$source_rows lists the rows of every key.
# Run from the package root: Rscript --vanilla mwtab_fc/merge_fc.R

out_dir <- "mwtab_fc"

files <- list.files(out_dir, "_fc\\.RData$", full.names = TRUE)
sources <- sub("_fc\\.RData$", "", basename(files))

read_fc <- function(file) {
  env <- new.env()
  load(file, envir = env)
  env$fc
}

# Geometric median of the values of one cell. NaN and NA carry no value.
collapse_values <- function(values) {
  if (length(values) == 1L) {
    return(values)
  }
  values <- values[!is.na(values)]
  if (!length(values)) {
    return(NA_real_)
  }
  exp(stats::median(log(values)))
}

unique_text <- function(values) {
  values <- unique(values[!is.na(values)])
  if (length(values)) paste(values, collapse = " | ") else NA_character_
}

source_rows <- list()
col_data <- list()
fold_change <- list()
n_values <- list()

for (i in seq_along(files)) {
  fc <- read_fc(files[i])
  row_data <- SummarizedExperiment::rowData(fc)
  keep <- grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", row_data$opt_global_inchi_key)
  if (!any(keep)) {
    next
  }
  source <- sources[i]
  columns <- sprintf("%s:%s", source, colnames(fc))

  names_of_rows <- row_data$chemical_name
  no_name <- is.na(names_of_rows)
  names_of_rows[no_name] <- row_data$opt_global_metabolite_name[no_name]
  source_rows[[source]] <- data.frame(
    opt_global_inchi_key = row_data$opt_global_inchi_key[keep],
    source = source,
    SML_ID = row_data$SML_ID[keep],
    name = names_of_rows[keep],
    uri = row_data$uri[keep],
    stringsAsFactors = FALSE
  )

  values <- SummarizedExperiment::assay(fc, "fold_change")[keep, , drop = FALSE]
  groups <- split(seq_len(nrow(values)), row_data$opt_global_inchi_key[keep])
  fold_change[[source]] <- do.call(rbind, lapply(groups, function(rows) {
    apply(values[rows, , drop = FALSE], 2L, collapse_values)
  }))
  n_values[[source]] <- do.call(rbind, lapply(groups, function(rows) {
    colSums(!is.na(values[rows, , drop = FALSE]))
  }))
  dimnames(fold_change[[source]]) <- list(names(groups), columns)
  dimnames(n_values[[source]]) <- list(names(groups), columns)

  column_data <- as.data.frame(SummarizedExperiment::colData(fc),
                               stringsAsFactors = FALSE)
  rownames(column_data) <- columns
  col_data[[source]] <- cbind(source = source, column_data,
                              stringsAsFactors = FALSE)
}

source_rows <- do.call(rbind, unname(source_rows))
rownames(source_rows) <- NULL
col_data <- do.call(rbind, unname(col_data))
keys <- sort(unique(source_rows$opt_global_inchi_key), method = "radix")

# Places the matrices of all objects on the shared rows and columns.
assemble <- function(parts, type) {
  merged <- matrix(type, nrow = length(keys), ncol = nrow(col_data),
                   dimnames = list(keys, rownames(col_data)))
  for (part in parts) {
    merged[rownames(part), colnames(part)] <- part
  }
  merged
}

fold_change <- assemble(fold_change, NA_real_)
n_values <- assemble(n_values, 0L)
storage.mode(n_values) <- "integer"

rows_of_key <- split(seq_len(nrow(source_rows)),
                     source_rows$opt_global_inchi_key)[keys]
merged_fc <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fold_change = fold_change, n_values = n_values),
  rowData = data.frame(
    opt_global_inchi_key = keys,
    name = vapply(rows_of_key, function(rows) {
      unique_text(source_rows$name[rows])
    }, character(1)),
    uri = vapply(rows_of_key, function(rows) {
      unique_text(source_rows$uri[rows])
    }, character(1)),
    n_rows = lengths(rows_of_key),
    n_sources = vapply(rows_of_key, function(rows) {
      length(unique(source_rows$source[rows]))
    }, integer(1)),
    row.names = keys,
    stringsAsFactors = FALSE
  ),
  colData = col_data,
  metadata = list(source_rows = source_rows)
)

save(merged_fc, file = file.path(out_dir, "merged_fold_change.RData"))
message(sprintf("%d keys x %d pairs from %d of %d objects", nrow(merged_fc),
                ncol(merged_fc), length(unique(col_data$source)), length(files)))
