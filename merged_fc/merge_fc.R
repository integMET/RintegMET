# Joins mztabm_fc/merged_fold_change.RData and mwtab_fc/merged_fold_change.RData on
# their InChIKey row names, saved as `merged_fc` in
# merged_fc/merged_fold_change.RData. Every key of either object is kept, and a key
# that one object lacks is NA in its fold_change columns and 0 in its n_values
# columns. The values themselves are not changed.
# Run from the package root: Rscript --vanilla merged_fc/merge_fc.R

# rownames() of a SummarizedExperiment needs its methods
invisible(loadNamespace("SummarizedExperiment"))

read_merged <- function(file) {
  env <- new.env()
  load(file, envir = env)
  env$merged_fc
}

unique_text <- function(values) {
  values <- unique(values[!is.na(values)])
  if (length(values)) paste(values, collapse = " | ") else NA_character_
}

parts <- list(read_merged("mztabm_fc/merged_fold_change.RData"),
              read_merged("mwtab_fc/merged_fold_change.RData"))
keys <- sort(unique(unlist(lapply(parts, rownames))), method = "radix")

# Places one assay of both objects side by side on the shared rows.
join_assay <- function(name, type) {
  joined <- lapply(parts, function(part) {
    values <- SummarizedExperiment::assay(part, name)
    filled <- matrix(type, nrow = length(keys), ncol = ncol(values),
                     dimnames = list(keys, colnames(values)))
    filled[rownames(values), ] <- values
    filled
  })
  do.call(cbind, joined)
}

# Stacks data frames, with NA for the columns that one of them lacks.
stack_frames <- function(frames) {
  columns <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(frame) {
    frame[setdiff(columns, names(frame))] <- NA_character_
    frame[columns]
  }))
}

source_rows <- stack_frames(lapply(parts, function(part) {
  S4Vectors::metadata(part)$source_rows
}))
rownames(source_rows) <- NULL
col_data <- stack_frames(lapply(parts, function(part) {
  as.data.frame(SummarizedExperiment::colData(part), stringsAsFactors = FALSE)
}))

rows_of_key <- split(seq_len(nrow(source_rows)),
                     source_rows$opt_global_inchi_key)[keys]
merged_fc <- SummarizedExperiment::SummarizedExperiment(
  assays = list(fold_change = join_assay("fold_change", NA_real_),
                n_values = join_assay("n_values", 0L)),
  rowData = data.frame(
    opt_global_inchi_key = keys,
    name = vapply(rows_of_key, function(rows) {
      unique_text(source_rows$name[rows])
    }, character(1)),
    uri = vapply(rows_of_key, function(rows) {
      unique_text(source_rows$uri[rows])
    }, character(1)),
    database_identifier = vapply(rows_of_key, function(rows) {
      unique_text(source_rows$database_identifier[rows])
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

save(merged_fc, file = "merged_fc/merged_fold_change.RData")
message(sprintf("%d keys x %d pairs", nrow(merged_fc), ncol(merged_fc)))
