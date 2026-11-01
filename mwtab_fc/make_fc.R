# Creates one fold-change SummarizedExperiment, `fc`, per mzTab-M file of
# integMET/mwtab2mztabm, with its study, analysis and InChIKeys populated, saved
# as mwtab_fc/<file>_fc.RData. InChIKeys are looked up in PubChem, so this takes
# several minutes.
# Run from the package root: Rscript --vanilla mwtab_fc/make_fc.R

pkgload::load_all(".", quiet = TRUE)
options(warn = 1)

# Last commit on main that changed inst/extdata/mztab
ref <- "3c615f295a4ca550c3f1091f18529fb19b487d4e"
out_dir <- "mwtab_fc"

# ST000956 merges the CIDs of rows that combine several metabolites. These merged
# values happen to be valid CIDs of unrelated compounds, so their keys are wrong.
merged_uris <- sprintf("https://pubchem.ncbi.nlm.nih.gov/compound/%s",
                       c("10151015", "1375810", "74057405", "62516"))

req <- httr2::request("https://api.github.com/repos/integMET/mwtab2mztabm")
req <- httr2::req_url_path_append(req, "git/trees", ref)
req <- httr2::req_url_query(req, recursive = 1)
tree <- httr2::resp_body_json(httr2::req_perform(req))$tree
paths <- vapply(tree, function(entry) entry$path, character(1))
paths <- grep("^inst/extdata/mztab/.*\\.mztab$", paths, value = TRUE)

mztab_dir <- tempfile("mztab")
dir.create(mztab_dir)

for (path in paths) {
  file <- file.path(mztab_dir, basename(path))
  url <- sprintf("https://raw.githubusercontent.com/integMET/mwtab2mztabm/%s/%s",
                 ref, path)
  utils::download.file(url, file, quiet = TRUE, mode = "wb")

  message(basename(path))
  se <- study_variable_summarized_experiment(file)
  fc <- study_variable_fold_change(se)
  fc <- populate_compact_identifier(fc, file)
  fc <- populate_analysis_id(fc, file)
  fc <- populate_inchi_key(fc)
  if ("opt_global_cmpdid" %in% names(SummarizedExperiment::rowData(fc))) {
    fc <- populate_inchi_key_from_cmpdid(fc)
  }
  if (startsWith(basename(path), "ST000956_")) {
    row_data <- SummarizedExperiment::rowData(fc)
    row_data$opt_global_inchi_key[row_data$uri %in% merged_uris] <- NA_character_
    SummarizedExperiment::rowData(fc) <- row_data
  }
  save(fc, file = file.path(out_dir, sub("\\.mztab$", "_fc.RData", basename(path))))
}
