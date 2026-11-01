# Creates one fold-change SummarizedExperiment, `fc`, per mzTab-M file of
# integMET/mtbls2mztabm, saved as mztabm_fc/<file>_fc.RData.
# Run from the package root: Rscript --vanilla mztabm_fc/make_fc.R

pkgload::load_all(".", quiet = TRUE)

# Last commit on devel that changed inst/extdata/mztabm
ref <- "4c5bd3a971719bccd798d1ee0b638083acab9cc0"
out_dir <- "mztabm_fc"

req <- httr2::request("https://api.github.com/repos/integMET/mtbls2mztabm")
req <- httr2::req_url_path_append(req, "contents/inst/extdata/mztabm")
req <- httr2::req_url_query(req, ref = ref)
listing <- httr2::resp_body_json(httr2::req_perform(req))

mztab_dir <- tempfile("mztabm")
dir.create(mztab_dir)

for (entry in listing) {
  file <- file.path(mztab_dir, entry$name)
  utils::download.file(entry$download_url, file, quiet = TRUE, mode = "wb")

  se <- study_variable_summarized_experiment(file)
  fc <- study_variable_fold_change(se)
  fc <- populate_compact_identifier(fc, file)
  # The chemical names of MTBLS135 do not match its database identifiers, so a key
  # could not be trusted to describe the row it is on.
  if (!startsWith(entry$name, "MTBLS135_")) {
    fc <- populate_inchi_key_from_database_identifier(fc)
  }
  save(fc, file = file.path(out_dir, sub("\\.mztab$", "_fc.RData", entry$name)))

  values <- SummarizedExperiment::assay(fc)
  message(sprintf("%s: %d x %d, %d finite of %d", entry$name, nrow(fc), ncol(fc),
                  sum(is.finite(values)), length(values)))
}
