#' Populate the identifiers.org compact identifier of the study
#'
#' Adds a `compact_identifier` column to `colData(se)`, for example in the result
#' of [study_variable_fold_change()], with the
#' [identifiers.org](https://identifiers.org/) compact identifier of the study
#' that the mzTab-M file reports.
#'
#' The study is taken from the `uri[n]` fields of the metadata (`MTD`) section. A
#' Metabolomics Workbench study page, e.g.
#' `https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000010`,
#' gives `mw.study:ST000010`, which resolves at
#' `https://identifiers.org/mw.study:ST000010`. A MetaboLights study page, e.g.
#' `https://www.ebi.ac.uk/metabolights/MTBLS1066`, gives
#' `metabolights:MTBLS1066`, which resolves at
#' `https://identifiers.org/metabolights:MTBLS1066`. Every column of `se` gets the
#' same value because a file reports one study; the identifier refers to the whole
#' study rather than to the analysis in the file, which
#' [populate_analysis_id()] adds for Metabolomics Workbench.
#'
#' @param se A `SummarizedExperiment`.
#' @param x Path to the mzTab-M file that `se` was made from, or an `MzTabM`
#'   object already read with [RmzTabM::readMzTabM()].
#' @param validate Passed to [RmzTabM::readMzTabM()] when `x` is a path.
#'   Validation downloads and runs the jmzTab-M validator, so it is off by
#'   default.
#'
#' @return `se` with `compact_identifier` in its `colData`, replacing any existing
#'   values. It is `NA` with a warning if the `uri[n]` fields do not name a single
#'   study.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- study_variable_summarized_experiment(file)
#' fc <- study_variable_fold_change(se)
#' fc <- populate_compact_identifier(fc, file)
#' fc$compact_identifier[1:3]
populate_compact_identifier <- function(se, x, validate = FALSE) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  x <- as_mztab(x, validate)
  uris <- RmzTabM::getMtdField(x, "^uri\\[")

  identifier <- uri_compact_identifier(uris)
  if (is.na(identifier)) {
    warning("No single Metabolomics Workbench or MetaboLights study found in the ",
            "uri[n] fields of the metadata section.", call. = FALSE)
  }

  se$compact_identifier <- rep(identifier, ncol(se))
  se
}

# Compact identifier of the one study that the URIs name, otherwise NA.
uri_compact_identifier <- function(uris) {
  workbench <- paste0("^https?://(www\\.)?metabolomicsworkbench\\.org/",
                      ".*[?&]StudyID=(ST[0-9]{6})(&.*)?$")
  metabolights <- paste0("^https?://(www\\.)?ebi\\.ac\\.uk/metabolights/",
                         "(MTBLS[0-9]+)(/.*)?$")
  studies <- unique(c(
    sub(workbench, "mw.study:\\2", grep(workbench, uris, value = TRUE)),
    sub(metabolights, "metabolights:\\2", grep(metabolights, uris, value = TRUE))
  ))
  if (length(studies) != 1L) {
    return(NA_character_)
  }
  studies
}
