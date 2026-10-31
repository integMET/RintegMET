#' Populate the Metabolomics Workbench analysis identifier
#'
#' Adds an `analysis_id` column to `colData(se)`, for example in the result of
#' [study_variable_fold_change()], with the Metabolomics Workbench analysis that
#' the mzTab-M file reports.
#'
#' The analysis is taken from the `custom[n]` fields of the metadata (`MTD`)
#' section, from the parameter named `Metabolomics Workbench analysis ID`, e.g.
#' `[,, Metabolomics Workbench analysis ID, AN000025]`. Every column of `se` gets
#' the same value because a file reports one analysis. The value is a plain
#' identifier rather than a compact identifier, because identifiers.org has no
#' prefix for Metabolomics Workbench analyses; use
#' [populate_compact_identifier()] for the study.
#'
#' @param se A `SummarizedExperiment`.
#' @param x Path to the mzTab-M file that `se` was made from, or an `MzTabM`
#'   object already read with [RmzTabM::readMzTabM()].
#' @param validate Passed to [RmzTabM::readMzTabM()] when `x` is a path.
#'   Validation downloads and runs the jmzTab-M validator, so it is off by
#'   default.
#'
#' @return `se` with `analysis_id` in its `colData`, replacing any existing
#'   values. It is `NA` with a warning if the `custom[n]` fields do not name a
#'   single analysis.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- study_variable_summarized_experiment(file)
#' fc <- study_variable_fold_change(se)
#' fc <- populate_analysis_id(fc, file)
#' fc$analysis_id[1:3]
populate_analysis_id <- function(se, x, validate = FALSE) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  x <- as_mztab(x, validate)
  params <- RmzTabM::getMtdField(x, "^custom\\[")

  analysis <- param_analysis_id(params)
  if (is.na(analysis)) {
    warning("No single Metabolomics Workbench analysis found in the custom[n] ",
            "fields of the metadata section.", call. = FALSE)
  }

  se$analysis_id <- rep(analysis, ncol(se))
  se
}

# Analysis that the mzTab parameters name, NA unless there is exactly one.
param_analysis_id <- function(params) {
  analysis <- paste0("^\\[[^,]*,[^,]*,\\s*Metabolomics Workbench analysis ID",
                     "\\s*,\\s*(AN[0-9]{6})\\s*\\]$")
  analyses <- unique(sub(analysis, "\\1", grep(analysis, params, value = TRUE)))
  if (length(analyses) != 1L) {
    return(NA_character_)
  }
  analyses
}
