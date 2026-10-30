#' Summarize the SML section by study variable as a SummarizedExperiment
#'
#' Builds a `SummarizedExperiment` with one row per row of the small molecule
#' summary (`SML`) section of an mzTab-M file and one column per study variable
#' (level) of its study variable groups.
#'
#' Each assay value is the arithmetic mean of the `abundance_assay[n]` columns
#' of the assays in the study variable's `-assay_refs`, as in
#' [sml_fold_change()].
#'
#' A sample usually belongs to one level of every study variable group, so the
#' columns of different groups share samples and are not independent.
#'
#' @inheritParams sml_fold_change
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with:
#'   \describe{
#'     \item{assay}{One assay, `mean_abundance`.}
#'     \item{colData}{The data frame returned by [samples_per_group()].}
#'     \item{rowData}{The `SML` columns `SML_ID`, `database_identifier`,
#'       `chemical_name` and `uri`, with `null` values as `NA`.}
#'   }
#'   Row names are the `SML_ID` values and column names are the study variable
#'   identifiers, e.g. `"study_variable[10]"`, since level names need not be
#'   unique across groups.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- sml_summarized_experiment(file)
#' se
#' SummarizedExperiment::colData(se)
#' SummarizedExperiment::assay(se)[1:3, se$Group == "TGF"]
sml_summarized_experiment <- function(x, na.rm = TRUE, validate = FALSE) {
  x <- as_mztab(x, validate)
  sml <- mztab_sml(x)
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")
  groups <- list_study_variable_groups(x)
  # Same order as the rows of samples_per_group()
  ids <- as.character(unlist(groups$study_variable_refs))
  annotation <- sml_annotation(sml)

  means <- vapply(ids, function(id) {
    study_variable_mean(sml, variables, id, na.rm)
  }, numeric(nrow(sml)))

  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      mean_abundance = matrix(means, nrow = nrow(sml), ncol = length(ids),
                              dimnames = list(annotation$SML_ID, ids))
    ),
    rowData = annotation,
    colData = samples_per_group(x)
  )
}

#' Convert study variable abundances to a SummarizedExperiment
#'
#' Builds a `SummarizedExperiment` directly from the values reported in an
#' mzTab-M file, with one row per row of the small molecule summary (`SML`)
#' section and one column per `study_variable[n]` of the metadata (`MTD`)
#' section.
#'
#' Unlike [sml_summarized_experiment()], which averages the assays itself, this
#' function does not compute anything, so the values and which of them are
#' missing depend on how the file was written. In the example file they are the
#' means of the available assays; other writers may leave
#' `abundance_study_variable[n]` as `null` as soon as any one assay of the study
#' variable is missing.
#'
#' @inheritParams list_study_variable_groups
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with:
#'   \describe{
#'     \item{assay}{One numeric assay, `abundance`, with the
#'       `abundance_study_variable[n]` columns of `SML`.}
#'     \item{rowData}{All other `SML` columns except the `SMH` line prefix and
#'       the remaining `abundance_` columns (`abundance_assay[n]` and
#'       `abundance_variation_study_variable[n]`). Columns are kept as text.}
#'     \item{colData}{`name`, the value of `study_variable[n]`, and one column
#'       per `study_variable[n]-` field, e.g. `assay_refs` and
#'       `average_function`, with their values as written in the file. Fields
#'       that a study variable lacks are `NA`.}
#'   }
#'   `null` values are `NA`. Row names are the `SML_ID` values and column names
#'   are the study variable identifiers, e.g. `"study_variable[10]"`, ordered by
#'   index.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- study_variable_summarized_experiment(file)
#' se
#' SummarizedExperiment::colData(se)[, c("name", "description")]
#' SummarizedExperiment::rowData(se)[1:3, c("SML_ID", "chemical_name")]
study_variable_summarized_experiment <- function(x, validate = FALSE) {
  x <- as_mztab(x, validate)
  sml <- mztab_sml(x)
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")

  ids <- grep("^study_variable\\[[0-9]+\\]$", names(variables), value = TRUE)
  ids <- ids[order(mtd_index(ids))]
  columns <- sprintf("abundance_%s", ids)
  missing <- setdiff(columns, colnames(sml))
  if (length(missing)) {
    stop("SML section lacks the column(s) ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  features <- sml[, !grepl("^(SMH$|abundance_)", colnames(sml)), drop = FALSE]
  features[] <- null_to_na(features)
  abundance <- null_to_na(sml[, columns, drop = FALSE])

  col_data <- data.frame(name = unname(variables[ids]), row.names = ids,
                         stringsAsFactors = FALSE)
  fields <- setdiff(sub("^study_variable\\[[0-9]+\\]-?", "", names(variables)),
                    "")
  # sprintf() rather than paste0() so that no study variables gives zero keys
  col_data[fields] <- lapply(fields, function(field) {
    unname(variables[sprintf("%s-%s", ids, field)])
  })

  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      abundance = matrix(as.numeric(abundance), nrow = nrow(sml),
                         ncol = length(ids),
                         dimnames = list(features[, "SML_ID"], ids))
    ),
    rowData = as.data.frame(features, stringsAsFactors = FALSE),
    colData = col_data
  )
}
