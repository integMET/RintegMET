#' Fold-changes between study variables of the same variable
#'
#' Builds a `SummarizedExperiment` of fold-changes from the result of
#' [study_variable_summarized_experiment()], with one column per pair of study
#' variables that belong to the same variable.
#'
#' Study variables are grouped by `colData(se)$description`. A description of
#' the form `"Variable TGF, value NO_TGF"` is unique to its study variable, so
#' only its variable part, `"TGF"`, is compared. Any other description is
#' compared as a whole, and a study variable with an `NA` description is left
#' out.
#'
#' Every pair is listed once, and its fold-change is the study variable that
#' comes later in `se` over the one that comes earlier. Take the reciprocal of
#' the assay values for the opposite direction.
#'
#' @param se A `SummarizedExperiment` from
#'   [study_variable_summarized_experiment()].
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with:
#'   \describe{
#'     \item{assay}{One numeric assay, `fold_change`. A missing abundance gives
#'       `NA`, a zero denominator gives `Inf`, and zero over zero gives `NaN`.}
#'     \item{rowData}{The `rowData` of `se`.}
#'     \item{colData}{`variable` (the group compared), `numerator` and
#'       `denominator` (study variable names), and `numerator_id` and
#'       `denominator_id` (study variable identifiers).}
#'   }
#'   Row names are those of `se`, and column names are of the form
#'   `"study_variable[2]/study_variable[1]"`.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- study_variable_summarized_experiment(file)
#' fc <- study_variable_fold_change(se)
#' fc
#' SummarizedExperiment::colData(fc)
#' SummarizedExperiment::assay(fc)[1:3, fc$variable == "TGF", drop = FALSE]
study_variable_fold_change <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  col_data <- SummarizedExperiment::colData(se)
  missing <- setdiff(c("name", "description"), names(col_data))
  if (length(missing)) {
    stop("colData(se) lacks the column(s) ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  groups <- sub("^Variable (.*), value .*$", "\\1", col_data$description)
  members <- split(seq_along(groups), factor(groups, levels = unique(groups)))
  members <- members[lengths(members) > 1L]
  # Two rows, the earlier study variable of each pair and then the later one
  pairs <- lapply(members, utils::combn, 2L)
  pairs <- matrix(as.integer(unlist(pairs)), nrow = 2L)
  denominator <- pairs[1L, ]
  numerator <- pairs[2L, ]

  ids <- colnames(se)
  abundance <- SummarizedExperiment::assay(se, "abundance")
  fold_change <- abundance[, numerator, drop = FALSE] /
    abundance[, denominator, drop = FALSE]
  columns <- sprintf("%s/%s", ids[numerator], ids[denominator])
  dimnames(fold_change) <- list(rownames(se), columns)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(fold_change = fold_change),
    rowData = SummarizedExperiment::rowData(se),
    colData = data.frame(
      variable = as.character(rep(names(members),
                                  choose(lengths(members), 2L))),
      numerator = col_data$name[numerator],
      denominator = col_data$name[denominator],
      numerator_id = ids[numerator],
      denominator_id = ids[denominator],
      row.names = columns,
      stringsAsFactors = FALSE
    )
  )
}
