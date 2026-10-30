#' Fold-change between two study variables for each small molecule
#'
#' Computes, for every row of the small molecule summary (`SML`) section of an
#' mzTab-M file, the fold-change of the mean abundance in an experimental
#' group over the mean abundance in a control group.
#'
#' Groups are given as study variables, i.e. the levels of a study variable
#' group such as `"NO_TGF"` and `"TGF"` of the group `TGF`; see
#' [list_study_variable_groups()] for the study variables in a file. The mean
#' of a group is the arithmetic mean of the `abundance_assay[n]` columns of the
#' assays in its `-assay_refs`. The `abundance_study_variable[n]` columns are
#' not used, since how they treat missing assays depends on the file's writer.
#'
#' @inheritParams list_study_variable_groups
#' @param control,experimental Name of the control and of the experimental
#'   study variable, e.g. `"NO_TGF"`, or its identifier, e.g.
#'   `"study_variable[10]"` (needed when a name is used by more than one
#'   study variable).
#' @param na.rm Whether to ignore missing assay abundances when averaging.
#'   With `FALSE`, a group mean is `NA` if any of its assays is missing.
#'
#' @return A data frame with one row per `SML` row, in the same order, and
#'   the columns `SML_ID`, `database_identifier`, `chemical_name`, `uri`,
#'   `control_mean`, `experimental_mean`, `fold_change`
#'   (`experimental_mean / control_mean`) and `log2_fold_change`. `null`
#'   values are `NA`; a control mean of zero gives an infinite fold-change.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' fc <- sml_fold_change(file, control = "NO_TGF", experimental = "TGF")
#' head(fc)
sml_fold_change <- function(x, control, experimental, na.rm = TRUE,
                            validate = FALSE) {
  x <- as_mztab(x, validate)
  sml <- mztab_sml(x)
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")

  control_mean <- study_variable_mean(sml, variables, control, na.rm)
  experimental_mean <- study_variable_mean(sml, variables, experimental, na.rm)
  fold_change <- experimental_mean / control_mean

  data.frame(
    sml_annotation(sml),
    control_mean = control_mean,
    experimental_mean = experimental_mean,
    fold_change = fold_change,
    log2_fold_change = log2(fold_change),
    stringsAsFactors = FALSE
  )
}

# Mean abundance over the assays of one study variable, for each SML row.
study_variable_mean <- function(sml, variables, name, na.rm) {
  id <- study_variable_id(variables, name)
  assays <- split_refs(variables[paste0(id, "-assay_refs")])
  if (length(assays) == 0L) {
    stop(id, " has no assay_refs.", call. = FALSE)
  }
  columns <- paste0("abundance_", assays)
  missing <- setdiff(columns, colnames(sml))
  if (length(missing)) {
    stop("SML section lacks the column(s) ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  abundance <- sml[, columns, drop = FALSE]
  abundance <- matrix(as.numeric(null_to_na(abundance)), nrow = nrow(sml))
  means <- rowMeans(abundance, na.rm = na.rm)
  # rowMeans() gives NaN when every assay of a row is missing
  means[is.nan(means)] <- NA_real_
  means
}

# Resolves a study variable name, or a "study_variable[n]" id, to its id.
study_variable_id <- function(variables, name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("Study variables must be given as a single string.", call. = FALSE)
  }
  ids <- grep("^study_variable\\[[0-9]+\\]$", names(variables), value = TRUE)
  if (name %in% ids) {
    return(name)
  }

  id <- ids[which(variables[ids] == name)]
  if (length(id) == 0L) {
    stop("No study variable named \"", name, "\". Available: ",
         paste(variables[ids], collapse = ", "), call. = FALSE)
  }
  if (length(id) > 1L) {
    stop("Study variable name \"", name, "\" is used by ",
         paste(id, collapse = ", "), "; pass one of these ids instead.",
         call. = FALSE)
  }
  id
}

# mzTab writes missing values as "null".
null_to_na <- function(x) {
  x <- trimws(x)
  x[x == "null"] <- NA_character_
  x
}
