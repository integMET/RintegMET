#' List the study variable groups of an mzTab-M file
#'
#' Reads an mzTab-M file with [RmzTabM::readMzTabM()] and returns one row per
#' `study_variable_group[n]` entry of its metadata (`MTD`) section, i.e. per
#' experimental factor, together with the study variables (factor levels) that
#' each group refers to.
#'
#' @param x Path to an mzTab-M file, or an `MzTabM` object already read with
#'   [RmzTabM::readMzTabM()].
#' @param validate Passed to [RmzTabM::readMzTabM()] when `x` is a path.
#'   Validation downloads and runs the jmzTab-M validator, so it is off by
#'   default.
#'
#' @return A data frame with one row per study variable group, ordered by
#'   group index, and the columns:
#'   \describe{
#'     \item{id}{Group identifier, e.g. `"study_variable_group[1]"`.}
#'     \item{name}{Group name, taken from the name element of the group's
#'       parameter.}
#'     \item{description}{Value of `-description`.}
#'     \item{type}{Name of the `-type` parameter, e.g.
#'       `"categorical variable"`.}
#'     \item{datatype}{Value of `-datatype`, e.g. `"xsd:string"`.}
#'     \item{n_samples}{Number of distinct samples across the group's study
#'       variables, through `-assay_refs` and `assay[n]-sample_ref`. An assay
#'       without a `-sample_ref` counts as a sample of its own.}
#'     \item{study_variable_refs}{List column of the study variable
#'       identifiers in `-study_variable_refs`.}
#'     \item{study_variables}{List column of the names of those study
#'       variables.}
#'   }
#'   Fields missing from the file are `NA`.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' groups <- list_study_variable_groups(file)
#' groups[, c("id", "name", "type")]
#' groups$study_variables
list_study_variable_groups <- function(x, validate = FALSE) {
  x <- as_mztab(x, validate)
  groups <- RmzTabM::getMtdField(x, "study_variable_group")
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")
  assay_fields <- RmzTabM::getMtdField(x, "assay\\[")

  ids <- grep("^study_variable_group\\[[0-9]+\\]$", names(groups), value = TRUE)
  ids <- ids[order(mtd_index(ids))]
  # sprintf() rather than paste0() so that no groups gives zero keys
  field <- function(suffix) unname(groups[sprintf("%s-%s", ids, suffix)])
  refs <- lapply(field("study_variable_refs"), split_refs)

  data.frame(
    id = ids,
    name = cv_name(groups[ids]),
    description = field("description"),
    type = cv_name(field("type")),
    datatype = field("datatype"),
    n_samples = vapply(refs, function(r) {
      n_samples(assay_fields, variable_assays(variables, r))
    }, integer(1)),
    study_variable_refs = I(refs),
    study_variables = I(lapply(refs, function(r) unname(variables[r]))),
    stringsAsFactors = FALSE
  )
}

#' Count the samples of each level of the study variable groups
#'
#' Returns the number of distinct samples of each study variable (level) of
#' each `study_variable_group[n]` entry of an mzTab-M file. See the `n_samples`
#' column of [list_study_variable_groups()] for how samples are counted, and
#' for the number of samples of a whole group.
#'
#' @inheritParams list_study_variable_groups
#'
#' @return A data frame with one row per level, ordered by group index and then
#'   by the order of the group's `-study_variable_refs`, and the columns `Group`
#'   (group name), `Level` (study variable name) and `Samples`. A group without
#'   study variables has no rows.
#'
#' @export
#' @examples
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' samples_per_group(file)
samples_per_group <- function(x, validate = FALSE) {
  x <- as_mztab(x, validate)
  groups <- list_study_variable_groups(x)
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")
  assay_fields <- RmzTabM::getMtdField(x, "assay\\[")
  ids <- as.character(unlist(groups$study_variable_refs))

  data.frame(
    Group = rep(groups$name, lengths(groups$study_variable_refs)),
    Level = as.character(unlist(groups$study_variables)),
    Samples = vapply(ids, function(id) {
      n_samples(assay_fields, variable_assays(variables, id))
    }, integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

# Reads `x` unless it already is an MzTabM object.
as_mztab <- function(x, validate) {
  if (inherits(x, "MzTabM")) {
    return(x)
  }
  RmzTabM::readMzTabM(x, validate = validate)
}

# SML section of an MzTabM object as a character matrix.
mztab_sml <- function(x) {
  sml <- RmzTabM::sml(x)
  if (!"SML_ID" %in% colnames(sml)) {
    stop("No small molecule summary (SML) section found.", call. = FALSE)
  }
  sml
}

# Columns that identify the small molecule of each SML row.
sml_annotation <- function(sml) {
  data.frame(
    SML_ID = trimws(sml[, "SML_ID"]),
    database_identifier = null_to_na(sml[, "database_identifier"]),
    chemical_name = null_to_na(sml[, "chemical_name"]),
    uri = null_to_na(sml[, "uri"]),
    # A one-row matrix keeps the column name when indexed; don't use it
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# "study_variable_group[10]" -> 10
mtd_index <- function(ids) {
  as.integer(sub("^.*\\[([0-9]+)\\]$", "\\1", ids))
}

# "study_variable[1]|study_variable[2]" -> c("study_variable[1]", ...)
split_refs <- function(x) {
  if (is.na(x)) {
    return(character(0))
  }
  trimws(strsplit(x, "|", fixed = TRUE)[[1]])
}

# Distinct assays referred to by the given study variables.
variable_assays <- function(variables, ids) {
  unique(unlist(lapply(ids, function(id) {
    split_refs(variables[paste0(id, "-assay_refs")])
  })))
}

# Number of distinct samples behind the given assays. An assay without a
# -sample_ref counts as a sample of its own.
n_samples <- function(assay_fields, assays) {
  # sprintf() rather than paste0() so that no assays gives zero keys
  refs <- unname(assay_fields[sprintf("%s-sample_ref", assays)])
  sum(is.na(refs)) + length(unique(refs[!is.na(refs)]))
}

# Name element of CV parameters; parseCvParameter() rejects NA, so missing
# fields are kept as NA.
cv_name <- function(x) {
  out <- rep(NA_character_, length(x))
  found <- !is.na(x)
  out[found] <- RmzTabM::parseCvParameter(x[found], 3L)
  out
}
