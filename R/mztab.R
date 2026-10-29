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
  if (!inherits(x, "MzTabM")) {
    x <- RmzTabM::readMzTabM(x, validate = validate)
  }
  groups <- RmzTabM::getMtdField(x, "study_variable_group")
  variables <- RmzTabM::getMtdField(x, "study_variable\\[")

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
    study_variable_refs = I(refs),
    study_variables = I(lapply(refs, function(r) unname(variables[r]))),
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

# Name element of CV parameters; parseCvParameter() rejects NA, so missing
# fields are kept as NA.
cv_name <- function(x) {
  out <- rep(NA_character_, length(x))
  found <- !is.na(x)
  out[found] <- RmzTabM::parseCvParameter(x[found], 3L)
  out
}
