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

# "study_variable[10]" -> 10
mtd_index <- function(ids) {
  as.integer(sub("^.*\\[([0-9]+)\\]$", "\\1", ids))
}

# mzTab writes missing values as "null".
null_to_na <- function(x) {
  x <- trimws(x)
  x[x == "null"] <- NA_character_
  x
}
