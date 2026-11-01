#' Populate InChIKeys from the URIs of small molecules
#'
#' Fills the `opt_global_inchi_key` column of `rowData(se)` by looking up each
#' molecule's `uri` in [PubChem](https://pubchem.ncbi.nlm.nih.gov/), for example
#' in the result of [study_variable_summarized_experiment()].
#'
#' Two forms of URI are recognized. A PubChem compound page, e.g.
#' `https://pubchem.ncbi.nlm.nih.gov/compound/74300`, gives the compound directly.
#' A KEGG compound page, e.g. `https://www.kegg.jp/entry/C02774`, gives the
#' compound that PubChem standardized from the substance deposited by KEGG.
#' The key describes the compound the URI refers to, so rows that share a URI
#' get the same key, and a KEGG entry for a class of lipids may give the key of
#' one of its members or none.
#'
#' This function needs internet access. It sends the identifiers in the URIs to
#' PubChem, at no more than five requests per second as PubChem asks. Network and
#' server errors stop with an error rather than leaving values missing.
#'
#' @param se A `SummarizedExperiment` whose `rowData` has a `uri` column.
#'
#' @return `se` with `opt_global_inchi_key` in its `rowData`. Existing values are
#'   kept, and the column is added if it is absent. Rows without a `uri` stay
#'   `NA`, and URIs that give no single InChIKey stay `NA` with a warning.
#'
#' @export
#' @examples
#' \donttest{
#' file <- system.file("extdata", "ST000010_AN000025.mztab",
#'                     package = "RintegMET")
#' se <- study_variable_summarized_experiment(file)
#' se <- populate_inchi_key(se[1:5, ])
#' SummarizedExperiment::rowData(se)[, c("uri", "opt_global_inchi_key")]
#' }
populate_inchi_key <- function(se) {
  fill_inchi_keys(se, "uri", uri_cid, "URIs")
}

#' Populate InChIKeys from the compound identifiers of small molecules
#'
#' Fills the `opt_global_inchi_key` column of `rowData(se)` by looking up each
#' molecule's `opt_global_cmpdid` in
#' [PubChem](https://pubchem.ncbi.nlm.nih.gov/), for example in the result of
#' [study_variable_fold_change()]. Use it for files that report their molecules
#' in this column rather than in `uri`, which [populate_inchi_key()] reads.
#'
#' Two forms of identifier are recognized. A KEGG compound, e.g. `C00041`, and an
#' HMDB metabolite in its five or seven digit form, e.g. `HMDB00688` or
#' `HMDB0000688`, each give the compound that PubChem standardized from the
#' substance deposited by that database. Any other identifier, such as an in-house
#' one, is left `NA` without a request. Rows that share an identifier get the same
#' key.
#'
#' This function needs internet access. It sends the identifiers to PubChem, at no
#' more than five requests per second as PubChem asks. Network and server errors
#' stop with an error rather than leaving values missing.
#'
#' @param se A `SummarizedExperiment` whose `rowData` has an `opt_global_cmpdid`
#'   column.
#'
#' @return `se` with `opt_global_inchi_key` in its `rowData`. Existing values are
#'   kept, and the column is added if it is absent. Rows without an identifier stay
#'   `NA`, and identifiers that give no single InChIKey stay `NA` with a warning.
#'
#' @export
#' @examples
#' \donttest{
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   rowData = data.frame(opt_global_cmpdid = c("C00041", "HMDB00688", "BA109"))
#' )
#' se <- suppressWarnings(populate_inchi_key_from_cmpdid(se))
#' SummarizedExperiment::rowData(se)
#' }
populate_inchi_key_from_cmpdid <- function(se) {
  fill_inchi_keys(se, "opt_global_cmpdid", cmpdid_cid, "identifiers")
}

#' Populate InChIKeys from the database identifiers of small molecules
#'
#' Fills the `opt_global_inchi_key` column of `rowData(se)` by looking up each
#' molecule's `database_identifier` in
#' [PubChem](https://pubchem.ncbi.nlm.nih.gov/), for example in the result of
#' [study_variable_fold_change()] for a file converted from MetaboLights. Use it
#' for files that report their molecules in this column rather than in `uri`, which
#' [populate_inchi_key()] reads.
#'
#' These forms of identifier are recognized:
#'
#' * A ChEBI entity, e.g. `CHEBI:28300`, a KEGG compound, e.g. `C00041` or
#'   `KEGG:C00041`, and an HMDB metabolite, e.g. `HMDB00688` or `HMDB:HMDB00688`,
#'   each give the compound that PubChem standardized from the substance deposited
#'   by that database.
#' * A PubChem compound, e.g. `PubChem:5950`. A bare number, e.g. `5950`, is taken
#'   to be one too.
#' * A CAS Registry Number, e.g. `56-41-7`, gives the compound that PubChem finds
#'   for that name, if it finds only one.
#'
#' Any other value is left `NA` without a request. That includes several
#' identifiers in one value, e.g. `CHEBI:16108|CHEBI:57597`, because they may name
#' different compounds. Rows that share an identifier get the same key.
#'
#' This function needs internet access. It sends the identifiers to PubChem, at no
#' more than five requests per second as PubChem asks. Network and server errors
#' stop with an error rather than leaving values missing.
#'
#' @param se A `SummarizedExperiment` whose `rowData` has a `database_identifier`
#'   column.
#'
#' @return `se` with `opt_global_inchi_key` in its `rowData`. Existing values are
#'   kept, and the column is added if it is absent. Rows without an identifier stay
#'   `NA`, and identifiers that give no single InChIKey stay `NA` with a warning.
#'
#' @export
#' @examples
#' \donttest{
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   rowData = data.frame(database_identifier = c("CHEBI:28300", "HMDB00688"))
#' )
#' se <- populate_inchi_key_from_database_identifier(se)
#' SummarizedExperiment::rowData(se)
#' }
populate_inchi_key_from_database_identifier <- function(se) {
  fill_inchi_keys(se, "database_identifier", database_identifier_cid,
                  "identifiers")
}

# Fills the missing InChIKeys of `se` from the identifiers in `column`, which
# `to_cid` turns into PubChem compound identifiers. `label` names them in warnings.
fill_inchi_keys <- function(se, column, to_cid, label) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  row_data <- SummarizedExperiment::rowData(se)
  if (!column %in% names(row_data)) {
    stop("rowData(se) lacks the column ", column, call. = FALSE)
  }

  identifiers <- row_data[[column]]
  keys <- row_data$opt_global_inchi_key
  if (is.null(keys)) {
    keys <- rep(NA_character_, nrow(se))
  }
  pending <- is.na(keys) & !is.na(identifiers)
  targets <- unique(identifiers[pending])
  found <- lookup_inchi_keys(targets, to_cid)
  keys[pending] <- found[match(identifiers[pending], targets)]

  if (anyNA(found)) {
    unresolved <- targets[is.na(found)]
    warning("No single InChIKey found for ", length(unresolved), " of ",
            length(targets), " ", label, ": ",
            paste(utils::head(unresolved, 5), collapse = ", "),
            if (length(unresolved) > 5L) ", ...", call. = FALSE)
  }

  row_data$opt_global_inchi_key <- keys
  SummarizedExperiment::rowData(se) <- row_data
  se
}

# InChIKey of each identifier, NA for those that give none.
lookup_inchi_keys <- function(identifiers, to_cid) {
  cids <- vapply(identifiers, to_cid, character(1), USE.NAMES = FALSE)
  known <- unique(cids[!is.na(cids)])
  batches <- split(known, ceiling(seq_along(known) / 100))
  # c() so that no batches still gives a vector to index, with NA for every URI
  keys <- c(character(0), unlist(unname(lapply(batches, cid_inchi_keys))))
  unname(keys[cids])
}

# PubChem compound identifier of a URI, NA if it is unknown or not recognized.
uri_cid <- function(uri) {
  pubchem <- "^https?://pubchem\\.ncbi\\.nlm\\.nih\\.gov/compound/([0-9]+)/?$"
  if (grepl(pubchem, uri)) {
    return(sub(pubchem, "\\1", uri))
  }
  kegg <- "^https?://www\\.kegg\\.jp/entry/(C[0-9]{5})/?$"
  if (grepl(kegg, uri)) {
    return(kegg_cid(sub(kegg, "\\1", uri)))
  }
  NA_character_
}

# PubChem compound identifier of a KEGG or HMDB identifier, NA if it is unknown
# or not recognized.
cmpdid_cid <- function(id) {
  if (grepl("^C[0-9]{5}$", id)) {
    return(kegg_cid(id))
  }
  if (grepl("^HMDB([0-9]{5}|[0-9]{7})$", id)) {
    return(hmdb_cid(id))
  }
  NA_character_
}

# PubChem compound identifier of a database_identifier value, NA if it is unknown
# or not recognized.
database_identifier_cid <- function(id) {
  id <- sub("^(KEGG|HMDB):", "", id)
  if (grepl("^CHEBI:[0-9]+$", id)) {
    return(source_cid("ChEBI", id))
  }
  if (grepl("^(PubChem:)?[0-9]+$", id)) {
    return(sub("^PubChem:", "", id))
  }
  if (grepl("^[0-9]{2,7}-[0-9]{2}-[0-9]$", id)) {
    return(single_cid(pubchem_get(sprintf("compound/name/%s/cids/TXT", id))))
  }
  cmpdid_cid(id)
}

kegg_cid <- function(kegg_id) {
  source_cid("KEGG", kegg_id)
}

# PubChem only knows the seven digit form of HMDB identifiers.
hmdb_cid <- function(hmdb_id) {
  digits <- sub("^HMDB", "", hmdb_id)
  padded <- paste0("HMDB", strrep("0", 7L - nchar(digits)), digits)
  source_cid("Human Metabolome Database (HMDB)", padded)
}

# Compound standardized from the substance deposited by a source. A substance
# without a single compound gives NA rather than a guess.
source_cid <- function(source, id) {
  single_cid(pubchem_get(sprintf("substance/sourceid/%s/%s/cids/TXT",
                                 utils::URLencode(source, reserved = TRUE), id),
                         cids_type = "standardized"))
}

# The compound identifier in a PubChem response, NA unless there is exactly one.
single_cid <- function(body) {
  if (is.na(body)) {
    return(NA_character_)
  }
  cids <- strsplit(trimws(body), "\\s+")[[1]]
  if (length(cids) != 1L) {
    return(NA_character_)
  }
  cids
}

# InChIKeys named by compound identifier; compounds without one are left out.
# PubChem rejects a whole request for one invalid identifier, so a rejected batch
# is split until only the invalid identifiers are left out.
cid_inchi_keys <- function(cids) {
  body <- pubchem_get(sprintf("compound/cid/%s/property/InChIKey/CSV",
                              paste(cids, collapse = ",")))
  if (is.na(body)) {
    if (length(cids) == 1L) {
      return(character(0))
    }
    half <- seq_len(length(cids) %/% 2L)
    return(c(cid_inchi_keys(cids[half]), cid_inchi_keys(cids[-half])))
  }
  properties <- utils::read.csv(text = body, colClasses = "character")
  keys <- properties$InChIKey
  names(keys) <- properties$CID
  # A compound without an InChIKey can be listed with an empty value
  keys[nzchar(keys)]
}

# Body of a PubChem PUG REST response. PubChem answers an identifier it does not
# know with 400 or 404, which gives NA; any other failure is an error.
pubchem_get <- function(path, ...) {
  req <- httr2::request("https://pubchem.ncbi.nlm.nih.gov/rest/pug")
  req <- httr2::req_url_path_append(req, path)
  req <- httr2::req_url_query(req, ...)
  req <- httr2::req_throttle(req, rate = 5)
  req <- httr2::req_retry(req, max_tries = 5)
  req <- httr2::req_error(req, is_error = function(resp) {
    status <- httr2::resp_status(resp)
    status >= 400L && !status %in% c(400L, 404L)
  })
  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 400L) {
    return(NA_character_)
  }
  httr2::resp_body_string(resp)
}
