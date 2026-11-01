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
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  row_data <- SummarizedExperiment::rowData(se)
  if (!"uri" %in% names(row_data)) {
    stop("rowData(se) lacks the column uri", call. = FALSE)
  }

  uris <- row_data$uri
  keys <- row_data$opt_global_inchi_key
  if (is.null(keys)) {
    keys <- rep(NA_character_, nrow(se))
  }
  pending <- is.na(keys) & !is.na(uris)
  targets <- unique(uris[pending])
  found <- uri_inchi_keys(targets)
  keys[pending] <- found[match(uris[pending], targets)]

  if (anyNA(found)) {
    unresolved <- targets[is.na(found)]
    warning("No single InChIKey found for ", length(unresolved), " of ",
            length(targets), " URIs: ", paste(utils::head(unresolved, 5),
                                              collapse = ", "),
            if (length(unresolved) > 5L) ", ...", call. = FALSE)
  }

  row_data$opt_global_inchi_key <- keys
  SummarizedExperiment::rowData(se) <- row_data
  se
}

# InChIKey of each URI, NA for those that give none.
uri_inchi_keys <- function(uris) {
  cids <- vapply(uris, uri_cid, character(1), USE.NAMES = FALSE)
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

# Compound standardized from the substance deposited by KEGG. A substance
# without a single compound gives NA rather than a guess.
kegg_cid <- function(kegg_id) {
  body <- pubchem_get(sprintf("substance/sourceid/KEGG/%s/cids/TXT", kegg_id),
                      cids_type = "standardized")
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
