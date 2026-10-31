kegg_uri <- function(id) paste0("https://www.kegg.jp/entry/", id)
pubchem_uri <- function(cid) {
  paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cid)
}

make_se <- function(uri, ...) {
  rows <- as.character(seq_along(uri))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = matrix(1, nrow = length(uri), ncol = 1,
                                     dimnames = list(rows, "study_variable[1]"))),
    rowData = data.frame(SML_ID = rows, uri = uri, ..., row.names = rows,
                         stringsAsFactors = FALSE)
  )
}

inchi_keys <- function(se) SummarizedExperiment::rowData(se)$opt_global_inchi_key

# PubChem with two KEGG substances and three compounds, recording its requests.
fake_pubchem <- function() {
  requests <- character(0)
  list(
    get = function(path, ...) {
      requests <<- c(requests, path)
      if (grepl("/KEGG/C00001/", path)) return("11\n")
      if (grepl("/KEGG/C00002/", path)) return("22\n")
      if (grepl("/KEGG/C00003/", path)) return("33\n44\n")
      if (grepl("/KEGG/", path)) return(NA_character_)
      cids <- strsplit(sub("^compound/cid/([^/]+)/.*$", "\\1", path), ",")[[1]]
      known <- c("11" = "AAA-N", "22" = "BBB-N", "55" = "CCC-N")
      cids <- intersect(cids, names(known))
      paste(c("\"CID\",\"InChIKey\"",
              sprintf("%s,\"%s\"", cids, known[cids])), collapse = "\n")
    },
    requests = function() requests
  )
}

test_that("KEGG and PubChem URIs give InChIKeys", {
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  se <- populate_inchi_key(make_se(c(kegg_uri("C00001"), pubchem_uri(55), NA,
                                     kegg_uri("C00002")),
                                   opt_global_inchi_key = NA_character_))

  expect_equal(inchi_keys(se), c("AAA-N", "CCC-N", NA, "BBB-N"))
})

test_that("the column is added when absent and other columns are kept", {
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  original <- make_se(kegg_uri("C00001"))
  se <- populate_inchi_key(original)

  expect_equal(inchi_keys(se), "AAA-N")
  row_data <- SummarizedExperiment::rowData(se)
  expect_named(row_data, c("SML_ID", "uri", "opt_global_inchi_key"))
  expect_equal(SummarizedExperiment::assay(se),
               SummarizedExperiment::assay(original))
})

test_that("existing values are kept and not looked up", {
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  se <- populate_inchi_key(make_se(c(kegg_uri("C00001"), kegg_uri("C00002")),
                                   opt_global_inchi_key = c("KEEP-N", NA)))

  expect_equal(inchi_keys(se), c("KEEP-N", "BBB-N"))
  expect_false(any(grepl("C00001", pubchem$requests())))
})

test_that("a repeated URI is looked up once and fills every row", {
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  se <- populate_inchi_key(make_se(rep(kegg_uri("C00001"), 2)))

  expect_equal(inchi_keys(se), c("AAA-N", "AAA-N"))
  expect_equal(sum(grepl("C00001", pubchem$requests())), 1)
})

test_that("compounds are requested in batches of 100", {
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  expect_warning(populate_inchi_key(make_se(pubchem_uri(1:101))), "98 of 101")

  requests <- grep("^compound/cid/", pubchem$requests(), value = TRUE)
  expect_length(requests, 2)
  expect_equal(lengths(strsplit(requests, ",")), c(100, 1))
})

test_that("URIs without a single InChIKey stay NA with a warning", {
  uris <- c(kegg_uri("C00001"),  # found
            kegg_uri("C00003"),  # substance with two compounds
            kegg_uri("C00004"),  # unknown to PubChem
            pubchem_uri(99),     # compound without an InChIKey
            "https://example.org/compound/1",
            "https://www.kegg.jp/entry/R00001")
  pubchem <- fake_pubchem()
  local_mocked_bindings(pubchem_get = pubchem$get)

  expect_warning(se <- populate_inchi_key(make_se(uris)), "5 of 6 URIs")

  expect_equal(inchi_keys(se), c("AAA-N", rep(NA, 5)))
  expect_false(any(grepl("example.org|R00001", pubchem$requests())))
})

test_that("nothing to look up makes no requests and no warning", {
  local_mocked_bindings(pubchem_get = function(...) stop("unexpected request"))

  se <- make_se(c(NA_character_, "https://www.kegg.jp/entry/C00001"),
                opt_global_inchi_key = c(NA, "KEEP-N"))

  expect_no_warning(result <- populate_inchi_key(se))
  expect_equal(inchi_keys(result), c(NA, "KEEP-N"))
})

test_that("a batch that PubChem does not know leaves its URIs NA", {
  local_mocked_bindings(pubchem_get = function(...) NA_character_)

  expect_warning(se <- populate_inchi_key(make_se(pubchem_uri(11))), "1 of 1")

  expect_equal(inchi_keys(se), NA_character_)
})

test_that("input without the expected contents is an error", {
  expect_error(populate_inchi_key(data.frame()), "SummarizedExperiment")

  se <- make_se("https://www.kegg.jp/entry/C00001")
  SummarizedExperiment::rowData(se)$uri <- NULL
  expect_error(populate_inchi_key(se), "uri")
})

test_that("pubchem_get() returns the body of a successful response", {
  httr2::local_mocked_responses(function(req) {
    expect_match(req$url, "/rest/pug/compound/cid/1%2C2/property/InChIKey/CSV$")
    httr2::response(200L, body = charToRaw("found"))
  })

  expect_equal(pubchem_get("compound/cid/1,2/property/InChIKey/CSV"), "found")
})

test_that("pubchem_get() sends query parameters", {
  httr2::local_mocked_responses(function(req) {
    expect_match(req$url, "cids_type=standardized$")
    httr2::response(200L, body = charToRaw("found"))
  })

  expect_equal(pubchem_get("substance", cids_type = "standardized"), "found")
})

test_that("pubchem_get() gives NA for identifiers PubChem does not know", {
  for (status in c(400L, 404L)) {
    httr2::local_mocked_responses(function(req) httr2::response(status))

    expect_identical(pubchem_get("substance"), NA_character_)
  }
})

test_that("pubchem_get() fails on server errors instead of giving NA", {
  httr2::local_mocked_responses(function(req) httr2::response(500L))

  expect_error(pubchem_get("substance"), "500")
})
