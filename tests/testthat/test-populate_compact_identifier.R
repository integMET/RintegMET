example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

workbench_uri <- function(study) {
  paste0("https://www.metabolomicsworkbench.org/data/DRCCMetadata.php",
         "?Mode=Study&StudyID=", study)
}

# Copy of the example file with its uri[n] fields replaced.
file_with_uris <- function(uris) {
  lines <- readLines(example_file())
  fields <- grep("^MTD\turi\\[", lines)
  new_fields <- sprintf("MTD\turi[%d]\t%s", seq_along(uris), uris)
  lines <- append(lines[-fields], new_fields, after = fields[1L] - 1L)
  file <- tempfile(fileext = ".mztab")
  writeLines(lines, file)
  file
}

make_fc <- function() {
  study_variable_fold_change(study_variable_summarized_experiment(example_file()))
}

test_that("the example file gives its Metabolomics Workbench study", {
  fc <- make_fc()
  result <- populate_compact_identifier(fc, example_file())

  expect_equal(result$compact_identifier, rep("mw.study:ST000010", ncol(fc)))
})

test_that("only the column is added", {
  fc <- make_fc()
  result <- populate_compact_identifier(fc, example_file())

  col_data <- SummarizedExperiment::colData(result)
  expect_named(col_data, c(names(SummarizedExperiment::colData(fc)),
                           "compact_identifier"))
  expect_equal(colnames(result), colnames(fc))
  expect_equal(SummarizedExperiment::assay(result),
               SummarizedExperiment::assay(fc))
  expect_equal(SummarizedExperiment::rowData(result),
               SummarizedExperiment::rowData(fc))
})

test_that("an MzTabM object is accepted", {
  fc <- make_fc()
  x <- RmzTabM::readMzTabM(example_file(), validate = FALSE)

  expect_equal(populate_compact_identifier(fc, x)$compact_identifier,
               rep("mw.study:ST000010", ncol(fc)))
})

test_that("existing values are replaced", {
  fc <- make_fc()
  fc$compact_identifier <- "old"

  result <- populate_compact_identifier(fc, example_file())

  expect_equal(unique(result$compact_identifier), "mw.study:ST000010")
  expect_equal(sum(names(SummarizedExperiment::colData(result)) ==
                     "compact_identifier"), 1)
})

test_that("the study is found among other URIs", {
  fc <- make_fc()
  file <- file_with_uris(c("https://example.org/study/1",
                           workbench_uri("ST000123")))

  expect_equal(unique(populate_compact_identifier(fc, file)$compact_identifier),
               "mw.study:ST000123")
})

test_that("URIs without a single study give NA with a warning", {
  fc <- make_fc()
  cases <- list(
    "https://example.org/data?StudyID=ST000010",
    "https://www.metabolomicsworkbench.org/data/?StudyID=PR000010",
    c(workbench_uri("ST000010"), workbench_uri("ST000011"))
  )
  for (uris in cases) {
    expect_warning(result <- populate_compact_identifier(fc, file_with_uris(uris)),
                   "No single")
    expect_equal(result$compact_identifier, rep(NA_character_, ncol(fc)))
  }
})

test_that("a repeated URI still gives its study", {
  expect_equal(uri_compact_identifier(rep(workbench_uri("ST000010"), 2)),
               "mw.study:ST000010")
})

test_that("URI variants give the study and no URIs give NA", {
  expect_equal(
    uri_compact_identifier(
      "http://metabolomicsworkbench.org/data/x.php?StudyID=ST000010&Mode=Study"
    ),
    "mw.study:ST000010"
  )
  expect_identical(uri_compact_identifier(character(0)), NA_character_)
})

test_that("input without the expected contents is an error", {
  expect_error(populate_compact_identifier(data.frame(), example_file()),
               "SummarizedExperiment")
})
