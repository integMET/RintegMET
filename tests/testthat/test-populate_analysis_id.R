example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

analysis_param <- function(id) {
  sprintf("[,, Metabolomics Workbench analysis ID, %s]", id)
}

# Copy of the example file with its custom[n] fields replaced.
file_with_params <- function(params) {
  lines <- readLines(example_file())
  fields <- grep("^MTD\tcustom\\[", lines)
  new_fields <- sprintf("MTD\tcustom[%d]\t%s", seq_along(params), params)
  lines <- append(lines[-fields], new_fields, after = fields[1L] - 1L)
  file <- tempfile(fileext = ".mztab")
  writeLines(lines, file)
  file
}

make_fc <- function() {
  study_variable_fold_change(study_variable_summarized_experiment(example_file()))
}

test_that("the example file gives its analysis", {
  fc <- make_fc()
  result <- populate_analysis_id(fc, example_file())

  expect_equal(result$analysis_id, rep("AN000025", ncol(fc)))
})

test_that("only the column is added", {
  fc <- make_fc()
  result <- populate_analysis_id(fc, example_file())

  col_data <- SummarizedExperiment::colData(result)
  expect_named(col_data, c(names(SummarizedExperiment::colData(fc)),
                           "analysis_id"))
  expect_equal(colnames(result), colnames(fc))
  expect_equal(SummarizedExperiment::assay(result),
               SummarizedExperiment::assay(fc))
  expect_equal(SummarizedExperiment::rowData(result),
               SummarizedExperiment::rowData(fc))
})

test_that("it keeps the compact identifier of the study", {
  fc <- populate_compact_identifier(make_fc(), example_file())
  result <- populate_analysis_id(fc, example_file())

  expect_equal(result$compact_identifier, fc$compact_identifier)
  expect_equal(unique(result$analysis_id), "AN000025")
})

test_that("an MzTabM object is accepted", {
  fc <- make_fc()
  x <- RmzTabM::readMzTabM(example_file(), validate = FALSE)

  expect_equal(unique(populate_analysis_id(fc, x)$analysis_id), "AN000025")
})

test_that("existing values are replaced", {
  fc <- make_fc()
  fc$analysis_id <- "old"

  result <- populate_analysis_id(fc, example_file())

  expect_equal(unique(result$analysis_id), "AN000025")
  expect_equal(sum(names(SummarizedExperiment::colData(result)) ==
                     "analysis_id"), 1)
})

test_that("the analysis is found among other parameters", {
  fc <- make_fc()
  file <- file_with_params(c("[,, Other parameter, value]",
                             analysis_param("AN000123")))

  expect_equal(unique(populate_analysis_id(fc, file)$analysis_id), "AN000123")
})

test_that("parameters without a single analysis give NA with a warning", {
  fc <- make_fc()
  cases <- list(
    "[,, Other parameter, AN000025]",
    analysis_param("ST000010"),
    c(analysis_param("AN000025"), analysis_param("AN000026"))
  )
  for (params in cases) {
    expect_warning(
      result <- populate_analysis_id(fc, file_with_params(params)),
      "No single"
    )
    expect_equal(result$analysis_id, rep(NA_character_, ncol(fc)))
  }
})

test_that("a file without custom fields gives NA with a warning", {
  lines <- readLines(example_file())
  file <- tempfile(fileext = ".mztab")
  writeLines(lines[!grepl("^MTD\tcustom\\[", lines)], file)

  expect_warning(result <- populate_analysis_id(make_fc(), file), "No single")
  expect_true(all(is.na(result$analysis_id)))
})

test_that("a repeated parameter still gives its analysis", {
  expect_equal(param_analysis_id(rep(analysis_param("AN000025"), 2)),
               "AN000025")
})

test_that("input without the expected contents is an error", {
  expect_error(populate_analysis_id(data.frame(), example_file()),
               "SummarizedExperiment")
})
