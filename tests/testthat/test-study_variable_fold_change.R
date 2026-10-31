example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

# Two variables: `dose` with three levels and `time` with one, so that `time`
# has nothing to be compared with.
make_se <- function(description = c("Variable dose, value 0",
                                    "Variable time, value 1",
                                    "Variable dose, value 1",
                                    "Variable dose, value 2")) {
  ids <- sprintf("study_variable[%d]", seq_along(description))
  abundance <- matrix(c(2, 0, NA, 0,
                        9, 9, 9, 9,
                        4, 5, 1, 0,
                        8, 0, 3, 7),
                      nrow = 4, dimnames = list(as.character(1:4), ids))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = abundance),
    rowData = data.frame(SML_ID = as.character(1:4),
                         chemical_name = c("A", NA, "C", "D"),
                         row.names = as.character(1:4)),
    colData = data.frame(name = c("0", "1", "1", "2"),
                         description = description,
                         row.names = ids)
  )
}

test_that("every pair within a variable gets one column", {
  fc <- study_variable_fold_change(make_se())

  expect_s4_class(fc, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assayNames(fc), "fold_change")
  expect_equal(colnames(fc), c("study_variable[3]/study_variable[1]",
                               "study_variable[4]/study_variable[1]",
                               "study_variable[4]/study_variable[3]"))
})

test_that("fold-changes are the later study variable over the earlier one", {
  assay <- SummarizedExperiment::assay(study_variable_fold_change(make_se()))

  expect_equal(unname(assay[1, ]), c(4 / 2, 8 / 2, 8 / 4))
  expect_equal(rownames(assay), as.character(1:4))
})

test_that("missing abundances and zeros follow arithmetic", {
  assay <- SummarizedExperiment::assay(study_variable_fold_change(make_se()))

  # Second molecule: 0, 5, 0 across the three levels of `dose`
  expect_equal(unname(assay[2, ]), c(Inf, NaN, 0))
  # Third molecule: NA, 1, 3
  expect_equal(unname(assay[3, ]), c(NA, NA, 3))
  # Fourth molecule: 0, 0, 7
  expect_equal(unname(assay[4, ]), c(NaN, Inf, Inf))
})

test_that("colData describes each pair", {
  col_data <- SummarizedExperiment::colData(study_variable_fold_change(make_se()))

  expect_named(col_data, c("variable", "numerator", "denominator",
                           "numerator_id", "denominator_id"))
  expect_equal(col_data$variable, rep("dose", 3))
  expect_equal(col_data$numerator, c("1", "2", "2"))
  expect_equal(col_data$denominator, c("0", "0", "1"))
  expect_equal(col_data$numerator_id,
               sprintf("study_variable[%d]", c(3, 4, 4)))
  expect_equal(col_data$denominator_id,
               sprintf("study_variable[%d]", c(1, 1, 3)))
})

test_that("rowData is kept", {
  se <- make_se()

  expect_equal(SummarizedExperiment::rowData(study_variable_fold_change(se)),
               SummarizedExperiment::rowData(se))
})

test_that("other descriptions are compared as a whole, without NA", {
  fc <- study_variable_fold_change(make_se(c("treated", "treated", "other", NA)))

  expect_equal(colnames(fc), "study_variable[2]/study_variable[1]")
  expect_equal(fc$variable, "treated")
})

test_that("no comparable study variables gives zero columns", {
  fc <- study_variable_fold_change(make_se(paste("group", 1:4)))

  expect_equal(dim(fc), c(4, 0))
  expect_equal(SummarizedExperiment::assayNames(fc), "fold_change")
})

test_that("a single small molecule still gives a matrix", {
  fc <- study_variable_fold_change(make_se()[1, ])

  expect_equal(dim(fc), c(1, 3))
  expect_equal(rownames(fc), "1")
})

test_that("input without the expected contents is an error", {
  expect_error(study_variable_fold_change(data.frame()), "SummarizedExperiment")

  se <- make_se()
  se$description <- NULL
  expect_error(study_variable_fold_change(se), "description")

  se <- make_se()
  SummarizedExperiment::assayNames(se) <- "other"
  expect_error(study_variable_fold_change(se), "abundance")
})

test_that("the example file compares the levels of each variable", {
  se <- study_variable_summarized_experiment(example_file())
  abundance <- SummarizedExperiment::assay(se)

  fc <- study_variable_fold_change(se)

  expect_equal(dim(fc), c(403, 13))
  expect_equal(as.vector(table(fc$variable)[c("FCS", "Hours", "REFED", "TGF")]),
               c(1L, 10L, 1L, 1L))
  tgf <- fc[, fc$variable == "TGF"]
  expect_equal(tgf$numerator, "TGF")
  expect_equal(tgf$denominator, "NO_TGF")
  expect_equal(as.vector(SummarizedExperiment::assay(tgf)),
               unname(abundance[, "study_variable[11]"] /
                        abundance[, "study_variable[10]"]))
})
