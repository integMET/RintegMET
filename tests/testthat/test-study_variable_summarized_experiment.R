example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

# Two study variables, declared out of order to test the ordering by index and
# with `description` only on the first one. `sml_lines` are the small molecules,
# one per string.
write_mztab <- function(sml_lines = c("SML\t  1\tHMDB:0000001\tA\t1\t3\t2\t5\t0.1\t0.2",
                                      "SML\t  2\tnull\tnull\t1\t3\tnull\t5\tnull\t0.2")) {
  path <- tempfile(fileext = ".mztab")
  writeLines(c(
    "MTD\tmzTab-version\t2.1.0-M",
    "MTD\tstudy_variable[10]\ttrt",
    "MTD\tstudy_variable[10]-assay_refs\tassay[2]",
    "MTD\tstudy_variable[2]\tctrl",
    "MTD\tstudy_variable[2]-assay_refs\tassay[1]",
    "MTD\tstudy_variable[2]-description\tcontrol",
    paste("SMH\tSML_ID\tdatabase_identifier\tchemical_name",
          "abundance_assay[1]\tabundance_assay[2]",
          "abundance_study_variable[2]\tabundance_study_variable[10]",
          "abundance_variation_study_variable[2]",
          "abundance_variation_study_variable[10]", sep = "\t"),
    sml_lines
  ), path)
  path
}

test_that("assay holds the abundance_study_variable values", {
  se <- study_variable_summarized_experiment(write_mztab())

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assayNames(se), "abundance")
  expect_equal(
    SummarizedExperiment::assay(se),
    matrix(c(2, NA, 5, 5), nrow = 2,
           dimnames = list(c("1", "2"),
                           c("study_variable[2]", "study_variable[10]")))
  )
})

test_that("rowData holds the SML columns other than SMH and abundances", {
  row_data <- SummarizedExperiment::rowData(
    study_variable_summarized_experiment(write_mztab())
  )

  expect_named(row_data, c("SML_ID", "database_identifier", "chemical_name"))
  expect_equal(row_data$SML_ID, c("1", "2"))
  expect_equal(row_data$database_identifier, c("HMDB:0000001", NA))
  expect_equal(row_data$chemical_name, c("A", NA))
})

test_that("colData holds the study_variable fields of MTD", {
  col_data <- SummarizedExperiment::colData(
    study_variable_summarized_experiment(write_mztab())
  )

  expect_equal(rownames(col_data), c("study_variable[2]", "study_variable[10]"))
  expect_named(col_data, c("name", "assay_refs", "description"))
  expect_equal(col_data$name, c("ctrl", "trt"))
  expect_equal(col_data$assay_refs, c("assay[1]", "assay[2]"))
  expect_equal(col_data$description, c("control", NA))
})

test_that("a single small molecule still gives a matrix", {
  se <- study_variable_summarized_experiment(
    write_mztab("SML\t1\tnull\tA\t1\t3\t2\t5\t0.1\t0.2")
  )

  expect_equal(dim(se), c(1, 2))
  expect_equal(rownames(se), "1")
  expect_equal(SummarizedExperiment::rowData(se)$chemical_name, "A")
})

test_that("a file without study variables gives zero columns", {
  path <- tempfile()
  writeLines(c("MTD\tmzTab-version\t2.1.0-M",
               "SMH\tSML_ID\tchemical_name",
               "SML\t1\tA"), path)

  expect_equal(dim(study_variable_summarized_experiment(path)), c(1, 0))
})

test_that("a study variable without an abundance column is an error", {
  path <- tempfile()
  writeLines(c("MTD\tmzTab-version\t2.1.0-M",
               "MTD\tstudy_variable[1]\tctrl",
               "SMH\tSML_ID\tchemical_name",
               "SML\t1\tA"), path)

  expect_error(study_variable_summarized_experiment(path),
               "abundance_study_variable[1]", fixed = TRUE)
})

test_that("a file without an SML section is an error", {
  path <- tempfile()
  writeLines("MTD\tmzTab-version\t2.1.0-M", path)

  expect_error(study_variable_summarized_experiment(path), "SML")
})

test_that("the example file gives one column per study variable", {
  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)
  sml <- RmzTabM::sml(mztab)
  columns <- sprintf("abundance_study_variable[%d]", 1:11)
  reported <- matrix(as.numeric(null_to_na(sml[, columns])), nrow = nrow(sml))

  se <- study_variable_summarized_experiment(mztab)

  expect_equal(dim(se), c(403, 11))
  expect_equal(colnames(se), sprintf("study_variable[%d]", 1:11))
  expect_equal(unname(SummarizedExperiment::assay(se)), reported)
  expect_equal(se$name[10:11], c("NO_TGF", "TGF"))
  expect_named(SummarizedExperiment::colData(se),
               c("name", "assay_refs", "average_function",
                 "variation_function", "description"))
  expect_equal(ncol(SummarizedExperiment::rowData(se)), 22)
  expect_false(any(grepl("abundance|^SMH$",
                         names(SummarizedExperiment::rowData(se)))))
})

test_that("a path and an MzTabM object give the same result", {
  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)

  expect_equal(study_variable_summarized_experiment(example_file()),
               study_variable_summarized_experiment(mztab))
})
