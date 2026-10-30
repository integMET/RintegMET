example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

# One group with two levels, where the two assays of `ctrl` come from the same
# sample. `sml_lines` are the small molecules, one per string.
write_mztab <- function(sml_lines = c("SML\t1\tHMDB:0000001\tA\thttps://a\t1\t3\t4\t8",
                                      "SML\t2\tnull\tnull\tnull\t2\tnull\t6\tnull",
                                      "SML\t3\tnull\tC\tnull\t0\t0\t5\t5")) {
  path <- tempfile(fileext = ".mztab")
  writeLines(c(
    "MTD\tmzTab-version\t2.1.0-M",
    "MTD\tstudy_variable_group[1]\t[,,treatment,]",
    paste0("MTD\tstudy_variable_group[1]-study_variable_refs\t",
           "study_variable[1]|study_variable[2]"),
    "MTD\tstudy_variable[1]\tctrl",
    "MTD\tstudy_variable[1]-assay_refs\tassay[1]|assay[2]",
    "MTD\tstudy_variable[2]\ttrt",
    "MTD\tstudy_variable[2]-assay_refs\tassay[3]|assay[4]",
    "MTD\tassay[1]-sample_ref\tsample[1]",
    "MTD\tassay[2]-sample_ref\tsample[1]",
    "MTD\tassay[3]-sample_ref\tsample[2]",
    "MTD\tassay[4]-sample_ref\tsample[3]",
    paste("SMH\tSML_ID\tdatabase_identifier\tchemical_name\turi",
          "abundance_assay[1]\tabundance_assay[2]",
          "abundance_assay[3]\tabundance_assay[4]", sep = "\t"),
    sml_lines
  ), path)
  path
}

test_that("assay holds the mean abundance of each study variable", {
  se <- sml_summarized_experiment(write_mztab())

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assayNames(se), "mean_abundance")
  expect_equal(
    SummarizedExperiment::assay(se),
    matrix(c(2, 2, 0, 6, 6, 5), nrow = 3,
           dimnames = list(c("1", "2", "3"),
                           c("study_variable[1]", "study_variable[2]")))
  )
})

test_that("na.rm = FALSE makes means with missing assays NA", {
  se <- sml_summarized_experiment(write_mztab(), na.rm = FALSE)

  expect_equal(unname(SummarizedExperiment::assay(se)),
               matrix(c(2, NA, 0, 6, NA, 5), nrow = 3))
})

test_that("colData is the data frame from samples_per_group()", {
  path <- write_mztab()
  col_data <- SummarizedExperiment::colData(sml_summarized_experiment(path))

  expect_equal(rownames(col_data), c("study_variable[1]", "study_variable[2]"))
  rownames(col_data) <- NULL
  expect_equal(as.data.frame(col_data), samples_per_group(path))
  expect_equal(col_data$Samples, c(1L, 2L))
})

test_that("rowData holds the SML annotation columns", {
  row_data <- SummarizedExperiment::rowData(sml_summarized_experiment(write_mztab()))

  expect_named(row_data, c("SML_ID", "database_identifier", "chemical_name",
                           "uri"))
  expect_equal(row_data$SML_ID, c("1", "2", "3"))
  expect_equal(row_data$database_identifier, c("HMDB:0000001", NA, NA))
  expect_equal(row_data$chemical_name, c("A", NA, "C"))
  expect_equal(row_data$uri, c("https://a", NA, NA))
})

test_that("a single small molecule still gives a matrix", {
  path <- write_mztab("SML\t1\tnull\tA\tnull\t1\t3\t4\t8")

  se <- sml_summarized_experiment(path)

  expect_equal(dim(se), c(1, 2))
  expect_equal(rownames(se), "1")
  expect_equal(rownames(sml_fold_change(path, "ctrl", "trt")), "1")
})

test_that("a file without groups gives zero columns", {
  path <- tempfile()
  writeLines(c("MTD\tmzTab-version\t2.1.0-M",
               "SMH\tSML_ID\tdatabase_identifier\tchemical_name\turi",
               "SML\t1\tnull\tA\tnull"), path)

  expect_equal(dim(sml_summarized_experiment(path)), c(1, 0))
})

test_that("a file without an SML section is an error", {
  path <- tempfile()
  writeLines("MTD\tmzTab-version\t2.1.0-M", path)

  expect_error(sml_summarized_experiment(path), "SML")
})

test_that("the example file gives one column per level", {
  se <- sml_summarized_experiment(example_file())

  expect_equal(dim(se), c(403, 11))
  expect_equal(colnames(se), sprintf("study_variable[%d]", 1:11))
  expect_equal(se$Level[10:11], c("NO_TGF", "TGF"))
  expect_equal(se$Samples[10:11], c(21L, 18L))
})

test_that("assay values agree with sml_fold_change()", {
  se <- sml_summarized_experiment(example_file())
  fc <- sml_fold_change(example_file(), "NO_TGF", "TGF")
  means <- SummarizedExperiment::assay(se)

  expect_equal(unname(means[, "study_variable[10]"]), fc$control_mean)
  expect_equal(unname(means[, "study_variable[11]"]), fc$experimental_mean)
})

test_that("means match every study variable abundance in the file", {
  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)
  columns <- sprintf("abundance_study_variable[%d]", 1:11)
  reported <- RmzTabM::sml(mztab)[, columns]
  reported <- matrix(as.numeric(null_to_na(reported)), nrow = nrow(reported))

  se <- sml_summarized_experiment(mztab)

  # The file stores rounded means
  expect_equal(unname(SummarizedExperiment::assay(se)), reported,
               tolerance = 1e-6)
})
