example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

# Two study variables with two assays each and three small molecules: one
# fully quantified, one with missing assays and one absent from the control.
write_mztab <- function(extra_mtd = character()) {
  path <- tempfile(fileext = ".mztab")
  writeLines(c(
    "MTD\tmzTab-version\t2.1.0-M",
    "MTD\tstudy_variable[1]\tctrl",
    "MTD\tstudy_variable[1]-assay_refs\tassay[1]|assay[2]",
    "MTD\tstudy_variable[2]\ttrt",
    "MTD\tstudy_variable[2]-assay_refs\tassay[3]|assay[4]",
    extra_mtd,
    paste("SMH\tSML_ID\tdatabase_identifier\tchemical_name\turi",
          "abundance_assay[1]\tabundance_assay[2]",
          "abundance_assay[3]\tabundance_assay[4]", sep = "\t"),
    "SML\t1\tHMDB:0000001\tA\thttps://hmdb.ca/metabolites/HMDB0000001\t1\t3\t4\t8",
    "SML\t2\tnull\tnull\tnull\t2\tnull\t6\tnull",
    "SML\t3\tnull\tC\tnull\t0\t0\t5\t5"
  ), path)
  path
}

test_that("fold-change is the ratio of the group means", {
  fc <- sml_fold_change(write_mztab(), control = "ctrl", experimental = "trt")

  expect_equal(fc$SML_ID, c("1", "2", "3"))
  expect_equal(fc$database_identifier, c("HMDB:0000001", NA, NA))
  expect_equal(fc$chemical_name, c("A", NA, "C"))
  expect_equal(fc$uri, c("https://hmdb.ca/metabolites/HMDB0000001", NA, NA))
  expect_equal(fc$control_mean, c(2, 2, 0))
  expect_equal(fc$experimental_mean, c(6, 6, 5))
  expect_equal(fc$fold_change, c(3, 3, Inf))
  expect_equal(fc$log2_fold_change, log2(c(3, 3, Inf)))
})

test_that("na.rm = FALSE makes means with missing assays NA", {
  fc <- sml_fold_change(write_mztab(), "ctrl", "trt", na.rm = FALSE)

  expect_equal(fc$control_mean, c(2, NA, 0))
  expect_equal(fc$fold_change, c(3, NA, Inf))
})

test_that("study variables can be given by id", {
  expect_equal(
    sml_fold_change(write_mztab(), "study_variable[1]", "study_variable[2]"),
    sml_fold_change(write_mztab(), "ctrl", "trt")
  )
})

test_that("unknown and ambiguous study variables are errors", {
  expect_error(sml_fold_change(write_mztab(), "none", "trt"),
               "No study variable named \"none\". Available: ctrl, trt")

  path <- write_mztab(c("MTD\tstudy_variable[3]\tctrl",
                        "MTD\tstudy_variable[3]-assay_refs\tassay[4]"))
  expect_error(sml_fold_change(path, "ctrl", "trt"),
               "study_variable\\[1\\], study_variable\\[3\\]")
})

test_that("a file without an SML section is an error", {
  path <- tempfile()
  writeLines("MTD\tmzTab-version\t2.1.0-M", path)

  expect_error(sml_fold_change(path, "ctrl", "trt"), "SML")
})

test_that("means match the study variable abundances in the file", {
  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)
  reported <- RmzTabM::sml(mztab)[, "abundance_study_variable[10]"]
  reported <- as.numeric(null_to_na(reported))

  fc <- sml_fold_change(mztab, "NO_TGF", "TGF")

  expect_equal(nrow(fc), 403)
  expect_equal(fc$control_mean, reported)
})

test_that("missing assays are ignored by default in the example file", {
  fc <- sml_fold_change(example_file(), "NO_TGF", "TGF")
  strict <- sml_fold_change(example_file(), "NO_TGF", "TGF", na.rm = FALSE)

  expect_gt(sum(!is.na(fc$fold_change)), sum(!is.na(strict$fold_change)))
  expect_equal(fc$fold_change, fc$experimental_mean / fc$control_mean)
})
