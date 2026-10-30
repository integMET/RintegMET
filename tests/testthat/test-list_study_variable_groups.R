example_file <- function() {
  system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
}

write_mtd <- function(lines) {
  path <- tempfile(fileext = ".mztab")
  writeLines(c("MTD\tmzTab-version\t2.1.0-M", lines), path)
  path
}

test_that("groups in the example file are listed", {
  groups <- list_study_variable_groups(example_file())

  expect_s3_class(groups, "data.frame")
  expect_equal(groups$id, sprintf("study_variable_group[%d]", 1:4))
  expect_equal(groups$name, c("FCS", "Hours", "REFED", "TGF"))
  expect_equal(groups$description[1], "mwTab sample factor FCS")
  expect_equal(groups$type, c("categorical variable", "continuous variable",
                              "categorical variable", "categorical variable"))
  expect_equal(groups$datatype, c("xsd:string", "xsd:decimal",
                                  "xsd:string", "xsd:string"))
})

test_that("study variable references are resolved to names", {
  groups <- list_study_variable_groups(example_file())

  expect_equal(groups$study_variable_refs[[1]],
               c("study_variable[1]", "study_variable[2]"))
  expect_equal(groups$study_variables[[1]], c("FCS", "NO_FCS"))
  expect_equal(groups$study_variables[[2]], c("72", "0", "12", "24", "2"))
  expect_equal(lengths(groups$study_variables), c(2, 5, 2, 2))
})

test_that("samples are counted across the study variables of a group", {
  expect_equal(list_study_variable_groups(example_file())$n_samples,
               c(39, 39, 39, 39))

  variable_lines <- c(
    "MTD\tstudy_variable_group[1]\t[,,dose,]",
    paste0("MTD\tstudy_variable_group[1]-study_variable_refs\t",
           "study_variable[1]|study_variable[2]"),
    "MTD\tstudy_variable[1]\tlow",
    "MTD\tstudy_variable[1]-assay_refs\tassay[1]|assay[2]",
    "MTD\tstudy_variable[2]\thigh",
    "MTD\tstudy_variable[2]-assay_refs\tassay[2]|assay[3]"
  )

  # Without sample_refs each distinct assay counts as one sample
  expect_equal(list_study_variable_groups(write_mtd(variable_lines))$n_samples,
               3)

  # Assays of the same sample count once
  path <- write_mtd(c(variable_lines,
                      "MTD\tassay[1]-sample_ref\tsample[1]",
                      "MTD\tassay[2]-sample_ref\tsample[2]",
                      "MTD\tassay[3]-sample_ref\tsample[2]"))
  expect_equal(list_study_variable_groups(path)$n_samples, 2)
})

test_that("groups are ordered by numeric index", {
  path <- write_mtd(c(
    "MTD\tstudy_variable_group[10]\t[,,ten,]",
    "MTD\tstudy_variable_group[2]\t[,,two,]"
  ))

  groups <- list_study_variable_groups(path)

  expect_equal(groups$name, c("two", "ten"))
})

test_that("missing fields become NA or empty", {
  path <- write_mtd("MTD\tstudy_variable_group[1]\t[,,dose,]")

  groups <- list_study_variable_groups(path)

  expect_equal(nrow(groups), 1)
  expect_true(is.na(groups$description))
  expect_true(is.na(groups$type))
  expect_true(is.na(groups$datatype))
  expect_equal(groups$n_samples, 0)
  expect_equal(groups$study_variable_refs[[1]], character(0))
})

test_that("an MzTabM object can be passed instead of a path", {
  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)

  expect_equal(list_study_variable_groups(mztab),
               list_study_variable_groups(example_file()))
})

test_that("a file without groups gives a zero-row data frame", {
  groups <- list_study_variable_groups(write_mtd("MTD\ttitle\tnone"))

  expect_equal(nrow(groups), 0)
  expect_named(groups, c("id", "name", "description", "type", "datatype",
                         "n_samples", "study_variable_refs",
                         "study_variables"))
})

test_that("samples_per_group() returns the samples of each level", {
  expected <- data.frame(
    Group = rep(c("FCS", "Hours", "REFED", "TGF"), c(2, 5, 2, 2)),
    Level = c("FCS", "NO_FCS", "72", "0", "12", "24", "2", "NO_REFED",
              "REFED_48", "NO_TGF", "TGF"),
    Samples = c(6L, 33L, 18L, 3L, 6L, 6L, 6L, 33L, 6L, 21L, 18L)
  )

  expect_equal(samples_per_group(example_file()), expected)

  mztab <- RmzTabM::readMzTabM(example_file(), validate = FALSE)
  expect_equal(samples_per_group(mztab), expected)
})

test_that("samples_per_group() counts assays of the same sample once", {
  path <- write_mtd(c(
    "MTD\tstudy_variable_group[1]\t[,,dose,]",
    paste0("MTD\tstudy_variable_group[1]-study_variable_refs\t",
           "study_variable[1]|study_variable[2]"),
    "MTD\tstudy_variable[1]\tlow",
    "MTD\tstudy_variable[1]-assay_refs\tassay[1]|assay[2]",
    "MTD\tstudy_variable[2]\thigh",
    "MTD\tstudy_variable[2]-assay_refs\tassay[3]",
    "MTD\tassay[1]-sample_ref\tsample[1]",
    "MTD\tassay[2]-sample_ref\tsample[1]",
    "MTD\tassay[3]-sample_ref\tsample[2]"
  ))

  expect_equal(samples_per_group(path),
               data.frame(Group = c("dose", "dose"), Level = c("low", "high"),
                          Samples = c(1L, 1L)))
})

test_that("samples_per_group() gives zero rows for a file without groups", {
  counts <- samples_per_group(write_mtd("MTD\ttitle\tnone"))

  expect_equal(nrow(counts), 0)
  expect_named(counts, c("Group", "Level", "Samples"))
})

test_that("unreadable input is an error", {
  expect_error(list_study_variable_groups(tempfile()), "does not exist")

  path <- tempfile()
  writeLines("SMH\tSML_ID", path)
  expect_error(list_study_variable_groups(path), "MTD")
})
