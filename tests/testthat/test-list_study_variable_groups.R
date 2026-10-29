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
                         "study_variable_refs", "study_variables"))
})

test_that("unreadable input is an error", {
  expect_error(list_study_variable_groups(tempfile()), "does not exist")

  path <- tempfile()
  writeLines("SMH\tSML_ID", path)
  expect_error(list_study_variable_groups(path), "MTD")
})
