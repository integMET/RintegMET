# RintegMET
Create iDMET network from mzTab-M files

## Installation

RintegMET parses mzTab-M files with
[RmzTabM](https://github.com/RforMassSpectrometry/RmzTabM), which depends on
Bioconductor packages, so install it with BiocManager:

```r
# install.packages("BiocManager")
BiocManager::install("RforMassSpectrometry/RmzTabM")
BiocManager::install("integMET/RintegMET")
```

## Usage

`list_study_variable_groups()` lists the `study_variable_group` entries
(experimental factors) in the metadata section of an mzTab-M file, together
with the study variables (factor levels) each group refers to. It takes a
file path or an `MzTabM` object from `RmzTabM::readMzTabM()`. File validation
with the jmzTab-M validator is off by default; turn it on with
`validate = TRUE`.

```r
library(RintegMET)

file <- system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
groups <- list_study_variable_groups(file)

groups[, c("id", "name", "type", "datatype")]
#>                        id  name                 type    datatype
#> 1 study_variable_group[1]   FCS categorical variable  xsd:string
#> 2 study_variable_group[2] Hours  continuous variable xsd:decimal
#> 3 study_variable_group[3] REFED categorical variable  xsd:string
#> 4 study_variable_group[4]   TGF categorical variable  xsd:string

groups$study_variables
#> [[1]]
#> [1] "FCS"    "NO_FCS"
#>
#> [[2]]
#> [1] "72" "0"  "12" "24" "2"
#>
#> [[3]]
#> [1] "NO_REFED" "REFED_48"
#>
#> [[4]]
#> [1] "NO_TGF" "TGF"
```

`samples_per_group()` counts the distinct samples of each level of each group.

```r
samples_per_group(file)
#>    Group    Level Samples
#> 1    FCS      FCS       6
#> 2    FCS   NO_FCS      33
#> 3  Hours       72      18
#> 4  Hours        0       3
#> 5  Hours       12       6
#> 6  Hours       24       6
#> 7  Hours        2       6
#> 8  REFED NO_REFED      33
#> 9  REFED REFED_48       6
#> 10   TGF   NO_TGF      21
#> 11   TGF      TGF      18
```

`sml_fold_change()` compares two of these study variables, a control and an
experimental group, and returns the fold-change of their mean abundances for
every row of the SML section. Means are taken over the assays of each study
variable, ignoring missing values unless `na.rm = FALSE`.

```r
fc <- sml_fold_change(file, control = "NO_TGF", experimental = "TGF")
```

`sml_summarized_experiment()` returns a `SummarizedExperiment` with one row per
SML row and one column per level. The `mean_abundance` assay holds the mean
abundance of each level, `colData` is the `samples_per_group()` table, and
`rowData` holds `SML_ID`, `database_identifier`, `chemical_name` and `uri`.

```r
se <- sml_summarized_experiment(file)
se
#> class: SummarizedExperiment
#> dim: 403 11
#> metadata(0):
#> assays(1): mean_abundance
#> rownames(403): 1 2 ... 402 403
#> rowData names(4): SML_ID database_identifier chemical_name uri
#> colnames(11): study_variable[1] study_variable[2] ... study_variable[10]
#>   study_variable[11]
#> colData names(3): Group Level Samples
```

`study_variable_summarized_experiment()` instead keeps the values as reported in
the file. The `abundance` assay holds the `abundance_study_variable[n]` columns,
`rowData` holds the other SML columns, and `colData` holds the `study_variable[n]`
fields of the metadata section. In the example file these values are the means
of the available assays, but other files may report them differently.

```r
se <- study_variable_summarized_experiment(file)
```
