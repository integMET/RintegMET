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
