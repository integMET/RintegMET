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

`study_variable_summarized_experiment()` converts an mzTab-M file to a
`SummarizedExperiment`, with one row per SML row and one column per study
variable. The `abundance` assay holds the `abundance_study_variable[n]` columns,
`rowData` holds the other SML columns, and `colData` holds the `study_variable[n]`
fields of the metadata section. In the example file these values are the means
of the available assays, but other files may report them differently.

It takes a file path or an `MzTabM` object from `RmzTabM::readMzTabM()`. File
validation with the jmzTab-M validator is off by default; turn it on with
`validate = TRUE`.

```r
library(RintegMET)

file <- system.file("extdata", "ST000010_AN000025.mztab", package = "RintegMET")
se <- study_variable_summarized_experiment(file)
se
#> class: SummarizedExperiment
#> dim: 403 11
#> metadata(0):
#> assays(1): abundance
#> rownames(403): 1 2 ... 402 403
#> rowData names(22): SML_ID SMF_ID_REFS ... opt_global_other_id
#>   opt_global_other_id_type
#> colnames(11): study_variable[1] study_variable[2] ...
#>   study_variable[10] study_variable[11]
#> colData names(5): name assay_refs average_function variation_function
#>   description
```

`populate_inchi_key()` fills the `opt_global_inchi_key` column of `rowData` by
looking up each molecule's `uri`, a KEGG or PubChem compound page, in PubChem. It
needs internet access, keeps existing values, and warns about URIs that give no
single InChIKey.

```r
se <- populate_inchi_key(se)
```
