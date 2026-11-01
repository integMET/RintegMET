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

`study_variable_fold_change()` compares the study variables of each variable,
taken from `colData(se)$description`, and returns a `SummarizedExperiment` with
one column per pair. The `fold_change` assay holds the study variable that comes
later in `se` over the one that comes earlier, `rowData` is that of `se`, and
`colData` holds the variable and the two study variables of each pair.

```r
fc <- study_variable_fold_change(se)
```

`populate_compact_identifier()` adds a `compact_identifier` column to `colData`
with the [identifiers.org](https://identifiers.org/) compact identifier of the
study, read from the `uri[n]` fields of the file's metadata section. It recognizes
Metabolomics Workbench studies, such as `mw.study:ST000010`, and MetaboLights
studies, such as `metabolights:MTBLS1066`, and warns and gives `NA` if the fields
do not name a single study.

```r
fc <- populate_compact_identifier(fc, file)
```

`populate_analysis_id()` adds an `analysis_id` column to `colData` with the
Metabolomics Workbench analysis, read from the `custom[n]` fields of the file's
metadata section, such as `AN000025`. It warns and gives `NA` if the fields do
not name a single analysis.

```r
fc <- populate_analysis_id(fc, file)
```

`populate_inchi_key()` fills the `opt_global_inchi_key` column of `rowData` by
looking up each molecule's `uri`, a KEGG or PubChem compound page, in PubChem. It
needs internet access, keeps existing values, and warns about URIs that give no
single InChIKey.

```r
se <- populate_inchi_key(se)
```

`populate_inchi_key_from_cmpdid()` does the same from the `opt_global_cmpdid`
column of `rowData`, for files that report their molecules there. It recognizes
KEGG compounds, such as `C00041`, and HMDB metabolites, such as `HMDB00688`, and
leaves any other identifier `NA` with a warning.

```r
fc <- populate_inchi_key_from_cmpdid(fc)
```

`populate_inchi_key_from_database_identifier()` does the same from the
`database_identifier` column of `rowData`, as in files converted from MetaboLights.
It recognizes ChEBI entities, such as `CHEBI:28300`, KEGG compounds, HMDB
metabolites, PubChem compounds and CAS Registry Numbers, and leaves any other value
`NA` with a warning, including several identifiers in one value.

```r
fc <- populate_inchi_key_from_database_identifier(fc)
```
