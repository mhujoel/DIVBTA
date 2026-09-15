# Diabetes trials in Cochrane reviews

Data and analysis code for a case study of clinical diabetes trials included in a sample of Cochrane reviews.

## Contents

| File | Description |
| --- | --- |
| [github.sas7bdat](github.sas7bdat) | SAS dataset used in the main analyses |
| [github-ref.sas7bdat](github-ref.sas7bdat) | SAS dataset documenting the source of each abstracted standard deviation |
| [README_DIABETES.docx](README_DIABETES.docx) | Variable definitions |
| [Diabetes_DIVBTA.R](Diabetes_DIVBTA.R) | R script for the key DIVBTA analyses |
| [Diabetes_rho.R](Diabetes_rho.R) | R script that produces the publication tables |

## Requirements

- SAS, or an R package that can read `.sas7bdat` files (for example [`haven`](https://haven.tidyverse.org/))
- R, plus the packages listed at the top of each script

## Usage

1. Download or clone this repository.
2. See [README_DIABETES.docx](README_DIABETES.docx) for variable definitions. GitHub will not preview this Word file; open it locally.
3. Load `github.sas7bdat` for the main analyses and `github-ref.sas7bdat` for standard-deviation source documentation.
4. Run `Diabetes_DIVBTA.R` and `Diabetes_rho.R`.

