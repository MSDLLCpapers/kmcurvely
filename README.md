
# kmcurvely

<!-- badges: start -->
<!-- badges: end -->

The goal of kmcurvely is to ...

## Installation

You can install the development version of kmcurvely from [GitHub](https://github.com/) with:

```r
# install.packages("remotes")
remotes::install_github("MSDLLCPapers/kmcurvely")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Example usage of kmcurvely function

meta <- meta_tte_example()
kmcurvely(meta = meta,
           population = "apat",
           observation = "efficacy_population",
           endpoint = "pfs;os",
           subgroup = "male;female")
```

