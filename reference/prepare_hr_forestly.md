# Prepare the dataset for the HR forest plot, alongwith its KM plots

Prepare the dataset for the HR forest plot, alongwith its KM plots

## Usage

``` r
prepare_hr_forestly(
  meta = NULL,
  population = NULL,
  observation = NULL,
  endpoint = NULL,
  subgroup = NULL,
  km_curves = NULL,
  arm_levels = NULL
)
```

## Arguments

- meta:

  A metadata object that contains the mappings for the variables of
  interest.

- population:

  A string indicating the population to use (e.g., "apat").

- observation:

  A string indicating the observation population (e.g.,
  "efficacy_population").

- endpoint:

  A semicolon-separated string of endpoints to be analyzed (e.g.,
  "pfs;os").

- subgroup:

  A semicolon-separated string of subgroups to filter by (e.g.,
  "age;gender"). The parameters should be defined at variable levels.
  Use "all" to include the result for overall population (e.g.,
  "all;age;gender").

- km_curves:

  A semicolon-separated string of subgroups for which KM curves should
  be plotted. The parameters should be defined at value levels for
  subgroups.

- arm_levels:

  A vector of character specifying the levels of arms, starting with the
  reference arm.

## Value

An metadata with HR per subgroup along with its KM plotting data

## Examples

``` r
prepare_hr_forestly(
  meta = meta_tte_example_new(),
  population = "apat",
  observation = "efficacy_population",
  endpoint = "pfs;os",
  subgroup = "age;gender",
  km_curves = "female;age65-80",
  arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
)
#> List of 17
#>  $ meta           :List of 7
#>  $ population     : chr "apat"
#>  $ observation    : chr "efficacy_population"
#>  $ parameter      : NULL
#>  $ n              :'data.frame': 10 obs. of  5 variables:
#>  $ order          : NULL
#>  $ group          : chr [1:3] "Placebo" "Xanomeline Low Dose" "Xanomeline High Dose"
#>  $ reference_group: num 1
#>  $ endpoint       : chr [1:2] "pfs" "os"
#>  $ subgroup       : chr [1:2] "age" "gender"
#>  $ kmcurves       : chr [1:2] "Female" "65-80"
#>  $ arm_comparison : chr [1:2] "Xanomeline Low Dose vs. Placebo" "Xanomeline High Dose vs. Placebo"
#>  $ event          :'data.frame': 10 obs. of  5 variables:
#>  $ hr_est         :'data.frame': 10 obs. of  5 variables:
#>  $ hr_ci_lower    :'data.frame': 10 obs. of  5 variables:
#>  $ hr_ci_upper    :'data.frame': 10 obs. of  5 variables:
#>  $ km_data        :'data.frame': 1590 obs. of  11 variables:
```
