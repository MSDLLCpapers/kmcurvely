# Format the outdata for the HR forest plot, along with its KM plots

Format the outdata for the HR forest plot, along with its KM plots

## Usage

``` r
format_hr_forestly(
  outdata,
  display = c("n", "event", "fig_hr"),
  digits_hr = 2,
  width_subgroup = 50,
  width_fig = 360,
  width_n = 40,
  width_event = 40,
  width_hr = 40,
  digits_surv = 2,
  footer_space = 150,
  hr_range = NULL,
  color = NULL,
  hr_label = "Treatment <- Favor -> Placebo"
)
```

## Arguments

- outdata:

  An `outdata` object created by
  [`prepare_hr_forestly()`](https://merck.github.io/kmcurvely/reference/prepare_hr_forestly.md).

- display:

  A character vector specifying which columns to display in the forest
  plot. Options include "n", "event", "hr", and "fig_hr". Default is
  c("n", "event", "fig_hr").

  - `n`: Number of participants in a comparison.

  - `event`: Number of events in a comparison.

  - `hr`: Hazard ratio estimate in a comparison.

  - `fig_hr`: Hazard ratio figure in a comparison.

- digits_hr:

  A numeric value specifying the number of digits to display for hazard
  ratios and confidence intervals. Default is 2.

- width_subgroup:

  A numeric value specifying the width of the subgroup column in pixels.

- width_fig:

  A numeric value specifying the width of the hazard ratio figure column
  in pixels.

- width_n:

  A numeric value specifying the width of the "n" column in pixels.

- width_event:

  A numeric value specifying the width of the "event" column in pixels.

- width_hr:

  A numeric value specifying the width of the hazard ratio column in
  pixels.

- digits_surv:

  A numeric value specifying the number of digits to display for the
  survival probability in KM curves.

- footer_space:

  A numeric value specifying the space for the footer in pixels.

- hr_range:

  A numeric vector of lower and upper limit of x-axis for the hazard
  ratio figure.

- color:

  A character vector of colors to use for the hazard ratio figures.
  Defaulte value supports up to 4 groups.

- hr_label:

  A character string specifying the label for the hazard ratio axis.

## Value

An `outdata` object.

## Examples

``` r
prepare_hr_forestly(
  meta = meta_tte_example_new(),
  population = "apat",
  observation = "efficacy_population",
  endpoint = "pfs;os",
  subgroup = "gender",
  km_curves = "female;male",
  arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
) |>
  format_hr_forestly()
#> List of 23
#>  $ meta                   :List of 7
#>  $ population             : chr "apat"
#>  $ observation            : chr "efficacy_population"
#>  $ parameter              : NULL
#>  $ n                      :'data.frame': 4 obs. of  5 variables:
#>  $ order                  : NULL
#>  $ group                  : chr [1:3] "Placebo" "Xanomeline Low Dose" "Xanomeline High Dose"
#>  $ reference_group        : num 1
#>  $ endpoint               : chr [1:2] "pfs" "os"
#>  $ subgroup               : chr "gender"
#>  $ kmcurves               : chr [1:2] "Female" "Male"
#>  $ arm_comparison         : chr [1:2] "Xanomeline Low Dose vs. Placebo" "Xanomeline High Dose vs. Placebo"
#>  $ event                  :'data.frame': 4 obs. of  5 variables:
#>  $ hr_est                 :'data.frame': 4 obs. of  5 variables:
#>  $ hr_ci_lower            :'data.frame': 4 obs. of  5 variables:
#>  $ hr_ci_upper            :'data.frame': 4 obs. of  5 variables:
#>  $ km_data                :'data.frame': 1416 obs. of  12 variables:
#>  $ tbl                    :'data.frame': 4 obs. of  14 variables:
#>  $ reactable_columns      :List of 14
#>  $ reactable_columns_group:List of 2
#>  $ display                : chr [1:3] "n" "event" "fig_hr"
#>  $ fig_hr_color           : chr [1:2] "#6ECEB2" "#00857C"
#>  $ color                  : chr [1:3] "#66203A" "#6ECEB2" "#00857C"
```
