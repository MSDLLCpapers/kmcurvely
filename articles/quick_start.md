# Interactive KM plots

``` r

library(dplyr)
library(survival)
library(ggplot2)
library(plotly)
library(r2rtf)
library(kmcurvely)
```

``` r

kmcurvely(
  meta = meta_tte_example(),
  population = "apat",
  observation = "efficacy_population",
  endpoint = "pfs;os",
  subgroup = "male;female",
  time_unit = c("days", "weeks", "months", "years")
)
```
