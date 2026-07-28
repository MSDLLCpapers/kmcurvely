# Generate Static KM plot by endpoints and subgroups and save as RTF/PNG files

Generate Static KM plot by endpoints and subgroups and save as RTF/PNG
files

## Usage

``` r
kmcurvely_static(
  surv_shared,
  tbl_at_risk,
  x_label = "Time in Weeks",
  y_label = "Survival rate",
  color,
  xlimit,
  break_x_by,
  population,
  legned_pos = c(0.1, 0.12),
  folder_path = tempdir(),
  zipname = "all_static_plots.zip"
)
```

## Arguments

- surv_shared:

  A data frame of survival data per endpoint per subgroup

- tbl_at_risk:

  A data frame containing the number of subjects at risk

- x_label:

  The label for the x-axis

- y_label:

  The label for the y-axis

- color:

  The color palette for different treament groups

- xlimit:

  The maximum value for the x-axis

- break_x_by:

  The interval for the x-axis breaks

- population:

  The population specify in kmcurvely.r

- legned_pos:

  The position of the legend

- folder_path:

  The path to save the RTF files

- zipname:

  The name of the zip file for all static KM plots

## Value

A zip file containing RTF files of the static KM plots

## Examples

``` r
key_data_objects <- attr(kmcurvely(), "key_data_objects")
#> Joining with `by = join_by(USUBJID, TRT01P, SEX)`
#> Joining with `by = join_by(time)`
#> Joining with `by = join_by(time)`
#> Joining with `by = join_by(time)`
#> Joining with `by = join_by(time)`
#> Joining with `by = join_by(time)`
#> Joining with `by = join_by(time)`
#> Warning: Ignoring unknown aesthetics: text

kmcurvely_static(
  surv_shared = key_data_objects$tbl_surv,
  tbl_at_risk = key_data_objects$tbl_at_risk,
  x_label = key_data_objects$x_label,
  y_label = key_data_objects$y_label,
  color = key_data_objects$color,
  xlimit = key_data_objects$xlimit,
  break_x_by = key_data_objects$break_x_by,
  population = key_data_objects$population
)
#> [1] "/tmp/RtmpFRR6jC/all_static_plots.zip"
```
