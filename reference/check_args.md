# Check Argument Types, Length or Dimension

Check Argument Types, Length or Dimension

## Usage

``` r
check_args(arg, type, length = NULL, dim = NULL)
```

## Arguments

- arg:

  an argument to be checked.

- type:

  a character vector of candidate argument type.

- length:

  a numeric value of argument length or NULL

- dim:

  a numeric vector of argument dimension or NULL.

## Value

Check failure detailed error message

## Details

if `type`, `length` or `dim` is NULL, the corresponding check will not
be executed.

## Specification

The contents of this section are shown in PDF user manual only.

## Examples

``` r
if (FALSE) { # \dontrun{
tbl < -as.data.frame(matrix(1:9, nrow = 3))
check_args(arg = tbl, type = c("data.frame"))

vec <- c("a", "b", "c")
check_args(arg = vec, type = c("character"), length = 3)
} # }
```
