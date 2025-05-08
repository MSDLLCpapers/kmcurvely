# Read in ADTTE dataset and create .rda file
kmcurvely_adtte <- haven::read_xpt("data-raw/adtte.xpt")
usethis::use_data(kmcurvely_adtte, overwrite = TRUE)
# Read in ADSL dataset and create .rda file
kmcurvely_adsl <- haven::read_xpt("data-raw/adsl.xpt")
usethis::use_data(kmcurvely_adsl, overwrite = TRUE)