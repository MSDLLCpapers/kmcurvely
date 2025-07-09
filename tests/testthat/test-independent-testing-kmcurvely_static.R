# Example usage of kmcurvely function
meta <- meta_tte_example()

# extract specific population
kmcur_data <- kmcurvely(
  meta = meta,
  population = "apat",
  observation = "efficacy_population",
  endpoint = "pfs",
  subgroup = "male;female"
)

# Extract key data objects from the kmcur_data
key_data_objects <- attr(kmcur_data, "key_data_objects")

out_kmplot <- kmcurvely_static(
  surv_shared = key_data_objects$tbl_surv,
  tbl_at_risk = key_data_objects$tbl_at_risk,
  x_label = key_data_objects$x_label,
  y_label = key_data_objects$y_label,
  color = key_data_objects$color,
  xlimit = key_data_objects$xlimit,
  break_x_by = key_data_objects$break_x_by,
  population = key_data_objects$population
)

plot_zipfile <- basename(out_kmplot)

# Get project root directory (assuming getwd() returns it)
project_dir <- getwd()

# Construct destination directory path
dest_dir <- file.path(project_dir, "tests", "testthat", "_snaps")

# Copy the file
success <- file.copy(from = out_kmplot, to = dest_dir, overwrite = TRUE)
if (!success) {
  stop("Failed to copy the file to the destination directory.")
}


print(paste("The plots are all wrapped in this .zip file: ", plot_zipfile, "at", dest_dir))

# Check if the zip file exists
if (file.exists(dest_dir)) {
  cat(plot_zipfile, "file exists at", dest_dir, "\n")
} else {
  cat(plot_zipfile, "file does not exist at", dest_dir, "\n")
}
