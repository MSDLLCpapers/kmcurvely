test_that("kmcurvely produces expected zip file output", {
  # Load example data
  meta <- meta_tte_example()

  # Run kmcurvely function with specified parameters
  kmcur_data <- kmcurvely(
    meta = meta,
    population = "apat",
    observation = "efficacy_population",
    endpoint = "pfs",
    subgroup = "male;female"
  )

  # Extract key data objects
  key_data_objects <- attr(kmcur_data, "key_data_objects")

  # Generate static KM plot zip file
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

  # Define destination directory for snapshots
  dest_dir <- file.path(getwd(), "tests", "testthat", "_snaps")

  # Ensure destination directory exists
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }

  # Copy the output zip file to the destination directory
  success <- file.copy(from = out_kmplot, to = dest_dir, overwrite = TRUE)
  expect_true(success, info = "Failed to copy the zip file to the destination directory")

  # Check that the zip file exists in the destination directory
  plot_zipfile <- basename(out_kmplot)
  expect_true(file.exists(file.path(dest_dir, plot_zipfile)),
    info = paste(plot_zipfile, "does not exist in", dest_dir)
  )
})
