test_that("kmcurvely reports parameters not defined in meta", {
  analysis_plan <- metalite::plan(
    analysis = "km_curvely",
    population = "apat",
    observation = "wk12",
    parameter = "ttde;male"
  )

  adsl <- kmcurvely_adtte
  adtte <- kmcurvely_adtte

  meta <- metalite::meta_adam(observation = adtte, population = adsl) |>
    metalite::define_plan(analysis_plan) |>
    metalite::define_population(
      name = "apat",
      var = c("USUBJID", "SAFFL", "TRTP", "SITEID", "SEX", "RACE", "AGE"),
      group = "TRTP",
      subset = SAFFL == "Y",
      label = "All Participants as Treated"
    ) |>
    metalite::define_observation(
      name = "wk12",
      var = c("USUBJID", "SAFFL", "TRTP", "SEX", "PARAM", "PARAMCD", "AVAL", "CNSR"),
      group = "TRTP",
      subset = SAFFL == "Y",
      label = "Weeks 0 to 12"
    ) |>
    metalite::define_parameter(
      name = "ttde",
      subset = PARAMCD == "TTDE",
      label = "Time to First Dermatologic Event"
    ) |>
    metalite::define_parameter(
      name = "male",
      subset = SEX == "M",
      label = "Male"
    ) |>
    metalite::define_analysis(
      name = "km_curvely",
      title = "KM Plots of All Participants and Subgroups",
      label = "KM curves"
    ) |>
    metalite::meta_build()

  expect_error(
    kmcurvely(
      meta = meta,
      population = "apat",
      observation = "wk12",
      endpoint = "ttde;undefined_endpoint",
      subgroup = "male"
    ),
    paste(
      "Endpoint(s) not defined in `meta`: undefined_endpoint",
      "Define each value with `metalite::define_parameter()` before calling `kmcurvely()`.",
      sep = "\n"
    ),
    fixed = TRUE
  )

  expect_error(
    kmcurvely(
      meta = meta,
      population = "apat",
      observation = "wk12",
      endpoint = "ttde",
      subgroup = "male;undefined_subgroup"
    ),
    paste(
      "Subgroup(s) not defined in `meta`: undefined_subgroup",
      "Define each value with `metalite::define_parameter()` before calling `kmcurvely()`.",
      sep = "\n"
    ),
    fixed = TRUE
  )

  expect_error(
    kmcurvely(
      meta = meta,
      population = "apat",
      observation = "wk12",
      endpoint = "undefined_endpoint",
      subgroup = "undefined_subgroup"
    ),
    paste(
      "Endpoint(s) not defined in `meta`: undefined_endpoint",
      "Subgroup(s) not defined in `meta`: undefined_subgroup",
      "Define each value with `metalite::define_parameter()` before calling `kmcurvely()`.",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
