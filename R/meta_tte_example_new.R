# Copyright (c) 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
# All rights reserved.
#
# This file is part of the kmcurvely program.
#
# kmcurvely is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Example Function for Time-to-Event Analysis
#'
#' This function prepares the metadata required for conducting time-to-event
#' analysis using Kaplan-Meier curves in the context of a clinical trial.
#' It cleans and processes the datasets, defines the analysis plan, and
#' structures the metadata necessary for further analysis.
#'
#' @return A list containing metadata structured for time-to-event analysis.
#'
#' @importFrom dplyr union_all rename mutate %>%
#' @importFrom metalite meta_adam define_plan define_population define_observation define_parameter define_analysis meta_build
#' @importFrom utils data
#' @export
#'
#' @examples
#' meta_data <- meta_tte_example_new()
#' # Use the returned meta_data in subsequent analysis functions.
meta_tte_example_new <- function() {
  # --------------------------------
  #     data cleaning
  # --------------------------------
  # read adsl dataset
  adsl <- kmcurvely_adsl

  # read adtte dataset
  adtte <- kmcurvely_adtte
  set.seed(2025)

  # manually add a new endpoint, i.e., PFS
  adtte_single_endpt <- adtte |> dplyr::filter(PARAMCD == unique(adtte$PARAMCD)[1])
  adtte_pfs <- adtte_single_endpt
  adtte_pfs$fail_time <- -100
  adtte_pfs$fail_time[which(adtte_pfs$TRTP == "Placebo")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Placebo") |> nrow(), log(2) / 6)
  adtte_pfs$fail_time[which(adtte_pfs$TRTP == "Xanomeline Low Dose")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Xanomeline Low Dose") |> nrow(), log(2) / (6 / 0.8))
  adtte_pfs$fail_time[which(adtte_pfs$TRTP == "Xanomeline High Dose")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Xanomeline High Dose") |> nrow(), log(2) / (6 / 0.6))

  adtte_pfs <- adtte_pfs |>
    dplyr::select(-c(AVAL, CNSR)) |>
    dplyr::mutate(
      PARAMCD = "PFS",
      enroll_time = simtrial::rpwexp_enroll(nrow(adtte_single_endpt), enroll_rate = data.frame(duration = 6, rate = nrow(adtte_single_endpt) / 6)),
      dropout_time = rexp(nrow(adtte_single_endpt), -log(0.98) / 12),
      cte = pmin(dropout_time, fail_time) + enroll_time,
      fail = (fail_time <= dropout_time) * 1,
      tte = pmin(cte, 36) - enroll_time,
      event = fail * (cte <= 36)
    ) |>
    dplyr::rename(AVAL = tte) |>
    dplyr::mutate(CNSR = 1 - event) |>
    dplyr::select(-c(enroll_time, dropout_time, cte, fail, event))

  # manually add a new endpoint, i.e., OS
  adtte_single_endpt <- adtte |> dplyr::filter(PARAMCD == unique(adtte$PARAMCD)[1])
  adtte_os <- adtte_single_endpt
  adtte_os$fail_time <- -100
  adtte_os$fail_time[which(adtte_os$TRTP == "Placebo")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Placebo") |> nrow(), log(2) / 20)
  adtte_os$fail_time[which(adtte_os$TRTP == "Xanomeline Low Dose")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Xanomeline Low Dose") |> nrow(), log(2) / (20 / 0.8))
  adtte_os$fail_time[which(adtte_os$TRTP == "Xanomeline High Dose")] <- rexp(adtte_single_endpt |> dplyr::filter(TRTP == "Xanomeline High Dose") |> nrow(), log(2) / (20 / 0.6))

  adtte_os <- adtte_os |>
    dplyr::select(-c(AVAL, CNSR)) |>
    dplyr::mutate(
      PARAMCD = "OS",
      enroll_time = simtrial::rpwexp_enroll(nrow(adtte_single_endpt), enroll_rate = data.frame(duration = 6, rate = nrow(adtte_single_endpt) / 6)),
      dropout_time = rexp(nrow(adtte_single_endpt), -log(0.98) / 12),
      cte = pmin(dropout_time, fail_time) + enroll_time,
      fail = (fail_time <= dropout_time) * 1,
      tte = pmin(cte, 36) - enroll_time,
      event = fail * (cte <= 36)
    ) |>
    dplyr::rename(AVAL = tte) |>
    dplyr::mutate(CNSR = 1 - event) |>
    dplyr::select(-c(enroll_time, dropout_time, cte, fail, event))

  adtte <- rbind(adtte_pfs, adtte_os)


  # define analysis plan
  plan <- metalite::plan(
    analysis = "interactive_km_curve", population = "apat",
    observation = "efficacy_population", parameter = "pfs;os;male;female"
  ) |>
    metalite::add_plan(
      analysis = "hr_forestly", population = "apat",
      observation = "efficacy_population", parameter = "pfs;os;male;female"
    )

  # define metadata
  meta <- metalite::meta_adam(
    population = adsl,
    observation = adtte
  ) |>
    metalite::define_plan(plan = plan) |>
    metalite::define_population(
      name = "apat",
      group = "TRT01P",
      subset = quote(EFFFL == "Y"),
      var = c("USUBJID", "TRT01P", "SEX", "AGEGR1")
    ) |>
    metalite::define_observation(
      name = "efficacy_population",
      group = "TRTP",
      subset = SAFFL == "Y",
      var = c("USUBJID", "TRTP", "SEX", "PARAMCD", "AVAL", "CNSR"),
      label = "Efficacy population"
    ) |>
    metalite::define_parameter(
      name = "pfs",
      subset = PARAMCD == "PFS",
      label = "PFS"
    ) |>
    metalite::define_parameter(
      name = "os",
      subset = PARAMCD == "OS",
      label = "Overall Survival"
    ) |>
    metalite::define_parameter(
      name = "male",
      subset = SEX == "M",
      label = "Male"
    ) |>
    metalite::define_parameter(
      name = "female",
      subset = SEX == "F",
      label = "Female"
    ) |>
    metalite::define_parameter(
      name = "age<65",
      subset = AGEGR1 == "<65",
      label = "Age < 65"
    ) |>
    metalite::define_parameter(
      name = "age65-80",
      subset = AGEGR1 == "65-80",
      label = "Age 65-80"
    ) |>
    metalite::define_parameter(
      name = "age>80",
      subset = AGEGR1 == ">80",
      label = "Age > 80"
    ) |>
    metalite::define_analysis(
      name = "interactive_km_curve",
      title = "KM curves",
      label = "km curve"
    ) |>
    metalite::define_analysis(
      name = "hr_forestly",
      title = "HR forest plot",
      label = "HR forest plot"
    ) |>
    metalite::meta_build()
}
