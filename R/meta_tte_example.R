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
#' meta_data <- meta_tte_example()
#'
#' # Use the returned meta_data in subsequent analysis functions.


meta_tte_example <- function() {
  # --------------------------------
  #     data cleaning
  # --------------------------------
  # read adsl dataset
  adsl <- kmcurvely_adsl
  adsl$TRT01P <- factor(
    adsl$TRT01P,
    levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"),
    labels = c("Placebo", "Low Dose", "High Dose"))

  # read adtte dataset
  adtte<-kmcurvely_adtte

  adtte <- adtte %>% dplyr::rename(TRT01P = TRTP)
  adtte$TRT01P <- factor(adtte$TRT01P,
    levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"),
    labels = c("Placebo", "Low Dose", "High Dose"))
  # manually add a new endpoint, i.e., PFS
  adtte <- adtte %>%
    dplyr::union_all(adtte %>% mutate(PARAMCD = "PFS", AVAL = AVAL - 3))

  # define analysis plan
  plan <- metalite::plan(
    analysis = "interactive_km_curve", population = "apat",
    observation = "efficacy_population", parameter = "pfs;ttde;male;female"
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
      var = c("USUBJID", "TRT01P", "SEX")
    ) |>
    metalite::define_observation(
      name = "efficacy_population",
      group = "TRT01P",
      subset = SAFFL == "Y",
      var = c("USUBJID", "TRT01P", "SEX", "PARAMCD", "AVAL", "CNSR"),
      label = "Efficacy population"
    ) |>
    metalite::define_parameter(
      name = "pfs",
      subset = PARAMCD == "PFS",
      label = "PFS"
    ) |>
    metalite::define_parameter(
      name = "ttde",
      subset = PARAMCD == "TTDE",
      label = "TTDE"
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
    metalite::define_analysis(
      name = "interactive_km_curve",
      title = "KM curves of PFS",
      label = "km curve"
    ) |>
    metalite::meta_build()
}
