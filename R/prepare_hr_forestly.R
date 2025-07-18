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

#' Prepare the dataset for the HR forest plot, alongwith its KM plots
#'
#' @param meta A metadata object that contains the mappings for the variables of interest.
#' @param population A string indicating the population to use (e.g., "apat").
#' @param observation A string indicating the observation population (e.g., "efficacy_population").
#' @param endpoint A semicolon-separated string of endpoints to be analyzed (e.g., "pfs;os").
#' @param subgroup A semicolon-separated string of subgroups to filter by (e.g., "male;female").
#' @param arm_levels A vector of character specifying the levels of arms, starting with the reference arm.
#' @return An metadata with HR per subgroup along with its KM plotting data
#'
#' @export
#'
#' @examples
#' prepare_hr_forestly(
#'   meta = meta_tte_example_new(),
#'   population = "apat",
#'   observation = "efficacy_population",
#'   endpoint = "pfs;os",
#'   subgroup = "male;female",
#'   arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
#' )
prepare_hr_forestly <- function(meta = NULL,
                                population = NULL,
                                observation = NULL,
                                endpoint = NULL,
                                subgroup = NULL,
                                arm_levels = NULL) {
  # obtain population/observation variables
  pop_var <- metalite::collect_adam_mapping(meta, population)$var
  obs_var <- metalite::collect_adam_mapping(meta, observation)$var

  # obtain the population/observation grouping variable
  pop_group <- metalite::collect_adam_mapping(meta, population)$group
  obs_group <- metalite::collect_adam_mapping(meta, observation)$group

  # obtain the population/observation subject identification number
  pop_id <- metalite::collect_adam_mapping(meta, population)$id
  obs_id <- metalite::collect_adam_mapping(meta, observation)$id

  # obtain the population/observation data
  pop <- metalite::collect_population_record(meta, population, var = pop_var)
  obs <- metalite::collect_observation_record(meta, population, observation, parameter = "", var = obs_var)

  # merge population and observation
  rename_lookup <- c(time = "AVAL", event = "CNSR", treatment = obs_group)

  surv_data <- pop |>
    dplyr::left_join(obs) |>
    dplyr::rename(dplyr::any_of(rename_lookup))

  # ----------------------------------------------- #
  #               run cox model or KM plots
  # ----------------------------------------------- #
  endpoint <- unlist(strsplit(endpoint, ";"))
  subgroup <- unlist(strsplit(subgroup, ";"))
  n_endpoint <- length(endpoint)
  n_subgroup <- length(subgroup)
  n_group <- length(unique(obs[[obs_group]]))

  # tables to store sample size, events, HR and KM plotting data
  n_tbl <- NULL
  event_tbl <- NULL
  hr_est_tbl <- NULL
  hr_ci_lower_tbl <- NULL
  hr_ci_upper_tbl <- NULL
  km_tbl <- NULL

  # loop for all possible endpoints
  for (endpt in endpoint) {
    # loop for all possible comparisons between arms
    for (comp in 2:n_group) {
      # loop for all possible subgroups
      for (subgrp in c("all", subgroup)) {
        endpt_filter <- metalite::collect_adam_mapping(meta, endpt)$subset
        endpt_label <- metalite::collect_adam_mapping(meta, endpt)$label
        arm_selected <- arm_levels[c(1, comp)]

        # filter the TTE data given the endpoint, arm comparison and subgroup
        if (subgrp == "all") {
          data_sub <- surv_data |> filter(!!endpt_filter, treatment %in% arm_selected)
          subgroup_label <- "All"
        } else {
          subgrp_filter <- metalite::collect_adam_mapping(meta, subgrp)$subset
          subgroup_label <- metalite::collect_adam_mapping(meta, subgrp)$label
          data_sub <- surv_data |> dplyr::filter(!!endpt_filter & !!subgrp_filter, treatment %in% arm_selected)
        }

        # set reference arms
        data_sub$treatment <- factor(data_sub$treatment, levels = arm_levels[c(1, comp)])

        # run cox model and get HR estimate of the given endpoint, given arm comparison and given subgroup
        cox_model <- survival::coxph(survival::Surv(time, event) ~ treatment, data = data_sub)
        cox_summary <- summary(cox_model)

        #  get the sample size
        n_tbl_new <- data.frame(
          n = cox_summary$n,
          endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. ")
        )
        n_tbl <- rbind(n_tbl, n_tbl_new)
        #  get the events
        event_tbl_new <- data.frame(
          event = cox_summary$nevent,
          endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. ")
        )
        event_tbl <- rbind(event_tbl, event_tbl_new)
        #  get the HR point estimate
        hr_est_tbl_new <- data.frame(
          hr_est = cox_summary$coef[2],
          endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. ")
        )
        hr_est_tbl <- rbind(hr_est_tbl, hr_est_tbl_new)
        #  get the HR lower CI
        hr_ci_lower_tbl_new <- data.frame(
          hr_ci_lower = cox_summary$conf.int[, "lower .95"],
          endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. ")
        )
        hr_ci_lower_tbl <- rbind(hr_ci_lower_tbl, hr_ci_lower_tbl_new)
        #  get the HR upper CI
        hr_ci_upper_tbl_new <- data.frame(
          hr_ci_upper = cox_summary$conf.int[, "upper .95"],
          endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. ")
        )
        hr_ci_upper_tbl <- rbind(hr_ci_upper_tbl, hr_ci_upper_tbl_new)


        # get the data for KM plotting
        km_tbl_new <- survival::survfit(survival::Surv(time, event) ~ treatment, data = data_sub, conf.type = "log-log") |>
          km_extract() |>
          dplyr::mutate(endpoint = endpt_label, subgroup = subgroup_label, arm_comparision = paste0(rev(arm_selected), collapse = " vs. "))
        km_tbl <- rbind(km_tbl, km_tbl_new)
      }
    }
  }

  # prepare outdata
  outdata <- metalite::outdata(
    meta = meta,
    population = population,
    observation = observation,
    parameter = NULL,
    endpoint = endpoint,
    subgroup = subgroup,
    order = NULL,
    group = NULL,
    reference_group = NULL,
    n = n_tbl,
    event = event_tbl,
    hr_est = hr_est_tbl,
    hr_ci_lower = hr_ci_lower_tbl,
    hr_ci_upper = hr_ci_upper_tbl,
    km_data = km_tbl
  )

  return(outdata)
}
