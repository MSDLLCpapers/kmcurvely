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
#' @param subgroup A semicolon-separated string of subgroups to filter by (e.g., "age;gender").
#'    The parameters should be defined at variable levels.
#'    Use "all" to include the result for overall population (e.g., "all;age;gender").
#' @param km_curves A semicolon-separated string of subgroups for which KM curves should be plotted.
#'    The parameters should be defined at value levels for subgroups.
#' @param digits_km_curves A numeric value specifying the number of digits
#'    for the survival probability in KM curves.
#' @param arm_levels A vector of character specifying the levels of arms, starting with the reference arm.
#'
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
#'   subgroup = "age;gender",
#'   km_curves = "female;age65-80",
#'   arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
#' )
prepare_hr_forestly <- function(meta = NULL,
                                population = NULL,
                                observation = NULL,
                                endpoint = NULL,
                                subgroup = NULL,
                                km_curves = NULL,
                                digits_km_curves = 2,
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
  rename_lookup <- c(time = "AVAL", treatment = obs_group)

  # Check arm levels
  if (!all(arm_levels %in% unique(obs[[obs_group]]))) {
    stop("`arm_levels` must be a character vector of treatment arms present in the input data.")
  }
  if (length(arm_levels) < 2) {
    stop("`arm_levels` must contain at least two treatment arms, including the reference arm.")
  }

  surv_data <- pop |>
    dplyr::left_join(obs, by = "USUBJID", suffix = c("", ".pop")) |>
    dplyr::mutate(
      event = 1 - CNSR
    ) |>
    dplyr::select(-dplyr::ends_with(".pop")) |>
    dplyr::rename(dplyr::any_of(rename_lookup)) |>
    dplyr::filter(treatment %in% arm_levels)

  surv_data$treatment <- factor(surv_data$treatment, levels = arm_levels)

  # ----------------------------------------------- #
  #               run cox model or KM plots
  # ----------------------------------------------- #
  endpoint <- unlist(strsplit(endpoint, ";"))
  subgroup <- unlist(strsplit(subgroup, ";"))
  if (!is.null(km_curves)) {
    km_curves <- unlist(strsplit(km_curves, ";"))
  }
  n_endpoint <- length(endpoint)
  n_subgroup <- length(subgroup)
  n_group <- length(unique(surv_data$treatment))
  arm_comparison <- sapply(arm_levels[2:n_group], function(x) {
    paste0(rev(c(arm_levels[1], x)), collapse = " vs. ")
  }, USE.NAMES = FALSE)

  # tables to store sample size, events, HR and KM plotting data
  res <- NULL
  res_name <- c("n_tbl", "event_tbl", "hr_est_tbl", "hr_ci_lower_tbl", "hr_ci_upper_tbl")
  km_tbl <- NULL

  # loop for all possible endpoints
  for (endpt in endpoint) {
    endpt_filter <- metalite::collect_adam_mapping(meta, endpt)$subset
    endpt_label <- metalite::collect_adam_mapping(meta, endpt)$label

    res_new <- lapply(subgroup, function(par_subgrp) {
      if (toupper(par_subgrp) == "ALL") {
        subgroup_var <- "ALL"
        subgroup_section <- "All"
        subgroup_val <- "All"
      } else {
        subgroup_var <- metalite::collect_adam_mapping(meta, par_subgrp)$var
        subgroup_section <- metalite::collect_adam_mapping(meta, par_subgrp)$label
        subgroup_val <- unique(surv_data[[subgroup_var]])
      }

      lapply(subgroup_val, function(subgrp) {
        # filter the TTE data given the endpoint, arm comparison and subgroup
        if (toupper(subgrp) == "ALL") {
          data_sub <- surv_data |> filter(!!endpt_filter)
        } else {
          subgrp_filter <- paste0(subgroup_var, " == '", subgrp, "'")
          data_sub <- surv_data |> dplyr::filter(!!endpt_filter & eval(parse(text = subgrp_filter)))
        }
        subgroup_label <- subgrp

        # run cox model and get HR estimate of the given endpoint, given arm comparison and given subgroup
        cox_model <- lapply(2:n_group, function(x) {
          data_sub_trt <- data_sub |> dplyr::filter(treatment %in% arm_levels[c(1, x)])
          # set reference arms
          data_sub_trt$treatment <- factor(data_sub_trt$treatment, levels = arm_levels[c(1, x)])
          survival::coxph(survival::Surv(time, event) ~ treatment, data = data_sub_trt)
        })

        cox_summary_list <- lapply(cox_model, summary)

        #  get the sample size
        n_tbl_new <- lapply(cox_summary_list, function(x) x$n)
        n_tbl_new <- do.call("cbind", n_tbl_new) |> data.frame(endpt_label, subgroup_label, subgroup_section)
        names(n_tbl_new) <- c(paste0("n_", 2:n_group), "endpoint", "subgroup", "subgroup_section")

        #  get the events
        event_tbl_new <- lapply(cox_summary_list, function(x) x$nevent)
        event_tbl_new <- do.call("cbind", event_tbl_new) |> data.frame(endpt_label, subgroup_label, subgroup_section)
        names(event_tbl_new) <- c(paste0("event_", 2:n_group), "endpoint", "subgroup", "subgroup_section")

        #  get the HR point estimate
        hr_est_tbl_new <- lapply(cox_summary_list, function(x) x$coef[2])
        hr_est_tbl_new <- do.call("cbind", hr_est_tbl_new) |> data.frame(endpt_label, subgroup_label, subgroup_section)
        names(hr_est_tbl_new) <- c(paste0("hr_est_", 2:n_group), "endpoint", "subgroup", "subgroup_section")

        #  get the HR lower CI
        hr_ci_lower_tbl_new <- lapply(cox_summary_list, function(x) x$conf.int[, "lower .95"])
        hr_ci_lower_tbl_new <- do.call("cbind", hr_ci_lower_tbl_new) |> data.frame(endpt_label, subgroup_label, subgroup_section)
        names(hr_ci_lower_tbl_new) <- c(paste0("hr_ci_lower_", 2:n_group), "endpoint", "subgroup", "subgroup_section")

        #  get the HR upper CI
        hr_ci_upper_tbl_new <- lapply(cox_summary_list, function(x) x$conf.int[, "upper .95"])
        hr_ci_upper_tbl_new <- do.call("cbind", hr_ci_upper_tbl_new) |> data.frame(endpt_label, subgroup_label, subgroup_section)
        names(hr_ci_upper_tbl_new) <- c(paste0("hr_ci_upper_", 2:n_group), "endpoint", "subgroup", "subgroup_section")

        list(
          n_tbl = n_tbl_new,
          event_tbl = event_tbl_new,
          hr_est_tbl = hr_est_tbl_new,
          hr_ci_lower_tbl = hr_ci_lower_tbl_new,
          hr_ci_upper_tbl = hr_ci_upper_tbl_new
        )
      })
    })
    res_new <- lapply(c("n_tbl", "event_tbl", "hr_est_tbl", "hr_ci_lower_tbl", "hr_ci_upper_tbl"), function(name) {
      do.call(rbind, lapply(res_new, function(res1) {
        do.call(rbind, lapply(res1, function(x) x[[name]]))
      }))
    })
    names(res_new) <- res_name
    if (!is.null(res)) {
      res <- lapply(res_name, function(x) {
        rbind(res[[x]], res_new[[x]])
      })
    } else {
      res <- res_new
    }
    names(res) <- res_name

    # get the data for KM plotting
    if (!is.null(km_curves)) {
      km_tbl_new <- lapply(km_curves, function(km_curve) {
        if (toupper(km_curve) == "ALL") {
          data_sub <- surv_data |> dplyr::filter(!!endpt_filter)
          km_curve_val <- "All"
        } else {
          km_curve_var <- metalite::collect_adam_mapping(meta, km_curve)$var
          km_curve_filter <- metalite::collect_adam_mapping(meta, km_curve)$subset
          data_sub <- surv_data |> dplyr::filter(!!endpt_filter & !!km_curve_filter)
          km_curve_val <- unique(data_sub[[km_curve_var]])
        }
        survival::survfit(survival::Surv(time, event) ~ treatment, data = data_sub, conf.type = "log-log") |>
          km_extract() |>
          dplyr::mutate(
            surv = round(surv, digits_km_curves),
            endpoint = endpt_label,
            subgroup = km_curve_val,
            text = paste0(endpoint, ": ", surv, "\n", "Number of participants at risk: ", n.risk)
          )
      })
      km_tbl_new <- do.call(rbind, km_tbl_new)
      km_tbl <- rbind(km_tbl, km_tbl_new)
    }
  }
  kmcurves <- as.character(unique(km_tbl$subgroup))

  # prepare outdata
  outdata <- metalite::outdata(
    meta = meta,
    population = population,
    observation = observation,
    parameter = NULL,
    endpoint = endpoint,
    subgroup = subgroup,
    kmcurves = kmcurves,
    order = NULL,
    group = arm_levels,
    reference_group = 1,
    arm_comparison = arm_comparison,
    n = res$n_tbl,
    event = res$event_tbl,
    hr_est = res$hr_est_tbl,
    hr_ci_lower = res$hr_ci_lower_tbl,
    hr_ci_upper = res$hr_ci_upper_tbl,
    km_data = km_tbl
  )

  return(outdata)
}
