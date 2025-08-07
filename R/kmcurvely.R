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

#' Generate interactive KM plot by endpoints and subgroups
#'
#' @name kmcurvely
#' @param meta A metadata object that contains the mappings for the variables of interest.
#' @param population A string indicating the population to use (e.g., "apat").
#' @param observation A string indicating the observation population (e.g., "efficacy_population").
#' @param endpoint A semicolon-separated string of endpoints to be analyzed (e.g., "pfs;ttde").
#' @param subgroup A semicolon-separated string of subgroups to filter by (e.g., "male;female").
#' @param x_label A string for the x-axis label (default is "Time in days/weeks/months/years").
#' @param y_label A string for the y-axis label (default is "Survival rate").
#' @param color A vector of colors for the plot (if NULL, defaults to a preset color palette).
#' @param time_unit A string specifying the unit of time; must be one of "days", "weeks", "months", or "years". Default is "days".
#'
#' @return An interactive KM plot
#'
#' @importFrom ggplot2 ggplot aes geom_step scale_color_manual ylab xlab theme_bw
#' @importFrom survival survfit Surv strata
#' @importFrom dplyr %>% filter left_join mutate rename row_number select starts_with any_of group_by everything union_all
#' @importFrom crosstalk SharedData filter_select
#' @importFrom uuid UUIDgenerate
#' @importFrom htmltools tagList browsable htmlDependency
#' @importFrom DT datatable
#' @importFrom plotly ggplotly highlight add_trace layout
#' @importFrom utils data
#' @import metalite
#' @export
#'

#' @examples
#' # Example usage of kmcurvely function
#'
#' meta <- meta_tte_example()
#' kmcurvely(
#'   meta = meta,
#'   population = "apat",
#'   observation = "efficacy_population",
#'   endpoint = "pfs;ttde",
#'   subgroup = "male;female"
#' )
#'
kmcurvely <- function(meta = meta_tte_example(),
                      population = "apat",
                      observation = "efficacy_population",
                      endpoint = "pfs;ttde",
                      subgroup = "male;female",
                      x_label = NULL,
                      y_label = "Survival rate",
                      color = NULL,
                      time_unit = c("days", "weeks", "months", "years")) {
  time_unit <- match.arg(time_unit, choices = time_unit)
  if (is.null(x_label)) {
    x_label <- paste("Time in", time_unit)
  }

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

  surv_data <- pop %>%
    dplyr::left_join(obs) %>%
    dplyr::mutate(
      event = 1 - CNSR
    ) %>%
    dplyr::rename(dplyr::any_of(rename_lookup))

  # ----------------------------------------------- #
  #   calculate survival fit of all subjects
  # ----------------------------------------------- #
  endpoint <- unlist(strsplit(endpoint, ";"))
  subgroup <- unlist(strsplit(subgroup, ";"))
  n_endpoint <- length(endpoint)
  n_subgroup <- length(subgroup)
  n_group <- length(unique(obs[[obs_group]]))

  tbl_surv <- NULL
  # obtain the survival per endpoint per subgroup
  for (endpt in endpoint) {
    for (subgrp in c("all", subgroup)) {
      endpt_filter <- metalite::collect_adam_mapping(meta, endpt)$subset
      endpt_label <- metalite::collect_adam_mapping(meta, endpt)$label

      if (subgrp == "all") {
        tbl_surv_new <- survival::survfit(survival::Surv(time, event) ~ treatment,
          data = surv_data %>% filter(!!endpt_filter),
          conf.type = "log-log"
        ) |>
          km_extract() |>
          dplyr::mutate(endpoint = endpt_label, subgroup = "All")
      } else {
        subgrp_filter <- metalite::collect_adam_mapping(meta, subgrp)$subset
        subgrp_label <- metalite::collect_adam_mapping(meta, subgrp)$label
        tbl_surv_new <- survival::survfit(survival::Surv(time, event) ~ treatment,
          data = surv_data %>% filter(!!endpt_filter, !!subgrp_filter),
          conf.type = "log-log"
        ) |>
          km_extract() |>
          dplyr::mutate(endpoint = endpt_label, subgroup = subgrp_label)
      }
      tbl_surv <- rbind(tbl_surv, tbl_surv_new)
    }
  }
  tbl_surv <- tbl_surv %>%
    dplyr::rename(group = strata) %>%
    dplyr::mutate(
      surv = round(surv, 2),
      text = paste0(endpoint, ": ", surv, "\n", "Number of participants at risk: ", n.risk),
      shared_id = paste0(endpoint, "-", subgroup)
    )

  # ----------------------------------------------- #
  #    count the number of subjects at risk
  #    by endpoint and subgroup
  # ----------------------------------------------- #
  tbl_at_risk <- NULL

  for (ept in endpoint) {
    for (subgrp in c("all", subgroup)) {
      endpt_label <- metalite::collect_adam_mapping(meta, ept)$label
      subgrp_label <- ifelse(subgrp == "all", "All", metalite::collect_adam_mapping(meta, subgrp)$label)

      # filter the endpoint and subgroup
      tbl_surv_sub <- tbl_surv %>% filter(endpoint == endpt_label, subgroup == subgrp_label)

      # range of the x-axis
      xlimit <- switch(time_unit,
        "years"  = ceiling(max(tbl_surv_sub$time)),
        "months" = round((max(tbl_surv_sub$time) + 3), -1),
        "weeks"  = round((max(tbl_surv_sub$time) + 20) / 20) * 20,
        "days"   = round((max(tbl_surv_sub$time) + 60) / 60) * 60
      )

      # intervals to break the x-axis
      break_x_by <- switch(time_unit,
        "years"  = max(ceiling(max(tbl_surv_sub$time) / 7), 1),
        "months" = max(round((max(tbl_surv_sub$time) / 7) / 3) * 3, 3),
        "weeks"  = max(round((max(tbl_surv_sub$time) / 7) / 20) * 20, 20),
        "days"   = max(round((max(tbl_surv_sub$time) / 7) / 60) * 60, 60)
      )

      # time points where number of risk is reported
      break_x <- seq(0, xlimit, 1)
      selected_time <- seq(0, xlimit, break_x_by)

      # a table summarizing the number of patients at risk
      tbl_at_risk_new <- lapply(break_x, function(tau) {
        sapply(split(tbl_surv_sub, tbl_surv_sub$group), function(x) subset(x, time > tau)[1, "n.risk"])
      })

      tbl_at_risk_new <- do.call(cbind, tbl_at_risk_new)
      tbl_at_risk_new[is.na(tbl_at_risk_new)] <- 0

      # a wide table with n-group rows and n-time columns
      tbl_at_risk_new <- tibble::as_tibble(tbl_at_risk_new) %>%
        dplyr::mutate(
          group = row.names(tbl_at_risk_new),
          endpoint = endpt_label, subgroup = subgrp_label
        )

      # a wide table with n-group rows and n-selected time columns
      tbl_at_risk_new_wide <- tbl_at_risk_new[, colnames(tbl_at_risk_new) %in% c("endpoint", "subgroup", "group", paste0("V", selected_time + 1))] %>%
        select("endpoint", "subgroup", "group", everything())

      colnames(tbl_at_risk_new_wide) <- c("endpoint", "subgroup", "group", selected_time)

      # a long table with each row for the survival prob at 1 endpoint 1 group 1 time point
      tbl_at_risk_new <- tbl_at_risk_new %>%
        tidyr::pivot_longer(
          cols = starts_with("V"),
          names_to = "time",
          values_to = "n_at_risk"
        ) %>%
        left_join(tibble::tibble(time = paste0("V", 1:length(break_x)), time_new = break_x)) %>%
        select(-time) %>%
        dplyr::rename(time = time_new)

      tbl_at_risk <- rbind(tbl_at_risk, tbl_at_risk_new)
    }
  }

  tbl_at_risk <- tbl_at_risk %>% dplyr::mutate(shared_id = paste0(endpoint, "-", subgroup))

  tbl_at_risk_select <- tbl_at_risk %>%
    group_by(endpoint, subgroup, group) %>%
    filter(row_number() %in% (selected_time + 1)) %>%
    tidyr::pivot_wider(names_from = time, values_from = n_at_risk)

  # ----------------------------------------------- #
  #    build shared data
  # ----------------------------------------------- #
  shareddata_id <- uuid::UUIDgenerate()

  surv_shared <- crosstalk::SharedData$new(tbl_surv,
    group = shareddata_id,
    key = ~shared_id
  )

  cnr_shared <- crosstalk::SharedData$new(subset(tbl_surv, n.censor >= 1),
    group = shareddata_id,
    key = ~shared_id
  )

  risk_hover_shared <- crosstalk::SharedData$new(tbl_at_risk,
    group = shareddata_id,
    key = ~shared_id
  )

  risk_table_shared <- crosstalk::SharedData$new(tbl_at_risk_select,
    group = shareddata_id,
    key = ~shared_id
  )

  # ----------------------------------------------- #
  #    build static KM plot
  # ----------------------------------------------- #
  # implement color
  if (is.null(color)) {
    color_pal <- c("#00857C", "#6ECEB2", "#BFED33", "#FFF063", "#0C2340", "#5450E4")
    color <- c("#66203A", rep(color_pal, length.out = n_group - 1))
  } else {
    color <- rep(color, length.out = n_group)
  }

  # build plot
  km_plot <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = surv_shared,
      aes(x = .data$time, y = .data$surv, group = .data$group, colour = .data$group, text = .data$text),
      direction = "hv",
      size = 0.7
    ) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::ylab(y_label) +
    ggplot2::xlab(x_label) +
    ggplot2::theme_bw()

  # ----------------------------------------------- #
  #    build selectors to filter endpoints and subgroup
  # ----------------------------------------------- #
  # Get first parameter name, create random id, and create filter for each
  default_param <- as.character(unique(tbl_surv$endpoint)[1])
  param_random_id <- paste0("filter_default_param_", uuid::UUIDgenerate(), "|", default_param)
  select_endpoint <- crosstalk::filter_select(
    id = param_random_id, sharedData = surv_shared,
    multiple = FALSE, group = ~endpoint, label = "Endpoint"
  )

  default_subgroup <- as.character(unique(tbl_surv$subgroup)[1])
  subgroup_random_id <- paste0("filter_default_subgroup_", uuid::UUIDgenerate(), "|", default_subgroup)
  select_subgroup <- crosstalk::filter_select(
    id = subgroup_random_id, sharedData = surv_shared,
    multiple = FALSE, group = ~subgroup, label = "Subgroup"
  )
  # ----------------------------------------------- #
  #    combine everything
  # ----------------------------------------------- #

  plot <- crosstalk::bscols(
    widths = c(3, 9),
    # 2 sets of filter to select endpoint and subgroup
    list(select_endpoint, select_subgroup),
    # interactive KM curve
    list(
      ggplotly(km_plot, tooltip = c("text"), dynamicTicks = TRUE) %>%
        highlight(on = "plotly_click", off = "plotly_doubleclick") %>%
        add_trace(
          data = cnr_shared,
          x = ~time,
          y = ~surv,
          color = ~group,
          colors = color,
          group_by = ~group,
          type = "scatter",
          mode = "markers",
          marker = list(symbol = "cross"),
          hoverinfo = "none",
          showlegend = FALSE,
          opacity = 1
        ) %>%
        layout(
          hovermode = "x unified",
          legend = list(orientation = "h", y = -0.2),
          title = list(text = "Kaplan-Meier curves")
        ),
      # at-risk table
      "Number of participants at risk",
      DT::datatable(risk_table_shared,
        options = list(columnDefs = list(list(
          visible = FALSE,
          targets = c(0)
        ))),
        rownames = FALSE
      )
    )
  )

  brew::brew(system.file("js/filter_default.js", package = "kmcurvely"),
    output = file.path(tempdir(), "filter_default.js")
  )
  # Extract key data and objects
  key_data_objects <- list(
    tbl_surv = tbl_surv,
    tbl_at_risk = tbl_at_risk,
    color = color,
    x_label = x_label,
    y_label = y_label,
    time_unit = time_unit,
    endpoint = endpoint,
    xlimit = xlimit,
    break_x_by = break_x_by,
    population = population
  )


  zip_path <- kmcurvely_static(
    surv_shared = tbl_surv,
    tbl_at_risk = tbl_at_risk,
    x_label = x_label,
    y_label = y_label,
    color = color,
    xlimit = xlimit,
    break_x_by = break_x_by,
    population = population
  )


  html_content <- htmltools::browsable(htmltools::tagList(
    htmltools::htmlDependency(
      # name and version: doesn't really matter, can be anything
      "default-parameter", "0.1.0",
      # src: path to directory containing the script (e.g., the current directory '.')
      src = tempdir(),
      # script: filename of script to include
      script = "filter_default.js",
      # Exclude all other files in src directory
      all_files = FALSE
    ),
    plot
  ))


  # HTML widget with download link
  download_link <- htmltools::a(
    href = paste0("data:application/zip;base64,", base64enc::base64encode(zip_path)),
    download = "all_static_plots.zip",
    target = "_blank",
    "Download All Plots"
  )

  # Combine plot and download link
  final_plot <- htmltools::browsable(htmltools::tagList(
    html_content,
    download_link
  ))

  # Attach the data objects as an attribute to the final plot
  attr(final_plot, "key_data_objects") <- key_data_objects

  return(final_plot)
}
