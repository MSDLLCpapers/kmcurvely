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

#' Display the HR forest plot, along with its KM plots
#'
#' @param outdata outdata An `outdata` object created by [format_hr_forestly()].
#' @param time_unit A string specifying the unit of time in x-axis of the KM plot.
#'    must be one of "days", "weeks", "months", or "years". Default is "days".
#' @param width A numeric value of width of the entire plot in pixels. Default is 1000.
#' @param height_kmoplot A numeric value of height of the KM plot in pixels. Default is 400.
#' @param height_at_risk A numeric value of height of the at-risk table in pixels.
#'    Default is 200.
#' @param max_page A numeric value of maximum number of subgroups to display per page.
#' @return An HR forest plot with subgroup KM plot saved as a `shiny.tag.list` object.
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
#' ) |>
#'   format_hr_forestly() |>
#'   hr_forestly()
hr_forestly <- function(outdata,
                        time_unit = c("days", "weeks", "months", "years"),
                        width = 1000,
                        height_kmoplot = 400,
                        height_at_risk = 200,
                        max_page = NULL) {
  time_unit <- match.arg(time_unit, choices = time_unit)
  x_label <- paste("Time in", time_unit)

  if (is.null(max_page)) {
    max_page <- if (length(unique(outdata$tbl$subgroup)) <= 100) c(10, 25, 50, 100) else c(10, 25, 50, 100, ceiling(length(unique(outdata$tbl$subgroup)) / 100) * 100)
  } else {
    max_page <- if (max_page <= 100) c(10, 25, 50, 100) else c(10, 25, 50, 100, max_page)
  }

  # Create SharedData for forest plot data
  tbl <- crosstalk::SharedData$new(outdata$tbl)

  # Filter of Endpoint
  default_endpoint <- as.character(unique(outdata$tbl$endpoint)[1])
  endpoint_id <- paste0("filter_default_endpoint_", uuid::UUIDgenerate(), "|", default_endpoint)
  filter_endpoint <- crosstalk::filter_select(
    id = endpoint_id,
    label = "Endpoint",
    sharedData = tbl,
    group = ~endpoint,
    multiple = FALSE
  )

  p_reactable <-
    reactable::reactable(
      data = tbl,
      columns = outdata$reactable_columns,
      columnGroups = outdata$reactable_columns_group,
      width = width,
      pageSizeOptions = max_page,
      details = function(index) {
        t_row <- outdata$tbl$subgroup[index]
        t_endpoint <- outdata$tbl$endpoint[index]

        t_details <- subset(
          outdata$km_data,
          (toupper(outdata$km_data$subgroup) %in% toupper(t_row)) &
            (toupper(outdata$km_data$endpoint) %in% toupper(t_endpoint))
        )
        cnr_details <- subset(t_details, n.censor >= 1)

        t_end <- subset(
          outdata$km_data,
          (toupper(outdata$km_data$endpoint) %in% toupper(t_endpoint))
        )
        xlimit <- ceiling(max(t_end$time))
        break_x_by <- max(ceiling(max(t_end$time) / 7), 1)
        break_x <- seq(0, xlimit, 1)
        selected_time <- seq(0, xlimit, break_x_by)

        # a table summarizing the number of patients at risk
        tbl_at_risk <- lapply(break_x, function(tau) {
          sapply(split(t_details, t_details$strata), function(x) subset(x, time > tau)[1, "n.risk"])
        })

        tbl_at_risk <- do.call(cbind, tbl_at_risk)
        tbl_at_risk[is.na(tbl_at_risk)] <- 0

        # a wide table with n-group rows and n-time columns
        tbl_at_risk <- tibble::as_tibble(tbl_at_risk) %>%
          dplyr::mutate(
            group = row.names(tbl_at_risk),
            endpoint = t_endpoint,
            subgroup = t_row
          )

        # a wide table with n-group rows and n-selected time columns
        tbl_at_risk_new_wide <- tbl_at_risk[, colnames(tbl_at_risk) %in% c("endpoint", "subgroup", "group", paste0("V", selected_time + 1))] %>%
          select("endpoint", "subgroup", "group", everything())

        colnames(tbl_at_risk_new_wide) <- c("endpoint", "subgroup", "group", selected_time)

        # a long table with each row for the survival prob at 1 endpoint 1 group 1 time point
        tbl_at_risk <- tbl_at_risk %>%
          tidyr::pivot_longer(
            cols = starts_with("V"),
            names_to = "time",
            values_to = "n_at_risk"
          ) %>%
          left_join(tibble::tibble(time = paste0("V", 1:length(break_x)), time_new = break_x)) %>%
          select(-time) %>%
          dplyr::rename(time = time_new)

        tbl_at_risk_select <- tbl_at_risk %>%
          dplyr::group_by(endpoint, subgroup, group) %>%
          dplyr::filter(row_number() %in% (selected_time + 1)) %>%
          tidyr::pivot_wider(names_from = time, values_from = n_at_risk)

        # Create kmplot
        km_plot <- ggplot2::ggplot() +
          ggplot2::geom_step(
            data = t_details,
            aes(
              x = .data$time,
              y = .data$surv,
              group = .data$strata,
              colour = .data$strata,
              text = .data$text
            ),
            direction = "hv",
            size = 0.7
          ) +
          ggplot2::scale_color_manual(values = outdata$color) +
          ggplot2::ylab(t_endpoint) +
          ggplot2::xlab(x_label) +
          ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1)) +
          ggplot2::scale_x_continuous(breaks = selected_time) +
          # ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = as.character(seq(0, 1, 0.1))) +
          # ggplot2::scale_x_continuous(limits = c(0, xlimit), breaks = selected_time, labels = as.character(selected_time)) +
          ggplot2::theme_bw()

        htmltools::div(
          htmltools::tagList(
            ggplotly(km_plot, tooltip = c("text"), dynamicTicks = FALSE, height = height_kmoplot) %>%
              highlight(on = "plotly_click", off = "plotly_doubleclick") %>%
              add_trace(
                data = cnr_details,
                x = ~time,
                y = ~surv,
                color = ~strata,
                colors = outdata$color,
                group = ~strata,
                type = "scatter",
                mode = "markers",
                marker = list(symbol = "cross"),
                hoverinfo = "none",
                showlegend = FALSE,
                opacity = 1
              ) %>%
              layout(
                hovermode = "x unified",
                legend = list(title = "", orientation = "h", y = -0.2)
              ),
            # at-risk table
            "Number of participants at risk",
            DT::datatable(
              tbl_at_risk_select,
              colnames = c(" " = 1),
              height = height_at_risk,
              options = list(
                dom = "t",
                ordering = FALSE,
                columnDefs = list(list(
                  visible = FALSE,
                  targets = c(1, 2)
                ))
              ),
              rownames = FALSE
            )
          )
        )
      }
    )

  p <- suppressWarnings(
    crosstalk::bscols(
      # Width of the select list and reactable
      widths = c(3, 12, 0),
      filter_endpoint,
      p_reactable
    )
  )

  # remove (All)
  brew::brew(
    system.file("js/filter_default.js", package = "kmcurvely"),
    output = file.path(tempdir(), "filter_default.js")
  )

  html_content <-
    htmltools::htmlDependency(
      name = "default-parameter",
      version = "0.1.0",
      src = tempdir(),
      script = c("filter_default.js"),
      all_files = FALSE
    )

  htmltools::browsable(
    htmltools::tagList(
      html_content,
      reactR::html_dependency_react(TRUE),
      forestly:::html_dependency_plotly(TRUE),
      forestly:::html_dependency_react_plotly(TRUE),
      p
    )
  )
}
