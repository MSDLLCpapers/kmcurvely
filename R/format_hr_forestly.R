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

#' Format the outdata for the HR forest plot, along with its KM plots
#'
#' @param outdata An `outdata` object created by [prepare_hr_forestly()].
#' @param display A character vector specifying which columns to display in the forest plot.
#'   Options include "n", "event", "hr", and "fig_hr".
#'   Default is c("n", "event", "fig_hr").
#'   - `n`: Number of participants in a comparison.
#'   - `event`: Number of events in a comparison.
#'   - `hr`: Hazard ratio estimate in a comparison.
#'   - `fig_hr`: Hazard ratio figure in a comparison.
#' @param digits A numeric value specifying the number of digits to
#'    display for hazard ratios and confidence intervals. Default is 2.
#' @param width_subgroup A numeric value specifying the width of
#'    the subgroup column in pixels.
#' @param width_fig A numeric value specifying the width of
#'    the hazard ratio figure column in pixels.
#' @param width_n A numeric value specifying the width of the "n" column in pixels.
#' @param width_event A numeric value specifying the width of
#'    the "event" column in pixels.
#' @param width_hr A numeric value specifying the width of the hazard ratio column in pixels.
#' @param footer_space A numeric value specifying the space for the footer in pixels.
#' @param hr_range A numeric vector of lower and upper limit of x-axis
#'    for the hazard ratio figure.
#' @param color A character vector of colors to use for the hazard ratio figures.
#'    Defaulte value supports up to 4 groups.
#' @param hr_label A character string specifying the label for the hazard ratio axis.
#'
#' @return An `outdata` object.
#'
#' @export
#'
#' @examples
#' prepare_hr_forestly(
#'   meta = meta_tte_example_new(),
#'   population = "apat",
#'   observation = "efficacy_population",
#'   endpoint = "pfs;os",
#'   subgroup = "gender",
#'   km_curves = "female;male",
#'   arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
#' ) |>
#'   format_hr_forestly()
format_hr_forestly <- function(outdata,
                               display = c("n", "event", "fig_hr"),
                               digits = 2,
                               width_subgroup = 50,
                               width_fig = 360,
                               width_n = 40,
                               width_event = 40,
                               width_hr = 40,
                               footer_space = 150,
                               hr_range = NULL,
                               color = NULL,
                               hr_label = "Treatment <- Favor -> Placebo") {
  display <- match.arg(
    display,
    c("n", "event", "hr", "fig_hr"),
    several.ok = TRUE
  )

  # display_total <- "total" %in% display
  display_hr <- "hr" %in% display

  # n_group_total <- index_total
  n_group <- length(outdata$group)
  n_group1 <- n_group - 1

  name_n <- names(outdata$n)[1:n_group1]
  name_event <- names(outdata$event)[1:n_group1]
  name_hr_est <- names(outdata$hr_est)[1:n_group1]

  # Input checking
  if (is.null(color)) {
    if (n_group <= 2) {
      color <- c("#66203A", "#00857C")
    } else {
      if (n_group1 > 3) stop("Please define color to display groups")
      color <- c("#66203A", rev(c("#00857C", "#6ECEB2", "#BFED33")[1:n_group1]))
    }
  }

  if (length(color) < n_group) {
    stop("Please define more color to display groups")
  }

  # Define table data
  tbl <- data.frame(
    outdata$n,
    outdata$event[1:n_group1],
    round(outdata$hr_est[1:n_group1], digits = digits),
    round(outdata$hr_ci_lower[1:n_group1], digits = digits),
    round(outdata$hr_ci_upper[1:n_group1], digits = digits),
    hr_fig = NA
  )
  col_names <- sapply(2:n_group, function(x) {
    names(tbl)[grep(paste0(".+_", x), names(tbl))]
  }, simplify = FALSE) |> unlist()
  tbl <- tbl[c("endpoint", "subgroup_section", "subgroup", col_names, "hr_fig")]
  rownames(tbl) <- NULL

  # Function to create hazard ratio figure
  tbl_hr <- data.frame(
    outdata$hr_est[1:n_group1],
    outdata$hr_ci_lower[1:n_group1],
    outdata$hr_ci_upper[1:n_group1]
  )

  if (is.null(hr_range)) {
    fig_hr_range <- range(tbl_hr, na.rm = TRUE)
    fig_hr_range[1] <- max(floor(fig_hr_range[1]), 0)
    fig_hr_range[2] <- if ((ceiling(fig_hr_range[2]) - fig_hr_range[2]) > 0.2) ceiling(fig_hr_range[2]) else ceiling(fig_hr_range[2]) + 0.5
  } else {
    if (hr_range[1] > range(tbl_hr, na.rm = TRUE)[1] |
      hr_range[2] < range(tbl_hr, na.rm = TRUE)[2]) {
      warning("There are data points outside the specified range for hazard ratio.")
    }
    fig_hr_range <- hr_range
  }
  fig_hr_color <- color[2:n_group]

  iter <- 1:ncol(outdata$hr_est[1:n_group1]) - 1
  text <- glue::glue("x[{iter}] + '(' + x_lower[{iter}] + ', ' + x_upper[{iter}] + ')'")
  js_hr_fig_cell <- sparkline_point_js(
    tbl = tbl,
    type = "cell",
    x = names(outdata$hr_est)[1:n_group1],
    x_lower = names(outdata$hr_ci_lower)[1:n_group1],
    x_upper = names(outdata$hr_ci_upper)[1:n_group1],
    y = 1:n_group1,
    xlim = fig_hr_range,
    color = fig_hr_color,
    width = width_fig,
    text = text,
    margin = c(0, 20, 0, 0, 0)
  )

  # Function to create Axis
  js_hr_fig_footer <- sparkline_point_js(
    tbl = tbl,
    x = names(outdata$hr_est)[1:n_group1],
    # tbl = data.frame(x = 1),
    # x = "x",
    y = -1,
    type = "footer",
    xlab = hr_label,
    xlim = fig_hr_range,
    height = footer_space,
    width = width_fig,
    # legend = FALSE,
    color = fig_hr_color,
    legend = TRUE,
    legend_label = outdata$group[2:n_group],
    legend_title = "",
    legend_position = -2,
    legend_type = "point",
    margin = c(footer_space - 20, 20, 0, 0, 0)
  )

  # Column Group information ----
  columnGroups <- list()
  name_group <- NULL
  for (i in 1:n_group1) {
    name_group <- c(name_n[i], name_event[i])
    if (display_hr) {
      name_group <- c(name_group, name_hr_est[i])
    }
    columnGroups[[i]] <- reactable::colGroup(
      name = gsub(" vs. ", " <br> vs. ", outdata$arm_comparison[i]),
      html = TRUE,
      columns = name_group
    )
  }

  # Column Definition ----

  # Format variables for group
  col_var <- list(
    subgroup_section = reactable::colDef(
      header = "Subgroup",
      minWidth = width_subgroup,
      align = "right"
    ),
    subgroup = reactable::colDef(
      header = "",
      minWidth = width_subgroup,
      align = "right"
    ),
    endpoint = reactable::colDef(
      header = "Endpoint",
      show = FALSE
    )
  )

  # n column format
  col_n <- lapply(name_n, function(x) {
    reactable::colDef(
      header = "n",
      minWidth = width_n,
      align = "center"
    )
  })
  names(col_n) <- name_n

  # event column format
  col_event <- lapply(name_event, function(x) {
    reactable::colDef(
      header = "#Events",
      minWidth = width_event,
      align = "center"
    )
  })
  names(col_event) <- name_event

  # Define diff column
  hr_name <- c(names(outdata$hr_est)[1:n_group1])
  col_hr <- lapply(
    hr_name,
    function(x) {
      i <- as.numeric(gsub("hr_est_", "", x, fixed = TRUE))
      reactable::colDef(
        header = "HR",
        minWidth = width_hr,
        show = display_hr,
        format = reactable::colFormat(digits = digits)
      )
    }
  )
  names(col_hr) <- hr_name

  # Define ci columns
  ci_name <- c(
    names(outdata$hr_ci_lower)[1:n_group1],
    names(outdata$hr_ci_upper)[1:n_group1]
  )
  col_ci <- lapply(
    ci_name,
    function(x) {
      reactable::colDef(show = FALSE)
    }
  )
  names(col_ci) <- ci_name

  # difference format
  col_hr_fig <- list(hr_fig = reactable::colDef(
    header = paste0("Hazard Ratio + 95% CI <br> vs. ", outdata$group[1]),
    defaultSortOrder = "desc",
    width = ifelse("fig_hr" %in% display, width_fig, 0),
    align = "center",
    sortable = FALSE,
    filterable = FALSE,
    cell = reactable::JS(js_hr_fig_cell),
    footer = reactable::JS(js_hr_fig_footer),
    html = TRUE,
    style = "font-size: 0px; padding: 0px; margin: 0px;",
    footerStyle = "font-size: 0px; padding: 0px; margin: 0px;"
  ))

  # Combine column definition
  columns <- c(
    col_var, col_n, col_event,
    col_hr, col_ci, col_hr_fig
  )

  # column hidden
  columns <- lapply(columns, function(x) {
    if (!"show" %in% names(x)) {
      x$show <- TRUE
    }
    return(x)
  })

  # Create outdata
  outdata$tbl <- tbl
  outdata$reactable_columns <- columns
  outdata$reactable_columns_group <- columnGroups
  outdata$display <- display
  outdata$fig_hr_color <- fig_hr_color
  outdata$color <- color

  outdata
}
