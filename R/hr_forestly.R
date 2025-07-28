#' @examples
#' prepare_hr_forestly(
#'   meta = meta_tte_example_new(),
#'   population = "apat",
#'   observation = "efficacy_population",
#'   endpoint = "pfs;os",
#'   subgroup = "male;female",
#'   arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
#' ) |>
#'   format_hr_forestly(
#'     hr_range = c(0, 3)
#'   ) |>
#'   hr_forestly(
#'     width = 1400,
#'     max_page = NULL
#'   )
#' Placeholder
hr_forestly <- function(outdata,
                        time_unit = c("days", "weeks", "months", "years"),
                        x_label = NULL,
                        width = 1400,
                        max_page = NULL) {
  time_unit <- match.arg(time_unit, choices = time_unit)
  if (is.null(x_label)) {
    x_label <- paste("Time in", time_unit)
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
      details = function(index) {
        t_row <- outdata$tbl$subgroup[index]
        t_endpoint <- outdata$tbl$endpoint[index]

        t_details <- subset(
          outdata$km_data,
          (toupper(outdata$km_data$subgroup) %in% toupper(t_row)) &
            (toupper(outdata$km_data$endpoint) %in% toupper(t_endpoint))
        )
        cnr_details <- subset(t_details, n.censor >= 1)

        # at risk
        # filter the endpoint and subgroup
        # tbl_surv_sub <- tbl_surv %>% filter(endpoint == endpt_label, subgroup == subgrp_label)

        # range of the x-axis
        xlimit <- switch(time_unit,
          "years"  = ceiling(max(t_details$time)),
          "months" = round((max(t_details$time) + 3), -1),
          "weeks"  = round((max(t_details$time) + 20) / 20) * 20,
          "days"   = round((max(t_details$time) + 60) / 60) * 60
        )

        # intervals to break the x-axis
        break_x_by <- switch(time_unit,
          "years"  = max(ceiling(max(t_details$time) / 7), 1),
          "months" = max(round((max(t_details$time) / 7) / 3) * 3, 3),
          "weeks"  = max(round((max(t_details$time) / 7) / 20) * 20, 20),
          "days"   = max(round((max(t_details$time) / 7) / 60) * 60, 60)
        )

        # time points where number of risk is reported
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
          group_by(endpoint, subgroup, group) %>%
          filter(row_number() %in% (selected_time + 1)) %>%
          tidyr::pivot_wider(names_from = time, values_from = n_at_risk)

        # Create kmplot
        km_plot <- ggplot2::ggplot() +
          ggplot2::geom_step(
            data = t_details,
            aes(
              x = .data$time, y = .data$surv, group = .data$strata,
              colour = .data$strata,
              text = .data$text
            ),
            direction = "hv",
            size = 0.7
          ) +
          ggplot2::scale_color_manual(values = outdata$color) +
          ggplot2::ylab(t_endpoint) +
          ggplot2::xlab(x_label) +
          ggplot2::theme_bw()

        htmltools::div(
          htmltools::tagList(
            ggplotly(km_plot, tooltip = c("text"), dynamicTicks = TRUE) %>%
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
                legend = list(orientation = "h", y = -0.2)
                # title = list(text = "Kaplan-Meier curves")
              ),
            # at-risk table
            "Number of participants at risk",
            DT::datatable(
              tbl_at_risk_select,
              options = list(columnDefs = list(list(
                visible = FALSE,
                targets = c(0)
              ))),
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
