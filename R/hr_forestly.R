#' @examples
#' prepare_hr_forestly(
#'   meta = meta_tte_example_new(),
#'   population = "apat",
#'   observation = "efficacy_population",
#'   endpoint = "pfs;os",
#'   subgroup = "male;female",
#'   arm_levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
#' ) |>
#' format_hr_forestly(
#'   hr_range = c(0, 3)
#' ) |>
#' hr_forestly(
#'   width = 1400,
#'   max_page = NULL
#' )
#' Placeholder
hr_forestly <- function(outdata,
                        width = 1400,
                        max_page = NULL) {

  # Create SharedData for forest plot data
  tbl <- crosstalk::SharedData$new(outdata$tbl)

  # Filter of Endpoint
  endpoint_id <- paste0("filter_endpoint_", uuid::UUIDgenerate())
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

        km_plot <- ggplot2::ggplot() +
          ggplot2::geom_step(data = t_details,
                             aes(x = .data$time, y = .data$surv, group = .data$strata,
                                 colour =.data$strata,
                                 text = .data$text
                                 ),
                             direction = "hv",
                             size = 0.7) +
          ggplot2::scale_color_manual(values = outdata$color) +
          ggplot2::ylab(t_endpoint) +
          ggplot2::theme_bw()

        htmltools::div(
          htmltools::tagList(

            ggplotly(km_plot, tooltip = c("text"), dynamicTicks = TRUE) %>%
              highlight(on = "plotly_click", off = "plotly_doubleclick") %>%
              add_trace(data = cnr_details,
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
                        opacity = 1) %>%
              layout(
                hovermode = "x unified",
                legend = list(orientation = "h", y = -0.2)
                # title = list(text = "Kaplan-Meier curves")
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

  htmltools::browsable(
    htmltools::tagList(
      reactR::html_dependency_react(TRUE),
      forestly:::html_dependency_plotly(TRUE),
      forestly:::html_dependency_react_plotly(TRUE),
      p
    )
  )
}