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
#' Generate Static KM plot by endpoints and subgroups and save as RTF/PNG files
#'
#' @name kmcurvely_static
#' @param surv_shared A data frame of survival data per endpoint per subgroup
#' @param tbl_at_risk A data frame containing the number of subjects at risk
#' @param x_label The label for the x-axis
#' @param y_label The label for the y-axis
#' @param color The color palette for different treament groups
#' @param xlimit The maximum value for the x-axis
#' @param break_x_by The interval for the x-axis breaks
#' @param population The population specify in kmcurvely.r
#' @param legned_pos  The position of the legend
#' @param folder_path The path to save the RTF files
#' @param zipname The name of the zip file for all static KM plots
#'
#' @return A zip file containing RTF files of the static KM plots
#' @import ggplot2
#' @import r2rtf
#' @export
#'
#' @examples
#' key_data_objects <- attr(kmcurvely(), "key_data_objects")
#'
#' kmcurvely_static(
#'   surv_shared = key_data_objects$tbl_surv,
#'   tbl_at_risk = key_data_objects$tbl_at_risk,
#'   x_label = key_data_objects$x_label,
#'   y_label = key_data_objects$y_label,
#'   color = key_data_objects$color,
#'   xlimit = key_data_objects$xlimit,
#'   break_x_by = key_data_objects$break_x_by,
#'   population = key_data_objects$population
#' )



kmcurvely_static <- function( surv_shared,
                              tbl_at_risk,
                              x_label = "Time in Weeks",
                              y_label = "Survival rate",
                              color,
                              xlimit,
                              break_x_by,
                              population,
                              legned_pos= c(0.1, 0.12),
                              folder_path=tempdir(),
                              zipname="all_static_plots.zip") {

  unique_shared_ids <- unique(surv_shared$shared_id)

  # Create an empty list to store the plots & tables
  plots <- list()
  par<- meta_tte_example()[["parameter"]]
  # Iterate over unique shared_id values
  for (shared_id in unique_shared_ids) {

    # Subset the data for the current shared_id
    subset_data <- filter(surv_shared, shared_id == .env$shared_id)
    #print(subset_data)
    # shared_id<-"PFS-ALL"

    endpoint_ <- tolower(unlist(strsplit(shared_id, "-"))[1])
    subgroup_ <- tolower(unlist(strsplit(shared_id, "-"))[2])

    # Extracting and reshaping the data
    extracted_data <- lapply(par, function(x) {
      data.frame(
        name = tolower(unlist(x$name)),
        label = unlist(x$label),
        stringsAsFactors = FALSE
      )
    })

    df <- do.call(rbind, extracted_data)


    endpoint_idx <- grep(endpoint_, df$name)
    if (length(endpoint_idx) > 0) {
      title_end <- df$label[endpoint_idx][1]
    }

    subgroup_idx <- grep(subgroup_, df$name)
    if (length(subgroup_idx) > 0) {
      title_sub <- df$label[subgroup_idx][1]
    } else{
      title_sub <- NULL
    }

    # Define the title of the plot

    title1 <- paste0("Kaplan-Meier Plot of ", title_end)
    title2 <- ifelse (is.null(title_sub), "", paste0(" (", title_sub, ")"))


    # Population Title
    title3 <- paste0( "(", toupper(population), " Population)")

    if(title2 != "") {
      plot_title <- paste(title1, title2, sep = "\n")
    } else {
      plot_title <- title1
    }

    legend_title <- "||| Censored"

    # filter missing data for plotting
    kmplot1 <-
      ggplot2::ggplot(data =  dplyr::filter(subset_data, !is.na(surv)),
                      ggplot2::aes(group = group, colour = group)) +
      ggplot2::geom_step(ggplot2::aes(x = time,
                                      y = surv,
                                      linetype = group,
                                      colour = group),
                         direction = "hv",
                         linewidth = 1,
                         show.legend = TRUE) +
      ggplot2::geom_point(data = subset(subset_data, n.censor >= 1),
                          ggplot2::aes(x = time, y = surv),
                          shape = 3,
                          size = 1,
                          show.legend = FALSE) +
      # Set up the theme for the plot
      ggplot2::theme_bw()+
      ggplot2::theme(
        legend.position = legned_pos,
        legend.title.align = 0.5,
        legend.background=element_blank(),
        panel.grid = ggplot2::element_line(color = "grey85",
                                           linewidth = 0.5,
                                           linetype = 1),
        axis.title.x = ggplot2::element_text(colour = "black"),
        # axis.title.y = ggplot2::element_text(colour = "black"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))  # Adjust the right margin of the y-axis label
      ) +
      # Define the legend for the plot
      ggplot2::guides(
        group = ggplot2::guide_legend(title = legend_title),
        colour = ggplot2::guide_legend(title = legend_title),
        linetype = ggplot2::guide_legend(title = legend_title),
        x = guide_axis(minor.ticks = TRUE),
        y = guide_axis(minor.ticks = TRUE)
      ) +
      ggplot2::scale_color_manual(values = color) +
      ggplot2::xlab(tools::toTitleCase(x_label)) +
      ggplot2::ylab(tools::toTitleCase(y_label)) +
      ggplot2::coord_cartesian(clip = "off") +
      # Control x-axis and y-axis breaks
      ggplot2::scale_x_continuous(breaks = seq(0, xlimit, by = break_x_by),
                                  limits = c(0, xlimit),
                                  expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1),
                                  limits = c(0, 1),
                                  expand = c(0.02, 0))


    tbl_at_risk_plot<- tbl_at_risk %>%
      group_by(endpoint, subgroup, group) %>%
      filter(row_number() %in% (seq(0, xlimit, by = break_x_by) + 1)) %>%
      filter(shared_id == .env$shared_id)

    # Create count table to append
    append_table <-
      ggplot2::ggplot(data = tbl_at_risk_plot,
                      ggplot2::aes(x = time,
                                   y = group,
                                   label = as.character(n_at_risk))
      ) +
      ggplot2::geom_text(
        lineheight = 1,
        size = 8 / ggplot2::.pt,
        show.legend = FALSE) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::ggtitle("Number of participants at risk \n") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        text = element_text(size = 8),
        axis.text = element_text(size = 8 * 1),
        plot.caption = element_text(size = 8 * 1),
        plot.title = element_text(size = 8 * 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()
      )+
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 1, unit = "lines"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(colour = "black")
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, xlimit, by = break_x_by),
        limits = c(0, xlimit),
        expand = c(0.05, 0)
      )

    kmplot_table <-
      patchwork::wrap_plots(kmplot1, append_table,
                            ncol = 1,
                            heights = c(8.5, 1.5))


    # Store the plot in the list
    plots[[shared_id]] <- kmplot_table

    # Save the plot as a PNG file
    png_filename <- file.path(folder_path,paste0("plot_", shared_id, ".png"))
    ggplot2::ggsave(png_filename, plot = kmplot_table, width = 8, height = 5)


    # Define the path of figure
    filename <- png_filename

    filename %>%
      r2rtf::rtf_read_figure() %>% # read PNG files from the file path
      r2rtf::rtf_title(plot_title, title3) %>% # add title or subtitle
      r2rtf::rtf_footnote("") %>% # add footnote
      r2rtf::rtf_source("") %>% # add data source
      r2rtf::rtf_figure(fig_width = 6, fig_height = 4) %>% # default setting of page and figure
      r2rtf::rtf_encode(doc_type = "figure") %>% # encode rtf as figure
      r2rtf::write_rtf(file = paste0(folder_path,"/","plot_", shared_id, ".rtf")) # write RTF to a file

    rtf_files <- list.files(folder_path, full.names = TRUE, pattern = "\\.rtf$")
    zip_path <- file.path(folder_path, zipname)
    # Ensure no `./` prefix in the filenames
    zip::zip(zip_path, files = rtf_files, mode = "cherry-pick")

  }

  return(zip_path)
}