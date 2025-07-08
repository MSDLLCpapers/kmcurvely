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

#' Extract Survival Estimates
#'
#' This function extracts the km estimates from survfit.
#'
#' @section Specification:
#' \if{latex}{
#'
#' Derivation logic in \code{km_extract}:
#'
#'  \itemize{
#'  \item Check if f.survfit argument type is valid.
#'  \item Summarize km estimates.
#'  \item Manipulate and define strata name.
#'  \item Define output data frame with:
#'  time, timepoint from km raw estimates;
#'  n.risk, number at risk at the specific timepoint;
#'  n.event, number of events at the specific timepoint;
#'  n.censor, number of censored patients at the specific timepoint;
#'  surv, survival probability from km estimates at the time point;
#'  upper/lower, upper and lower limit of the km estimates (survival rate);
#'  strata, the strata of the population (from the survfit), if no just use the 'All';
#'  strata_name, name of strata variable to be displayed.
#'  }
#'  }
#'
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @param f.survfit Survival fit object from survfit.
#'
#' @return Data frame of timepoint estimates from km methods.
#' @keywords internal
#'
km_extract <- function(f.survfit) {
  # input type check
  check_args(f.survfit, type = c("survfit"))

  # Summarize KM estimator
  x <- summary(f.survfit, times = c(0, unique(f.survfit$time)), extend = TRUE)

  # Define Strata
  if ("strata" %in% names(x)) {
    s <- data.frame(do.call(rbind, strsplit(as.character(x$strata), ", ")), stringsAsFactors = FALSE)
    strata <- apply(s, 2, function(x) gsub(".*=", "", x))
    strata_name <- apply(s, 2, function(x) gsub("=.*", "", x))
    if (ncol(strata) > 1) {
      strata <- apply(strata, 1, paste, collapse = ", ")
      strata_name <- apply(strata_name, 1, paste, collapse = ", ")
    }

    strata <- factor(strata, levels = unique(strata))
    strata_name <- tools::toTitleCase(strata_name)
  } else {
    strata <- factor("All")
    strata_name <- ""
  }

  # Define output data frame
  data.frame(
    time = x$time,
    n.risk = x$n.risk,
    n.event = x$n.event,
    n.censor = x$n.censor,
    surv = x$surv,
    lower = x$lower,
    upper = x$upper,
    strata = strata,
    strata_name = strata_name
  )
}
