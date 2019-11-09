#' @title plot_j0j0
#'
#' @description helper function to plot the output of j0j0r::j0j0()
#'
#' @param j0j0_df \code{data.frame} output of j0j0r::j0j0() containing j0j0 elements
#'
#' @param wrap_by \code{charater} the plot will facet wrap by this variable in j0j0_df
#'
#' @param color_by \code{charater} lines will be colored according to this variable in j0j0_df
#'
#' @return \code{ggplot}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang sym !!
#'
#' @export
plot_j0j0 <- function(j0j0_df, wrap_by, color_by) {

  scaleFUN <- function(x) sprintf("%.2e", x)

  j0j0_df %>%
    dplyr::mutate(
      real = Re(j0j0),
      imaginary = Im(j0j0)
    ) %>%
    tidyr::gather("real", "imaginary", key = "component", value = "j0j0") %>%
    dplyr::mutate(component = forcats::fct_rev(as.factor(component)),) %>%
    dplyr::mutate(
      frequency = frequency/1e6,
      element =
        paste0(
          directions,
          ", ",
          component,
          " part"
        )
    ) %>%
    dplyr::group_by(element) %>%
    dplyr::filter(!(all(j0j0 == 0))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      !!color_by := forcats::fct_rev(as.factor(.data[[color_by]])),
      !!wrap_by := as.factor(.data[[wrap_by]])
    ) %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes_string(
        x = "frequency",
        y = "j0j0",
        color = color_by
      )
    ) +
    ggplot2::geom_line(size = 1.5) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::facet_wrap( facets = sym(wrap_by)) +
    ggplot2::xlab("Frequency in MHz") +
    ggplot2::theme(
      legend.position = "top",
      text = ggplot2::element_text(size = 25)
    ) +
    ggplot2::scale_y_continuous(labels = scaleFUN)

}


