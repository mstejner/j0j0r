#' @title calculate_distribution_data_frame
#'
#' @description function to evaluate velocity distributions and return in a data
#'   frame
#'
#' @param v_perp \code{numeric} value of perpendicular velocity
#' @param v_par \code{numeric} value of parallel velocity
#' @param particles \code{list} list of particles with distribution setups.
#'
#' @return \code{data.frame} with evaluated distribution in column
#'   distribution_value
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @export
calculate_distribution_data_frame <- function(particles, v_par, v_perp){

  if (!rlang::is_named(particles)) {
    particles <-
      rlang::set_names(particles, unlist(purrr::transpose(particles)[["name"]]))
  }

  dist_df <- expand.grid(
    name = names(particles),
    v_par = v_par,
    v_perp = v_perp
  )

  dist_df[["p_scale"]] <-
    sapply(
      dist_df$name,
      function(x){particles[[x]][["distribution"]][["p_scale"]]}
    )

  dist_df[["A"]] <-
    sapply(dist_df$name, function(x){particles[[x]][["A"]]})


  dist_df <- dist_df %>%
    dplyr::mutate(
      p_par = .data[["v_par"]] * .data[["A"]] * const[["amu"]],
      p_perp = .data[["v_perp"]] * .data[["A"]] * const[["amu"]],
      distribution =
        purrr::transpose(particles)[["distribution"]][.data[["name"]]]
    )

  dist_df[["value"]] <-
    furrr::future_pmap_dbl(
      .l = dist_df %>%
        dplyr::select(
          .data[["distribution"]],
          .data[["p_perp"]],
          .data[["p_par"]]
          ),
      .f = eval_distribution
    )

  dist_df[["distribution"]] <- NULL

  dist_df %>% as.data.frame()
}


#' @title plot_dist
#'
#' @description helper function to plot momentum distributions
#'
#' @param dist_df \code{data.frame} output of
#'   calculate_distribution_data_frame() containing the evaluated distribution
#'
#' @param velocity_dist \code{logical} plot velocity distribution (TRUE) or
#'   momentum distribution (FALSE)
#' @return \code{ggplot}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
plot_dist <- function(dist_df, velocity_dist = TRUE){

  v_par_varies <- length(unique(dist_df[["v_par"]])) != 1
  v_perp_varies <- length(unique(dist_df[["v_perp"]])) != 1
  both_vary = v_par_varies & v_perp_varies

  if (velocity_dist) {
    dist_df[["value"]] <-
      dist_df[["value"]] * (dist_df[["A"]] * const[["amu"]])^3
    ylab <- "Density, s^3/m^6"
  } else {
    ylab <- "Density, s^3/(kg m^6)"
  }


  if (both_vary) {
    dist_df %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes_string(y = "v_par", x = "v_perp", z = "value")
      ) +
      ggplot2::geom_raster(
        mapping = ggplot2::aes_string(fill = "value"),
        interpolate = TRUE
      ) +
      ggplot2::geom_contour(colour = "white", alpha = 0.2) +
      viridis::scale_fill_viridis(option = "plasma") +
      ggplot2::ylab(unname(latex2exp::TeX("$v_{par}$ in m/s^2"))) +
      ggplot2::xlab(latex2exp::TeX("$v_{perp}$ in m/s^2")) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(title = "Density")) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(angle = -45, size = 17)
      ) +
      ggplot2::scale_x_continuous(labels = scales::scientific) +
      ggplot2::coord_fixed(ratio = 1, expand = TRUE)
  } else {
    x <- ifelse(
      test = v_par_varies,
      yes = "v_par",
      no = "v_perp"
    )

    xlab <- ifelse(
      test = v_par_varies,
      yes = "$v_{par}$ in m/s^2",
      no = "$v_{perp}$ in m/s^2"
    )

    dist_df %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes_string(x = x, y = "value", color = "name")
      ) +
      ggplot2::geom_line(size = 1.3) +
      ggplot2::theme(legend.position = "top") +
      ggplot2::ylab(latex2exp::TeX(ylab)) +
      ggplot2::xlab(latex2exp::TeX(xlab)) +
      ggplot2::theme(text = ggplot2::element_text(size = 17))
  }
}
