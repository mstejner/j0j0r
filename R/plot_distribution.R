#' @title plot_distribution
#'
#' @description function to plot velocity distributions
#'
#' @param v_perp \code{numeric} value of perpendicular velocity
#' @param v_par \code{numeric} value of parallel velocity
#' @param particles \code{list} list of particles with distribution setups.
#' @param logscale \code{logical} if TRUE sets log scale on y-axis
#'
#' @return \code{ggplot} 1D plot of distributions
#'
plot_distribution <- function(particles, v_par, v_perp, logscale = TRUE){

  dist_df <- expand.grid(
    distname = names(particles),
    v_par = v_par,
    v_perp = v_perp
  )

  dist_df[["p_scale"]] <-
    sapply(
      dist_df$distname,
      function(x){particles[[x]][["distribution"]][["p_scale"]]}
    )

  dist_df[["A"]] <-
    sapply(dist_df$distname, function(x){particles[[x]][["A"]]})


  dist_df <- dist_df %>%
    dplyr::mutate(
      p_par = v_par * A * const[["amu"]],
      p_perp = v_perp * A * const[["amu"]],
      distribution = purrr::transpose(particles)[["distribution"]][distname]
    )

  dist_df[["distval"]] <-
    purrr::pmap_dbl(
      .l = dist_df %>% dplyr::select(distribution, p_perp, p_par),
      .f = eval_distribution
    )

  plot <- dist_df %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(x = v_par, y = distval, color = distname)
    ) +
    ggplot2::geom_line()

  if (logscale == TRUE) {plot <- plot + ggplot2::scale_y_log10()}

  plot
}
