#' @title calculate_distribution_data_frame
#'
#' @description function to evaluate velocity distributions and return in a data
#'   frame
#'
#' @param v_perp \code{numeric} value of perpendicular velocity
#' @param v_par \code{numeric} value of parallel velocity
#' @param particles \code{list} list of particles with distribution setups.
#' @param logscale \code{logical} if TRUE sets log scale on y-axis
#'
#' @return \code{data.frame} with evaluated distribution in column
#'   distribution_value
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
      p_par = v_par * A * const[["amu"]],
      p_perp = v_perp * A * const[["amu"]],
      distribution = purrr::transpose(particles)[["distribution"]][name]
    )

  dist_df[["value"]] <-
    furrr::future_pmap_dbl(
      .l = dist_df %>% dplyr::select(distribution, p_perp, p_par),
      .f = eval_distribution
    )

  dist_df[["distribution"]] <- NULL

  dist_df %>% as.data.frame()
}
