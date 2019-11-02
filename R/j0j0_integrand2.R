#' @title j0j0_integrand
#'
#' @description Function to calculate the integrand in eq. 25 of Bindslev 1996,
#'   JATP 58, 983.
#'
#' @param p_perp \code{numeric} component of momentum vector perpendicular to
#'   the magnetic field divided by p_scale. Should be near 1.
#' @param directions \code{character} vector with two spatial directions:
#'   acombination of x, y, and z. e.g  xx, xy or yz.
#' @param k_perp \code{numeric} length of component of wavevector perpendicular
#'   to the magnetic field.
#' @param k_par \code{numeric} length of component of wavevector parallel to the
#'   magnetic field.
#' @param omega \code{numeric} angular cyclotron frequency of current
#'   fluctuations.
#' @param omega_c \code{numeric} particle angular cyclotron frequency (omega_c =
#'   qB/m).
#' @param mass \code{numeric} particle mass in kg.
#' @param distribution \code{list} with the velocity distribution.
#'
#' @return \code{numeric}
#'
#' @export
j0j0_integrand2 <- function(
  p_perp,
  directions,
  k_perp,
  k_par,
  omega,
  omega_c,
  mass,
  distribution
){

  # v_par_max <-
  #   optimize_distribution_perp(distribution, p_perp)[["maximum"]] *
  #   distribution[["p_scale"]] / mass

  p_perp <- p_perp * distribution[["p_scale"]]

  v_scale <- distribution[["p_scale"]] / mass

  l_0 <- round((omega - 0 * k_par) / omega_c)
  l_scale <- max(c(10,  abs(round((omega - v_scale * k_par) / omega_c))))

  go_up_or_down <- list("up" = "up", "down" = "down")

  j0j0_sum <- lapply(
    X = go_up_or_down,
    FUN = j0j0r:::j0j0_integrand_sum,
    l_0 = l_0,
    l_scale = l_scale,
    p_perp = p_perp,
    directions = directions,
    k_perp = k_perp,
    k_par = k_par,
    omega = omega,
    omega_c = omega_c,
    mass = mass,
    distribution = distribution
    )

  2 * pi * p_perp * (j0j0_sum[["up"]] + j0j0_sum[["down"]])
}

#' @title j0j0_integrand_sum
#'
#' @description Function to calculate the sum in the j0j0 integrand.
#'
#' @param p_perp \code{numeric} component of momentum vector perpendicular to
#'   the magnetic field divided by p_scale. Should be near 1.
#' @param directions \code{character} vector with two spatial directions:
#'   acombination of x, y, and z. e.g  xx, xy or yz.
#' @param k_perp \code{numeric} length of component of wavevector perpendicular
#'   to the magnetic field.
#' @param k_par \code{numeric} length of component of wavevector parallel to the
#'   magnetic field.
#' @param omega \code{numeric} angular cyclotron frequency of current
#'   fluctuations.
#' @param omega_c \code{numeric} particle angular cyclotron frequency (omega_c =
#'   qB/m).
#' @param mass \code{numeric} particle mass in kg.
#'
#' @param distribution \code{list} with the velocity distribution.
#'
#' @param l_0 \code{integer} Starting point for sum over l values.
#'
#' @param l_scale \code{integer} Relevant range of l values
#'
#' @param go_up_or_down \code{character} starting from l_0, go "up" or go "down"
#'
#' @return \code{numeric}
#'
#' @export
j0j0_integrand_sum <- function(
  go_up_or_down,
  l_0,
  l_scale,
  p_perp,
  directions,
  k_perp,
  k_par,
  omega,
  omega_c,
  mass,
  distribution
){

  niter <- 0
  converged <- FALSE
  j0j0_sum <- 0
  while (!converged) {

    l_values <- l_0 +
      switch(
        go_up_or_down,
        up =  seq(niter * l_scale, (niter + 1) * l_scale - 1, by = 1),
        down = seq(-niter * l_scale - 1, -(niter + 1) * l_scale, by = -1)
      )

    sum_terms <- j0j0r:::j0j0_sum_terms(
      l_values,
      p_perp,
      directions,
      k_perp,
      k_par,
      omega,
      omega_c,
      mass,
      distribution
    )

    sum_sum_terms <- sum(sum_terms)

    j0j0_sum <- j0j0_sum + sum_sum_terms

    converged <- abs(j0j0_sum) > 0 & abs(sum_sum_terms / j0j0_sum) < 1e-4

    #if (j0j0_sum == 0 & max(abs(l_values - l_0)) > ) {converged <- TRUE}
    if (j0j0_sum == 0 & niter > 5) {converged <- TRUE}

    if (!converged & niter > 1000) {
      warning("could not converge integrand in 1000 iterations")
      break()
    }

    niter <- niter + 1
  }

  j0j0_sum
}
