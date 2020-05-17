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
j0j0_integrand <- function(p_perp,
                           directions,
                           k_perp,
                           k_par,
                           omega,
                           omega_c,
                           mass,
                           distribution
){

  p_perp <- p_perp * distribution[["p_scale"]]

  l_0 <- round((omega - 0 * k_par) / omega_c)
  l_scale <- 6

  sum_up <- j0j0_integrand_sum(
    go_up_or_down = "up",
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

  sum_down <- j0j0_integrand_sum(
    go_up_or_down = "down",
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

  2 * pi * p_perp * (sum_up + sum_down)
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
  continue <- TRUE
  j0j0_sum <- 0
  while (continue) {

    l_values <- l_0 +
      switch(
        go_up_or_down,
        up =  (niter * l_scale):((niter + 1) * l_scale - 1),
        down = (-niter * l_scale - 1):(-(niter + 1) * l_scale)
      )

    sum_terms <- j0j0_sum_terms(
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

    j0j0_sum <- j0j0_sum + sum(sum_terms)

    continue <- j0j0_sum == 0 || abs(sum_terms[[length(sum_terms)]]) > 1e-10 * abs(j0j0_sum)

    if ((j0j0_sum == 0 & niter > 3) || niter > 1000) break

    niter <- niter + 1
  }

  if (niter > 1000) {
    warning("could not converge integrand in 1000 iterations")
    break()
  }

  j0j0_sum
}

#' @title j0j0_sum_terms
#'
#' @description Function to calculate the terms of the sum in the integrand of
#'   eq. 25 of Bindslev 1996, JATP 58, 983.
#'
#' @param l_values \code{numeric} vector of l-values at which to evaluate terms.
#' @param p_perp \code{numeric} component of momentum vector perpendicular to
#'   the magnetic field.
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
#' @return \code{complex}
#'
#' @export
j0j0_sum_terms <- function(
  l_values,
  p_perp,
  directions,
  k_perp,
  k_par,
  omega,
  omega_c,
  mass,
  distribution
){

  v_perp <- p_perp / mass
  v_par  <- (omega - l_values * omega_c) / k_par
  p_par  <- v_par * mass

  if (identical(directions[[1]], directions[[2]])) {
    clcl <-
      cl(l_values, directions[[1]], k_perp, k_par, v_perp, v_par, omega_c)^2
  } else {
    clcl <-
      cl(l_values, directions[[1]], k_perp, k_par, v_perp, v_par, omega_c) *
      cl(l_values, directions[[2]], k_perp, k_par, v_perp, v_par, omega_c)
  }

  fl <- eval_distribution(distribution, p_perp, p_par)

  clcl * fl
}
