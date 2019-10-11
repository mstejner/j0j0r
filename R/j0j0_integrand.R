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
j0j0_integrand <- function(
  p_perp,
  directions,
  k_perp,
  k_par,
  omega,
  omega_c,
  mass,
  distribution
){

  # v_par_max <- optimize_distribution_perp(distribution, p_perp)[["maximum"]] *
  #   distribution[["p_scale"]] / m

  l_scale <- 3

  v_par_max <- 0
  p_perp <- p_perp * distribution[["p_scale"]]

  l_0 <- round((omega - v_par_max * k_par) / omega_c)
  l_values <- l_0 + seq(-l_scale, l_scale, 1)

  niter <- 0
  converged <- FALSE
  sum_terms <- numeric()
  sum_terms_all <- numeric()
  l_values_done <- numeric()
  while (!converged) {
    niter <- niter + 1

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

    l_values_done <- sort(c(l_values_done, l_values))
    sum_terms_all <- c(sum_terms_all, sum_terms)
    max_term <- max(abs(sum_terms_all))

    l_values <- new_l_values(max_term, l_values_done, sum_terms, nl = l_scale)

    converged <- length(l_values) == 0
    if (max_term == 0 & niter > 5) {converged <- TRUE}

    if (!converged & niter > 50) {
      warning("could not converge integrand in 50 iterations")
      break()
    }

    l_scale <- ceiling(l_scale * 1.2)

  }
  2 * pi * p_perp * sum(sum_terms_all)
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

#' @title new_l_values
#'
#' @description  construct a new vextor of l-values by expanding up and/or dowm
#'
#' @param max_term \code{numeric} maximum term encountered so far
#' @param l_values_done \code{numeric} vector of l-values already done.
#' @param sum_terms \code{complex} sum terms encountered so far.
#' @param nl \code{numeric} number of l-values to expand with.
#'
#' @return \code{numeric}
#'
new_l_values <- function(max_term, l_values_done, sum_terms, nl){
  l_values <- numeric()
  if (max_term == 0) {
    l_values <- expand_l_values(l_values, l_values_done, "up", nl)
    l_values <- expand_l_values(l_values, l_values_done, "down", nl)
  } else {
    if (abs(sum_terms[[1]]) / max_term > 1e-6) {
      l_values <- expand_l_values(l_values, l_values_done, "down", nl)
    }
    if (abs(utils::tail(sum_terms, 1)) /  max_term > 1e-6) {
    l_values <- expand_l_values(l_values, l_values_done, "up", nl)
    }
  }

  l_values
}

#' @title expand_l_values
#'
#' @description  expand a vector of l-values up or down by nl
#'
#' @param l_values \code{numeric} vector of l-values to expand.
#' @param l_values_done \code{numeric} vector of l-values already done.
#' @param updown \code{character} expand "up" or "down".
#' @param nl \code{numeric} number of l-values to expand with
#'
#' @return \code{numeric}
#'
expand_l_values <- function(l_values, l_values_done, updown, nl){
  switch(
    updown,
    up =  c(l_values, seq(max(l_values_done) + 1, max(l_values_done) + nl, 1)),
    down = c(seq(min(l_values_done) - nl, min(l_values_done) - 1, 1), l_values)
  )
}

