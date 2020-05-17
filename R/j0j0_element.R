#' @title j0j0_element
#'
#' @description Function to calculate terms in the  unscreened current
#'   fluctuation correlation tensor. I.e th integral in in eq. 25 of Bindslev
#'   1996, JATP 58, 983.
#'
#' @param directions \code{character} vector with two spatial directions: a
#'   combination of x, y, and z. e.g  xx, xy or yz.
#' @param k \code{numeric} length of the fluctuation wavevector.
#' @param phi \code{numeric} angle (in degrees) between the magnetic field and
#'   the fluctuation wavevector
#' @param frequency \code{numeric} fluctuation frequency in Hz.
#' @param B \code{numeric} strength of magnetic field in Tesla.
#' @param A \code{numeric} particle mass number.
#' @param Z \code{numeric} particle charge number.
#' @param distribution \code{list} with the velocity distribution and associated
#'   settings.
#' @param integration_method \code{character} method to use for integration. One
#'   of "stats", "hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas"
#' @return \code{complex}
#'
#' @export
j0j0_element <- function(directions,
                         k,
                         phi,
                         frequency,
                         B,
                         A,
                         Z,
                         distribution,
                         integration_method
){

  directions <- unlist(strsplit(directions, split = ""))

  mass <- A * const$amu
  charge <- Z * const$qe

  omega_c <- charge * B / mass
  omega = 2 * pi * frequency

  k_perp <- sin(phi * pi / 180) * k
  k_par <- cos(phi * pi / 180) * k

  frontfaktor <-  (2 * pi)^2 * mass * charge^2 / k_par
  signs <- c("x" = 1, "y" = -const$i, "z" = 1)

  integral_settings <- create_integral_settings(
    directions = directions,
    k_perp = k_perp,
    k_par = k_par,
    omega = omega,
    omega_c = omega_c,
    mass = mass,
    distribution = distribution,
    integration_method = integration_method
  )

  integral <- evaluate_integral(integral_settings, integration_method)

  frontfaktor * integral * distribution[["p_scale"]] * prod(signs[directions])
}

#' @title create_integral_settings
#' @description creates a list of arguments for the integral
#' @param integral_settings list of settings for the integral
#' @param integration_method name of integration method
#' @return \code{numeric} value of integral
evaluate_integral <- function(integral_settings, integration_method){

  if (integration_method == "stats") {
    integrator = stats::integrate
  } else {
    integrator = cubature::cubintegrate
  }

  integral <- do.call(
    what = integrator,
    args = integral_settings
  )

  if (integration_method == "stats") {
    integral <- integral[["value"]]
  } else {
    integral <- integral[["integral"]]
  }
  integral
}

#' @title create_integral_settings
#' @description creates a list of arguments for the integral
#' @inheritParams j0j0_element
#' @inheritParams j0j0_integrand
#' @return \code{list} with settinngs for the integral
create_integral_settings <- function(directions,
                                     k_perp,
                                     k_par,
                                     omega,
                                     omega_c,
                                     mass,
                                     distribution,
                                     integration_method
){

  integral_settings <- list(
    lower = 0,
    upper = Inf,
    directions = directions,
    k_perp = k_perp,
    k_par = k_par,
    omega = omega,
    omega_c = omega_c,
    mass = mass,
    distribution = distribution
  )

  if (integration_method == "stats") {
    integral_settings <- c(
      integral_settings,
      list(
        subdivisions = 1000L,
        f = vector_integrand
      )
    )
  } else {
    integral_settings <- c(
      integral_settings,
      list(
        nVec = 1024,
        f = matrix_integrand,
        method = integration_method
      )
    )
  }
}

#' @title vector_integrand
#' @description Vectorized wrapper for j0j0_integrand.
#' @inheritParams j0j0_integrand
#' @return \code{numeric}
#' @export
vector_integrand <- function(p_perp,
                             directions,
                             k_perp,
                             k_par,
                             omega,
                             omega_c,
                             mass,
                             distribution
) {
  purrr::map_dbl(
    .x = p_perp,
    .f = j0j0_integrand,
    directions = directions,
    k_perp = k_perp,
    k_par = k_par,
    omega = omega,
    omega_c = omega_c,
    mass = mass,
    distribution = distribution
  )
}

#' @title matrix_integrand
#' @description Vectorized wrapper for j0j0_integrand. Returns a matrix for use
#'   with the curbature package.
#' @inheritParams j0j0_integrand
#' @return \code{matrix}
#' @export
matrix_integrand <- function(p_perp,
                             directions,
                             k_perp,
                             k_par,
                             omega,
                             omega_c,
                             mass,
                             distribution
) {
  matrix(
    apply(
      X = p_perp,
      MARGIN = 2,
      FUN = j0j0_integrand,
      directions = directions,
      k_perp = k_perp,
      k_par = k_par,
      omega = omega,
      omega_c = omega_c,
      mass = mass,
      distribution = distribution
    ),
    ncol = ncol(p_perp)
  )
}
