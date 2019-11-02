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
j0j0_element <- function(
  directions,
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

  vector_integrand <- Vectorize(
    FUN = j0j0_integrand2,
    vectorize.args = "p_perp"
  )


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

  assertthat::assert_that(
    integration_method %in% c(
      "stats", "hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas"
      ),
    msg = "integration_method not recognized"
  )

  if (integration_method == "stats") {
    integral_settings <- c(
      integral_settings,
      list(subdivisions = 1000L,
           f = vector_integrand,
           integrator = stats::integrate
      )
    )
  } else {
    integral_settings <- c(
      integral_settings,
      list(nVec = 1024,
           f = matrix_integrand,
           method = integration_method,
           integrator = cubature::cubintegrate
      )
    )
  }


  integral <- do.call(
    what = integral_settings[["integrator"]],
    args = integral_settings[names(integral_settings) != "integrator"]
  )
  if (integration_method == "stats") {
    integral <- integral[["value"]]
  } else {
    integral <- integral[["integral"]]
  }

  # integral <- stats::integrate(
  #   f = vector_integrand,
  #   lower = 0,
  #   upper = Inf,
  #   subdivisions = 1000L,
  #   directions = directions,
  #   k_perp = k_perp,
  #   k_par = k_par,
  #   omega = 2 * pi * frequency,
  #   omega_c = omega_c,
  #   m = m,
  #   distribution = distribution
  # )[["value"]]

  frontfaktor * integral * distribution[["p_scale"]] * prod(signs[directions])
}


#' @title matrix_integrand
#'
#' @description Vectorized wrapper for j0j0_integrand. Returns a matrix for use
#'   with the curbature package.
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
#' @return \code{matrix}
#'
#' @export
matrix_integrand <- function(
  p_perp,
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
      FUN = j0j0_integrand2,
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