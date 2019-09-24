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
#'
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
  distribution
){

  directions <- unlist(strsplit(directions, split = ""))

  m <- A * const$amu
  q <- Z * const$qe

  omega_c <- q * B / m

  k_perp <- sin(phi * pi / 180) * k
  k_par <- cos(phi * pi / 180) * k

  frontfaktor <-  (2 * pi)^2 * m * q^2 / k_par
  signs <- c("x" = 1, "y" = -const$i, "z" = 1)

  vector_integrand <- Vectorize(
    FUN = j0j0_integrand,
    vectorize.args = "p_perp"
  )

  integral <- stats::integrate(
    f = vector_integrand,
    lower = 0,
    upper = Inf,
    directions = directions,
    k_perp = k_perp,
    k_par = k_par,
    omega = 2 * pi * frequency,
    omega_c = omega_c,
    m = m,
    distribution = distribution
  )

  frontfaktor * integral[["value"]] *
    distribution[["p_scale"]] * prod(signs[directions])
}