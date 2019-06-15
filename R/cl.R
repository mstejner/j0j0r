#' @title cl.
#'
#' @description Function to calculate cl terms, eq. 26 of Bindslev 1996, JATP
#'   58, 983.
#'
#' @param l \code{numeric} vector of l-values (Index of Bessel functions).
#' @param direction \code{character} spatial direction: x, y or z.
#' @param k_perp \code{numeric} length of component of wavevector perpendicular
#'   to the magnetic field.
#' @param k_par \code{numeric} length of component of wavevector parallel to the
#'   magnetic field.
#' @param v_perp \code{numeric} particle velocity component perpendicular to the
#'   magnetic field.
#' @param v_par \code{numeric} particle velocity component  parallel to the
#'   magnetic field.
#' @param omega_c \code{numeric} particle angular cyclotron frequency (omega_c =
#'   qB/m).
#'
#' @return \code{complex} value of cl term
#'
#' @export
cl <- function(l, direction, k_perp, k_par, v_perp, v_par, omega_c){
  rho <- v_perp / omega_c
  switch(
    direction,
    x = (l * omega_c / k_perp) * besselJ(k_perp * rho, l),
    y = const$i * v_perp * (besselJ(k_perp * rho, l - 1) + besselJ(k_perp * rho, l + 1)) / 2,
    z = v_par * besselJ(k_perp * rho, l)
  )
}