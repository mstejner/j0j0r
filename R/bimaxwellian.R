#' @title bimaxwellian_expr
#'
#' @description expression for a bimaxwellian momentum distribution in
#'   cylindrical coordinates
bimaxwellian_expr <- expression(
  (n / (sqrt(pi)^3 * p_tperp^2 * p_tpar)) *
    exp(-((p_par - p_drift)^2/p_tpar^2 + p_perp^2/p_tperp^2))
)

#' @title bimaxwellian_func
#'
#' @description function to evaluate a bimaxwellian momentum distribution
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_tperp \code{numeric} themal momentum perpendicular to magnetic field
#' @param p_tpar \code{numeric} themal momentum parallel to magnetic field
#' @param p_drift \code{numeric} drift momentum parallel to magnetic field
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
bimaxwellian_func <- function(p_perp, p_par, n, p_tperp, p_tpar, p_drift){
  eval(bimaxwellian_expr)
}

#' @title bimaxwellian_grad
#'
#' @description function to calculate the gradient of a  bimaxwellian momentum
#'   distribution with respect to parallel and perpencicular momentum
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_tperp \code{numeric} thermal momentum perpendicular to magnetic
#'   field
#' @param p_tpar \code{numeric} thermal momentum parallel to magnetic field
#' @param p_drift \code{numeric} drift momentum parallel to magnetic field
#'
#' @return \code{list}
#'
bimaxwellian_grad <- deriv(
  expr = bimaxwellian_expr,
  namevec = c("p_perp","p_par"),
  function.arg = c("p_perp", "p_par", "n", "p_tperp", "p_tpar", "p_drift")
)

#' @title bimaxwellian_setup
#'
#' @description Function to return a simple maxwellian momentum distribution
#'
#' @param n \code{numeric} Particle density.
#' @param T_eV_perp \code{numeric} Perpendicular temperaure in eV.
#' @param T_eV_par \code{numeric} Parallel temperaure in eV.
#' @param v_drift \code{numeric} drift velocity parallel to magnetic field in
#'   m/s.
#' @param A \code{numeric} Particle mass number
#' @param Z \code{numeric} Particle charge number
#' @param name \code{character} Name of distribution/particle
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
bimaxwellian_setup <- function(n, T_eV_perp, T_eV_par, v_drift, A, Z, name){

  distribution <- list(
    function_name = "bimaxwellian_func",
    gradient = "bimaxwellian_grad",
    distargs = list(
      n = n,
      p_tperp = find_p_term(T_eV_perp, A),
      p_tpar = find_p_term(T_eV_par, A),
      p_drift = A * const[["amu"]] * v_drift
    ),
    p_scale = find_p_term(T_eV_perp, A)
  )

  list(
    name = name,
    Z = Z,
    A = A,
    distribution = distribution
  )

}
