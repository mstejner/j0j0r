#Rules
# * in cylindrical coordinates: p = c(p_perp, p_par)
# * must be normalized so the integral is 1
# * p must be scaled so it is near 1


#' @title generalized_lorentzian_expr
#'
#' @description expression for a generalized Lorentzian (kappa) distribution,
#'   see
#'   https://www.spenvis.oma.be/help/background/distributions/distributions.html.
#'
generalized_lorentzian_expr <- expression(
  (2 * pi * n / (sqrt(pi) * theta)^3) *
    (gamma(kp + 1) / (kp^1.50 * gamma(kp - 0.5))) *
    (1 + (p_par^2 + p_perp^2) / (kp * theta^2))^(-kp - 1)
)

#' @title find_glz_theta
#'
#' @description Function to calculate the theta parameter in a generalized Lorentzian
#'
#' @param T_eV \code{numeric} Temperature in eV.
#' @param kp \code{numeric} kappa parameter (spectral index).
#' @param A \code{numeric} Particle mass number
#'
#' @return \code{numeric} theta parameter
#'
#' @export
find_glz_theta <- function(T_eV, kp, A) {
  kT <- T_eV * const[["qe"]]
  m = A * const[["amu"]]
  sqrt((2 * kp - 3) * kT * m / kp)
}


#' @title generalized_lorentzian_func
#'
#' @description function to evaluate a generalized Lorentzian distribution
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param theta \code{numeric} Particle themal momentum parameter
#' @param kp \code{numeric} kappa parameter (spectral index). (kappa is a
#'   function in base R so here we use kp for kappa)
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
generalized_lorentzian_func <- function(p_perp, p_par, n, theta, kp){
  eval(generalized_lorentzian_expr)
}

#' @title generalized_lorentzian_grad
#'
#' @description function to calculate the gradient of a generalized Lorentzian
#'   distribution momentum distribution with respect to parallel and
#'   perpendicular momentum
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param theta \code{numeric} Particle themal momentum parameter
#' @param kp \code{numeric} kappa parameter (spectral index).
#'
#' @return \code{list}
#'
generalized_lorentzian_grad <- deriv(
  expr = generalized_lorentzian_expr,
  namevec = c("p_perp","p_par"),
  function.arg = c("p_perp", "p_par", "n", "theta", "kp")
)

#' @title generalized_lorentzian_setup
#'
#' @description Function to return a simple maxwellian momentum distribution
#'
#' @param n \code{numeric} Particle density.
#' @param T_eV \code{numeric} Temperaure in eV.
#' @param kp \code{numeric} kappa parameter (spectral index).
#' @param A \code{numeric} Particle mass number
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
generalized_lorentzian_setup <- function(n, T_eV, kp, A){
  list(
    function_name = "generalized_lorentzian_func",
    gradient = "generalized_lorentzian_grad",
    distargs = list(
      n = n,
      kp = kp,
      theta = find_glz_theta(T_eV, kp, A)
    ),
    p_scale = find_glz_theta(T_eV, kp, A)
  )
}
