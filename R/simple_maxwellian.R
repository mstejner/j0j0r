#Rules
# * in cylindrical coordinates: p = c(p_perp, p_par)
# * must be normalized so the integral is 1
# * p must be scaled so it is near 1

#' @title maxwellian_expr
#'
#' @description expression for a simple maxwellian momentum distribution in
#'   cylindrical coordinates
maxwellian_expr <- expression(
  n / (sqrt(pi) * p_term)^3 * exp(-(p_par^2 + p_perp^2) / p_term^2)
)

#' @title maxwellian_func
#'
#' @description function to evaluate a simple maxwellian momentum distribution
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_term \code{numeric} Particle themal momentum
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
maxwellian_func <- function(p_perp, p_par, n, p_term){
  n / (sqrt(pi) * p_term)^3 * exp(-(p_par^2 + p_perp^2) / p_term^2)
}

#' @title maxwellian_grad
#'
#' @description function to calculate the gradient of a  simple maxwellian
#'   momentum distribution with respect to parallel and perpencicular momentum
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_term \code{numeric} Particle themal momentum
#'
#' @return \code{list}
#'
maxwellian_grad <- deriv(
  expr = maxwellian_expr,
  namevec = c("p_perp","p_par"),
  function.arg = c("p_perp", "p_par", "n", "p_term")
)



#' @title find_p_term
#'
#' @description Function to calculate thermal momentum for a simple maxwellian
#'
#' @param T_eV \code{numeric} Temperaure in eV.
#' @param A \code{numeric} Particle mass number
#'
#' @return \code{numeric} thermal momentum
#'
#' @export
find_p_term <- function(T_eV, A) {
  kT <- T_eV * const[["qe"]]
  m = A * const[["amu"]]
  sqrt(2 * kT * m )
}


#' @title maxwellian_setup
#'
#' @description Function to return a simple maxwellian momentum distribution
#'
#' @param n \code{numeric} Particle density.
#' @param T_eV \code{numeric} Temperaure in eV.
#' @param A \code{numeric} Particle mass number
#' @param Z \code{numeric} Particle charge number
#' @param name \code{character} Name of distribution/particle
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
maxwellian_setup <- function(n, T_eV, A, Z, name = "maxwellian"){

  distribution <- list(
    function_name = "maxwellian_func",
    gradient = "maxwellian_grad",
    distargs = list(
      n = n,
      p_term = find_p_term(T_eV, A)
    ),
    p_scale = find_p_term(T_eV, A)
  )

  list(
    name = name,
    Z = Z,
    A = A,
    distribution = distribution
  )

}

