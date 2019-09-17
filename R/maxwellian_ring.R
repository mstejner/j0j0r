#Rules
# * in cylindrical coordinates: p = c(p_perp, p_par)
# * must be normalized so the integral is 1
# * p must be scaled so it is near 1

#' @title maxwellian_ring_expr
#'
#' @description expression for a maxwellian ring distribution in
#'   cylindrical coordinates
maxwellian_ring_expr <- expression(
  2 * pi * n / (sqrt(pi) * p_width)^3 * K *
    exp(-(sqrt(p_par^2 + p_perp^2) - p_rad)^2 / p_width^2)
)

#' @title maxwellian_ring_func
#'
#' @description function to evaluate a simple maxwellian momentum distribution
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_width \code{numeric} Width of ring momentum distribution
#' @param p_term \code{numeric} Radius of ring momentum distribution
#' @param K \code{numeric} Integration constant
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
maxwellian_ring_func <- function(p_perp, p_par, n, p_width, p_rad, K){
  eval(maxwellian_ring_expr)
}

#' @title maxwellian_ring_grad
#'
#' @description function to calculate the gradient of a  maxwellian ring
#'   momentum distribution with respect to parallel and perpencicular momentum
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param p_width \code{numeric} Width of ring momentum distribution
#' @param p_term \code{numeric} Radius of ring momentum distribution
#' @param K \code{numeric} Integration constant
#'
#' @return \code{list}
#'
maxwellian_ring_grad <- deriv(
  expr = maxwellian_ring_expr,
  namevec = c("p_perp","p_par"),
  function.arg = c("p_perp", "p_par", "n", "p_width", "p_rad", "K")
)


#' @title maxwellian_ring_setup
#'
#' @description Function to return a maxwellian ring momentum distribution
#'
#' @param n \code{numeric} Particle density.
#' @param A \code{numeric} Particle mass number
#' @param v_width \code{numeric} Width of ring in terms of velocity.
#' @param v_rad \code{numeric} Radius of ring in terms of velocity.
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
maxwellian_ring_setup <- function(n, A, v_width, v_rad){
  dist_unnormalized <- list(
    function_name = "maxwellian_ring_func",
    gradient = "maxwellian_ring_grad",
    distargs = list(
      n = 1,
      p_width = v_width * A * const[["amu"]],
      p_rad = v_rad *  A * const[["amu"]],
      K = 1
    ),
    p_scale = v_rad *  A * const[["amu"]]
  )

  list(
    function_name = "maxwellian_ring_func",
    gradient = "maxwellian_ring_grad",
    distargs = list(
      n = n,
      p_width = v_width *  A *const[["amu"]],
      p_rad = v_rad *  A * const[["amu"]],
      K = 1 / integrate_distribution(dist_unnormalized)
    ),
    p_scale = v_rad *  A * const[["amu"]]
  )
}

