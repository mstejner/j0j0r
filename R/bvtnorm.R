#Rules
# * in cylindrical coordinates: p = c(p_perp, p_par)
# * must be normalized so the integral is 1
# * p must be scaled so it is near 1

#' @title bvtnorm_expr
#'
#' @description expression for a bivariate normal momentum distribution in
#'   cylindrical coordinates
bvtnorm_expr <- expression(
  n * exp(
    -(1 / 2) *
      (c(p_perp, p_par) - center) %*%
      (solve(covariance) %*%
      (c(p_perp, p_par) - center))
  ) /
    sqrt(det(covariance))
)

#' @title bvtnorm_func
#'
#' @description function to evaluate a bivariate normal momentum distribution
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param center \code{numeric} center of distribution (2d vector, c(p_perp,
#'   p_par))
#' @param covariance \code{numeric} covariance of distribution (2 by 2 matrix)
#' @param K \code{numeric} normalizations constant
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
bvtnorm_func <- function(p_perp, p_par, n, center, covariance, K){
  mapply(
    FUN = function(p_perp, p_par, n, center, covariance, K){
      K * n *
        mvtnorm::dmvnorm(
          x = c(p_perp, p_par),
          mean = center,
          sigma = covariance
        )
    },
    p_perp = p_perp,
    p_par = p_par,
    MoreArgs = list(
      n = n,
      center = center,
      covariance = covariance,
      K = K
    )
  )
}

#' #' @title bvtnorm_grad
#' #'
#' #' @description function to calculate the gradient of a a bivariate normal
#' #'   momentum distribution with respect to parallel and perpencicular momentum
#' #'
#' #' @param p_perp \code{numeric} value of perpendicular momentum component
#' #' @param p_par \code{numeric} value of parallel momentum component
#' #' @param n \code{numeric} particle density.
#' #' @param center \code{numeric} center of distribution (2d vector)
#' #' @param covariance \code{numeric} covariance of distribution (2 by 2 matrix)
#' #'
#' #' @return \code{list}
#' #'
#' bvtnorm_grad <- deriv(
#'   expr = bvtnorm_expr,
#'   namevec = c("p_perp","p_par"),
#'   function.arg = c("p_perp", "p_par", "n", "center", "covariance")
#' )


#' @title bvtnorm_setup
#'
#' @description Function to return a bivariate normal momentum distribution
#'
#' @param n \code{numeric} Particle density.'
#' @param center \code{numeric} center of velocity distribution (2d vector)
#' @param covariance \code{numeric} covariance of velocity distribution (2 by 2 matrix)
#' @param A \code{numeric} Particle mass number
#' @param Z \code{numeric} Particle charge number
#' @param name \code{character} Name of distribution/particle
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
bvtnorm_setup <- function(n, center, covariance, A, Z, name){

  unnormalized_dist <- list(
    function_name = "bvtnorm_func",
    gradient = "bvtnorm_grad",
    distargs = list(
      n = 1,
      center = center * A * const[["amu"]],
      covariance = covariance * (A * const[["amu"]])^2,
      K = 1
    ),
    p_scale = sqrt(max(covariance)) * A * const[["amu"]]
  )

  distribution <- list(
    function_name = "bvtnorm_func",
    gradient = "bvtnorm_grad",
    distargs = list(
      n = n,
      center = center * A * const[["amu"]],
      covariance = covariance * (A * const[["amu"]])^2,
      K = 1 / integrate_distribution(unnormalized_dist)
    ),
    p_scale = sqrt(max(covariance)) * A * const[["amu"]]
  )

  list(
    name = name,
    Z = Z,
    A = A,
    distribution = distribution
  )

}
