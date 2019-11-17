#' @title eval_distribution
#'
#' @description evaluate a distribution at specific values of p_perp and p_par
#'
#' @param distribution \code{list} with the momentum distribution.
#' @param p_perp \code{numeric} value of perpendicular momentum component.
#' @param p_par \code{numeric} value of parallel momentum component.
#'
#' @return \code{numeric} value of distribution at (p_perp, p_par).
#'
#' @export
eval_distribution <- function(distribution, p_perp, p_par){
  do.call(
    what = distribution[["function_name"]],
    args = c(
      list(
        p_perp = p_perp,
        p_par = p_par
      ),
      distribution[["distargs"]]
    )
  )
}

#' @title eval_distribution_list
#'
#' @description evaluate a list of distributions at specific values of p_perp
#'   and p_par
#'
#' @param distribution_list \code{list} with a list of momentum distributions.
#' @param p_perp \code{numeric} value of perpendicular momentum component.
#' @param p_par \code{numeric} value of parallel momentum component.
#'
#' @return \code{numeric} sum of values of distributions at (p_perp, p_par).
#'
#' @export
eval_distribution_list <- function(distribution_list, p_perp, p_par){
  rowSums(
    sapply(
      X = distribution_list,
      FUN = eval_distribution,
      p_perp = p_perp,
      p_par = p_par

    )
  )
}

#' @title integrand_distribution
#'
#' @description function to calculate the integrand in an integral of a
#'   distribution over all momentum space. The integral is carried out in
#'   cylindrical coordinates, so the volume element in the integral is p_perp *
#'   dp_perp * dp_par. Therefore the distribution is here multiplied by p_perp.
#'   Also note the input p must be of order 1 and is scaled in the function by
#'   distribution[["p_scale"]].
#'
#' @param distribution \code{list} with the momentum distribution.
#' @param p \code{numeric} momentum vector: p = c(p_perp, p_par). Of order 1.
#'
#' @return \code{numeric} value of integrand at p = c(p_perp, p_par)
#'
#' @export
integrand_distribution <- function(p, distribution){
  #TODO: kan vektoriseres?
  p <- p * distribution[["p_scale"]]
  2 * pi * p[[1]] *
    eval_distribution(
    distribution = distribution,
    p_perp = p[[1]],
    p_par = p[[2]]
  )
}

#' @title integrate_distribution
#'
#' @description function to calculate integral of a
#'   distribution over all momentum space. The integral is carried out in
#'   cylindrical coordinates.
#'
#' @param distribution \code{list} with the momentum distribution.
#'
#' @return \code{numeric} value of integral
#'
#' @export
integrate_distribution <- function(distribution){
  #TODO: kan bruge vektoriseret integrand_distribution
  cubature::adaptIntegrate(
    f = integrand_distribution,
    lowerLimit = c(0, -Inf),
    upperLimit = c(Inf, Inf),
    distribution = distribution
  )[["integral"]]*distribution[["p_scale"]]^2
}

#' @title integrate_distribution_list
#'
#' @description function to calculate integral of a summed list of
#'   distributions over all momentum space. The integral is carried out in
#'   cylindrical coordinates.
#'
#'
#' @param distribution_list \code{list} wit a list of momentum distributions.
#'
#' @return \code{numeric}
#'
#' @export
integrate_distribution_list <- function(distribution_list){
  #TODO: kan goeres smartere ....
  sum(
    sapply(
      X = distribution_list,
      FUN = integrate_distribution
    )
  )
}

#' @title optimize_distribution_perp
#'
#' @description function to find the optimum of a momentum distribution for
#'   constant perpendicular momentum.
#'
#'
#' @param distribution \code{list} with the momentum distributions.
#' @param p_perp \code{numeric} value of perpendicular momentum component.
#' @param maximum \code{logical} maximize (TRUE) or minimize (FALSE).
#'
#' @return \code{numeric} values of p_par and the distribution at the optimum
#'
#' @export
optimize_distribution_perp <- function(distribution, p_perp, maximum = TRUE){
  stats::optimize(
    f = function(p_par, p_perp, distribution){
      eval_distribution(
        distribution = distribution,
        p_perp = p_perp * distribution[["p_scale"]],
        p_par = p_par * distribution[["p_scale"]]
      )
    },
    distribution = distribution,
    p_perp = p_perp,
    lower = -100,
    upper = 100,
    maximum = maximum
  )
}


#' @title root_distribution_perp
#'
#' @description function to find a root (to within tol) of a momentum
#'   distribution for constant perpendicular momentum.
#'
#'
#' @param distribution \code{list} with the momentum distributions.
#' @param p_perp \code{numeric} value of perpendicular momentum component.
#' @param tol \code{numeric} tolerance.
#'
#' @return \code{numeric} values of p_par and the distribution at the root
#'
#' @export
root_distribution_perp <- function(distribution, p_perp, tol = 1e-10){
  stats::uniroot(
    f = function(p_par, p_perp, distribution){
      eval_distribution(
        distribution = distribution,
        p_perp = p_perp * distribution[["p_scale"]],
        p_par = p_par * distribution[["p_scale"]]
      ) /
        eval_distribution(
          distribution = distribution,
          p_perp = 0,
          p_par = 0
        ) - tol
    },
    distribution = distribution,
    p_perp = p_perp,
    interval = c(0, 1) * 100
  )
}


#' @title eval_homogeneous_distribution
#'
#' @description evaluate a homogeneous distribution at a specific values of p
#'
#' @param distribution \code{list} with the momentum distribution.
#' @param p \code{numeric} length of momentum vector.
#'
#' @return \code{numeric} value of distribution at p.
#'
#' @export
eval_homogeneous_distribution <- function(distribution, p){
  do.call(
    what = distribution[["function_name"]],
    args = c(
      list(
        p = p
      ),
      distribution[["distargs"]]
    )
  )
}

#' @title integrand_homogeneous_distribution
#'
#' @description function to calculate the integrand in an integral of a
#'   distribution over all momentum space. The integral is carried out in
#'   spherical coordinates. Therefore the distribution is here multiplied by
#'   p^2. Also note the input p must be of order 1 and is scaled in the function
#'   by distribution[["p_scale"]].
#'
#' @param distribution \code{list} with the momentum distribution.
#' @param p \code{numeric} length ofmomentum vector: p = norm(c(p_perp, p_par)).
#'   Of order 1.
#'
#' @return \code{numeric} value of integrand at p
#'
#' @export
integrand_homogeneous_distribution <- function(p, distribution){
  #TODO: kan vektoriseres?
  p <- p * distribution[["p_scale"]]
  2 * pi * p^2 *
    eval_homogeneous_distribution(
      distribution = distribution,
      p = p
    )
}

#' @title integrate_homogeneous_distribution
#'
#' @description function to calculate integral of a homogeneous
#'   distribution over all momentum space. The integral is carried out in
#'   spherical coordinates.
#'
#' @param distribution \code{list} with the momentum distribution.
#'
#' @param limits \code{list} integration limits.
#'
#' @return \code{numeric} value of integral
#'
#' @export
integrate_homogeneous_distribution <- function(distribution, limits = c(0, Inf)){
  #TODO: kan bruge vektoriseret integrand_distribution
  cubature::adaptIntegrate(
    f = integrand_homogeneous_distribution,
    lowerLimit = limits[[1]],
    upperLimit = limits[[2]],
    distribution = distribution
  )[["integral"]]*distribution[["p_scale"]]
}
