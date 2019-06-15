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
  eval_distribution(
    distribution = distribution,
    p_perp = p[[1]],
    p_par = p[[2]]
  ) * p[[1]]
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
#'
#' @return \code{numeric} values of p_par and the distribution at the optimum
#'
#' @export
optimize_distribution_perp <- function(distribution, p_perp){
  stats::optimise(
    f = function(p_par, p_perp, distribution){
      eval_distribution(
        distribution = distribution,
        p_perp = p_perp * distribution[["p_scale"]],
        p_par = p_par * distribution[["p_scale"]]
      )
    },
    distribution = distribution,
    p_perp = p_perp,
    interval = c(-1, 1) * 100,
    maximum = TRUE
  )
}
