#' @title wilkie_setup
#'
#' @description Function to return a setup for a wilkie slowdown momentum
#'   distribution, based on Wilkie 2018, https://arxiv.org/abs/1808.01934v2, eq.
#'   2.9.
#'
#' @param b \code{numeric} b-parameter in Wilkie 2018 eq. 2.9. Quantifies
#'   importance of transport. Suggested in the range 0 (no transport) to 10
#'   (significant effect).
#' @param n \code{numeric} Particle density.
#' @param A \code{numeric} Particle mass number
#' @param Z \code{numeric} Particle charge number
#' @param birth_energy \code{numeric} Particle birth energy in eV.
#' @param n_e \code{numeric} Electron density.
#' @param T_e_eV \code{numeric} Electron temperature in eV.
#' @param ions \code{data frame} with information on all ion species, containing
#'   columns "n" (ion density), "A" (ion mass number), and "Z" (ion charge
#'   number)
#'
#' @return \code{list} with momentum distribution setup
#'
#' @export
wilkie_setup <- function(b, n, A, Z, birth_energy, n_e, T_e_eV, ions){
  m <- A * const[["amu"]]

  p_c <- m * critical_velocity(n_e, T_e_eV, ions)
  p_b <- m * birth_velocity(birth_energy, m)

  unnormalized_dist <- list(
    function_name = "wilkie_func_p",
    gradient = "wilkie_grad",
    distargs = list(
      n = 1,
      p_c = p_c,
      p_b = p_b,
      b = b,
      K = 1
    ),
    p_scale = p_b
  )

  K <-  1 / integrate_homogeneous_distribution(unnormalized_dist)

  list(
    function_name = "wilkie_func",
    gradient = "wilkie_grad",
    distargs = list(
      n = n,
      p_c = p_c,
      p_b = p_b,
      b = b,
      K = K
    ),
    p_scale = p_b
  )
}


#' @title wilkie_expr
#'
#' @description expression for a Wilkie slowdown momentum distribution,
#'   https://arxiv.org/abs/1808.01934v2, eq. 2.9, in cylindrical coordinates
wilkie_expr <- expression(
  # (2 * pi * n * K) *
  #   (tau_s / 4 * pi) *
  #   (1 / (p_c^3 + sqrt(p_perp^2 + p_par^2)^3)) *
  #   (
  #     (sqrt(p_perp^2 + p_par^2)^3 / p_b^3) *
  #       ((p_b^3 + p_c^3) / (sqrt(p_perp^2 + p_par^2)^3 + p_b^3))
  #   )^(b / 3)
  #
  #   (2 * pi * n * K) *
  n * K * (1 / (p_c^3 + sqrt(p_perp^2 + p_par^2)^3)) *
    (
      (sqrt(p_perp^2 + p_par^2)^3 / p_b^3) *
        ((p_b^3 + p_c^3) / (sqrt(p_perp^2 + p_par^2)^3 + p_b^3))
    )^(b / 3)
)


#' @title wilkie_expr_p
#'
#' @description expression for a Wilkie slowdown momentum distribution,
#'   https://arxiv.org/abs/1808.01934v2, eq. 2.9, using only length of p
wilkie_expr_p <- expression(
    n * K * (1 / (p_c^3 + p^3)) *
    ((p^3 / p_b^3) * ((p_b^3 + p_c^3) / (p^3 + p_b^3)))^(b / 3)
)

#' @title wilkie_func
#'
#' @description function to evaluate a Wilkie slowdown momentum distribution,
#'   https://arxiv.org/abs/1808.01934v2, eq. 2.9
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param K \code{numeric} integration constant
#' @param p_c \code{numeric} critical momentum
#' @param p_b \code{numeric} birth momentum
#' @param b \code{numeric} transport parameter
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
wilkie_func <- function(p_perp, p_par, n, K, p_c, p_b, b){
  eval(wilkie_expr) * as.numeric(sqrt(p_perp^2 + p_par^2) < p_b)
}

#' @title wilkie_func_p
#'
#' @description function to evaluate a Wilkie slowdown momentum distribution,
#'   https://arxiv.org/abs/1808.01934v2, eq. 2.9, using only length pf p
#'
#' @param p \code{numeric} length of momentum vector, p = norm(c(p_perp, p_par))
#' @param n \code{numeric} particle density.
#' @param K \code{numeric} integration constant
#' @param p_c \code{numeric} critical momentum
#' @param p_b \code{numeric} birth momentum
#' @param b \code{numeric} transport parameter
#'
#' @return \code{numeric} value of momentum distribution at (p_perp, p_par)
#'
wilkie_func_p <- function(p, n, K, p_c, p_b, b){
  eval(wilkie_expr_p) * as.numeric(p < p_b)
}

#' @title wilkie_grad
#'
#' @description function to calculate the gradient of a Wilkie slowdown
#'   momentum distribution with respect to parallel and perpencicular momentum
#'
#' @param p_perp \code{numeric} value of perpendicular momentum component
#' @param p_par \code{numeric} value of parallel momentum component
#' @param n \code{numeric} particle density.
#' @param K \code{numeric} integrations constant
#' @param p_c \code{numeric} critical momentum
#' @param p_b \code{numeric} birth momentum
#' @param b \code{numeric} transport parameter
#'
#' @return \code{list}
#'
wilkie_grad <- deriv(
  expr = wilkie_expr,
  namevec = c("p_perp","p_par"),
  function.arg = c("p_perp", "p_par", "n", "K", "p_c", "p_b", "b")
)


#' @title fast_ion_slowdown_time
#'
#' @description Function to calculate fast ion slowdown time, eq. 2.2 of Wilkie
#'   2018, https://arxiv.org/abs/1808.01934v2
#'
#' @param n_e \code{numeric} electron density.
#' @param T_i_eV \code{numeric} Ion temperaure in eV.
#' @param T_e_eV \code{numeric} Electron temperaure in eV.
#' @param A \code{numeric} Ion mass number
#' @param Z \code{numeric} Ion charge number
#'
#' @return \code{numeric} plasma parameter
#'
#' @export
fast_ion_slowdown_time <- function(n_e, T_i_eV, T_e_eV, A, Z){
  m_i <- A * const[["amu"]]
  m_e <- const[["m_e"]]
  v_te <- thermal_velocity(T_e_eV, m_e)
  Lambda <- plasma_parameter(n_e, T_i_eV, T_e_eV)

  (3 / (16 * sqrt(pi))) *
    (m_i * m_e * v_te^3) /
    (Z^2 * exp(4) * n_e * log(Lambda))
}

#' @title plasma_parameter
#'
#' @description Function to calculate the plasma parameter (eq. 1.8 of Swanson
#'   2008)
#'
#' @param n_e \code{numeric} electron density.
#' @param T_i_eV \code{numeric} Ion temperaure in eV.
#' @param T_e_eV \code{numeric} Electron temperaure in eV.
#'
#' @return \code{numeric} plasma parameter
#'
#' @export
plasma_parameter <- function(n_e, T_i_eV, T_e_eV) {
  (4 * pi / 3) * n_e * debye_length(n_e, T_i_eV, T_e_eV)^3
}

#' @title debye_length
#'
#' @description Function to calculate the debye length (eq. 1.5 of Swanson 2008)
#'
#' @param n_e \code{numeric} electron density.
#' @param T_i_eV \code{numeric} Ion temperaure in eV.
#' @param T_e_eV \code{numeric} Electron temperaure in eV.
#'
#' @return \code{numeric} Debye length
#'
#' @export
debye_length <- function(n_e, T_i_eV, T_e_eV){
  T_i <- T_i_eV * const[["qe"]]
  T_e <- T_e_eV * const[["qe"]]
  e <- const[["qe"]]
  eps0 <- const[["epsilon_0"]]

  kd2 <- (n_e * e^2 / eps0) * (1 / T_i + 1 / T_e)

  sqrt(1 / kd2)
}

#' @title thermal_velocity
#'
#' @description Function to calculate a particles thermal velocity
#'
#' @param T_eV \code{numeric} Particle temperaure in eV.
#' @param m \code{numeric} particle mass in kg.
#'
#' @return \code{numeric} thermal velocity
#'
#' @export
thermal_velocity <- function(T_eV, m){
  sqrt(2 * T_eV * const[["qe"]] / m)
}

#' @title critical_velocity
#'
#' @description Function to calculate the criticl velocity, eq. 2.1 of Wilkie
#'   2018, https://arxiv.org/abs/1808.01934v2.
#'
#' @param n_e \code{numeric} Electron density.
#' @param T_e_eV \code{numeric} Electron temperaure in eV.
#' @param ions \code{data frame} with information on all ion species,
#'   containing columns "n" (ion density), "A" (ion mass number), and "Z"
#'   (ion charge number)
#'
#' @return \code{numeric} critical velocity
#'
#' @export
critical_velocity <- function(n_e, T_e_eV, ions){
  m_e <- const[["m_e"]]
  n_i <- ions[["n"]]
  m_i <- ions[["A"]] * const[["amu"]]
  Z_i <- ions[["Z"]]
  v_te <- thermal_velocity(T_e_eV, m_e)

  v_te * ((3 * sqrt(pi) / 4) * sum(n_i * m_e * Z_i^2 / (n_e * m_i)))^(1 / 3)
}

#' @title birth_velocity
#'
#' @description Function to calculate the birth velocity of a perticle given its
#'   mass and birth energy.
#'
#' @param birth_energy \code{numeric} birth energy in eV.
#' @param m \code{numeric} particle mass in kg.
#'
#' @return \code{numeric} birth velocity
#'
#' @export
birth_velocity <- function(birth_energy, m){
  sqrt(2 * birth_energy * const[["qe"]] / m)
}
