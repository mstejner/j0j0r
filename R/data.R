#' Examples of results for j0j0 for a sweep of frequency, directions and a
#' number of different velocity distributions:
#'
#'   * maxwellian
#'   * bimaxwellian
#'   * lorentzian
#'   * ring
#'   * bvtnorm
#'   * slowdown_b0
#'   * slowdown_b5
#'
#'
#' @format A tibble with 12642 rows and 10 variables:
#' \describe{
#'   \item{k}{length of wave vector in 1/m}
#'   \item{phi}{angle between wave vector and magnetic field in degrees}
#'   \item{frequency }{in Hz}
#'   \item{B}{Magnetic field strength in Tesla}
#'   \item{particle}{particle name}
#'   \item{integration_method}{name of integration method}
#'   \item{directions}{string with direction names}
#'   \item{A}{mass number}
#'   \item{Z}{charge number}
#'   \item{j0j0}{result for j0j0}
#'   }
#'
#' @source File /data/raw/j0j0_examples
"j0j0_examples"


#' Examples of results for j0j0 for a sweep of frequency, directions and two
#' values of phi for a maxwellian distribution.
#'
#'
#' @format A tibble with 3618 rows and 10 variables:
#' \describe{
#'   \item{k}{length of wave vector in 1/m}
#'   \item{phi}{angle between wave vector and magnetic field in degrees}
#'   \item{frequency }{in Hz}
#'   \item{B}{Magnetic field strength in Tesla}
#'   \item{particle}{particle name}
#'   \item{integration_method}{name of integration method}
#'   \item{directions}{string with direction names}
#'   \item{A}{mass number}
#'   \item{Z}{charge number}
#'   \item{j0j0}{result for j0j0}
#'   }
#'
#' @source File /data/raw/j0j0_examples
"maxwellian_example"


#'  examples of results for j0j0 for a sweep of frequency and phi angle forthe
#'  zz-component for three distributions:
#'
#'  * maxwellian
#'  * ring
#'  * slowdown_b5
#'
#' @format A tibble with 12,642 rows and 10 variables:
#' \describe{
#'   \item{k}{length of wave vector in 1/m}
#'   \item{phi}{angle between wave vector and magnetic field in degrees}
#'   \item{frequency }{in Hz}
#'   \item{B}{Magnetic field strength in Tesla}
#'   \item{particle}{particle name}
#'   \item{integration_method}{name of integration method}
#'   \item{directions}{string with direction names}
#'   \item{A}{mass number}
#'   \item{Z}{charge number}
#'   \item{j0j0}{result for j0j0}
#'   }
#'
#' @source File /data/raw/j0j0_examples
"phi_sweep"


#'  Examples 1D velocity distributions at v_perp = 0
#'
#'  * maxwellian
#'  * ring
#'  * slowdown_b5
#'
#' @format A data frame with 7000 rows and 8 variables:
#' \describe{
#'   \item{name}{Distribution name}
#'   \item{v_par}{Parallel velocity}
#'   \item{v_perp}{Perpendiculer velocity}
#'   \item{p_scale}{Momentum scale - used to scale momentum in integrals, so it is near 1}
#'   \item{A}{mass number}
#'   \item{p_par}{parallel momentum}
#'   \item{p_perp}{perpendicular momentum}
#'   \item{value}{momentum distribution value (i.e. density)}
#'   }
#'
#' @source File /data/raw/j0j0_examples
"distribution_examples"

#'  Examples of 2D velocity distributions with variation in both v_par and
#'  v_perp for a number of velocity distributions:
#'
#'  * maxwellian
#'  * bimaxwellian
#'  * lorentzian
#'  * ring
#'  * bvtnorm
#'  * slowdown_b5
#'  * slowdown_b0
#'
#' @format A data frame with 630000 rows and 8 variables:
#' \describe{
#'   \item{name}{Distribution name}
#'   \item{v_par}{Parallel velocity}
#'   \item{v_perp}{Perpendiculer velocity}
#'   \item{p_scale}{Momentum scale - used to scale momentum in integrals, so it is near 1}
#'   \item{A}{mass number}
#'   \item{p_par}{parallel momentum}
#'   \item{p_perp}{perpendicular momentum}
#'   \item{value}{momentum distribution value (i.e. density)}
#'   }
#'
#' @source File /data/raw/j0j0_examples
"dist_2D_examples"

