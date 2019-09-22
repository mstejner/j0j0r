#' @title j0j0_spectrum
#'
#' @description Wrapper function for j0j0 to calculate a spectrum of j0j0 values
#'   for multiple frequencies and directions. Parallelized with future/furrr.
#'
#' @param directions \code{character} one or more spatial directions ("x", "y"
#'   or "z").
#' @param k \code{numeric} length of the fluctuation wavevector.
#' @param phi \code{numeric} angle (in degrees) between the magnetic field and
#'   the fluctuation wavevector
#' @param frequencies \code{numeric} fluctuation frequencies in Hz.
#' @param B \code{numeric} strength of magnetic field in Tesla.
#' @param particles \code{list} with mass, charge, and momentum distribution
#'
#' @return \code{data.frame}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
j0j0_spectrum <- function(
  k,
  phi,
  frequencies,
  directions = c("x", "y", "z"),
  B,
  particles
) {

  spec <-
    purrr::cross_df(
      .l = list(
        k = k,
        phi = phi,
        frequency = frequencies,
        d1 = directions,
        d2 = directions,
        B = B,
        particle = names(particles)
      ),
      .filter = function(k, phi, frequency, d1, d2, B, particle) { d1 > d2 }
    ) %>%
    dplyr::mutate(
      directions = paste0(.data[["d1"]], .data[["d2"]]),
      A = sapply(.data[["particle"]], function(x){particles[[x]][["A"]]}),
      Z = sapply(.data[["particle"]], function(x){particles[[x]][["Z"]]}),
      distribution = lapply(.data[["particle"]], function(x){particles[[x]][["distribution"]]})
    ) %>%
    dplyr::select(-.data[["d1"]], -.data[["d2"]])

  spec %>%
    dplyr::mutate(
      j0j0 =
        furrr::future_pmap(
          .l = spec %>% dplyr::select(-.data[["particle"]]),
          .f = j0j0r::j0j0,
          .progress = TRUE
        ) %>%
        unlist()
    ) %>%
    dplyr::select(-.data[["distribution"]])

}
