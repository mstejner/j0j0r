#' @title j0j0_spectrum
#'
#' @description Wrapper function for j0j0 to calculate a spectrum of j0j0 values
#'   for multiple frequencies and directions
#'
#' @param directions \code{character} one or more spatial directions.
#' @param k \code{numeric} length of the fluctuation wavevector.
#' @param phi \code{numeric} angle (in degrees) between the magnetic field and
#'   the fluctuation wavevector
#' @param frequencies \code{numeric} fluctuation frequencies in Hz.
#' @param B \code{numeric} strength of magnetic field in Tesla.
#' @param particle \code{list} with mass, charge, and momentum distribution
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
  particle
){

  vector_j0j0 <- Vectorize(
    FUN = j0j0,
    vectorize.args = c("frequency", "directions")
  )

  df <- expand.grid(
    frequency = frequencies,
    direction1 = directions,
    direction2 = directions,
    stringsAsFactors = FALSE
    ) %>%
    dplyr::filter(.data[["direction1"]] >= .data[["direction2"]]) %>%
    dplyr::mutate(
      directions = paste0(.data[["direction1"]], .data[["direction2"]])
    )

  df %>%
    dplyr::mutate(
      j0j0 = vector_j0j0(
        directions = .data[["directions"]],
        k = k,
        phi = phi,
        frequency = .data[["frequency"]],
        B = B,
        A = particle[["A"]],
        Z = particle[["Z"]],
        distribution = particle[["distribution"]]
      )
    )
}