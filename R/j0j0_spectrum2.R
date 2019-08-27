#' @title j0j0_spectrum2
#'
#' @description Wrapper function for j0j0 to calculate a spectrum of j0j0 values
#'   for multiple frequencies and directions.
#'
#' @param directions \code{character} one or more spatial directions.
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
j0j0_spectrum2 <- function(
  k,
  phi,
  frequencies,
  directions = c("x", "y", "z"),
  B,
  particles
){

  speclist = purrr::cross(
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
  )

  specdf <- lapply(
    X = speclist,
    FUN = function(x){as.data.frame(x, stringsAsFactors = FALSE)}
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(directions = paste0(.data[["d1"]], .data[["d2"]])) %>%
    dplyr::select(-.data[["d1"]], -.data[["d2"]])

  specdf %>%
    dplyr::mutate(
      j0j0 = furrr::future_map(
        .x = speclist,
        .f = function(x, particles) {
          j0j0r::j0j0(
            directions = paste0(x[["d1"]], x[["d2"]]),
            k = x[["k"]],
            phi = x[["phi"]],
            frequency = x[["frequency"]],
            B = x[["B"]],
            A = particles[[x[["particle"]]]][["A"]],
            Z = particles[[x[["particle"]]]][["Z"]],
            distribution = particles[[x[["particle"]]]][["distribution"]]
          )
        },
        particles = particles
      ) %>%
        unlist()
    )
}
