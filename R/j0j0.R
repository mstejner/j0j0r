#' @title j0j0
#'
#' @description Wrapper function for j0j0_element to calculate a spectrum of
#'   j0j0 values for multiple frequencies and directions. Parallelized with
#'   future/furrr.
#'
#' @param directions \code{character} one or more spatial directions ("x", "y"
#'   or "z").
#' @param k \code{numeric} length of the fluctuation wavevector.
#' @param phi \code{numeric} angle (in degrees) between the magnetic field and
#'   the fluctuation wavevector
#' @param frequencies \code{numeric} fluctuation frequencies in Hz.
#' @param B \code{numeric} strength of magnetic field in Tesla.
#' @param particles \code{list} with mass, charge, and momentum distribution.
#' @param integration_method \code{character} method to use for integration. One
#'   of "stats", "hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas"
#'
#' @return \code{data.frame}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
j0j0 <- function(k,
                 phi,
                 frequencies,
                 directions = c("xx", "xy", "xz", "yy", "yz", "zz"),
                 B,
                 particles,
                 integration_method
) {

  spec <- initiate_spec(
    k = k,
    phi = phi,
    frequencies = frequencies,
    directions = directions,
    B = B,
    particles = particles,
    integration_method = integration_method
  )


  spec[["j0j0"]] <- furrr::future_pmap(
    .l = dplyr::select(spec, -"particle"),
    .f = j0j0r::j0j0_element,
    .progress = TRUE
  ) %>%
    unlist()

  dplyr::select(spec, -"distribution")
}

#' @title initiate_spec
#'
#' @description creates a \code{data.frame}, spec, with all combinations of the
#'   values in the inpurt parameters
#'
#' @inheritParams j0j0
#'
#' @return \code{data.frame}
#'
#' @export
#'
initiate_spec <- function(k,
                          phi,
                          frequencies,
                          directions,
                          B,
                          particles,
                          integration_method
) {

  assertthat::assert_that(
    all(directions %in% c("xx", "xy", "xz", "yy", "yz", "zz")),
    msg = "direction not recognized"
  )

  assertthat::assert_that(
    all(
      integration_method %in% c(
        "stats", "hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas"
      )
    ),
    msg = "integration_method not recognized"
  )

  settings <- list(
    k = k,
    phi = phi,
    frequency = frequencies,
    directions = directions,
    B = B,
    particle = names(particles),
    integration_method = integration_method
  )

  spec <- purrr::cross_df(settings)

  particles <- purrr::transpose(particles)

  spec[["A"]] <- unlist(particles$A[spec$particle])
  spec[["Z"]] <- unlist(particles$Z[spec$particle])
  spec[["distribution"]] <- particles$distribution[spec$particle]

  spec
}