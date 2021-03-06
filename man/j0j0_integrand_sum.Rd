% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/j0j0_integrand.R
\name{j0j0_integrand_sum}
\alias{j0j0_integrand_sum}
\title{j0j0_integrand_sum}
\usage{
j0j0_integrand_sum(
  go_up_or_down,
  l_0,
  l_scale,
  p_perp,
  directions,
  k_perp,
  k_par,
  omega,
  omega_c,
  mass,
  distribution
)
}
\arguments{
\item{go_up_or_down}{\code{character} starting from l_0, go "up" or go "down"}

\item{l_0}{\code{integer} Starting point for sum over l values.}

\item{l_scale}{\code{integer} Relevant range of l values}

\item{p_perp}{\code{numeric} component of momentum vector perpendicular to
the magnetic field divided by p_scale. Should be near 1.}

\item{directions}{\code{character} vector with two spatial directions:
acombination of x, y, and z. e.g  xx, xy or yz.}

\item{k_perp}{\code{numeric} length of component of wavevector perpendicular
to the magnetic field.}

\item{k_par}{\code{numeric} length of component of wavevector parallel to the
magnetic field.}

\item{omega}{\code{numeric} angular cyclotron frequency of current
fluctuations.}

\item{omega_c}{\code{numeric} particle angular cyclotron frequency (omega_c =
qB/m).}

\item{mass}{\code{numeric} particle mass in kg.}

\item{distribution}{\code{list} with the velocity distribution.}
}
\value{
\code{numeric}
}
\description{
Function to calculate the sum in the j0j0 integrand.
}
