% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxwellian_ring.R
\name{maxwellian_ring_grad}
\alias{maxwellian_ring_grad}
\title{maxwellian_ring_grad}
\usage{
maxwellian_ring_grad(p_perp, p_par, n, p_width, p_rad, K)
}
\arguments{
\item{p_perp}{\code{numeric} value of perpendicular momentum component}

\item{p_par}{\code{numeric} value of parallel momentum component}

\item{n}{\code{numeric} particle density.}

\item{p_width}{\code{numeric} Width of ring momentum distribution}

\item{p_rad}{\code{numeric} Radius of ring momentum distribution}

\item{K}{\code{numeric} Integration constant}
}
\value{
\code{list}
}
\description{
function to calculate the gradient of a  maxwellian ring
  momentum distribution with respect to parallel and perpencicular momentum
}
