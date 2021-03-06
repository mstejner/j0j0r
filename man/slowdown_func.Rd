% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wilkie_sd.R
\name{slowdown_func}
\alias{slowdown_func}
\title{slowdown_func}
\usage{
slowdown_func(p_perp, p_par, n, K, p_c, p_b, b)
}
\arguments{
\item{p_perp}{\code{numeric} value of perpendicular momentum component}

\item{p_par}{\code{numeric} value of parallel momentum component}

\item{n}{\code{numeric} particle density.}

\item{K}{\code{numeric} integration constant}

\item{p_c}{\code{numeric} critical momentum}

\item{p_b}{\code{numeric} birth momentum}

\item{b}{\code{numeric} transport parameter}
}
\value{
\code{numeric} value of momentum distribution at (p_perp, p_par)
}
\description{
function to evaluate a slowdown momentum distribution,
  https://arxiv.org/abs/1808.01934v2, eq. 2.9
}
