
<!-- README.md is generated from README.Rmd. Please edit that file -->

# j0j0r

<!-- badges: start -->

<!-- badges: end -->

The goal of `j0j0r` is to calculate the unscreened current correlation
tensor of the plasma fluctuation model of [Bindslev 1996, Journal of
Atmospheric and Terrestial Physics, **58**,
983](https://www.sciencedirect.com/science/article/pii/0021916995001298).

The code is meant to be distribution-agnostic, i.e. it should be able to
handle any reasonable form of the momentum distribution. Integrals over
the momentum distribution are calculated numerically. Analytic
solutions, possible only in special cases, are not used/supported.

The package contains functions to set up a number of momentum
distributions:

  - An isotropic Maxwellian distribution.

  - A bi-Maxwellian distribution, with drift along the magnetic field
    (z-direction).

  - A [generalized Lorentzian /
    Kappa](https://www.spenvis.oma.be/help/background/distributions/distributions.html)
    distribution.

  - A ring-distribution.

  - A bivariate normal distribution.

  - A slowdown distribution: an isotropic fast-ion slowdown with
    transport modifications as formulated by George J. Wilkie . See
    <https://arxiv.org/abs/1808.01934v2>.

Apart from those, users are free to input new distribution functions of
their own design.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mstejner/j0j0r")
```

The package is not yet on CRAN.

## Vignette

The example below shows the basic work flow for a Maxwellian
distribution. The package vignette gives a more thorough introductiion
and a discussion of the effects of strongly non-Maxwwellian
distributions. It can be viewed with:

``` r
devtools::build_vignettes()
vignette("j0j0r_vignette")
```

## Example

To run the example below, first attach the `j0j0r` and `magrittr`
packages:

``` r
library(j0j0r)
library(magrittr)
```

The code is parallelized using the
[future](https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html)
and [furrr](https://cran.r-project.org/web/packages/furrr/index.html)
packages. To make use of the parallelization, first set a plan for the
`future` package. Different plans will be appropriate for different
operating systems. On Windows it could be:

``` r
future::plan(strategy = "multisession")
```

A Maxwellian distribution can be set up with:

``` r
maxwellian_deuterium <- maxwellian_setup(
  n = 4e19,
  T_eV = 2000,
  A = 2,
  Z = 1,
  name = "maxwellian"
)
```

It can evaluated and be plotted with:

``` r
calculate_distribution_data_frame(
  particles = list(
    maxwellian_deuterium = maxwellian_deuterium
    ), 
  v_par = seq(-1e6, 1e6, length.out = 300), 
  v_perp = seq(0, 1e6, length.out = 300) 
  ) %>% 
  plot_dist() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))
```

<img src="man/figures/README-maxwdist-1.png" width="100%" style="display: block; margin: auto;" />

Assuming a magnetic field of 2.5 T, a wave vector length of 2000/m, and
two values of the resolved angle, the elements of the current
correlation tensor can now be calculated with:

``` r
maxwellian_example <- j0j0(
  k = 2000,
  phi = c(60, 86),
  frequencies = seq(0, 400e6, by = 2e6),
  directions = c("x", "y", "z"),
  B = 2.5,
  particles = list(
    maxwellian = maxwellian_deuterium
  ),
  integration_method = "stats"
)
```

And the results can be plotted with:

``` r
plot_j0j0(maxwellian_example, wrap_by = "element", color_by = "phi")
```

<img src="man/figures/README-maxwwj0j0-1.png" width="100%" style="display: block; margin: auto;" />
