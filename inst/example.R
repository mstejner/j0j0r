library(magrittr)
library(tictoc)

deuterium_maxw = list(
  name = "deuterium_maxw",
  Z = 1,
  A = 2,
  distribution = j0j0r::maxwellian_setup(n = 4e19, T_eV = 2000, A = 2)
)

deuterium_bimaxw = list(
  name = "deuterium_bimaxw",
  Z = 1,
  A = 2,
  distribution = j0j0r::bimaxwellian_setup(
    n = 4e19,
    T_eV_perp = 2000,
    T_eV_par = 2000,
    v_drift = 2e5,
    A = 2
    )
)


phi <- 86
k <- 2 * pi / (j0j0r::const$c / 100e9)
B = 2.5

future::plan(strategy = "multisession", .skip = TRUE)

tic()
specdf <- j0j0r::j0j0_spectrum2(
  k = k,
  phi = c(86),
  frequencies = seq(0, 400e6, length.out = 401),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    deuterium_maxw = deuterium_maxw,
    deuterium_bimaxw = deuterium_bimaxw
    )
)
toc()

specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = as.factor(phi))) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)


specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = directions, linetype = as.factor(phi))) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0)


specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = particle)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)
