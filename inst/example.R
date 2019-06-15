
deuterium = list(
  name = "deuterium",
  Z = 1,
  A = 2,
  distribution = j0j0r::maxwellian_setup(n = 4e19, T_eV = 2000, A = 2)
)


phi <- 86
k <- 2 * pi / (const$c / 100e9)
B = 2.5

specdf <- j0j0r::j0j0_spectrum(
  k = k,
  phi = phi,
  frequencies = seq(0, 500e6, length.out = 101),
  directions = c("x", "y", "z"),
  B = B,
  particle = deuterium
)

specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = abs(j0j0), color = directions)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)
