devtools::install("C:/ctsr/j0j0r")
devtools::load_all("C:/ctsr/j0j0r")
library(magrittr)
library(tictoc)

phi <- 86
k <- 2 * pi / (j0j0r::const$c / 100e9)
B = 2.5

future::plan(strategy = "multisession", .skip = TRUE)


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




tic()
specdf <- j0j0r::j0j0_spectrum2(
  k = k,
  phi = c(86),
  frequencies = seq(0, 400e6, length.out = 401),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    deuterium_maxw = deuterium_maxw
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



deuterium_kp2 = list(
  name = "deuterium_kp2",
  Z = 1,
  A = 2,
  distribution = j0j0r::generalized_lorentzian_setup(n = 4e19, T_eV = 2000, kp = 2, A = 2)
)

deuterium_kp15 = list(
  name = "deuterium_kp15",
  Z = 1,
  A = 2,
  distribution = j0j0r::generalized_lorentzian_setup(n = 4e19, T_eV = 2000, kp = 5, A = 2)
)

deuterium_kp50 = list(
  name = "deuterium_kp50",
  Z = 1,
  A = 2,
  distribution = j0j0r::generalized_lorentzian_setup(n = 4e19, T_eV = 2000, kp = 50, A = 2)
)

phi <- 86
k <- 2 * pi / (j0j0r::const$c / 100e9)
B = 2.5

future::plan(strategy = "multisession", .skip = TRUE)

tic()
specdf <- j0j0r::j0j0_spectrum2(
  k = k,
  phi = c(86),
  frequencies = seq(0, 400e6, length.out = 40),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    deuterium_kp15 = deuterium_kp15
  )
)
toc()

specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = particle)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)

deuterium_ring = list(
  name = "deuterium_ring",
  Z = 1,
  A = 2,
  distribution = maxwellian_ring_setup(n = 4e19, A = 2, v_width = 1e5, v_rad = 0.5e6)
)

distributions = list(
  deuterium_maxw = deuterium_maxw,
  deuterium_bimaxw = deuterium_bimaxw,
  deuterium_kp50 = deuterium_kp50,
  deuterium_kp15 = deuterium_kp15,
  deuterium_kp2 = deuterium_kp2,
  deuterium_ring = deuterium_ring
)
v_par <- seq(-2e6, 2e6, length.out = 1000)
v_perp <- 0
plot_distribution(distributions, v_par, v_perp) + ggplot2::scale_y_log10(limits = c(1e80, 1e83))


n <- 4e19
A <- 2
center <- c(0, 1e6)
T_eV <- 2e3
v_term <- find_p_term(T_eV, A) / (A * const[["amu"]])
covariance <- rbind(c(v_term^2, 0), c(0, (v_term)^2))

bvtnorm <- bvtnorm_setup(n, A, center, covariance)
bvtnorm
integrate_distribution(bvtnorm)

v_par <- seq(-2e6, 2e6, length.out = 1000)
v_perp <- 0
plot_distribution(particles = list(bvtnorm = list(A = 2, distribution = bvtnorm)), v_par, v_perp)


deuterium_wilkie= list(
  name = "deuterium_wilkie",
  Z = 1,
  A = 2,
  distribution = wilkie_setup(
    b = 3,
    n = 4e19,
    A = 2,
    Z = 1,
    birth_energy = 50e3,#3.52e6,
    n_e = 4e19,
    T_e_eV = 2e3,
    ions = data.frame(
      Z = c(1),
      A = c(2),
      n = c(4e19)
    )
  )
)

v_par <- seq(-5e6, 5e6, length.out = 10000)
v_perp <- 0
plot_distribution(particles = list(wilkie = deuterium_wilkie, maxw = deuterium_maxw), v_par, v_perp, logscale = FALSE)


tic()
specdf <- j0j0r::j0j0_spectrum2(
  k = k,
  phi = c(86),
  frequencies = seq(0, 100e7, 1e6),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    dwk = deuterium_wilkie
  )
)
toc()
specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = particle)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)


specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = directions)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0)
