devtools::install("C:/ctsr/j0j0r")
devtools::load_all("C:/ctsr/j0j0r")
library(magrittr)
library(tictoc)

phi <- 86
k <- 2 * pi / (j0j0r::const$c / 100e9)
B = 2.5
n = 4e19
A = 2
Z = 1

future::plan(strategy = "multisession", .skip = TRUE)

maxwellian = j0j0r::maxwellian_setup(
  n = n,
  T_eV = 2000,
  A = A,
  Z = Z,
  name = "maxwellian"
)

bimaxwellian <-  j0j0r::bimaxwellian_setup(
  n = n,
  T_eV_perp = 2000,
  T_eV_par = 2000,
  v_drift = 2e5,
  A = A,
  Z = Z,
  name =  "bimaxwellian"
)

lorentzian <- j0j0r::generalized_lorentzian_setup(
  n = n,
  T_eV = 2000,
  kp = 15,
  A = A,
  Z = Z,
  name = "gen_lorentzian"
)

ring = j0j0r::maxwellian_ring_setup(
  n = n,
  v_width = 1e5,
  v_rad = 0.5e6,
  A = A,
  Z = Z,
  name = "ring"
)


center <- c(0, 0.5e6)
T_eV <- 0.1e3
v_term <- find_p_term(T_eV, A) / (A * const[["amu"]])
covariance <- rbind(c(v_term^2, 0), c(0, (v_term)^2))

bvtnorm <- j0j0r::bvtnorm_setup(
  n = 4e19,
  center = center,
  covariance = covariance,
  A = A,
  Z = Z,
  name = "bivariate_normal"
)

wilkie = j0j0r::wilkie_setup(
  b = 3,
  n = n,
  A = A,
  Z = Z,
  birth_energy = 50e3,#3.52e6,
  n_e = n,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = c(1),
    A = c(2),
    n = n
  ),
  name = "deuterium_wilkie"
)

dist_plot <- plot_distribution(
  particles = list(
    maxwellian = maxwellian,
    bimaxwellian = bimaxwellian,
    lorentzian = lorentzian,
    ring = ring,
    bvtnorm = bvtnorm,
    wilkie = wilkie
  ),
  v_par = seq(-5e6, 5e6, length.out = 10000),
  v_perp = 0,
  logscale = TRUE
)
dist_plot + ggplot2::scale_y_log10(limits = c(1e78, 1e84))


tic()
specdf <- j0j0r::j0j0_spectrum(
  k = k,
  phi = phi,
  frequencies = seq(0, 400e6, length.out = 41),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    maxwellian = maxwellian,
    bimaxwellian = bimaxwellian,
    lorentzian = lorentzian,
    ring = ring,
    bvtnorm = bvtnorm,
    wilkie = wilkie
  )
)
toc()

specdf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = frequency, y = Re(j0j0) + Im(j0j0), color = particle)) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::facet_wrap( ~ directions)
