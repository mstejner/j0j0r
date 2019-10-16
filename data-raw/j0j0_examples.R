#script to generate j0j0_examples, distribution_examples

phi <- 70
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
  T_eV_par = 1000,
  v_drift = 2e5,
  A = A,
  Z = Z,
  name =  "bimaxwellian"
)

lorentzian <- j0j0r::generalized_lorentzian_setup(
  n = n,
  T_eV = 2000,
  kp = 10,
  A = A,
  Z = Z,
  name = "lorentzian"
)

ring = j0j0r::maxwellian_ring_setup(
  n = n,
  v_width = 1.5e5,
  v_rad = 0.5e6,
  A = A,
  Z = Z,
  name = "ring"
)


center <- c(5e5, 3e5)
T_eV <- 2e3
v_term <- j0j0r::find_p_term(T_eV, A) / (A * j0j0r::const[["amu"]])
covariance <- rbind(c(v_term^2/4, v_term^2/3), c(v_term^2/5, (v_term)^2))

bvtnorm <- j0j0r::bvtnorm_setup(
  n = 4e19,
  center = center,
  covariance = covariance,
  A = A,
  Z = Z,
  name = "bivariate_normal"
)


j0j0r::calculate_distribution_data_frame(
  particles = list(
    bvtnorm = bvtnorm
  ),
  v_par = seq(-1e6, 1e6, length.out = 300),
  v_perp = seq(0, 1e6, length.out = 300)
) %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(y = v_par, x = v_perp, z = value)
  ) +
  ggplot2::geom_raster(ggplot2::aes(fill = value)) +
  ggplot2::geom_contour(colour = "white", alpha = 0.2) +
  viridis::scale_fill_viridis(option = "plasma") +
  ggplot2::ylab(unname(latex2exp::TeX("$v_{par}$ in m/s^2"))) +
  ggplot2::xlab(latex2exp::TeX("$v_{perp}$ in m/s^2")) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(title = "Density")) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = -45)
  ) +
  ggplot2::scale_x_continuous(labels = scales::scientific)

wilkie_b5 = j0j0r::wilkie_setup(
  b = 5,
  n = n,
  A = A,
  Z = Z,
  birth_energy = 20e3,
  n_e = n,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = c(1),
    A = c(2),
    n = n
  ),
  name = "deuterium_wilkie"
)

wilkie_b0 = j0j0r::wilkie_setup(
  b = 0,
  n = n,
  A = A,
  Z = Z,
  birth_energy = 20e3,
  n_e = n,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = c(1),
    A = c(2),
    n = n
  ),
  name = "deuterium_wilkie"
)

distribution_examples <- j0j0r::calculate_distribution_data_frame(
  particles = list(
    maxwellian = maxwellian,
    bimaxwellian = bimaxwellian,
    lorentzian = lorentzian,
    ring = ring,
    bvtnorm = bvtnorm,
    wilkie_b5 = wilkie_b5,
    wilkie_b0 = wilkie_b0
  ),
  v_par = seq(-2.5e6, 2.5e6, length.out = 1000),
  v_perp = 0
)

ggplot2::ggplot( data = distribution_examples,
  mapping = ggplot2::aes(x = v_par, y = value, color = name)
) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::ylab("Density") +
  ggplot2::xlab(latex2exp::TeX("$v_{par}$ in m/s^2")) +
  ggplot2::theme(text = ggplot2::element_text(size = 17))



j0j0_examples <- j0j0r::j0j0(
  k = k,
  phi = phi,
  frequencies = seq(0, 600e6, by = 2e6),
  directions = c("x", "y", "z"),
  B = B,
  particles = list(
    maxwellian = maxwellian,
    bimaxwellian = bimaxwellian,
    lorentzian = lorentzian,
    ring = ring,
    bvtnorm = bvtnorm,
    wilkie_b0 = wilkie_b0,
    wilkie_b5 = wilkie_b5
  ),
  integration_method = "hcubature"
)



j0j0r::j0j0_element(
  directions = "xx",
  k = k,
  phi = phi,
  frequency = 600e6,
  B = B,
  A = A,
  Z = Z,
  distribution = lorentzian$distribution,
  integration_method = "hcubature"
)

usethis::use_data(j0j0_examples, distribution_examples, overwrite = TRUE)
