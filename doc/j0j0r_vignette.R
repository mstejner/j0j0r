## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = TRUE, results = FALSE, warning = FALSE, message = FALSE-----
library(j0j0r)

## ---- echo = TRUE, results = FALSE, eval = FALSE-------------------------
#  future::plan(strategy = "multisession")

## ---- echo = TRUE, results = FALSE---------------------------------------
maxwellian_deuterium <- maxwellian_setup(
  n = 4e19,
  T_eV = 2000,
  A = 2,
  Z = 1,
  name = "maxwellian"
)

slowdown_deuterium <- slowdown_setup(
  b = 3,
  n = 4e19,
  A = 2,
  Z = 1,
  birth_energy = 10e3,
  n_e = 4e19,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = 1,
    A = 2,
    n = 4e19
  ),
  name = "slowdown"
)

## ----fig1, echo = TRUE, fig.align = "center"-----------------------------
dist_examples <- calculate_distribution_data_frame(
  particles = list(
    maxwellian = maxwellian_deuterium, 
    slowdown = slowdown_deuterium
    ), 
  v_par = seq(-1.5e6, 1.5e6, length.out = 1000), 
  v_perp = 0
) 

plot_dist(dist_examples)

## ---- echo = TRUE, fig.align = "center", eval = TRUE---------------------
calculate_distribution_data_frame(
  particles = list(
    maxwellian = maxwellian_deuterium, 
    slowdown = slowdown_deuterium
    ), 
  v_par = seq(-1e6, 1e6, length.out = 300), 
  v_perp = seq(0, 1e6, length.out = 300) 
  ) %>% 
  plot_dist() +
  ggplot2::facet_wrap( ~ name)

## ---- echo = TRUE, results = "asis"--------------------------------------
integrate_distribution(distribution = maxwellian_deuterium[["distribution"]])

## ---- echo = TRUE, results = FALSE, eval = FALSE-------------------------
#  future::plan(strategy = "multisession")

## ---- echo = TRUE, results = FALSE, eval = FALSE-------------------------
#  
#  maxwellian_example <-j0j0(
#    k = 2 * pi / (j0j0r::const$c / 100e9),
#    phi = c(60, 86),
#    frequencies = seq(0, 400e6, by = 2e6),
#    directions = c("x", "y", "z"),
#    B = 2.5,
#    particles = list(
#      maxwellian = maxwellian_deuterium
#    ),
#    integration_method = "stats"
#  )

## ----fig3, echo = TRUE, fig.width = 16, fig.height = 8, results = FALSE, fig.align = "center"----
plot_j0j0(maxwellian_example, wrap_by = "element", color_by = "phi")

## ---- echo = TRUE--------------------------------------------------------
bimaxwellian <-  bimaxwellian_setup(
  n = 4e9,
  T_eV_perp = 2000,
  T_eV_par = 1000,
  v_drift = 2e5,
  A = 2,
  Z = 1,
  name =  "bimaxwellian"
)

## ----fig4, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
pp1 <- j0j0r::distribution_examples %>% 
  dplyr::filter(name %in% c("maxwellian", "bimaxwellian")) %>% 
  plot_dist()

pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("bimaxwellian")) %>% 
  plot_dist() +
  ggplot2::ggtitle("bimaxwellian")

gridExtra::grid.arrange(pp1, pp2, ncol = 2)


## ----fig5, echo = FALSE, fig.width = 16, fig.height = 8, results = FALSE, fig.align = "center"----
j0j0r::j0j0_examples %>%
  dplyr::filter(particle %in% c("maxwellian", "bimaxwellian")) %>% 
  plot_j0j0(wrap_by = "element", color_by = "particle")

## ---- echo = TRUE--------------------------------------------------------
lorentzian <- generalized_lorentzian_setup(
  n = 4e19,
  T_eV = 2000,
  kp = 10,
  A = 2,
  Z = 1,
  name = "lorentzian"
)

## ----fig6, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
pp1 <- j0j0r::distribution_examples %>% 
  dplyr::filter(name %in% c("maxwellian", "lorentzian")) %>% 
  plot_dist() + ggplot2::scale_y_log10()

pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("lorentzian")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Generalized Lorentzian")

gridExtra::grid.arrange(pp1, pp2, ncol = 2)


## ----fig7, echo = FALSE, fig.width = 16, fig.height = 8, results = FALSE, fig.align = "center"----
j0j0r::j0j0_examples %>% 
  dplyr::filter(particle %in% c("maxwellian", "lorentzian")) %>% 
  plot_j0j0(wrap_by = "element", color_by = "particle")

## ---- echo = TRUE--------------------------------------------------------
ring = maxwellian_ring_setup(
  n = 4e19,
  v_width = 1.5e5,
  v_rad = 0.5e6,
  A = 2,
  Z = 1,
  name = "ring"
)

## ----fig8, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
pp1 <- j0j0r::distribution_examples %>% 
  dplyr::filter(name %in% c("maxwellian", "ring")) %>% 
  plot_dist()

pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("ring")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Ring distribution")

gridExtra::grid.arrange(pp1, pp2, ncol = 2)


## ----fig9, echo = FALSE, fig.width = 16, fig.height = 8, results = FALSE, fig.align = "center"----
j0j0r::j0j0_examples %>% 
  dplyr::filter(particle %in% c("maxwellian", "ring")) %>% 
  plot_j0j0(wrap_by = "element", color_by = "particle")

## ---- echo = TRUE--------------------------------------------------------
center <- c(5e5, 3e5)
T_eV <- 2e3
A <- 2
v_term <- find_p_term(T_eV, A) / (A * j0j0r::const[["amu"]])
covariance <- rbind(c(v_term^2/4, v_term^2/3), c(v_term^2/5, (v_term)^2))

bvtnorm <- bvtnorm_setup(
  n = 4e19,
  center = center,
  covariance = covariance,
  A = 2,
  Z = 1,
  name = "bivariate_normal"
)

## ----fig10, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
pp1 <- j0j0r::distribution_examples %>% 
  dplyr::filter(name %in% c("maxwellian", "bvtnorm")) %>% 
  plot_dist() + 
  ggplot2::scale_y_log10()

pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("bvtnorm")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Bivariate normal distribution")

gridExtra::grid.arrange(pp1, pp2, ncol = 2)


## ----fig11, echo = FALSE, fig.width = 16, fig.height = 8, results = FALSE,fig.align = "center"----
j0j0r::j0j0_examples %>% 
  dplyr::filter(particle %in% c("maxwellian", "bvtnorm")) %>% 
  plot_j0j0(wrap_by = "element", color_by = "particle")

## ---- echo = TRUE--------------------------------------------------------
n <- 4e19
A <- 2
Z <- 1

slowdown_b5 = slowdown_setup(
  b = 5,
  n = n,
  A = A,
  Z = Z,
  birth_energy = 20e3,
  n_e = n,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = Z,
    A = A,
    n = n
  ),
  name = "deuterium_slowdown_b5"
)

slowdown_b0 = slowdown_setup(
  b = 0,
  n = n,
  A = A,
  Z = Z,
  birth_energy = 20e3,
  n_e = n,
  T_e_eV = 2e3,
  ions = data.frame(
    Z = Z,
    A = A,
    n = n
  ),
  name = "deuterium_slowdown_b0"
)

## ----fig12, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
j0j0r::distribution_examples %>% 
  dplyr::filter(name %in% c("maxwellian", "slowdown_b0", "slowdown_b5")) %>% 
  plot_dist() 

pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("slowdown_b0")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Slowdown distribution, b = 0")

pp3 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("slowdown_b5")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Slowdown distribution, b = 5")

gridExtra::grid.arrange(pp2, pp3, ncol = 2)


## ----fig13, echo = FALSE, fig.width = 16, fig.height = 8, results = FALSE, fig.align = "center"----
j0j0r::j0j0_examples %>% 
  dplyr::filter(particle %in% c("maxwellian", "slowdown_b0", "slowdown_b5")) %>% 
  plot_j0j0(wrap_by = "element", color_by = "particle")

## ----echo = FALSE--------------------------------------------------------
A <- 2
Z <- 1
mass <- A * j0j0r::const$amu
charge <- Z * j0j0r::const$qe
B <- 2.5
omega_c <- charge * B / mass

omega = 2 * pi * 115e6
phi = 86
k <- 2 * pi / (j0j0r::const$c / 100e9)
k_par <- cos(phi * pi / 180) * k

l_values <- seq(-20, 20)

v_parallel_86 <- (omega - l_values * omega_c) / k_par

phi = 80
k_par <- cos(phi * pi / 180) * k
v_parallel_80 <- (omega - l_values * omega_c) / k_par


## ----fig14, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
pp2 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("maxwellian")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Maxwellian") +
  ggplot2::geom_hline(
    yintercept = v_parallel_86[abs(v_parallel_86) < 1e6], 
    linetype = "dashed", 
    color = "black",
    size = 2.2
    ) +
    ggplot2::geom_hline(
    yintercept = v_parallel_80[abs(v_parallel_80) < 1e6], 
    linetype = "dashed", 
    color = "white",
    size = 1.5
    )

pp3 <- dist_2D_examples %>% 
  dplyr::filter(name %in% c("slowdown_b5")) %>% 
  plot_dist() +
  ggplot2::ggtitle("Slowdown distribution, b = 5") +
    ggplot2::geom_hline(
    yintercept = v_parallel_86[abs(v_parallel_86) < 1.5e6], 
    linetype = "dashed", 
    color = "black",
    size = 2.2
    ) +
    ggplot2::geom_hline(
    yintercept = v_parallel_80[abs(v_parallel_80) < 1.5e6], 
    linetype = "dashed", 
    color = "white",
    size = 1.5
    )

gridExtra::grid.arrange(pp2, pp3, ncol = 2)

## ----fig15, echo = FALSE, fig.width = 16, fig.height = 8, fig.align = "center"----
mass <- 2 * j0j0r::const$amu
charge <- 1 * j0j0r::const$qe

omega_c <- charge * B / mass

omega = 2 * pi * 115e6
phi = 86
k_par <- cos(phi * pi / 180) * k
k_perp <- sin(phi * pi / 180) * k

df <- expand.grid(
  l = seq(5,7), 
  v_perp = seq(1e4, 1e6, 1e2), 
  direction = c("x", "y", "z"), 
  stringsAsFactors = FALSE
  )

df <- df %>% dplyr::mutate(
  k_perp = k_perp, 
  k_par = k_par,
  omega_c = omega_c,
  v_par = (omega - l * omega_c) / k_par
)


df[["clcl"]] <- unlist(purrr::pmap(.l = df, .f = cl))^2

df[["factor"]] <- df[["clcl"]] * df[["v_perp"]] * mass
df[["l"]] <- as.factor(df[["l"]])
df[["direction"]] <-  paste0(df[["direction"]], df[["direction"]])

ggplot2::ggplot(
  data = df, 
  mapping = ggplot2::aes(x = v_perp, y = factor, color = l)
  ) +
  ggplot2::geom_line(size = 1.5) +
  ggplot2::scale_y_log10() +
  ggplot2::facet_wrap( ~ direction) +
  ggplot2::theme(
      legend.position = "top",
      text = ggplot2::element_text(size = 25),
      axis.text.x = ggplot2::element_text(angle = -45, size = 20)
      ) +
  ggplot2::ylab(unname(latex2exp::TeX("abs(p_{perp}  c_l c_l^*)"))) +
  ggplot2::xlab(latex2exp::TeX("$v_{perp}$ in m/s^2")) +
  ggplot2::scale_x_continuous(labels = scales::scientific)

## ----fig16, echo = FALSE, fig.width = 16, fig.height = 16, fig.align = "center"----
phi_sweep %>% 
  plot_j0j0(wrap_by = "particle", color_by = "phi") + 
  ggplot2::facet_wrap( ~particle, ncol = 1, scales = "free") +
  ggplot2::coord_cartesian(xlim = c(0,500)) + 
  ggplot2::labs(color = expression(phi))

