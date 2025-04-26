library(ggplot2)
library(RColorBrewer)
remotes::install_github("jonocarroll/ggeasy")
library(ggeasy)

## Define parameters:
alpha <- 1 # permeability parameter
delta <- 0.05 # blob width parameter

## Create plotting grid:
grid_size <- 300
x <- seq(-0.015, 0.015, length.out = grid_size)
y <- seq(-0.015, 0.015, length.out = grid_size)
grid <- expand.grid(x = x, y = y)
grid$r <- sqrt(grid$x^2 + grid$y^2)
r <- grid$r

## Calculate exact Laplacian of H2:
laplacian_H2_exact <- function(r, alpha, delta) {
  R2 <- r^2 + delta^2
  sqrt_R2 <- sqrt(R2)
  exp_term <- exp(-alpha * sqrt_R2)
  one_minus_exp <- 1 - exp_term
  pi4 <- 4 * pi
  
  term1 <- 3 * alpha * exp_term / (pi4 * R2^2)
  term2 <- -alpha^2 * r^2 * exp_term / (pi4 * R2^(5/2))
  term3 <- 9 * exp_term / (2 * pi * R2^(5/2))
  term4 <- -5 * alpha * r^2 * exp_term / (2 * pi * R2^3)
  term5 <- 45 * exp_term / (pi4 * alpha * R2^3)
  term6 <- -45 * r^2 * exp_term / (pi4 * R2^(7/2))
  term7 <- -105 * r^2 * exp_term / (pi4 * alpha * R2^4)
  term8 <- -45 * one_minus_exp / (pi4 * alpha^2 * R2^(7/2))
  term9 <- 105 * r^2 * one_minus_exp / (pi4 * alpha^2 * R2^(9/2))
  
  return(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9)
}

## Define blob function:
phi_delta <- function(r, delta, alpha) {
  R2 <- r^2 + delta^2
  blob <- (3 * delta^2) / (4 * pi * R2^(5/2))
  lapH2 <- laplacian_H2_exact(r, alpha, delta)
  return(blob - delta^2 * lapH2)
}

## Calculate blob function for specified parameters:
grid$phi <- phi_delta(grid$r, delta, alpha)

## Plot as a circular heat map:
ggplot(grid, aes(x = x, y = y, fill = phi)) +
  geom_raster(interpolate = TRUE) +
  coord_fixed() +
  scale_fill_distiller(palette = "PuRd", direction = 1, name = expression(phi[delta](r))) +
  labs(
    title = expression("Blob function for" ~delta~ 0.05),
    x = "x", y = "y"
  ) +
  theme_minimal(base_size = 14) +
  ggeasy::easy_center_title() 



