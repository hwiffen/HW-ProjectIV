library(ggplot2)
library(dplyr)
remotes::install_github("jonocarroll/ggeasy")
library(ggeasy)

#######################################################################################

### Local linear regression estimator:

gauss.kernel<-function(u){
  return(dnorm(u))
}

s_nj <- function(x.data, x, h, j){
  return(sum(gauss.kernel((x.data - x) / h)*(x.data - x)^j))
}

v_i <- function(x.data, x, h){
  s_n1 <- s_nj(x.data, x, h, 1)
  s_n2 <- s_nj(x.data, x, h, 2)
  return(gauss.kernel((x.data - x) / h)*(s_n2 - (x.data - x)*s_n1))
}

m_ll <- function(x.data, y.data, x.grid=x.data, h){
  p <- length(x.grid)
  m.est <- rep(0, p)
  for (i in 1:p){
    m.est[i] <- sum(y.data*v_i(x.data, x.grid[i], h)) / sum(v_i(x.data, x.grid[i], h))
  }
  return(as.data.frame(cbind(x.grid, m.est)))
}
#######################################################################################

## Read in and prepare raw particle spread data:
df <- read.csv("Stokes_2000_1000_particle_widths.csv") %>%
  mutate(avg_width = rowMeans(select(., x_width, y_width))) %>%
  arrange(z)

x = df$avg_width
y = df$z

m.hat_ll <- m_ll(x, y, h=0.1)

df$label <- "Raw data"

## Run local linear regression estimator on raw data:
df_fit <- m.hat_ll %>%
  rename(avg_width = x.grid, z = m.est) %>%
  mutate(type = "Local linear regression estimator")

## Create plot with raw data, LLRE and splitting point:
p <- ggplot() +
  geom_path(data = df, aes(x = avg_width, y = z, linetype = label, color = label), linewidth = 0.5) +
  geom_line(data = df_fit, aes(x = avg_width, y = z, linetype = type, color = type), linewidth = 0.5) +
  geom_hline(aes(yintercept = -420, linetype = "Splitting point", color = "Splitting point"), linewidth = 0.5) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Raw data" = "#00FFFF",
      "Local linear regression estimator" = "#FF00FF",
      "Splitting point" = "black"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Raw data" = "solid",
      "Local linear regression estimator" = "longdash",
      "Splitting point" = "dotted"
    )
  ) +
  scale_y_continuous(trans = "identity") +
  labs(
    title = bquote("Particle spread vs vertical height for " ~ N[0] == 2000),
    x = "Average width",
    y = "z"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  ggeasy::easy_center_title()

## Show plot:
print(p)


