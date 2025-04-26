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

## Define colour map to align with Python plots:
autumn.colours <- colorRampPalette(c("#FF0000", "#FFFF00")) 

## Define spread analysis .csv files to read in:
files <- c(
  "Magnetic_0.01_5_500_300_particle_widths.csv",
  "Magnetic_0.1_5_500_300_particle_widths.csv",
  "Magnetic_1_5_500_300_particle_widths.csv",
  "Magnetic_5_5_500_300_particle_widths.csv",
  "Magnetic_3_5_500_300_particle_widths.csv"
)
labels <- c("0.01", "0.1", "1", "3", "5")

## Define custom locations to add plot labels:
label_data_magnetic <- data.frame(
  label = factor(labels[1:5], levels = labels),
  x = c(59, 58, 52.5, 37.5, 31.8),
  y = c(-78, -90, -89, -83, -84)
)

colours <- autumn.colours(length(files))

## Define function to add local linear regression estimates to plot:
add_ll_fit_to_plot <- function(filepath, h = 0.1, color = NULL, label = NULL, base_plot = NULL) {
  df <- read.csv(filepath) %>%
    mutate(avg_width = rowMeans(select(., x_width, y_width))) %>%
    arrange(z)
  
  x <- df$avg_width
  y <- df$z
  m.hat_ll <- m_ll(x, y, h = h)
  df_fit <- data.frame(
    avg_width = df$avg_width,
    z = m.hat_ll
  )
  df_fit$label <- factor(label, levels = labels)
  
  if (is.null(base_plot)) {
    base_plot <- ggplot() +
      labs(
        title = bquote("Particle spread in magentic field simulations for " ~ alpha == 5),
        x = "Average width",
        y = "z",
        color = expression(beta)  
      ) +
      scale_y_continuous(trans = "identity") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top") +
      ggeasy::easy_center_title()
  }
  
  base_plot +
    geom_line(data = df_fit, aes(x = avg_width, y = z.m.est, color=label) , linetype = 5, linewidth = 0.5)
}

## Create plot using first .csv file:
p <- add_ll_fit_to_plot(files[1], h = 1, label = labels[1])

## Overlap other .csv file results onto plot:
for (i in 2:length(files)) {
  p <- add_ll_fit_to_plot(files[i], h = 1, label = labels[i], base_plot = p)
}

## Customise plot appearance:
p <- p + scale_color_manual(values = setNames(colors, labels))

p <- p +
  scale_x_continuous(limits = c(1, 60))  +
  scale_y_continuous(limits = c(-90, 0)) +
  geom_text(data = label_data_magnetic,
            aes(x = x, y = y, label = label, color = label),
            hjust = 0,
            size = 4,
            show.legend = FALSE)

## Show the final plot:
print(p)
