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
winter.colours <- colorRampPalette(c("#0000FF", "#00FFFF", "#00FF80"))

## Define spread analysis .csv files to read in:
files <- c(
  "BM_0.01_2000_500_particle_widths.csv",
  "BM_0.1_2000_500_particle_widths.csv",
  "BM_0.3_2000_500_particle_widths.csv",
  "BM_0.5_2000_500_particle_widths.csv",
  "BM_1_2000_500_particle_widths.csv",
  "BM_3_2000_500_particle_widths.csv",
  "BM_5_2000_500_particle_widths.csv",
  "BM_10_2000_500_particle_widths.csv",
  "BM_15_2000_500_particle_widths.csv"
)
labels <- c("0.01", "0.1", "0.3", "0.5", "1", "3", "5", "10", "15")

## Define custom locations to add plot labels:
label_data_small <- data.frame(
  label = factor(labels[1:5], levels = labels),
  x = c(1.4, 4.1, 7.9, 8.8, 9.5),
  y = c(-470, -450, -425, -329, -222)  
)

label_data_large <- data.frame(
  label = factor(labels[6:9], levels = labels),
  x = c(10, 45, 80, 110),
  y = c(-85, -83, -81, -80)
)
colours <- winter.colours(length(labels))

## Define functions to add local linear regression estimates to plot (small and large alpha):
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
        title = bquote("Particle spread in high"~alpha~"Brinkman simulations"),
        x = "Average width",
        y = "z",
        color = expression(alpha)  
      ) +
      scale_y_continuous(trans = "identity") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top") +
      ggeasy::easy_center_title()
  }
  base_plot +
    geom_line(data = df_fit, aes(x = avg_width, y = z.m.est, color=label) , linetype = 5, linewidth = 0.5)
}

add_ll_fit_to_plot_small <- function(filepath, h = 0.1, color = NULL, label = NULL, base_plot = NULL) {
  df <- read.csv(filepath) %>%
    filter(z >= -500) %>%
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
        title = bquote("Particle spread in low"~alpha~"Brinkman simulations"),
        x = "Average width",
        y = "z",
        color = expression(alpha)  
      ) +
      scale_y_continuous(trans = "identity") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top") +
      ggeasy::easy_center_title()
  }
  base_plot +
    geom_line(data = df_fit, aes(x = avg_width, y = z.m.est, color=label) , linetype = 5, linewidth = 0.5)
}


## Plot for small alpha .csv files:
small_indices <- 1:5
p1 <- add_ll_fit_to_plot_small(files[small_indices[1]], h = 0.8, label = labels[small_indices[1]])
for (i in small_indices[-1]) {
  p1 <- add_ll_fit_to_plot_small(files[i], h = 0.8, label = labels[i], base_plot = p1)
}
p1 <- p1 + 
  scale_color_manual(values = setNames(colours, labels)) +
  scale_x_continuous(limits = c(1, 10.5)) +
  scale_y_continuous(limits = c(-500, 0)) +
  geom_text(data = label_data_small,
            aes(x = x, y = y, label = label, color = label),
            hjust = 0,
            size = 4,
            show.legend = FALSE)

## Plot for large alpha .csv files:
large_indices <- 6:9
p2 <- add_ll_fit_to_plot(files[large_indices[1]], h = 0.5, label = labels[large_indices[1]])
for (i in large_indices[-1]) {
  p2 <- add_ll_fit_to_plot(files[i], h = 0.5, label = labels[i], base_plot = p2)
}
p2 <- p2 + 
  scale_color_manual(values = setNames(colours, labels)) +
  scale_x_continuous(limits = c(1, 127)) +
  scale_y_continuous(limits = c(-90, 0)) +
  geom_text(data = label_data_large,
            aes(x = x, y = y, label = label, color = label),
            hjust = 0,
            size = 4,
            show.legend = FALSE)

## Show the final plots:
p1
p2