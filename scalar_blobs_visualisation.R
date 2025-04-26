library(ggplot2)
library(dplyr)
library(ggeasy)
remotes::install_github("jonocarroll/ggeasy")

## Create a dataframe of Gaussian blobs:
x_vals <- seq(-0.1, 0.1, length.out = 1000)
deltas <- c(0.01, 0.02, 0.05)
label_data_blob <- data.frame(
  label = factor(deltas[1:3], levels = deltas),
  x = c(0.012, 0.027, 0.085), 
  y = c(27, 12, 3.5)   
)

df <- expand.grid(x = x_vals, delta = deltas) %>%
  mutate(gaussian = (1 / (sqrt(2 * pi) * delta)) * exp(-x^2 / (2 * delta^2)),
         delta = factor(delta))

## Base plot:
p <- ggplot(df, aes(x = x, y = gaussian, color = delta)) +
  geom_line(size = 1.2) +
  scale_color_brewer(type="seq",palette = "PuRd", name = expression(delta), direction=-1) +
  labs(title = "Scalar blob functions with the Dirac delta function",
       x = "x", y = expression(phi[delta](x))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  ggeasy::easy_center_title() +
  coord_cartesian(ylim = c(0, 45))

## Add Dirac delta measure for comparison:
p + geom_segment(aes(x = 0, xend = 0, y = 0, yend = 45),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "red", linetype = "dashed", size = 1.2) +
  annotate("text", x = -0.01, y = 43, label = expression(delta(x)), , size = 5, color='red')+
  geom_text(data = label_data_blob,
            aes(x = x, y = y, label = label, color = label),
            hjust = 0,
            size = 4,
            show.legend = FALSE)



