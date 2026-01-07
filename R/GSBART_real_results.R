#### Reproducible code for "GS-BART: Bayesian Additive Regression Trees with Graph-split Decision Rules 
#### for Generalized Spatial Nonparametric Regressions"
#### Code for generating figures and tables in Real data analysis (Section 4)

library(ggplot2)
library(ggpubr)
library(sf)

### GS-BART model results ======================================================

load('data/real_continuous_res.Rdata')

load('data/real_count_res.Rdata')

load('data/real_classification_res.Rdata')

table_4 = real_continuous_res

table_5 = real_count_res

table_6 = real_classification_res

load('data/KingHousePDP.Rdata')
figure_4 =  ggplot(sf_bnd) + 
  geom_sf(color = 'black', fill = NA, lwd = 0.75) + 
  geom_sf(data = mesh_reg$mesh, aes(fill=y.test.unstandardized), color = NA, alpha =0.8) +
  scale_fill_viridis_c(option = "H", direction = 1, name = '')+
  theme_bw() + xlim(c(122.44,121.82)) + ylim(c(47.23,47.78))


p5 <-ggplot(sf_bnd) + 
  geom_sf(color = 'black', fill = NA, lwd = 0.8) + 
  geom_point(data = coords_ho, aes(x = long, y = lat, color = group), inherit.aes = FALSE, size = 3) +
  scale_colour_hue(guide = "none") + 
  theme_bw() + theme(
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "Longitude", y = "Latitude")


p6 <- ggplot(data = df_new, aes(x = x, y = y, group = group, fill = group, color = group, ymin=CI_low, ymax=CI_high)) +
  geom_line() + scale_colour_hue() + scale_fill_hue() + 
  geom_ribbon(alpha=0.3, linewidth = 0) +
  theme_bw() + theme(
    legend.position="none",
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "living area (square feet)", y = "Predicted log sale price")

figure_5 = ggarrange(p5, p6, nrow = 1, ncol = 2, widths = c(5,5), heights = c(7))

