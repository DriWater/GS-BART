#### Reproducible code for "GS-BART: Bayesian Additive Regression Trees with Graph-split Decision Rules 
#### for Generalized Spatial Nonparametric Regressions"
#### Code for generating figures and tables in Simulation Studies (Section 4)

library(ggplot2)
library(ggpubr)
library(sf)

### GS-BART model results ======================================================

load('data/sim_continuous_res.Rdata')

load('data/sim_count_res.Rdata')

load('data/sim_classification_res.Rdata')

table_1 = data.frame('U-SHAPE' = c(Ushape_continuous_summary$MSPE_mean, Ushape_continuous_summary$MSPE_sd),
                     'TORUS' = c(Torus_continuous_summary$MSPE_mean, Torus_continuous_summary$MSPE_sd),
                     'FRIEDMAN' = c(Friedman_continuous_summary$MSPE_mean, Friedman_continuous_summary$MSPE_sd))

rownames(table_1) <- c('MSE', 'sd')

table_2 = data.frame('U-SHAPE' = c(Ushape_count_summary$RMSPE_mean, Ushape_count_summary$RMSPE_sd),
                     'TORUS' = c(Torus_count_summary$RMSPE_mean, Torus_count_summary$RMSPE_sd))

rownames(table_2) <- c('RMSE', 'sd')

table_3 = data.frame('U-SHAPE' = c(Ushape_class_summary$ACC_mean, Ushape_class_summary$ACC_sd),
                     'TORUS' = c(Torus_class_summary$ACC_mean, Torus_class_summary$ACC_sd))

rownames(table_3) <- c('ACC', 'sd')

load('data/UshapeVis.Rdata')
range_all <- range(c(sim$f_ho_true, GSBart_fit$yhat.test.unstandardized))
range_len = range_all[2] - range_all[1]
true.range.lower = (range(sim$f_ho_true)[1] - range_all[1])/range_len
true.range.upper = (range(sim$f_ho_true)[2] - range_all[1])/range_len
GSBart.range.lower = (range(GSBart_fit$yhat.test.unstandardized)[1] - range_all[1])/range_len
GSBart.range.upper = (range(GSBart_fit$yhat.test.unstandardized)[2] - range_all[1])/range_len

p1 <-  ggplot(data = sim$sf_cluster[3,2]) +
  geom_sf(data = mesh_ho$mesh, aes(fill=sim$f_ho_true),inherit.aes = FALSE, 
          color = NA,  lwd = 0.01) +
  scale_fill_viridis_c(option = "H", direction = 1, name = "",
                       begin = true.range.lower, end= true.range.upper)+
  geom_sf(color = '#CC79A7', fill = NA, lwd = .1) +
  geom_sf(data = sf_bnd, color = 'black', fill = NA, lwd = 1) +
  theme_bw() + theme(
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  labs(x = bquote(S[h]), y = bquote(S[v]), title = "True f")

p2 <- ggplot(mesh_ho$mesh) + 
  geom_sf(aes(fill=GSBart_fit$yhat.test.unstandardized),inherit.aes = FALSE, 
          color = NA,  lwd = 0.01) +
  scale_fill_viridis_c(option = "H", direction = 1, name = "",
                       begin = GSBart.range.lower, end= GSBart.range.upper)+
  geom_sf(data = sf_bnd, color = 'black', fill = NA, lwd = 1) +
  theme_bw() + theme(
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  labs(x = bquote(S[h]), y = bquote(S[v]), title = "GS-BART")

figure_1 = ggarrange(p1, p2, nrow = 1, ncol = 2, widths = c(5,5), heights = c(6), common.legend = TRUE, legend="right")

thres = .9
GSBart_poster_sd_tmp <- GSBart_poster_sd
GSBart_poster_sd_tmp[GSBart_poster_sd>thres] <- NA
figure_2 = ggplot(sim$sf_cluster[3,2]) + 
  geom_sf(data = mesh_ho$mesh, aes(fill=GSBart_poster_sd_tmp),inherit.aes = FALSE, 
          color = NA,  lwd = 0.02) +
  scale_fill_gradientn(colours = cm.colors(10),na.value = "grey50", name = "") +
  geom_sf(color = '#6495ED', fill = NA, lwd = .3) +
  geom_sf(data = sf_bnd, color = 'black', fill = NA, lwd = 1) +
  theme_bw() + theme(
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  labs(x = bquote(S[h]), y = bquote(S[v]), title = "GS-BART")

GSBart_abs_err <- GSBart_fit$abs_err
GSBart_abs_err_tmp <- GSBart_abs_err
GSBart_abs_err_tmp[GSBart_abs_err > 1] <- NA
figure_3 = ggplot(sim$sf_cluster[3,2]) + 
  geom_sf(data = mesh_ho$mesh, aes(fill=GSBart_abs_err_tmp),inherit.aes = FALSE, 
          color = NA,  lwd = 0.01) +
  scale_fill_gradientn(colours = cm.colors(10),na.value = "grey50", name = "") +
  geom_sf(color = '#6495ED', fill = NA, lwd = .1) +
  geom_sf(data = sf_bnd, color = 'black', fill = NA, lwd = 1) +
  theme_bw() + theme(
    axis.text.x = element_text(hjust = 0.5, size  = 7),
    axis.text.y = element_text(hjust = 0.5, size  = 7),
    plot.title = element_text(hjust = 0.5, size  = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  labs(x = bquote(S[h]), y = bquote(S[v]), title = "GS-BART")