rm(list=ls())
library(GSBart)

# NYC Education
load("data/NYEducation.Rdata")
NYEdu_Time=Sys.time()
NYEdu_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, 
                   hyperpar = hyperpar, verbose = F, nthreads = 1, seed = 1234) 
NYEdu_Time=difftime(Sys.time(), NYEdu_Time, units = "secs")
NYEdu_MSPE = mean((NYEdu_GSB$yhat.test.mean - y.test.unstandardized)^2)
NYEdu_MAPE = mean(abs(NYEdu_GSB$yhat.test.mean - y.test.unstandardized))


# King County House
load("data/KingHouse.Rdata")
KingHouse_Time=Sys.time()
KingHouse_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, verbose = F, 
                       nthreads = 1, seed = 1234)
KingHouse_Time=difftime(Sys.time(), KingHouse_Time, units = "secs")
KingHouse_MSPE = mean((KingHouse_GSB$yhat.test.mean - y.test.unstandardized)^2)
KingHouse_MAPE = mean(abs(KingHouse_GSB$yhat.test.mean - y.test.unstandardized))

load('data/KingHousePDP.Rdata')
p4 <- ggplot(sf_bnd) + 
  geom_sf(color = 'black', fill = NA, lwd = 0.75) + 
  geom_sf(data = mesh_reg$mesh, aes(fill=y.test.unstandardized), color = NA, alpha =0.8) +
  scale_fill_viridis_c(option = "H", direction = 1, name = 'log mean house price')+
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


# US Election
load("data/USelection.Rdata")
USElection_Time=Sys.time()
USElection_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, verbose = T, 
                       nthreads = 1, seed = 1234)
USElection_Time=difftime(Sys.time(), USElection_Time, units = "secs")
USElection_MSPE = mean((USElection_GSB$yhat.test.mean - y.test.unstandardized)^2)
USElection_MAPE = mean(abs(USElection_GSB$yhat.test.mean - y.test.unstandardized))

real_continuous_res = data.frame(MSPE = c(NYEdu_MSPE, KingHouse_MSPE, USElection_MSPE),
                                 MAPE = c(NYEdu_MAPE, KingHouse_MAPE, USElection_MAPE))
rownames(real_continuous_res) <- c('NYC Education','KingHouse', 'US Elecation')

save(p4, p5, p6, real_continuous_res, file = 'data/real_continuous_res.RData', compress = 'xz')
