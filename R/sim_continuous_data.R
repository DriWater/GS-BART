rm(list=ls())
library(bayestestR)
library(GSBart)
library(purrr)

load("data/sim_input.RData")

repetitions = 50

## Ushape Simulation
Ushape_continuous_list = vector(50, mode="list")
sigma = 0.15
n = nrow(sim_Ushape$X); n_ho = nrow(sim_Ushape$X_ho)

for(repetition in 1:repetitions){
  set.seed(1234+repetition)
  Y0 = sim_Ushape$f_true + rnorm(n, 0, sigma)
  Y0_ho = sim_Ushape$f_ho_true + rnorm(n_ho, 0, sigma)
  GSBART_Fit = GSBart::gsbart(Y0, sim_Ushape$Graphs, 200, 15, 200, sim_Ushape$graphs_weight, nthreads = 7, verbose = T, seed = 1234)
  GSBart_MSPE = mean((Y0_ho - GSBART_Fit$yhat.test.mean)^2)
  GSBart_MAPE = mean(abs(Y0_ho - GSBART_Fit$yhat.test.mean))
  GSBart.upper.lvl = NULL
  GSBart.lower.lvl = NULL
  for(i in 1:ncol(GSBres$yhat.test)){
    tmp <- ci(GSBres$yhat.test[,i], method = "HDI")
    GSBart.upper.lvl <- append(GSBart.upper.lvl, tmp$CI_high)
    GSBart.lower.lvl <- append(GSBart.lower.lvl, tmp$CI_low)
  }
  GSBart_coverage = mean((Y0_ho<GSBart.upper.lvl)&(Y0_ho>GSBart.lower.lvl))
  GSBart_HDI_len = mean(GSBart.upper.lvl - GSBart.lower.lvl)
  GSBart_poster_sd = mean(apply(GSBres$yhat.test, 2, sd))
  Ushape_continuous_list[[i]] = data.frame(
    models = c("GSBART"),
    MSPE = c(GSBart_MSPE),
    MAPE = c(GSBart_MAPE),
    times = c(GSBart_Time),
    HDI_coverage = c(GSBart_coverage),
    HDI_len = c(GSBart_HDI_len),
    poster_sd = c(GSBart_poster_sd),
    sim = c(repetition)
  )
}

Ushape_continuous_summary = Ushape_continuous_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(MSPE, MAPE, HDI_coverage, HDI_len),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(poster_sd, times),
           mean,
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))

Ushape_continuous_summary

load('data/UshapeVis.RData')
range_all <- range(c(sim$f_ho_true, GSBart_fit$yhat.test.unstandardized))
range_len = range_all[2] - range_all[1]
true.range.lower = (range(f_ho_true)[1] - range_all[1])/range_len
true.range.upper = (range(f_ho_true)[2] - range_all[1])/range_len
GSBart.range.lower = (range(GSBart_fit$yhat.test.unstandardized)[1] - range_all[1])/range_len
GSBart.range.upper = (range(GSBart_fit$yhat.test.unstandardized)[2] - range_all[1])/range_len

p1 <-  ggplot(data = sf_cluster[3,2]) + 
  geom_sf(data = mesh_ho$mesh, aes(fill=f_ho_true),inherit.aes = FALSE, 
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

GSBart_abs_err <- GSBart_fit$abs_err
range_all <- range(c(GSBart_abs_err))
range_len = range_all[2] - range_all[1]
GSBart.range.lower = 1 - (range(GSBart_abs_err)[1] - range_all[1])/range_len 
GSBart.range.upper = 1- (range(GSBart_abs_err)[2] - range_all[1])/range_len 

GSBart_abs_err_tmp <- GSBart_abs_err
GSBart_abs_err_tmp[GSBart_abs_err > 1] <- NA
p3 <- ggplot(sf_cluster[3,2]) + 
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

# Torus Simulation
Torus_continuous_list = vector(50, mode="list")
sigma = 0.1
n = nrow(sim_Torus$X); n_ho = nrow(sim_Torus$X_ho)

for(repetition in 1:repetitions){
  set.seed(1234+repetition)
  Y0 = sim_Torus$f_true + rnorm(n, 0, sigma)
  Y0_ho = sim_Torus$f_ho_true + rnorm(n_ho, 0, sigma)
  GSBART_Fit = GSBart::gsbart(Y0, sim_Torus$Graphs, 200, 15, 200, sim_Torus$graphs_weight, nthreads = 7, verbose = T, seed = 1234)
  GSBart_MSPE = mean((Y0_ho - GSBART_Fit$yhat.test.mean)^2)
  GSBart_MAPE = mean(abs(Y0_ho - GSBART_Fit$yhat.test.mean))
  GSBart.upper.lvl = NULL
  GSBart.lower.lvl = NULL
  for(i in 1:ncol(GSBres$yhat.test)){
    tmp <- ci(GSBres$yhat.test[,i], method = "HDI")
    GSBart.upper.lvl <- append(GSBart.upper.lvl, tmp$CI_high)
    GSBart.lower.lvl <- append(GSBart.lower.lvl, tmp$CI_low)
  }
  GSBart_coverage = mean((Y0_ho<GSBart.upper.lvl)&(Y0_ho>GSBart.lower.lvl))
  GSBart_HDI_len = mean(GSBart.upper.lvl - GSBart.lower.lvl)
  GSBart_poster_sd = mean(apply(GSBres$yhat.test, 2, sd))
  Torus_continuous_list[[i]] = data.frame(
    models = c("GSBART"),
    MSPE = c(GSBart_MSPE),
    MAPE = c(GSBart_MAPE),
    times = c(GSBart_Time),
    HDI_coverage = c(GSBart_coverage),
    HDI_len = c(GSBart_HDI_len),
    poster_sd = c(GSBart_poster_sd),
    sim = c(repetition)
  )
}

Torus_continuous_summary = Torus_continuous_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(MSPE, MAPE, HDI_coverage, HDI_len),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(poster_sd, times),
           mean,
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))

Torus_continuous_summary


# Torus Simulation
Friedman_continuous_list = vector(50, mode="list")
sigma = 1; n = nrow(sim_Friedman$X); n_ho = nrow(sim_Friedman$X_ho)

for(repetition in 1:repetitions){
  set.seed(1234+repetition)
  Y0 = sim_Friedman$f_true + rnorm(n, 0, sigma)
  Y0_ho = sim_Friedman$f_ho_true + rnorm(n_ho, 0, sigma)
  GSBART_Fit = GSBart::gsbart(Y0, sim_Friedman$Graphs, 200, 15, 200, sim_Friedman$graphs_weight,
                              nthreads = 1, verbose = T, seed = 1234)
  GSBart_MSPE = mean((Y0_ho - GSBART_Fit$yhat.test.mean)^2)
  GSBart_MAPE = mean(abs(Y0_ho - GSBART_Fit$yhat.test.mean))
  GSBart.upper.lvl = NULL
  GSBart.lower.lvl = NULL
  for(i in 1:ncol(GSBres$yhat.test)){
    tmp <- ci(GSBres$yhat.test[,i], method = "HDI")
    GSBart.upper.lvl <- append(GSBart.upper.lvl, tmp$CI_high)
    GSBart.lower.lvl <- append(GSBart.lower.lvl, tmp$CI_low)
  }
  GSBart_coverage = mean((Y0_ho<GSBart.upper.lvl)&(Y0_ho>GSBart.lower.lvl))
  GSBart_HDI_len = mean(GSBart.upper.lvl - GSBart.lower.lvl)
  GSBart_poster_sd = mean(apply(GSBres$yhat.test, 2, sd))
  Friedman_continuous_list[[i]] = data.frame(
    models = c("GSBART"),
    MSPE = c(GSBart_MSPE),
    MAPE = c(GSBart_MAPE),
    times = c(GSBart_Time),
    HDI_coverage = c(GSBart_coverage),
    HDI_len = c(GSBart_HDI_len),
    poster_sd = c(GSBart_poster_sd),
    sim = c(repetition)
  )
}

Friedman_continuous_summary = Friedman_continuous_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(MSPE, MAPE, HDI_coverage, HDI_len),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(poster_sd, times),
           mean,
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))

save(Ushape_continuous_summary, Torus_continuous_summary, Friedman_continuous_summary, 
     p1, p2, p3, file = 'data/sim_continuous.RData', compress = 'xz')