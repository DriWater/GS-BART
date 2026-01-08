rm(list=ls())
library(G2SBart)
library(purrr)
library(dplyr)

load("data/sim_input.Rdata")

repetitions = 50

## U-SHAPE Simulation
Ushape_count_list = vector(50, mode="list")
n = nrow(sim_Ushape$X); n_ho = nrow(sim_Ushape$X_ho)

for(repetition in seq_len(repetitions)){
  set.seed(1234+repetition)
  Y0 = exp((sim_Ushape$f_true+7)/3)
  Y0_ho = exp((sim_Ushape$f_ho_true+7)/3)
  Y <- rpois(n, Y0)
  GSBart_Time = Sys.time() 
  GSBART_Fit <- g2sbart(Y, sim_Ushape$Graphs, 35, 15, 35, sim_Ushape$graphs_weight, 
                        family = 'poisson', nthreads = 6, sparse = F, verbose = F, seed = 1234)
  GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
  GSBart_RMSPE = sqrt(mean((Y0_ho - exp(colMeans(GSBART_Fit$phi.test)))^2))
  GSBart_MSPE = mean((Y0_ho - exp(colMeans(GSBART_Fit$phi.test)))^2)
  Ushape_count_list[[repetition]] = data.frame(
    models = c("GSBART"),
    RMSPE = c(GSBart_RMSPE),
    MSPE = c(GSBart_MSPE),
    times = c(GSBart_Time),
    sim = c(repetition)
  )
}

Ushape_count_summary = Ushape_count_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(RMSPE, MSPE),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(times), mean,
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))


# Torus Simulation
Torus_count_list = vector(50, mode="list")
n = nrow(sim_Torus$X); n_ho = nrow(sim_Torus$X_ho)

for(repetition in seq_len(repetitions)){
  set.seed(1234+repetition)
  Y0 <- exp((sim_Torus$f_true + 17)/8)
  Y0_ho <- exp((sim_Torus$f_ho_true + 17)/8)
  Y <- rpois(n, Y0)
  GSBart_Time = Sys.time() 
  GSBART_Fit <- g2sbart(Y, sim_Torus$Graphs, 35, 15, 35, sim_Torus$graphs_weight,
                        family = 'poisson', nthreads = 7, sparse = F, verbose = F, seed = 1234)
  GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
  GSBart_RMSPE = sqrt(mean((Y0_ho - exp(colMeans(GSBART_Fit$phi.test)))^2))
  GSBart_MSPE = mean((Y0_ho - exp(colMeans(GSBART_Fit$phi.test)))^2)
  Torus_count_list[[repetition]] = data.frame(
    models = c("GSBART"),
    RMSPE = c(GSBart_RMSPE),
    MSPE = c(GSBart_MSPE),
    times = c(GSBart_Time),
    sim = c(repetition)
  )
}

Torus_count_summary = Torus_count_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(RMSPE, MSPE),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(times), mean, .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))

save(Ushape_count_summary, Torus_count_summary, file = 'data/sim_count_res.Rdata', compress = 'xz')
