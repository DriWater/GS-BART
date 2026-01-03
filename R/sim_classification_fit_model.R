rm(list=ls())
library(bayestestR)
library(G2SBart)
library(purrr)
library(caret)

load("data/sim_input.RData")

repetitions = 50

## Ushape Simulation
Ushape_class_list = vector(50, mode="list")
sigma = 0.15
n = nrow(sim_Ushape$X); n_ho = nrow(sim_Ushape$X_ho)

for(repetition in 1:repetitions){
  set.seed(1234+repetition)
  message("Sampling Y")
  Y0 = sim_Ushape$f_true + rnorm(n, 0, sigma)
  Y0_ho = sim_Ushape$f_ho_true + rnorm(n_ho, 0, sigma)
  
  multinom_resp = as.numeric(cut(Y0,c(-Inf,quantile(Y0,c(.2,.3,.65,.8)),Inf),labels=1:5))-1
  multinom_resp_test =as.numeric(cut(Y0_ho,c(-Inf,quantile(Y0,c(.2,.3, .65,.8)),Inf),labels=1:5))-1
  K = length(unique(multinom_resp))
  message("Fitting GSBART")
  MGSB_Time = Sys.time()
  MGSB_Fit = G2SBart::mgsbart2(multinom_resp, sim_Ushape$Graphs, 200, 15, 200, 
                               sim_Ushape$graphs_weight, seed = 1234)
  MGSB_Time = Sys.time() - MGSB_Time
  MGSB_y_test = MGSB_Fit$yhat.test
  MGSB_ACC = mean(MGSB_y_test == multinom_resp_test)
  MGSB_CM = confusionMatrix(factor(MGSB_y_test, levels = (0:(K - 1))),
                            factor(multinom_resp_test, levels = (0:(K - 1))))
  Ushape_class_list[[repetition]] = data.frame(
    models = c("GSBART"), 
    ACC = c(MGSB_ACC),
    times = c(MGSB_Time), 
    sim = c(repetition)
  )
  sim_results_df
}

Ushape_class_summary = Ushape_class_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(ACC),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(times), mean,
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))


# Torus Simulation
Torus_class_list = vector(50, mode="list")
sigma = 0.1
n = nrow(sim_Torus$X); n_ho = nrow(sim_Torus$X_ho)

for(repetition in 1:repetitions){
  set.seed(1234+repetition)
  message("Sampling Y")
  Y0 = sim$f_true + rnorm(n, 0, sigma) 
  Y0_ho = sim$f_ho_true + rnorm(n_ho, 0, sigma)
  
  multinom_resp = as.numeric(cut(Y0,c(-Inf,quantile(sim$Y,c(.2,.3,.65,.8)),Inf),labels=1:5))-1
  multinom_resp_test =as.numeric(cut(Y0_ho,c(-Inf,quantile(sim$Y,c(.2,.3, .65,.8)),Inf),labels=1:5))-1
  K = length(unique(multinom_resp))
  message("Fitting GSBART")
  MGSB_Time = Sys.time()
  MGSB_Fit = G2SBart::mgsbart2(multinom_resp, sim_Torus$Graphs, 200, 15, 200, 
                               sim_Torus$graphs_weight, seed = 1234)
  MGSB_Time = Sys.time() - MGSB_Time
  MGSB_y_test = MGSB_Fit$yhat.test
  MGSB_ACC = mean(MGSB_y_test == multinom_resp_test)
  MGSB_CM = confusionMatrix(factor(MGSB_y_test, levels = (0:(K - 1))),
                            factor(multinom_resp_test, levels = (0:(K - 1))))
  Torus_class_list[[repetition]] = data.frame(
    models = c("GSBART"), 
    ACC = c(MGSB_ACC),
    times = c(MGSB_Time), 
    sim = c(repetition)
  )
  sim_results_df
}

Torus_class_summary = Torus_class_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(
    across(c(MGSB_ACC),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    across(c(times), mean, .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%  
  arrange(match(models, c("GSBART"), desc(models)))

save(Ushape_class_summary, Torus_class_summary, file = 'data/sim_class.RData', compress = 'xz')
