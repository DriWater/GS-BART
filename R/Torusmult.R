rm(list=ls())

library(igraph)
library(dplyr)
library(caret)
library(BART)
library(data.tree)
library(matrixStats)
library(Matrix)
library(purrr)
library(BART)
library(parallel)
library(doParallel)
library(xgboost)


load("~/Downloads/code/test_code/testToruscntorg.Rdata")
hyperpar['mu0'] <- NULL
n_rounds = 50
hyperpar["b"] = 0.5 * 6 / (n_rounds);
modelpar['tausq'] = (3/(2*sqrt(length(InputGraphs))))^2

# load('/Users/shurenhe/Downloads/code/GSBart_test/Torus/TorusRegressionfG.Rdata')

# load("~/Downloads/code/test_code/testTorusmult.Rdata")

# load("~/Downloads/GSBart_rebuttal/Torus/GSBTorusmultres.Rdata")

# tmp = c()
# 
# for(i in 1:length(sim.list)){
#   if(is.na(sim.list[[i]]$ACC)){
#     tmp = append(tmp, i)
#   }
# }
# 
# sim.summary = sim.list[-tmp]%>%
#   list_rbind() %>%
#   group_by(models) %>%
#   summarise(across(c('ACC', 'times'),
#                    list("mean" = mean, "sd" = sd),
#                    .names = "{.col}_{.fn}")) %>%
#   arrange(match(models, c("GSBART"), desc(models)))
# 
# print(sim.summary)


# save(sim, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/testTorusmultorg.Rdata")
# 
# load("~/Downloads/GSBart_rebuttal/Torus/GSBTorusmultres.Rdata")
# load("~/Downloads/GSBart_rebuttal/Torus/GSBTorusmultorg.Rdata")
nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)

sigmasq_y = 0.1
n_train = nrow(sim$X)
n_test = nrow(sim$X_ho)
n_rounds = length(InputGraphs)

MCMC = 6000
BURNIN = 1000
THIN = 5
ndpost = (MCMC - BURNIN)/THIN
nskip = BURNIN

xgb_grid = expand.grid(
  nrounds = n_rounds,
  max_depth = 20,
  eta = 5 / seq(10, 100, length = 10),
  gamma = seq(0, 1, length = 10),
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

Standardize = function(Y, std_par = NULL) {
  if(is.null(std_par)) {
    ymean = mean(Y)
    yscale = 2 * max(abs(Y))
  } else {
    ymean = std_par['mean']
    yscale = std_par['scale']
  }
  Y = (Y - ymean) / yscale
  std_par = c('mean' = ymean, 'scale' = yscale)
  return(list(Y = Y, std_par = std_par))
}

repetitions = 50
sim.list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'data.tree','matrixStats', 'xgboost', 'BART', 'caret')) %dopar% {
  set.seed(1234+repetition)
  Y0 = sim$f_true + rnorm(n_train, 0, sigmasq_y) 
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, sigmasq_y)
  
  ## standardize data
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  # b = 0.5*var(sim$Y)/length(M)
  # hyperpar_BGG['b'] = b
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  
  multinom_resp = as.numeric(cut(sim$Y,c(-Inf,quantile(sim$Y,c(.2,.3,.65,.8)),Inf),labels=1:5))-1
  multinom_resp_test =as.numeric(cut(sim$Y_ho,c(-Inf,quantile(sim$Y,c(.2,.3, .65,.8)),Inf),labels=1:5))-1
  # multinom_resp = as.numeric(cut(sim$Y,c(-Inf,quantile(sim$Y,c(.2,.4,.6,.8)),Inf),labels=1:5))-1 
  # multinom_resp_test =as.numeric(cut(sim$Y_ho,c(-Inf,quantile(sim$Y,c(.2,.4, .6,.8)),Inf),labels=1:5))-1 
  
  K = length(unique(multinom_resp)) 
  
  # message("Fitting BART")
  # Fit BART
  BART_Time = Sys.time()
  capture.output(BART_Fit <- mbart2(sim$train_X,
                                  factor(multinom_resp+1,levels=1:K),
                                     x.test=sim$test_X, type = "pbart", ntree = n_rounds,
                                     ndpost = ndpost, nskip = nskip, keepevery = THIN),
               file = 'NUL')
  BART_Time = Sys.time() - BART_Time
  #
  BART_y_test = BART_Fit$prob.test.mean %>%
   cbind(class = ., i = rep(1:n_test, each = BART_Fit$K)) %>%
   aggregate(class ~ i, data = ., FUN = "which.max") %>%
   `[`(., ,2) %>%
   `-`(1)
  BART_ACC = mean(BART_y_test == multinom_resp_test)
  BART_CM = confusionMatrix(factor(BART_y_test, levels = (0:(K - 1))),
                             factor(multinom_resp_test, levels = (0:(K - 1))))
  
  # Fit XGBoost
  XGBoost_Time = Sys.time()
  XGBoost_Fit = caret::train(x = sim$train_X,
                               y = factor(multinom_resp, levels = (0:(K - 1))),
                               #trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = 'xgbTree',
                               nthread = 1L,
                               verbose = FALSE, verbosity = 0)
  XGBoost_Time = Sys.time() - XGBoost_Time
  
  XGBoost_y_test = predict(XGBoost_Fit, sim$test_X)
  XGBoost_ACC = mean(XGBoost_y_test == multinom_resp_test)
  XGBoost_CM = confusionMatrix(factor(XGBoost_y_test, levels = (0:(K - 1))),
                                 factor(multinom_resp_test, levels = (0:(K - 1))))
  
  
  # Fit Random Forest
  # message("Fitting RF")
  # RF_Time = Sys.time()
  # RF_Fit = train(x = sim$train_X,
  #                 y = factor(multinom_resp, levels = (0:(K - 1))),
  #                 trControl = rf_trcontrol,
  #                 method = 'parRF', verbose = FALSE, verbosity = 0)
  # RF_Time = Sys.time() - RF_Time
  # 
  # RF_y_test = predict(RF_Fit, sim$test_X)
  # RF_ACC = mean(RF_y_test == multinom_resp_test)
  # RF_CM = confusionMatrix(factor(RF_y_test, levels = (0:(K - 1))),
  #                          factor(multinom_resp_test, levels = (0:(K - 1))))
  
  # sim_results_df = data.frame(
  # models = c("GGBoost", "BART", "XGBoost", "RF"),
  # ACC = c(GGBoost_ACC, BART_ACC, XGBoost_ACC, RF_ACC),
  # sim = rep(repetition, 4)
  # )
  sim_results_df = data.frame(
    models = c("XGBoost", "BART"),
    ACC = c(XGBoost_ACC, BART_ACC),
    times = c(XGBoost_Time, BART_Time),
    sim = rep(repetition, 2)
  )
  
  sim_results_df
  
}

sim.summary = sim.list%>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('ACC', 'times'), 
                   list("mean" = mean, "sd" = sd), 
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("XGBoost", "BART"), desc(models)))

print(sim.summary)