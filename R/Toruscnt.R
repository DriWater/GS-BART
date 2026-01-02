rm(list=ls())
library(doParallel)
library(purrr)
library(dplyr)
library(G2SBart)

setwd("/Users/shurenhe/Downloads/GSBart_rebuttal/Torus/")
load("torus_dist_mat_grid.RData")

# source('genTorusExample.R')
source('~/Downloads/code/GSBart_test/VAR_Model/genTorusExample.R')

n_bin = 100
n_rounds = 50
nref = 100
num_ost = 5
hyperpar = c("max_splits" = 20, 
             "lr" = 6/n_rounds, 
             "n_rounds" = n_rounds, 
             "lambda" = 0, 
             "alpha" = 0, 
             "nref" = nref, 
             "n_bin" = n_bin,
             "num_ost" = num_ost)

p = 5

sim = genTorusExample (p, hyperpar)

# source("/Users/shurenhe/Downloads/code/GSBart_test/Multi_Model/Generate_ChainfromBART.R")
# cutpoints_lst = get_cutpoint_from_bart(sim$train_X, sim$test_X, sim$f_true)
# Graphs = sim$Graphs
# for(i in 1:length(Graphs)){
#   Graphs_tmp <- Generate_GraphCandidates_graphbin(g0=NULL, sim$train_X, X_ho = sim$test_X, n_bin = NULL, cutpoints_lst = cutpoints_lst)
#   Graphs[[i]][(num_ost+1):(num_ost + ncol(sim$train_X))] <- Graphs_tmp
# }
# sim$Graphs = Graphs

n_train = nrow(sim$X)
n_test = nrow(sim$X_ho)
sigmasq_y = 0.1

alpha=0.95
beta=2
rb = .9
a = 3

hyperpar=c( 
  alpha=alpha,
  beta=beta,
  rb=rb, 
  split_eps=0,
  a = a,
  max_depth = 6,
  tree_iter = 12
)

p = p + 1

modelpar=c(sigmasq=var(sim$f_true+.01), tausq = (6/(2*sqrt(n_rounds)))^2)
hyperpar <- as.list(hyperpar); modelpar <- as.list(modelpar)

graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/(num_ost + 2), (num_ost + 2)), rep(min(p / (p + 2), 0.85)/p,p))

InputGraphs = lapply(sim$Graphs, function(x) {
  lapply(x, function(y) {
    stopifnot(length(igraph::graph.attributes(y$g)$tr2mesh) == n_train,
              length(igraph::graph.attributes(y$g)$ho2mesh) == n_test,
              length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2tr))) == n_train,
              length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2ho))) == n_test)
    return(list(root = y$root,
                children = y$children,
                mesh2tr = igraph::vertex.attributes(y$g)$mesh2tr,
                mesh2ho = igraph::vertex.attributes(y$g)$mesh2ho,
                tr2mesh = igraph::graph.attributes(y$g)$tr2mesh,
                ho2mesh = igraph::graph.attributes(y$g)$ho2mesh))
  }) 
})
# Reduce every index by 1 to make everything 0 indexed
InputGraphs = lapply(InputGraphs, function(x) {
  lapply(x, function(y) {
    list(root = y$root - 1,
         children = unname(lapply(y$children, function(z){unname(z) - 1})),
         mesh2tr = unname(lapply(y$mesh2tr, function(z){unname(z) - 1})),
         mesh2ho = unname(lapply(y$mesh2ho, function(z){unname(z) - 1})),
         tr2mesh = unname(y$tr2mesh - 1),
         ho2mesh = unname(y$ho2mesh - 1))
  }) 
})

save(sim, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/TorusSim.Rdata")

load('~/Downloads/code/test_code/TorusSim.Rdata')


nworkers = 8
cl <- makeCluster(nworkers)
registerDoParallel(cl)

n_rounds = length(InputGraphs)
n_train = nrow(sim$X)
n_test = nrow(sim$X_ho)
repetitions = 50

sim.list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'bayestestR','G2SBart', 
                                                             'caret')) %dopar% {
                                                               set.seed(1234+repetition)
                                                               message("Sampling Y")
                                                               sim$Y0 <- exp((sim$f_true + 17)/8)
                                                               sim$Y0_ho <- exp((sim$f_ho_true + 17)/8)
                                                               sim$Y <- rpois(n_train, sim$Y0)
                                                               sim$Y_ho <- rpois(n_test, sim$Y0_ho)
                                                               hyperpar['b'] = 0.5*var(log(sim$Y+.01))/n_rounds
                                                               nu = 3; q = 0.9;  quant = qchisq(1-q, nu)
                                                               hyperpar['lambda'] = (var(log(sim$Y+.01))*quant)/nu
                                                               hyperpar['mu0'] = mean(log(sim$Y+0.01))/n_rounds
                                                               GSBart_Time = Sys.time()
                                                               GSBres = g2sbart(sim$Y, InputGraphs, 200, 15, 200, family = "poisson", graphs_weight = graphs_weight, 
                                                                                nthreads = 7, hyperpar = hyperpar, modelpar = modelpar, seed = 1234)
                                                               GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
                                                               GSBart_RMSE = sqrt(mean((sim$Y0_ho - exp(colMeans(GSBres$phi.test)))^2))
                                                               GSBart_MSE = mean((sim$Y0_ho - exp(colMeans(GSBres$phi.test)))^2)
                                                               
                                                               sim_results_df = data.frame(
                                                                 models = c("GSBart"),
                                                                 RMSPE = c(GSBart_RMSE),
                                                                 MSPE = c(GSBart_MSE),
                                                                 times = c(GSBart_Time),
                                                                 sim = c(repetition)
                                                               )
                                                               sim_results_df
                                                             }



# XGBoost hyperparameters
xgb_grid = expand.grid(
  # nrounds = n_rounds,
  max_depth = 20,
  eta = 5 / seq(10, 100, length = 10),
  gamma = seq(0, 1, length = 10),
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

n_train = nrow(sim$X); n_test = nrow(sim$X_ho)
n_rounds = length(InputGraphs)
repetitions = c(1:50)
sim.list = vector(mode = "list", length = length(repetitions))

MCMC = 5000
BURNIN = 1000
THIN = 5
ndpost = (MCMC - BURNIN)
nskip = BURNIN

nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)
repetitions = c(1:50)
sim.list <- foreach(i = 1:length(repetitions), .packages= c('BART','xgboost', 'caret' )) %dopar% {
  repetition = repetitions[i]
  set.seed(1234+repetition)
  sim$Y0 = exp((sim$f_true+17)/8)
  sim$Y0_ho = exp((sim$f_ho_true+17)/8)
  sim$Y <- rpois(n_train, sim$Y0)
  sim$Y_ho <- rpois(n_test, sim$Y0_ho)
  
  XGB_Time = Sys.time()
  XGB_RMSE_CV = Inf
  for(i in 1:nrow(xgb_grid)){
    cv <- xgb.cv(data = sim$train_X, label = sim$Y, nrounds = 50, params = as.list(xgb_grid[i,]),
                 nfold = 5, metrics = "rmse", objective='count:poisson', verbose = F)
    RMSE <- mean(cv$evaluation_log$test_rmse_mean[20:50])
    if(RMSE < XGB_RMSE_CV){
      XGB_RMSE_CV <- RMSE
      para_select <- as.list(xgb_grid[i,])
    }
  }
  XGB_Fit<- xgboost(data = sim$train_X, label = sim$Y, nrounds = 50, params = para_select,
                    metrics = "rmse", objective='count:poisson')
  XGB_Time = Sys.time() - XGB_Time
  b = predict(XGB_Fit, sim$test_X)
  XGB_RMSE = sqrt(mean((sim$Y0_ho - b)^2))
  XGB_MSPE = mean((sim$Y0_ho - b)^2)
  
  # d <- data.frame(Y = sim$Y)
  # d <- cbind(d, sim$train_X)
  # options(mc.cores=5)
  # GLMM_Time = Sys.time()
  # m_spatial <- glmmfields(Y ~ ., family = poisson(link = "log"),
  #                         data = d, lat = "x", lon = "y", nknots = 60, iter = 2000, chains = 5,
  #                         prior_intercept = student_t(3, 0, 12),
  #                         prior_beta = student_t(3, 0, 3),
  #                         prior_sigma = half_t(3, 0, 6),
  #                         prior_gp_theta = half_t(3, 0, 10),
  #                         prior_gp_sigma = half_t(3, 0, 6),
  #                         seed = 1234 # passed to rstan::sampling()
  # )
  # Y_pred <- predict(
  #   m_spatial,
  #   newdata = test_X,
  #   type = "response"
  # )$estimate
  # GLMM_Time = difftime(Sys.time(), GLMM_Time, units = "secs")
  # GLMM_RMSE = sqrt(mean((sim$Y_ho - Y_pred)^2))
  
  message("Fitting BART")
  # Fit BART
  BART_Time = Sys.time()
  BART_Fit = wbart(sim$train_X, sim$Y, sim$test_X,
                   ntree = n_rounds, ndpost = ndpost, nskip = nskip)
  BART_Time = difftime(Sys.time(), BART_Time, units = "secs")
  BART_y_test = colMeans(BART_Fit$yhat.test[seq(1, ndpost, THIN), ])
  BART_RMSE = sqrt(mean((sim$Y0_ho - BART_y_test)^2))
  BART_MSPE = mean((sim$Y0_ho - BART_y_test)^2)
  
  
  # sim_1_results_df = data.frame(
  # models = c("GGBoost", "BART", "XGBoost", "RF"),
  # ACC = c(GGBoost_ACC_1, BART_ACC_1, XGBoost_ACC_1, RF_ACC_1),
  # sim = rep(repetition, 4)
  # )
  sim_results_df =  data.frame(
    models = c("BART", "XGBoost"),
    RMSE = c(BART_RMSE, XGB_RMSE),
    MSPE = c(BART_MSPE, XGB_MSPE),
    times = c(BART_Time, XGB_Time),
    sim = rep(repetition, 2)
  )
  sim_results_df
}

sim.summary = sim.list %>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('RMSE','MSPE','times'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("BART", "XGBoost"), desc(models)))

print(sim.summary)

save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Torus/Toruscntpartial.Rdata")

XGB_Time = Sys.time()
XGB_RMSE_CV = Inf
for(i in 1:nrow(xgb_grid)){
  cv <- xgb.cv(data = sim$train_X, label = sim$Y, nrounds = 50, params = as.list(xgb_grid[i,]),
               nfold = 5, metrics = "rmse", objective='count:poisson', verbose = F)
  RMSE <- mean(cv$evaluation_log$test_rmse_mean[20:50])
  if(RMSE < XGB_RMSE_CV){
    XGB_RMSE_CV <- RMSE
    para_select <- as.list(xgb_grid[i,])
  }
}
XGB_Fit<- xgboost(data = sim$train_X, label = sim$Y, nrounds = 50, params = para_select,
                  metrics = "rmse", objective='count:poisson')
XGB_Time = Sys.time() - XGB_Time
b = predict(XGB_Fit, sim$test_X)
XGB_RMSE = sqrt(mean((sim$Y0_ho - b)^2))
print(XGB_RMSE)

# load("~/Downloads/GSBart_rebuttal/Torus/GSBToruscntres.Rdata")
# 
# sim.summary = sim.list %>%
#   list_rbind() %>%
#   group_by(models) %>%
#   summarise(across(c('RMSE', 'MSE', 'times'),
#                    list("mean" = mean, "sd" = sd),
#                    .names = "{.col}_{.fn}")) %>%
#   arrange(match(models, c("GSBart"), desc(models)))
# 
# print(sim.summary)
# 
# save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Torus/AddiVTorusreg.Rdata")
# 
# load("/Users/shurenhe/Downloads/GSBart_rebuttal/Torus/AddiVTorusreg.Rdata")
# 
# sim.list1 <- sim.list
# 
# load("~/Downloads/GSBart_rebuttal/Torus/AddiVTorusreg1.Rdata")
# 
# sim.list.comb <- vector(mode = 'list', length = 50)
# sim.list.comb[1:20] <- sim.list
# sim.list.comb[21:50] <- sim.list1
# 
# sim.list <- sim.list.comb

load("~/Downloads/GSBart_rebuttal/Torus/GSBToruscntorg.Rdata")
library(dplyr)
library(purrr)

sim.summary = sim.list%>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'RMSE', 'times'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("GSBART"), desc(models)))

print(sim.summary)

save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Torus/GSBToruscntorg.Rdata")
