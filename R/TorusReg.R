rm(list=ls())

# library(dplyr)
# library(caret)
# library(purrr)
# library(testthat)
# library(parallel)
# library(sf)
# library(doParallel)
# library(GGally)
# library(igraph)
# library(nngeo)


load("~/Downloads/code/test_code/TorusSim.Rdata")
library(GSBart)

set.seed(1236)
n_train = length(sim$f_true)
n_test = length(sim$f_ho_true)
noise = 0.1
Y0 = sim$f_true + rnorm(n_train, 0, noise)
Y0_ho = sim$f_ho_true + rnorm(n_test, 0, noise)

GSBart_Time=Sys.time()
GSBART_Fit = GSBart::gsbart(Y0, Graphs = InputGraphs, 200, 15, 200, graphs_weight, nthreads = 7, seed = 1234)
GSBart_Time=difftime(Sys.time(), GSBart_Time, units = "secs")
repetition = 1
set.seed(1234+repetition)                                                 
Y0 = sim$f_true + rnorm(n_train, 0, noise)
Y0_ho = sim$f_ho_true + rnorm(n_test, 0, noise)

GSBart_Time=Sys.time()
GSBART_Fit = gsbart(Y0, InputGraphs, 200, 15, 200, graphs_weight, nthreads = 7, seed = 1234)
GSBart_Time=difftime(Sys.time(), GSBart_Time, units = "secs")

# n_rounds = 50
# 
# modelpar['tausq'] =  (0.5/(2*sqrt(n_rounds)))^2
# hyperpar['mu0'] = 0
# 
# save(sim, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/testTorusReg.Rdata")

n_train = nrow(sim$X); n_test = nrow(sim$X_ho)
noise = 0.1; repetitions = 50

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

Unstandardize = function(x, std_par) {
  return(x * std_par['scale'] + std_par['mean']) 
} 

coords_all = rbind(sim$coords, sim$coords_ho)
df_coords <- data.frame(lat = coords_all[,1], long = coords_all[,2])
sf_coords <- st_as_sf(df_coords, coords = c(1:2))

nn.list = st_nn(sf_coords, sf_coords, k = 5)
g0 = graph.adjlist(nn.list) %>% as.undirected() %>% igraph::simplify()
E(g0)$weight = runif(ecount(g0), 0, 1)
components(g0)$membership
adj_matrix <- as_adjacency_matrix(g0, sparse = FALSE)

ggnet2(g0, mode = st_coordinates(sf_coords), node.size = 0, edge.size = 0.5, edge.color = "orange")

library(deldir)
tri <- deldir(df_coords[,1], df_coords[,2])
dis <- sqrt( (tri$delsgs$x1 - tri$delsgs$x2)^2 + (tri$delsgs$y1 - tri$delsgs$y2)^2 )
edges_df <- tri$delsgs[, c("ind1", "ind2")]
edges_df <- edges_df[dis<1.28,]
# Step 3: Create igraph object from edges
g_tri <- graph_from_edgelist(as.matrix(edges_df), directed = FALSE)
E(g_tri)$weight = runif(ecount(g_tri), 0, 1)
# components(g_tri)$membership
adj_tri_matrix <- as_adjacency_matrix(g_tri, sparse = FALSE)
ggnet2(g_tri, mode = st_coordinates(sf_coords), node.size = 0, edge.size = 0.5, edge.color = "orange")

# voronoi <- st_voronoi(st_union(sf_coords), sf_bnd) 
# sf_mesh = st_intersection(st_cast(voronoi), sf_bnd)
# neighboring = st_relate(sf_mesh, pattern = "****1****") ## check which polygons share boundaries
# g0 = graph.adjlist(neighboring) %>% as_undirected() %>% igraph::simplify()
# library(network)
# plotGraph(st_coordinates(st_centroid(sf_mesh)), g0)


X_all <- rbind(sim$train_X, sim$test_X)
vertex_id_train = c(1:n_train)
vertex_id_test = c((n_train+1):(n_train+n_test))
X_all_std <- apply(X_all, 2, function(x){(x - min(x)) / (max(x) - min(x))})
X_train_std = X_all_std[vertex_id_train,]
X_test_std = X_all_std[vertex_id_test,]


nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)
repetitions = c(1:50)

nburnin = 15; ndpost = 200
sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'bayestestR','matrixStats','stochtree')) %dopar% {
  # sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'data.tree','matrixStats','FlexBart')) %dopar% {
  message("Sampling Y")
  set.seed(1234+repetition)
  sigma_y = noise
  Y0 = sim$f_true + rnorm(n_train, 0, noise)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, noise)
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  
  
  XBART_Time=Sys.time()
  XBART_Fit <- stochtree::bart(sim$train_X, sim$Y, X_test = sim$test_X, num_gfr = (nburnin+ndpost), num_burnin = 0, num_mcmc = 0, 
                               mean_forest_params = list(num_trees = 50), general_params = list(keep_gfr = T))
  XBART_Time = difftime(Sys.time(), XBART_Time, units = "secs")
  XBART_Fit$time = XBART_Time
  XBART_Y_out = t(apply(XBART_Fit$y_hat_test, 2, function(X){Unstandardize(X, stdY$std_par)}))
  XBART_pred_unstandardized = colMeans(XBART_Y_out[(nburnin+1):(ndpost+nburnin),])
  XBART_MSPE = mean((XBART_pred_unstandardized - Y0_ho)^2)
  XBART_MAPE = mean(abs(XBART_pred_unstandardized - Y0_ho))
  
  XBART.upper.lvl = NULL
  XBART.lower.lvl = NULL
  for(i in 1:ncol(XBART_Y_out)){
    tmp <- ci(XBART_Y_out[,i], method = "HDI")
    XBART.upper.lvl <- append(XBART.upper.lvl, tmp$CI_high)
    XBART.lower.lvl <- append(XBART.lower.lvl, tmp$CI_low)
  }
  XBART_coverage = mean((Y0_ho<XBART.upper.lvl)&(Y0_ho>XBART.lower.lvl))
  XBART_HDI_len = mean(XBART.upper.lvl - XBART.lower.lvl)
  XBART_poster_sd = mean(apply(XBART_Y_out, 2, sd))
  
  
  sim_results_df = data.frame(
    models = c("XBART"),
    MSPE = c(XBART_MSPE),
    MAPE = c(XBART_MAPE),
    # MSPE_std = c(Flex_BART_MSPE_std), 
    times = c(XBART_Time),  
    HDI_coverage = c(XBART_coverage),
    HDI_len = c(XBART_HDI_len),
    poster_sd = c(XBART_poster_sd),
    sim = c(repetition)
  )
  sim_results_df
}

sim_summary = sim_list%>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'times', 'HDI_coverage', 'HDI_len', 'poster_sd'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("XBART"), desc(models)))

print(sim_summary)

sim.list <- foreach(i = 1:length(repetitions), .packages= c('Matrix','igraph', 'data.tree','flexBART')) %dopar% {
  repetition = repetitions[i]
  set.seed(1234+repetition)
  Y0 = sim$f_true + rnorm(n_train, 0, noise)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, noise)
  ## standardize data
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  
  Flex_BART_Time=Sys.time()
  Flex_BART_Fit = flexBART::network_BART(sim$Y, vertex_id_train, X_train_std, vertex_id_test, verbose = T,
                                         X_test_std, A = adj_tri_matrix, M = 50, nd = 400, burn = 2000, thin = 10)
  # Flex_BART_Fit1 = flexBART(sim$Y, X_cont_train = X_train_std, X_cont_test =  X_test_std,
  #                          M = 100, nd = 400, burn = 2000, thin = 10)
  Flex_BART_Time = difftime(Sys.time(), Flex_BART_Time, units = "secs")
  Flex_BART_Y_out = t(apply(Flex_BART_Fit$yhat.test, 1, function(X){Unstandardize(X, stdY$std_par)}))
  Flex_BART_pred_unstandardized = colMeans(Flex_BART_Y_out)
  Flex_BART_MSPE = mean((Flex_BART_pred_unstandardized - Y0_ho)^2)
  Flex_BART_MAPE = mean(abs(Flex_BART_pred_unstandardized - Y0_ho))
 
  Flex.BART.upper.lvl = NULL
  Flex.BART.lower.lvl = NULL
  for(i in 1:ncol(Flex_BART_Y_out)){
    tmp <- ci(Flex_BART_Y_out[,i], method = "HDI")
    Flex.BART.upper.lvl <- append(Flex.BART.upper.lvl, tmp$CI_high)
    Flex.BART.lower.lvl <- append(Flex.BART.lower.lvl, tmp$CI_low)
  }
  Flex_BART_coverage = mean((Y0_ho<Flex.BART.upper.lvl)&(Y0_ho>Flex.BART.lower.lvl))
  Flex_BART_HDI_len = mean(Flex.BART.upper.lvl - Flex.BART.lower.lvl)
  Flex_BART_poster_sd = mean(apply(Flex_BART_Y_out, 2, sd))
  
  
  sim_results_df = data.frame(
    models = c("Flex-BART"),
    MSPE = c(Flex_BART_MSPE),
    MAPE = c(Flex_BART_MAPE),
    # MSPE_std = c(Flex_BART_MSPE_std), 
    times = c(Flex_BART_Time),  
    HDI_coverage = c(Flex_BART_coverage),
    HDI_len = c(Flex_BART_HDI_len),
    poster_sd = c(Flex_BART_poster_sd),
    sim = c(repetition)
  )
  sim_results_df
}

sim_summary = sim_list%>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'times', 'HDI_coverage', 'HDI_len', 'poster_sd'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("Flex-BART"), desc(models)))

print(sim_summary)

load("/Users/shurenhe/Downloads/GSBart_rebuttal/Torus/testFlexTorusReg.Rdata")

source("/Users/shurenhe/Downloads/code/GSBart/AddiVortes/AddiVortesMainCode.R")

# load("~/Downloads/code/test_code/testToruscnt1.Rdata")

sim.list <- foreach(i = 1:length(repetitions), .packages= c('Matrix','igraph', 'data.tree','matrixStats','FNN', 'invgamma',
                                                            'RandomForestsGLS', 'BRISC')) %dopar% {
  repetition = repetitions[i]
  set.seed(1234+repetition)
  sigma_y = 0.1
  Y0 = sim$f_true + rnorm(n_train, 0, sigma_y)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, sigma_y)
  ## standardize data
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  AddiV_Time = Sys.time()
  AddiV_Fit<- AddiVortes_Algorithm(sim$Y, sim$train_X, m = 50, max_iter = 4000, burn_in = 1000, sd = sd(sim$Y),
                                   YTest = sim$Y_ho, XTest = sim$test_X, IntialSigma = "Linear", thinning = 10)
  AddiV_Time = difftime(Sys.time(), AddiV_Time, units = "secs")
  AddiV_pred_unstandardized = Unstandardize(AddiV_Fit$yhat.test.mean, stdY$std_par)
  AddiV_MSPE = mean((AddiV_pred_unstandardized - Y0_ho)^2)
  AddiV_MAPE = mean(abs(AddiV_pred_unstandardized - Y0_ho))
  
  # RFGLS_Time = Sys.time()
  # est_known<-RFGLS_estimate_spatial(sim$coords[,1:2], sim$Y, sim$train_X, ntree=50, cov.model="exponential", nthsize=20, param_estimate=TRUE, h = 5)
  # RFGLS_Time = difftime(Sys.time(), RFGLS_Time, units = "secs")
  # RFGLS_y_test<-try(RFGLS_predict_spatial(est_known,sim$coords_ho[,1:2], sim$test_X)$prediction)
  # if (!("try-error" %in% class(RFGLS_y_test))) {
  #   RFGLS_pred_unstandardized <- Unstandardize(RFGLS_y_test, stdY$std_par)
  #   RFGLS_MSPE = mean((RFGLS_pred_unstandardized - Y0_ho)^2)
  #   RFGLS_MAPE = mean(abs(RFGLS_pred_unstandardized - Y0_ho))
  # }else{
  #   RFGLS_MSPE = NA
  #   RFGLS_MAPE = NA
  # }
  # 
  # BRISC_Time = Sys.time()
  # BRISC_Fit <- BRISC_estimation(sim$coords[,1:2], sim$Y, sim$X)
  # BRISC_Time = Sys.time() - BRISC_Time
  # BRISC_y_test <- BRISC_prediction(BRISC_Fit, sim$coords_ho[,1:2],sim$X_ho)$prediction
  # BRISC_pred_unstandardized = Unstandardize(BRISC_y_test, stdY$std_par)
  # BRISC_MSPE = mean((BRISC_pred_unstandardized - Y0_ho)^2)
  # BRISC_MAPE = mean(abs(BRISC_pred_unstandardized - Y0_ho))
  # 
  # sim_results_df = data.frame(
  #   models = c("AddiV", "RFGLS", "BRISC"),
  #   MSPE = c(AddiV_MSPE, RFGLS_MSPE, BRISC_MSPE),
  #   MAPE = c(AddiV_MAPE, RFGLS_MAPE, BRISC_MAPE),
  #   times = c(AddiV_Time, RFGLS_Time, BRISC_Time),
  #   sim = rep(repetition, 3)
  # )
  sim_results_df = data.frame(
    models = c("AddiV"),
    MSPE = c(AddiV_MSPE),
    MAPE = c(AddiV_MAPE),
    times = c(AddiV_Time),
    sim = rep(repetition, 1)
  )
  sim_results_df
}

sim.summary = sim.list %>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'times'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("AddiV"), desc(models)))

print(sim.summary)

save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Torus/TorusregAddiv.Rdata")

repetitions = c(1:50)

xgb_grid = expand.grid(
  nrounds = n_rounds,
  max_depth = 20,
  eta = 5 / seq(10, 100, length = 10),
  gamma = seq(0, 1, length = 10),
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

niter = 9000
nburnin = 4000
THIN = 5
ndpost = niter - nburnin
nskip = nburnin

sim.list <- foreach(i = 1:length(repetitions), .packages= c('xgboost', 'caret', 'BART')) %dopar% {
  repetition = repetitions[i]
  set.seed(1234+repetition)
  sigma_y = 0.1
  Y0 = sim$f_true + rnorm(n_train, 0, sigma_y)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, sigma_y)
  ## standardize data
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  # # Fit XGBoost
  XGBoost_Time = Sys.time()
  XGBoost_Fit = caret::train(x = sim$train_X,
                             y = sim$Y,
                             tuneGrid = xgb_grid,
                             method = 'xgbTree',
                             nthread = 1L,
                             verbose = FALSE, verbosity = 0)
  XGBoost_Time = Sys.time() - XGBoost_Time
  XGBoost_y_test = predict(XGBoost_Fit, sim$test_X)
  XGBoost_pred_unstandardized = Unstandardize(XGBoost_y_test, stdY$std_par)
  XGBoost_MSPE = mean((XGBoost_pred_unstandardized - Y0_ho)^2)
  XGBoost_MAPE = mean(abs(XGBoost_pred_unstandardized - Y0_ho))
  
  # Fit BART
  BART_Time = Sys.time()
  BART_Fit = wbart(sim$train_X, sim$Y, sim$test_X,
                   ntree = n_rounds, ndpost = ndpost, nskip = nskip)
  BART_Time = difftime(Sys.time(), BART_Time, units = "secs")
  BART_y_test = colMeans(BART_Fit$yhat.test[seq(1, ndpost, THIN), ])
  BART_pred_unstandardized = Unstandardize(BART_y_test, stdY$std_par)
  BART_MSPE = mean((BART_pred_unstandardized - Y0_ho)^2)
  BART_MAPE = mean(abs(BART_pred_unstandardized - Y0_ho))
  
  sim_results_df = data.frame(
    models = c("XGBoost", "BART"),
    MSPE = c(XGBoost_MSPE, BART_MSPE),
    MAPE = c(XGBoost_MAPE, BART_MAPE),
    times = c(XGBoost_Time, BART_Time),
    sim = rep(repetition, 2)
  )
  sim_results_df
}

sim.summary = sim.list %>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'times'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("XGBoost", "BART"), desc(models)))


save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Torus/Torusregpartial1.Rdata")


stop(cl)
  

library('RandomForestsGLS')
RFGLS_Time = Sys.time()
est_known<-RFGLS_estimate_spatial(sim$coords[,1:2], sim$Y, sim$train_X, ntree=50, cov.model="exponential", nthsize=20, param_estimate=TRUE, h = 5)
RFGLS_Time = difftime(Sys.time(), RFGLS_Time, units = "secs")
RFGLS_y_test<-try(RFGLS_predict_spatial(est_known,sim$coords_ho[,1:2], sim$test_X)$prediction)
if (!("try-error" %in% class(RFGLS_y_test))) {
  RFGLS_pred_unstandardized <- Unstandardize(RFGLS_y_test, stdY$std_par)
  RFGLS_MSPE = mean((RFGLS_pred_unstandardized - Y0_ho)^2)
  RFGLS_MAPE = mean(abs(RFGLS_pred_unstandardized - Y0_ho))
}else{
  RFGLS_MSPE = NA
  RFGLS_MAPE = NA
}
print(RFGLS_MSPE); print(RFGLS_MAPE)

library(BRISC)
# # Fit BRISC 
BRISC_Time = Sys.time()
BRISC_Fit <- BRISC_estimation(sim$coords[,1:2], sim$Y, sim$X)
BRISC_Time = Sys.time() - BRISC_Time
BRISC_y_test <- BRISC_prediction(BRISC_Fit, sim$coords_ho[,1:2],sim$X_ho)$prediction
BRISC_pred_unstandardized = Unstandardize(BRISC_y_test, stdY$std_par)
BRISC_MSPE = mean((BRISC_pred_unstandardized - Y0_ho)^2)
BRISC_MAPE = mean(abs(BRISC_pred_unstandardized - Y0_ho))

print(BRISC_MSPE); print(BRISC_MAPE)

load("~/Downloads/GSBart_rebuttal/Torus/GSBTorusorgReg.Rdata")

load("~/Downloads/GSBart_rebuttal/Friedman/GSBFriedmanReg.Rdata")
