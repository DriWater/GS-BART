rm(list=ls())

library(dplyr)
library(caret)
library(purrr)
library(testthat)
library(parallel)
library(doParallel)
library(flexBART)
library(sf)
library(igraph)
library(nngeo)
library(bayestestR)
library(GSBart)

source('/Users/shurenhe/Downloads/code/GSBart_test/Ushape/Ushape_GenData_fun.R')

n_bin = 100
n_rounds = 50
nref = 100
num_ost = 5

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

hyperpar = c("n_rounds" = n_rounds,
             "nref" = nref,
             "num_ost" = num_ost,
             "n_bin"=n_bin)

n_train = 800; n_test = 200; p = 5
sim = genUshape_fG(n_train, n_test, p = p, hyperpar, seed=1234)
# 
# save(sim, file = "~/Downloads/code/test_code/Ushapesim.Rdata")

load("~/Downloads/code/test_code/Ushapesim.Rdata")
sim_Ushape = sim
sim_Ushape$train_X = NULL; sim_Ushape$test_X = NULL

n_train = length(sim$f_true); n_test = length(sim$f_ho_true)

X_all <- rbind(sim$train_X, sim$test_X)
vertex_id_train = c(1:n_train)
vertex_id_test = c((n_train+1):(n_train+n_test))
X_all_std <- apply(X_all, 2, function(x){(x - min(x)) / (max(x) - min(x))})
X_train_std = X_all_std[vertex_id_train,]
X_test_std = X_all_std[vertex_id_test,]

num_ost = 5; n_train = 800; n_test = 200; p = 5
graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/(num_ost + 2), (num_ost + 2)), rep(min(p / (p + 2), 0.85)/p,p))
noise = 0.15; repetitions = 50

nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)
set.seed(1234)

repetition = 10
set.seed(1234+repetition)                                                 
Y0 = sim$f_true + rnorm(n_train, 0, noise)
Y0_ho = sim$f_ho_true + rnorm(n_test, 0, noise)


df_coords_all = rbind(sim$coords, sim$coords_ho)
df_coords <- data.frame(lat = df_coords_all[,1], lon = df_coords_all[,2])
sf_coords <- st_as_sf(df_coords, coords = c(1:2))

library(deldir)
tri <- deldir(df_coords[,1], df_coords[,2])
dis <- sqrt( (tri$delsgs$x1 - tri$delsgs$x2)^2 + (tri$delsgs$y1 - tri$delsgs$y2)^2 )
edges_df <- tri$delsgs[, c("ind1", "ind2")]
edges_df <- edges_df[dis<0.2,]
# Step 3: Create igraph object from edges
g_tri <- graph_from_edgelist(as.matrix(edges_df), directed = FALSE)
E(g_tri)$weight = runif(ecount(g_tri), 0, 1)
components(g_tri)$membership
adj_tri_matrix <- as_adjacency_matrix(g_tri, sparse = FALSE)
vertex_id_train = c(1:n_train)

ggnet2(g_tri, mode = st_coordinates(sf_coords), node.size = 0, edge.size = 0.5, edge.color = "orange")

sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'bayestestR','matrixStats','flexBART', 'RandomForestsGLS', 'caret', 'BRISC','xgboost')) %dopar% {
  # sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'data.tree','matrixStats','FlexBart')) %dopar% {
  message("Sampling Y")
  repetition = 3
  set.seed(1234+repetition)
  sigma_y = noise
  Y0 = sim$f_true + rnorm(n_train, 0, sigma_y)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, sigma_y)

  GSBART_Fit = GSBart::gsbart(Y0, sim$Graphs, 200, 35, 200, graphs_weight, nthreads = 6, seed = 514)
  GSBart_MSPE = mean((Y0_ho - GSBART_Fit$yhat.test.mean)^2)
  GSBart_MAPE = mean(abs(Y0_ho - GSBART_Fit$yhat.test.mean))
  
  print(c(GSBart_MSPE,GSBart_MAPE))
  save(GSBART_Fit, file =  "~/Downloads/GSBart_rebuttal/Ushape/GSBres10.Rdata")
  stdY = Standardize(Y0)
  sim$Y = stdY$Y
  sim$Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']

  Flex_BART_Time=Sys.time()
  Flex_BART_Fit = flexBART::network_BART(sim$Y, vertex_id_train, X_train_std, vertex_id_test,
                          X_test_std, A = adj_tri_matrix, M = 50, nd = 600, burn = 0, thin = 10)
  # Flex_BART_Fit1 = flexBART(sim$Y, X_cont_train = X_train_std, X_cont_test =  X_test_std,
  #                          M = 100, nd = 400, burn = 2000, thin = 10)
  Flex_BART_Time = difftime(Sys.time(), Flex_BART_Time, units = "secs")
  Flex_BART_Y_out = t(apply(Flex_BART_Fit$yhat.test, 1, function(X){Unstandardize(X, stdY$std_par)}))
  Flex_BART_pred_unstandardized = colMeans(Flex_BART_Y_out)
  Flex_BART_MSPE = mean((Flex_BART_pred_unstandardized - Y0_ho)^2)
  # Flex_BART_MSPE_std = Flex_BART_MSPE/mean(Y0)
  Flex_BART_MAPE = mean(abs(Flex_BART_pred_unstandardized - Y0_ho))
  save(Flex_BART_Fit, Flex_BART_Y_out , Flex_BART_pred_unstandardized, Y0, Y0_ho, file = "~/Downloads/GSBart_rebuttal/Ushape/FLBres10.Rdata")
  
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


nburnin = 15; ndpost = 200
sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'bayestestR','matrixStats','stochtree')) %dopar% {
  # sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'data.tree','matrixStats','FlexBart')) %dopar% {
  message("Sampling Y")
  set.seed(1234+repetition)
  sigma_y = noise
  Y0 = sim$f_true + rnorm(n_train, 0, sigma_y)
  Y0_ho = sim$f_ho_true + rnorm(n_test, 0, sigma_y)
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

# 
# load('~/Downloads/code/test_code/testUshapeReg.Rdata')
# 
# # set.seed(1234)
# # Graph_len = length(InputGraphs[[1]])
# for(i in 1:length(InputGraphs)){
#   InputGraphs[[i]] = InputGraphs[[1]]
# }
# 
# for(i in 1:length(InputGraphs)){
#   for(j in 1:length(InputGraphs[[i]])){
#     InputGraphs[[i]][[j]]$mesh2tr <- InputGraphs[[i]][[j]]$mesh2trList
#     InputGraphs[[i]][[j]]$mesh2ho <- InputGraphs[[i]][[j]]$mesh2hoList
#     # InputGraphs[[i]][[j]]$mesh2ho <- NULL
#     InputGraphs[[i]][[j]]$mesh2trList <- NULL
#     InputGraphs[[i]][[j]]$mesh2hoList <- NULL
#     InputGraphs[[i]][[j]]$ho2mesh <- NULL
#   }
# }
# 
# hyperpar['mu0'] = NULL

# repetitions = c(1:50)
sigmasq_y = 0.15; n = nrow(sim$X); n_ho = nrow(sim$X_ho) # Ushape
p = 5; num_ost = 5
graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/(num_ost + 2), (num_ost + 2)), rep(min(p / (p + 2), 0.85)/p,p))
# n_rounds= 50
# modelpar['tausq'] = (0.5/(2*sqrt(n_rounds)))^2
# nu = 3; q = 0.9; quant = qchisq(1-q, nu)
# hyperpar['tree_iter'] = 12; hyperpar['max_depth'] = 6

# sim.list = vector(mode = "list", length = length(repetitions))

# nworkers <- detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
repetitions = 50

set.seed(1234)
# subInputGraphs = vector(mode = "list", length = length(sim$Graphs))
# Graph_len = length(sim$Graphs[[1]])
# for( i in 1:length(sim$Graphs)){
#   graphs_id = sample(1:Graph_len, 9, prob = graphs_weight)
#   graphs_id = graphs_id[order(graphs_id)]
#   subInputGraphs[[i]] = sim$Graphs[[i]][graphs_id]
# }
# 

sim.list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'bayestestR','matrixStats','GSBart', 
                                                             'caret')) %dopar% {
  set.seed(1234+repetition)
  message("Sampling Y")
  Y0 = sim$f_true + rnorm(n, 0, sigmasq_y)
  Y0_ho = sim$f_ho_true + rnorm(n_ho, 0, sigmasq_y)
  GSBart_Time = Sys.time() 
  GSBres <- gsbart(Y0, sim$Graphs, 25, 15, 25, nthreads = 1, sparse = T, verbose = F, seed = 1234)
  GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
  GSBart_MSPE = mean((Y0_ho - GSBres$yhat.test.mean)^2)
  GSBart_MAPE = mean(abs(Y0_ho - GSBres$yhat.test.mean))
  # vcntmean = colMeans(GSBres$varcount) 
  # spatial_prop = sum(vcntmean[1:7])/sum(vcntmean)
  # x1_prop = sum(vcntmean[8])/sum(vcntmean)
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

  sim_results_df = data.frame(
    models = c("GSBart"),
    MSPE = c(GSBart_MSPE),
    MAPE = c(GSBart_MAPE),
    times = c(GSBart_Time),
    # spatial_prop = spatial_prop,
    # x1_prop = x1_prop,
    HDI_coverage = c(GSBart_coverage),
    HDI_len = c(GSBart_HDI_len),
    poster_sd = c(GSBart_poster_sd),
    sim = c(repetition)
  )
  sim_results_df
}

stopCluster(cl)
# save(sim.list, file = "~/Downloads/GSBart_rebuttal/Ushape/GSBUshapereg.Rdata")

sim.summary = sim.list%>%
  list_rbind() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'times', 'HDI_coverage', 'HDI_len',
                     # 'spatial_prop', 'nonspatial_prop',
                     'poster_sd'),
  # summarise(across(c('MSPE', 'MAPE', 'times', 'spatial_prop', 'x1_prop'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("GSBART"), desc(models)))

print(sim.summary)

save(sim.list, sim.summary, file = "~/Downloads/GSBart_rebuttal/Ushape/GSBUshaperegres.Rdata")

sim.summary1
sim.summary

load("~/Downloads/GSBart_rebuttal/Ushape/GSBUshaperegnonsparse.Rdata")

# save(sim_list, sim_summary, file = "~/Downloads/code/test_code/testFlexUshapeReg.Rdata")
load("/Users/shurenhe/Downloads/GSBart_rebuttal/Ushape/testFlexUshapeReg.Rdata")
load("/Users/shurenhe/Downloads/GSBart_rebuttal/Ushape/GSBUshapereg.Rdata")

library(ggplot2)
load("~/Downloads/GSBart_rebuttal/Ushape/GSBres10.Rdata")
load("~/Downloads/GSBart_rebuttal/Ushape/FLBres10.Rdata")
GSB_trace <- apply(GSBART_Fit$yhat.test, 1, function(X){mean((Y0_ho - X)^2)})
FLB_trace <- apply(Flex_BART_Y_out[201:600,], 1, function(X){mean((Y0_ho - X)^2)})
GSB_est = colMeans(GSBART_Fit$yhat.test)
FLB_est = colMeans(Flex_BART_Y_out[201:600, ])
FLB_est_half = colMeans(Flex_BART_Y_out[401:600, ])
FL_iter = c(1:length(FLB_trace)*10)
GS_iter = c(1:length(GSB_trace))

GSB_MAE = mean(abs(Y0_ho - GSB_est))
FLB_MAE = mean(abs(Y0_ho - FLB_est))

GSB_MSE = mean((Y0_ho - GSBART_Fit$yhat.test.mean)^2)
FLB_MSE = mean((Y0_ho - FLB_est)^2)
FLB_MSE_half = mean((Y0_ho - FLB_est_half)^2)

mean(coda::effectiveSize(GSBART_Fit$yhat.test))/2 # 3:36.29  6: 29.59985 6_1: 33.28917 6_2: 31.34395 10:31.711
mean(coda::effectiveSize(Flex_BART_Y_out[401:600,]))/20 # 3: 1.30 3_half: 2.62 6: 0.7348807 ; 10: 3.92

df_mse <- data.frame(
  Iteration = c(GS_iter),
  MSPE = c(GSB_trace),
  Method = c(rep("bold(GS-BART)",length(GSB_trace)))
)

text_data <- data.frame(
  Method = c("bold(GS-BART)"),
  label = c("MSPE:0.080; ESS: 32.45")
)

pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/Ushape/Ushape_GSB10_traceplot.pdf", width  = 6, height = 4)

p1 <- ggplot(df_mse, aes(x = Iteration, y = MSPE, color = Method)) +
  geom_line(col = 'black', alpha = 1) +
  # geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
  #                              xintercept = c(14)),  aes(xintercept = xintercept),
  #            linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_label(data = text_data,
             aes(x = 35, y = 0.265, label = label),
             color = "black", size = 4,
             fill = "white", 
             label.size = 0.4, 
             inherit.aes = FALSE)+
  theme_minimal(base_size = 14) +
  # coord_cartesian(ylim = c(0.095, 0.3)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15)  # Adjust size as needed
  ) +
  ggtitle("Ushape (sim 1): GS-BART")
p1
dev.off()

df_mse <- data.frame(
  Iteration = c(FL_iter),
  MSPE = c(FLB_trace),
  Method = c(rep("bold(flexBART)",length(FLB_trace)))
)

text_data <- data.frame(
  Method = c("bold(flexBART)"),
  label = c("MSPE: 0.183; ESS: 1.56")
)

mean(coda::effectiveSize(GSBART_Fit$yhat.test))/2 # 29.59985
mean(coda::effectiveSize(Flex_BART_Y_out[201:600,]))/40 #  0.7348807

mean(coda::effectiveSize(Flex_BART_Y_out[401:600,]))/20 #  0.7348807

pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/Ushape/Ushape_FLB10_traceplot.pdf", width  = 6, height = 4)

p1 <- ggplot(df_mse, aes(x = Iteration, y = MSPE, color = Method)) +
  geom_line(col = 'black', alpha = 1) +
  # geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
  #                              xintercept = c(14)),  aes(xintercept = xintercept),
  #            linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_label(data = text_data,
             aes(x = 600, y = 0.265, label = label),
             color = "black", size = 4,
             fill = "white", 
             label.size = 0.4, 
             inherit.aes = FALSE)+
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0.095, 0.27)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15)  # Adjust size as needed
  ) +
  ggtitle("Ushape (sim 2): flexBART")
p1
dev.off()

# 
# df_mse <- data.frame(
#   Iteration = c(GS_iter, FL_iter),
#   MSE = c(GSB_trace, FLB_trace),
#   Method = c(rep("bold(GS-BART)",length(GSB_trace)),rep("bold(Flex-BART)", length(FLB_trace)))
# )
# 
# df_m <- data.frame(
#   Y = c(Y0_ho, Y0_ho),
#   Y_hat = c(GSB_est, FLB_est),
#   Method = c(rep("bold(GS-BART)",length(Y0_ho)),rep("bold(Flex-BART)", length(Y0_ho)))
# )
# library(ggpubr)
# 
# text_data <- data.frame(
#   Method = c("bold(GS-BART)", "bold(Flex-BART)"),
#   Iteration = c(350, 1000),
#   MSE = rep(0.5,2),
#   label = "burn-in"
# )
# 
# text_data2 <- data.frame(
#   Method = c("bold(GS-BART)",
#              "bold(Flex-BART)"),
#   Y = rep(-6, 2),
#   Y_hat = rep(4, 2),
#   label = c("MSE = 0.094",
#             "MSE = 0.107")
# )
# 
# text_data3 <- data.frame(
#   Method = c("bold(GS-BART)",
#              "bold(Flex-BART)"),
#   Y = rep(-6, 2),
#   Y_hat = rep(3, 2),
#   label = c("MAE = 0.215",
#             "MAE = 0.241")
# )

# pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/Ushape/Ushape_GSB_FLB3_compare_new.pdf", width  = 8, height = 7)
# p1 <- ggplot(df_mse, aes(x = Iteration, y = MSE, color = Method)) +
#   geom_line(alpha = 0.8) +
#   facet_wrap(~ Method, scales = "free_x", labeller = label_parsed) +
#   geom_vline(data = data.frame(Method = c("bold(GS-BART)", "bold(Flex-BART)"),
#              xintercept = c(700, 2000)),  aes(xintercept = xintercept),
#              linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
#   geom_text(data = text_data,
#             aes(x = Iteration, y = MSE, label = label),
#             color = "blue", size = 6,
#             inherit.aes = FALSE) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
# 
# p2 <- ggplot(df_m, aes(x = Y_hat, y = Y, color = Method)) +
#   geom_point(alpha = 0.8) +
#   facet_wrap(~ Method, scales = "fixed", labeller = label_parsed) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_text(data = text_data2,
#             aes(x = Y, y = Y_hat, label = label),
#             color = "black", size = 4,
#             inherit.aes = FALSE) +
#   geom_text(data = text_data3,
#             aes(x = Y, y = Y_hat, label = label),
#             color = "black", size = 4,
#             inherit.aes = FALSE) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none") +
#   labs(x = expression(hat(Y)), y = expression(Y))
# ggarrange(p1, p2, nrow = 2, ncol = 1, widths = c(3), heights = c(5, 5))
# dev.off()

# load("~/Downloads/code/test_code/Ushapesim.Rdata")
# n_rounds = 50
# nref = 100; num_ost = 7
# hyperpar = c("max_splits" = 20,
#              "lr" = 6/n_rounds,
#              "n_rounds" = n_rounds,
#              "lambda" = 0,
#              "alpha" = 0,
#              "nref" = nref,
#              "num_ost" = num_ost)
# 
# noise = 0.15
# n_train = 800
# Y0 = sim$f_true + rnorm(n_train, 0, noise)
# 
# Standardize = function(Y, std_par = NULL) {
#   if(is.null(std_par)) {
#     ymean = mean(Y)
#     yscale = 2 * max(abs(Y))
#   } else {
#     ymean = std_par['mean']
#     yscale = std_par['scale']
#   }
#   Y = (Y - ymean) / yscale
#   std_par = c('mean' = ymean, 'scale' = yscale)
#   return(list(Y = Y, std_par = std_par))
# }
# 
# Unstandardize = function(x, std_par) {
#   return(x * std_par['scale'] + std_par['mean']) 
# } 
# 
# stdY= Standardize(Y0)
# 
# for(i in 1:length(sim$Graphs)){
#   for(j in 1:length(sim$Graphs[[i]])){
#     sim$Graphs[[i]][[j]]$mesh2ho <- NULL
#   }
# }
# 
# sim$Graphs[[1]][[1]]$mesh2ho
# tmp <- GGBoost::GGBoostTrain(stdY$Y, sim$Graphs, hyperpar, loss = "mse", threads = 7, verbose = TRUE)
# 
