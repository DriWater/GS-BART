rm(list=ls())
library(BART) 
library(doRNG)
library(doParallel)
library(foreach)
library(purrr)
library(data.tree)
library(igraph)
library(matrixStats)
library(Matrix)
library(igraph)
library(dplyr)
library(xgboost)
library(BRISC)
library(caret)
library(sf)
library(G2SBart)

setwd("/Users/shurenhe/Downloads/code/GSBart_test/Cancer/")
load("US_CancerMortality.Rdata")
load("/Users/shurenhe/Downloads/code/GSBart_test/Cancer/USCancer_new.RData")

source('~/Downloads/code/GraphCon/GSBart_graph_fun.R')
X_all = rbind(X, X_ho)
X_all[in_train,]= X
X_all[-in_train,] = X_ho
n_rounds = 50
Graphs = Generate_GraphList_parallel(sf_coords= NULL, sf_bnd = NULL, g0 = g0, X = X_all, in_train = in_train, nrounds = n_rounds,
                            nost = 5, nbin = 100, nref = 100, nthreads = 5, chain.type = "stratified")

library(G2SBart)
# load('~/Downloads/code/test_code/USCancer.Rdata')
# 
# InputGraphs = lapply(Graphs, function(x) {
#   lapply(x, function(y) {
#     return(list(root = y$root,
#                 children = y$children,
#                 mesh2tr = igraph::vertex.attributes(y$g)$mesh2tr,
#                 mesh2ho = igraph::vertex.attributes(y$g)$mesh2ho,
#                 tr2mesh = igraph::graph.attributes(y$g)$tr2mesh))
#   }) 
# })
# 
# InputGraphs = lapply(InputGraphs, function(x) {
#   lapply(x, function(y) {
#     list(root = y$root - 1,
#          children = unname(lapply(y$children, function(z){unname(z) - 1})),
#          mesh2tr = unname(lapply(y$mesh2tr, function(z){unname(z) - 1})),
#          mesh2ho = unname(lapply(y$mesh2ho, function(z){unname(z) - 1})),
#          tr2mesh = unname(y$tr2mesh - 1))
#   }) 
# })

n_rounds = length(Graphs)
n_ost=7
n_train = length(Y)
n_test = length(Y_ho)
p = length(Graphs[[1]]) - n_ost


alpha=0.95
beta=2
nu = 3; q = 0.9
rb=.9
split_eps=0
quant = qchisq(1-q, nu)

hyperpar=c(
  alpha=alpha,
  beta=beta,
  nu = nu,
  rb=rb,
  split_eps=split_eps,
  max_depth = 6,
  tree_iter = 15,
  a = 3,
  mu0 = mean(log(Y/(5*exp(offset_train))+.01))/n_rounds,
  b = 0.5*var(log(Y/(5*exp(offset_train))+.01))/n_rounds
)

graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/n_ost,n_ost), rep(min(p / (p + 2), 0.85)/p,p))

modelpar=c(tausq = (6/(2*sqrt(n_rounds)))^2)

hyperpar <- as.list(hyperpar); modelpar <- as.list(modelpar)

hyperpar['tree_iter'] = 15; hyperpar['max_depth'] = 6
hyperpar['b'] = 0.5*var(log(Y/(5*exp(offset_train))+.01))/n_rounds
hyperpar['mu0'] = mean(log(Y/(5*exp(offset_train))+.01))/n_rounds
modelpar['tausq']  = (6/(2*sqrt(n_rounds)))^2

GSBart_Time = Sys.time()
GSBart_Fit <- g2sbart(Y, Graphs, 200, 15, 200, "poisson", graphs_weight, hyperpar, modelpar, offset_train = 5*exp(offset_train), 
                      offset_test = 5*exp(offset_test), nthreads = 9, sparse = F, verbose = T, seed = 1234)
GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
Y_pred_ho <- colMeans(t(t(exp(GSBart_Fit$phi.test))* 5 * exp(offset_test)))
GSBart_RMSE = sqrt(mean((Y_ho - Y_pred_ho)^2))
c(GSBart_RMSE)

GSBart_Fit$time <- as.numeric(GSBart_Time)
print(GSBart_Fit$time)
plot(GSBart_Fit$tau[15:215], type = 'l')

geweke_test=function(sigma_out){
  geweke_results=(coda::geweke.diag(sigma_out))
  # Extract z-scores
  zvals <- geweke_results$z
  
  # Compute two-sided p-values
  pvals <- 2 * (1 - pnorm(abs(zvals)))
  return(pvals)
}

geweke_test(GSBart_Fit$tau[15:215])

save(GSBart_Fit, file = "/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/GSBart_Fit_50rounds_5.RData")

load("/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/GSBart_Fit_50rounds_5.RData")
# thin_idx = seq(140, 21400, 100)
# thin_idx = c(1400:21400)


GSB_trace <- sqrt(colMeans((Y_ho - t(exp(GSBart_Fit$phi.test))*5*exp(offset_test))^2))
GS_iter = c(1:length(GSB_trace))

# GSB_trace[which.max(GSB_trace)] = GSB_trace[which.max(GSB_trace)] - 100

df_mse <- data.frame(
  Iteration = c(GS_iter),
  RMSPE = c(GSB_trace),
  Method = c(rep("bold(GS-BART)",length(GSB_trace)))
)

text_data <- data.frame(
  Method = c("bold(GS-BART)"),
  label = c("RMSPE:22.90; ESS:81.08")
)

mean(coda::effectiveSize(GSBart_Fit$phi.test))/2

pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/USCancer_traceplot1.pdf", width  = 6, height = 4)

p1 <- ggplot(df_mse, aes(x = Iteration, y = RMSPE, color = Method)) +
  geom_line(col = 'black', alpha = 1) +
  # geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
  #                              xintercept = c(14)),  aes(xintercept = xintercept),
  #            linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_label(data = text_data,
             aes(x = 35, y = 87, label = label),
             color = "black", size = 4,
             fill = "white", 
             label.size = 0.4, 
             inherit.aes = FALSE)+
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(20, 90)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15)  # Adjust size as needed
  ) +
  ggtitle("USCancer: GS-BART (count)")
p1
dev.off()




sqrt(mean((Y_pred_ho- Y_ho)^2))

n_ost=7
n_train = length(Y)
n_test = length(Y_ho)
p = length(Graphs[[1]]) - n_ost


alpha=0.95
beta=2
nu = 3; q = 0.9
rb=.9
split_eps=0
quant = qchisq(1-q, nu)


hyperpar=c(
  alpha=alpha,
  beta=beta,
  nu = nu,
  rb=rb,
  split_eps=split_eps,
  max_depth = 6,
  mu0 = 0,
  tree_iter = 15,
  a = 3,
  b = 0.5*var(log(Y/(2*exp(offset_train))+.01))/n_rounds
)

graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/n_ost,n_ost), rep(min(p / (p + 2), 0.85)/p,p))

modelpar=c(sigmasq = (2/(2*sqrt(n_rounds)))^2, tausq = (6/(2*sqrt(n_rounds)))^2)

save(Y, Y_ho, offset_train, offset_test, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/USCancer.Rdata")


nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)

par_grid = expand.grid(
  sigmasq = seq(1,6, 0.5),
  a = seq(2,6,1)
)

repetitions = nrow(par_grid)

sim.list = parallel::mclapply(1:repetitions, function(repetition){
  modelpar['sigmasq'] = par_grid[repetition,1]
  hyperpar['a'] = par_grid[repetition,2]
  GSBart_Time = Sys.time()
  GSBart_Fit <- GenGSBart(Y/(2*exp(offset_train)), 25, 15, 1250, InputGraphs, graphs_weight,
                          hyperpar, modelpar, 2, 1, verbose = F, seed = 1234)
  GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
  GSBart_RMSE = sqrt(mean((Y_ho - exp(colMeans(GSBart_Fit$phi.test))*2*exp(offset_test))^2))
  sim_results_df = data.frame(
    models = c("GSBart"),
    MSPE = c(GSBart_RMSE),
    times = c(GSBart_Time),
    sim = rep(repetition, 1)
  )
  return(sim_results_df)
}, mc.cores = 8)


load('~/Downloads/code/test_code/USCancer.Rdata')

var(log(Y/2*exp(offset_train) + 0.1))


save(sim.list, file = "~/Downloads/GSBart_rebuttal/USCancer/GSBUSCancerres.Rdata")

load("~/Downloads/GSBart_rebuttal/USCancer/GSBUSCancerres1.Rdata")
rmsevec = vector(mode = 'numeric', length = length(sim.list))
for(i in 1:length(sim.list)){
  rmsevec[i] = sim.list[[i]]$MSPE
}

stop(cl)

library(FlexBartVar)
MCMC = 3000
BURNIN = 1000
THIN = 5
ndpost = MCMC - BURNIN
nskip = BURNIN

X_all = rbind(X, X_ho)
X_all_std <- apply(X_all, 2, function(x){(x - min(x)) / (max(x) - min(x))})

Y_std <- Y/(5*exp(offset_train))
# sdy = sd(Y_std); meany = mean(Y_std)
# Y_std = (Y_std - meany)/sdy

tmp1 = lm(Y_std ~ X_all_std[1:2483, ])

FlexBart_Time=Sys.time()
FlexBart_Fit =  FlexBartVar::FlexBart(X_all_std[1:2483, ], Y/(5*exp(offset_train)), X_all_std[2484:3103, ], probs = Matrix(diag(ncol(X_all_std))), sigma_hat = summary(tmp1)$sigma, sigma_mu = sqrt((3/(2*sqrt(100)))^2), kappa = 1, 
                         update_sigma_mu = T, update_kappa = T, num_trees = 50, num_burn = 2000, num_thin = 1, num_save = 4000, num_print = 100)
FlexBart_Time=difftime(Sys.time(), FlexBart_Time, units = "secs")
FlexBart_Fit$time <- FlexBart_Time

FlexBart_pred_unstandardized = colMeans(exp(FlexBart_Fit$mu_test))
Flex_pred_ho = (FlexBart_pred_unstandardized) *(5*exp(offset_test))
FlexBart_RMSE = sqrt(mean((Flex_pred_ho - Y_ho)^2))

save(FlexBart_Fit, file = "/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/FlexBart_Fit50rounds_5.RData")

load("/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/FlexBart_Fit.RData")

mean(coda::effectiveSize(FlexBart_Fit$mu_test))

Flex_trace <- colMeans((Y_ho - t(exp(FlexBart_Fit$mu_test))*5*exp(offset_test))^2)
Flex_iter = c(1:length(Flex_trace))

Flex_df_mse <- data.frame(
  Iteration = c(Flex_iter),
  RMSE = c(sqrt(Flex_trace)),
  Method = c(rep("bold(BART (count))",length(Flex_iter)))
)

text_data <- data.frame(
  Method = c("bold(GS-BART)"),
  label = c("RMSPE:32.30; ESS:1.05")
)

pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/USCancer/USCancer_BART_poi_traceplot1.pdf", width  = 6, height = 4)

p1 <- ggplot(Flex_df_mse, aes(x = Iteration, y = RMSE, color = Method)) +
  geom_line(col = 'black', alpha = 1) +
  # geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
  #                              xintercept = c(14)),  aes(xintercept = xintercept),
  #            linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_label(data = text_data,
             aes(x = 650, y = 540, label = label),
             color = "black", size = 4,
             fill = "white", # 背景色
             label.size = 0.4, # 边框粗细，设为0可以去掉边框
             inherit.aes = FALSE)+
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 560)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15)  # Adjust size as needed
  ) +
  ggtitle("USCancer: BART (count)")
p1
dev.off()
