rm(list=ls())
library(dplyr)
library(caret)
library(purrr)
library(testthat)
library(parallel)
library(doParallel)
library(flexBART)
library(GSBart)

source('/Users/shurenhe/Downloads/code/GSBart/GSBart_graph_fun.R')

p = 5
f <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

set.seed(1234)
sigma <- 1
n <- 500

niter=6000
nskip=2000
ndpost = niter - nskip
thin=10
n_rounds = 50


x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", 1:p)
Ey <- f(x)

# load('/Users/shurenhe/Downloads/code/GSBart_test/Friedman Example/Friedman_test.Rdata')
in_train <- createDataPartition(Ey, p = .8, list = FALSE)
X = x[in_train,]; X_ho = x[-in_train,]
n_train = nrow(X); n_test = nrow(X_ho)

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
# function to unstandardize
Unstandardize = function(x, std_par) {
  return(x * std_par['scale'] + std_par['mean'])
}

Graphs = lapply(1:n_rounds, function(x){Generate_GraphCandidates_graphbin(X = X, X_ho = X_ho, n_bin = 100, chain.type = 'regular')})

alpha=0.95
beta=2
nu = 3; q = 0.9
rb=.5
quant = qchisq(1-q, nu)

hyperpar=c(
  alpha=0.95,
  beta=2,
  nu = 3,
  mu0 = 0,
  rb=0.9,
  split_eps=0,
  max_depth = 6,
  tree_iter = 12,
  a = 3
)

graphs_weight = rep(0.2,5)


modelpar=c(tausq = (0.5/(2*sqrt(n_rounds)))^2)
hyperpar <- as.list(hyperpar); modelpar <- as.list(modelpar)

InputGraphs = lapply(Graphs, function(x) {
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

save(Ey, X, X_ho, in_train, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/testFriedmanReg.Rdata")

load("~/Downloads/code/test_code/testFriedmanReg.Rdata")


nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)

repetitions = 50
sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'data.tree','matrixStats','BART', 'bayestestR')) %dopar% {
  message("Sampling Y")
  repetition = 7
  set.seed(1234+repetition)

  set.seed(1234)
  y <- rnorm(n, Ey, sigma)
  Y0=y[in_train]; Y0_ho=y[-in_train]
  
  GSBart_Time=Sys.time()
  # GSBART_Fit = gsbart(Y0, InputGraphs, 200, 15, 200, graphs_weight, nthreads = 5, seed = 1234, verbose = T)
  GSBART_Fit = gsbart(Y0, InputGraphs, 200, 15, 200, graphs_weight, nthreads = 5, seed = 1234, verbose = T)
  GSBart_Time=difftime(Sys.time(), GSBart_Time, units = "secs")

  ## standardize data
  stdY = Standardize(Y0)
  Y = stdY$Y
  Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']

  m = ncol(X)
  if(m < n_train) {
    df = data.frame(X ,Y)
    lmf = lm(Y~.,df)
    sigest = summary(lmf)$sigma
  } else {
    sigest = sd(Y)
  }

  hyperpar['lambda_s'] = (sigest*sigest*quant)/nu
  hyperpar['prob_split_by_x'] = 1  # probability for a split on x
  hyperpar['b'] = 0.5*sigest*sigest/n_rounds
  modelpar['sigmasq_y']=sigest*sigest
  ndpost=50; nskip=15; nkeep=50;
  nthreads=7

  GSB_fit = GSBart::gsbart(Y0, Graphs, ndpost=50, nskip=15, nkeep=50, nthreads =nthreads,
                           verbose=T, seed=1234)


  GSBart_MSPE=mean((Y0_ho-GSB_fit$yhat.test.mean)^2)
  GSBart_MAPE=mean(abs(GSB_fit$yhat.test.mean-Y0_ho))

  # Fit BART
  BART_Time=Sys.time()
  BART_Fit = wbart(X, Y, X_ho, sparse = F, augment = F, ntree = 100, ndpost = 4000, nskip = 2000)
  BART_Time = difftime(Sys.time(), BART_Time, units = "secs")
  BART_Fit$time = BART_Time
  BART_Y_out = t(apply(BART_Fit$yhat.test[seq(1, ndpost, thin), ], 1, function(X){Unstandardize(X, stdY$std_par)}))
  BART.pred.unstandardized = colMeans(BART_Fit$yhat.test[seq(1, 4000, 10), ])
  BART_MSPE = mean((BART.pred.unstandardized - Y_ho)^2)
  BART_MAPE = mean(abs(BART.pred.unstandardized - Y_ho))
  save(BART_Fit, BART_Y_out , BART.pred.unstandardized, Y0, Y0_ho, file = "~/Downloads/GSBart_rebuttal/Friedman/BARTres7.Rdata")


  # Fit GSBart Here
  # GSBart_Time=Sys.time()
  # GSB_res = suppressWarnings(GSBart(Graphs = Graphs, y.train = Y, modelpar = modelpar, hyperpar =  hyperpar, 115, 15, verbose = T, nthreads = 1))
  # GSBart_Time=difftime(Sys.time(), GSBart_Time, units = "secs")
  # GSB_res$time = GSBart_Time
  # GSBart_Y_out = t(apply(GSB_res$Y_new_out, 1, function(X){Unstandardize(X, stdY$std_par)}))
  # GSBart_pred_unstandardized = colMeans(GSBart_Y_out)
  # GSBart_MSPE = mean((GSBart_pred_unstandardized - Y0_ho)^2)
  # GSBart_MAPE = mean(abs(GSBart_pred_unstandardized - Y0_ho))

  # Fit XBART
  # XBART_Time = Sys.time()
  # XBART_Fit = XBART(matrix(Y), X, num_trees = 50, num_sweeps = 115, burnin = 15)
  # XBART_Time = difftime(Sys.time(), XBART_Time, units = "secs")
  # XBART_Y_pred = predict(XBART_Fit, as.matrix(X_ho))
  # XBART_Y_out = t(apply(XBART_Y_pred, 1, function(X){Unstandardize(X, stdY$std_par)}))
  # XBART_pred_unstandardized = colMeans(XBART_Y_out)
  # XBART_MSPE = mean((XBART_pred_unstandardized - Y0_ho)^2)
  # XBART_MAPE = mean(abs(XBART_pred_unstandardized - Y0_ho))


  BART.upper.lvl = NULL
  BART.lower.lvl = NULL
  for(i in 1:ncol(BART_Y_out)){
    tmp <- ci(BART_Y_out[,i], method = "HDI")
    BART.upper.lvl <- append(BART.upper.lvl, tmp$CI_high)
    BART.lower.lvl <- append(BART.lower.lvl, tmp$CI_low)
  }

  # GSBart.upper.lvl = NULL
  # GSBart.lower.lvl = NULL
  # for(i in 1:ncol(GSBart_Y_out)){
  #   tmp <- ci(GSBart_Y_out[,i], method = "HDI")
  #   GSBart.upper.lvl <- append(GSBart.upper.lvl, tmp$CI_high)
  #   GSBart.lower.lvl <- append(GSBart.lower.lvl, tmp$CI_low)
  # }

  # XBART.upper.lvl = NULL
  # XBART.lower.lvl = NULL
  # for(i in 1:ncol(XBART_Y_out)){
  #   tmp <- ci(XBART_Y_out[,i], method = "HDI")
  #   XBART.upper.lvl <- append(XBART.upper.lvl, tmp$CI_high)
  #   XBART.lower.lvl <- append(XBART.lower.lvl, tmp$CI_low)
  # }

  BART_coverage = mean((Y0_ho<BART.upper.lvl)&(Y0_ho>BART.lower.lvl))
  # GSBart_coverage = mean((Y0_ho<GSBart.upper.lvl)&(Y0_ho>GSBart.lower.lvl))
  # XBART_coverage = mean((Y0_ho<XBART.upper.lvl)&(Y0_ho>XBART.lower.lvl))

  BART_HDI_len = mean(BART.upper.lvl - BART.lower.lvl)
  # GSBart_HDI_len = mean(GSBart.upper.lvl - GSBart.lower.lvl)
  # XBART_HDI_len = mean(XBART.upper.lvl - XBART.lower.lvl)

  BART_poster_sd = mean(apply(BART_Y_out, 2, sd))
  # GSBart_poster_sd = mean(apply(GSBart_Y_out, 2, sd))
  # XBART_poster_sd = mean(apply(XBART_Fit$pred_res, 2, sd))


  # sim_results_df = data.frame(
  #   models = c("BART", "GSBart","XBART"),
  #   MSPE = c(BART_MSPE, GSBart_MSPE, XBART_MSPE),
  #   MAPE = c(BART_MAPE, GSBart_MAPE, XBART_MAPE),
  #   HDI_coverage = c(BART_coverage,GSBart_coverage,XBART_coverage),
  #   HDI_len = c(BART_HDI_len,GSBart_HDI_len,XBART_HDI_len),
  #   poster_sd = c(BART_poster_sd,GSBart_poster_sd,XBART_poster_sd),
  #   times = c(as.numeric(BART_Time), as.numeric(GSBart_Time), as.numeric(XBART_Time)),
  #   sim = rep(repetition, 3)
  # )

  # sim_results_df = data.frame(
  #   models = c("BART", "XBART"),
  #   MSPE = c(BART_MSPE, XBART_MSPE),
  #   MAPE = c(BART_MAPE, XBART_MAPE),
  #   HDI_coverage = c(BART_coverage, XBART_coverage),
  #   HDI_len = c(BART_HDI_len, XBART_HDI_len),
  #   poster_sd = c(BART_poster_sd, XBART_poster_sd),
  #   times = c(as.numeric(BART_Time), as.numeric(XBART_Time)),
  #   sim = rep(repetition, 2)
  # )

  sim_results_df = data.frame(
    models = c("BART"),
    MSPE = c(BART_MSPE),
    MAPE = c(BART_MAPE),
    HDI_coverage = c(BART_coverage),
    HDI_len = c(BART_HDI_len),
    poster_sd = c(BART_poster_sd),
    times = c(as.numeric(BART_Time)),
    sim = rep(repetition, 1)
  )

  sim_results_df
}

n_train = nrow(X); n_test = nrow(X_ho)
repetitions = 50
sim_list <- foreach(repetition = 1:repetitions, .packages= c('Matrix','igraph', 'flexBART','data.tree','matrixStats','BART', 'bayestestR')) %dopar% {
  message("Sampling Y")
  set.seed(1234+repetition)

  y <- rnorm(n, Ey, sigma)
  Y0=y[in_train]; Y0_ho=y[-in_train]

  ## standardize data
  stdY = Standardize(Y0)
  Y = stdY$Y
  Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']

  # Fit BART
  # BART_Time=Sys.time()
  # BART_Fit = wbart(X, Y, X_ho, sparse = F, augment = F, ntree = 50, ndpost = ndpost, nskip = nskip)
  # BART_Time = difftime(Sys.time(), BART_Time, units = "secs")
  # BART_Fit$time = BART_Time
  # BART_Y_out = t(apply(BART_Fit$yhat.test[seq(1, ndpost, thin), ], 1, function(X){Unstandardize(X, stdY$std_par)}))
  # BART.pred.unstandardized = colMeans(BART_Y_out)
  # BART_MSPE = mean((BART.pred.unstandardized - Y0_ho)^2)
  # BART_MAPE = mean(abs(BART.pred.unstandardized - Y0_ho))
  X_all <- rbind(X, X_ho)
  X_all_std <- apply(X_all, 2, function(x){(x - min(x)) / (max(x) - min(x))})
  X_train_std = X_all_std[1:n_train,]
  X_test_std = X_all_std[(n_train+1):(n_train+n_test),]

  # Flex-BART
  Flex_BART_Time=Sys.time()
  Flex_BART_Fit = flexBART(Y, X_cont_train = X_train_std, X_cont_test = X_test_std, M = 100, nd = 400, burn = 2000, thin = 10)
  Flex_BART_Time = difftime(Sys.time(), Flex_BART_Time, units = "secs")
  Flex_BART_Y_out = t(apply(Flex_BART_Fit$yhat.test, 1, function(X){Unstandardize(X, stdY$std_par)}))
  Flex_BART_pred_unstandardized = colMeans(Flex_BART_Y_out)
  Flex_BART_MSPE = mean((Flex_BART_pred_unstandardized - Y0_ho)^2)
  Flex_BART_MAPE = mean(abs(Flex_BART_pred_unstandardized - Y0_ho))
  # save(Flex_BART_Fit, Flex_BART_Y_out , Flex_BART_pred_unstandardized, Y0, Y0_ho, file = "~/Downloads/GSBart_rebuttal/Friedman/FLBres33.Rdata")

  # Fit GSBart Here
  # GSBart_Time=Sys.time()
  # GSB_res = suppressWarnings(GSBart(Graphs = Graphs, y.train = Y, modelpar = modelpar, hyperpar =  hyperpar, 115, 15, verbose = T, nthreads = 1))
  # GSBart_Time=difftime(Sys.time(), GSBart_Time, units = "secs")
  # GSB_res$time = GSBart_Time
  # GSBart_Y_out = t(apply(GSB_res$Y_new_out, 1, function(X){Unstandardize(X, stdY$std_par)}))
  # GSBart_pred_unstandardized = colMeans(GSBart_Y_out)
  # GSBart_MSPE = mean((GSBart_pred_unstandardized - Y0_ho)^2)
  # GSBart_MAPE = mean(abs(GSBart_pred_unstandardized - Y0_ho))

  # Fit XBART
  # XBART_Time = Sys.time()
  # XBART_Fit = XBART(matrix(Y), X, num_trees = 50, num_sweeps = 115, burnin = 15)
  # XBART_Time = difftime(Sys.time(), XBART_Time, units = "secs")
  # XBART_Y_pred = predict(XBART_Fit, as.matrix(X_ho))
  # XBART_Y_out = t(apply(XBART_Y_pred, 1, function(X){Unstandardize(X, stdY$std_par)}))
  # XBART_pred_unstandardized = colMeans(XBART_Y_out)
  # XBART_MSPE = mean((XBART_pred_unstandardized - Y0_ho)^2)
  # XBART_MAPE = mean(abs(XBART_pred_unstandardized - Y0_ho))


  Flex.BART.upper.lvl = NULL
  Flex.BART.lower.lvl = NULL
  for(i in 1:ncol(Flex_BART_Y_out)){
    tmp <- ci(Flex_BART_Y_out[,i], method = "HDI")
    Flex.BART.upper.lvl <- append(Flex.BART.upper.lvl, tmp$CI_high)
    Flex.BART.lower.lvl <- append(Flex.BART.lower.lvl, tmp$CI_low)
  }

  # GSBart.upper.lvl = NULL
  # GSBart.lower.lvl = NULL
  # for(i in 1:ncol(GSBart_Y_out)){
  #   tmp <- ci(GSBart_Y_out[,i], method = "HDI")
  #   GSBart.upper.lvl <- append(GSBart.upper.lvl, tmp$CI_high)
  #   GSBart.lower.lvl <- append(GSBart.lower.lvl, tmp$CI_low)
  # }

  # XBART.upper.lvl = NULL
  # XBART.lower.lvl = NULL
  # for(i in 1:ncol(XBART_Y_out)){
  #   tmp <- ci(XBART_Y_out[,i], method = "HDI")
  #   XBART.upper.lvl <- append(XBART.upper.lvl, tmp$CI_high)
  #   XBART.lower.lvl <- append(XBART.lower.lvl, tmp$CI_low)
  # }

  Flex_BART_coverage = mean((Y0_ho<Flex.BART.upper.lvl)&(Y0_ho>Flex.BART.lower.lvl))
  # GSBart_coverage = mean((Y0_ho<GSBart.upper.lvl)&(Y0_ho>GSBart.lower.lvl))
  # XBART_coverage = mean((Y0_ho<XBART.upper.lvl)&(Y0_ho>XBART.lower.lvl))

  Flex_BART_HDI_len = mean(Flex.BART.upper.lvl - Flex.BART.lower.lvl)
  # GSBart_HDI_len = mean(GSBart.upper.lvl - GSBart.lower.lvl)
  # XBART_HDI_len = mean(XBART.upper.lvl - XBART.lower.lvl)

  Flex_BART_poster_sd = mean(apply(Flex_BART_Y_out, 2, sd))
  # GSBart_poster_sd = mean(apply(GSBart_Y_out, 2, sd))
  # XBART_poster_sd = mean(apply(XBART_Fit$pred_res, 2, sd))


  # sim_results_df = data.frame(
  #   models = c("BART", "GSBart","XBART"),
  #   MSPE = c(BART_MSPE, GSBart_MSPE, XBART_MSPE),
  #   MAPE = c(BART_MAPE, GSBart_MAPE, XBART_MAPE),
  #   HDI_coverage = c(BART_coverage,GSBart_coverage,XBART_coverage),
  #   HDI_len = c(BART_HDI_len,GSBart_HDI_len,XBART_HDI_len),
  #   poster_sd = c(BART_poster_sd,GSBart_poster_sd,XBART_poster_sd),
  #   times = c(as.numeric(BART_Time), as.numeric(GSBart_Time), as.numeric(XBART_Time)),
  #   sim = rep(repetition, 3)
  # )

  # sim_results_df = data.frame(
  #   models = c("BART", "XBART"),
  #   MSPE = c(BART_MSPE, XBART_MSPE),
  #   MAPE = c(BART_MAPE, XBART_MAPE),
  #   HDI_coverage = c(BART_coverage, XBART_coverage),
  #   HDI_len = c(BART_HDI_len, XBART_HDI_len),
  #   poster_sd = c(BART_poster_sd, XBART_poster_sd),
  #   times = c(as.numeric(BART_Time), as.numeric(XBART_Time)),
  #   sim = rep(repetition, 2)
  # )

  sim_results_df = data.frame(
    models = c("Flex-BART"),
    MSPE = c(Flex_BART_MSPE),
    MAPE = c(Flex_BART_MAPE),
    HDI_coverage = c(Flex_BART_coverage),
    HDI_len = c(Flex_BART_HDI_len),
    poster_sd = c(Flex_BART_poster_sd),
    times = c(as.numeric(Flex_BART_Time)),
    sim = rep(repetition, 1)
  )

  sim_results_df
}



sim_summary = sim_list %>%
  list_rbind() %>%
  na.omit() %>%
  group_by(models) %>%
  summarise(across(c('MSPE', 'MAPE', 'HDI_coverage', 'HDI_len', 'poster_sd', 'times'),
                   list("mean" = mean, "sd" = sd),
                   .names = "{.col}_{.fn}")) %>%
  arrange(match(models, c("Flex-BART"), desc(models)))
sim_summary

load("~/Downloads/GSBart_rebuttal/Friedman/GSBFriedmanReg.Rdata")
