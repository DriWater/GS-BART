rm(list=ls())
library(G2SBart)

# US Cancer
load("data/USCancer.Rdata")
USCancer_Time = Sys.time()
USCancer_Fit <- g2sbart(Y, Graphs, 200, 15, 200, "poisson", graphs_weight, hyperpar, modelpar, offset_train = 5*exp(offset_train), 
                      offset_test = 5*exp(offset_test), nthreads = 9, verbose = T, seed = 1234)
USCancer_Time = difftime(Sys.time(), USCancer_Time, units = "secs")
Y_pred_ho <- colMeans(t(t(exp(USCancer_Fit$phi.test))* 5 * exp(offset_test)))
USCancer_RMSE = sqrt(mean((Y_ho - Y_pred_ho)^2))



