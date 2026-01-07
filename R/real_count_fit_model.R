rm(list=ls())
library(G2SBart)

# US Cancer
load("data/USCancer.Rdata")
USCancer_Time = Sys.time()
USCancer_Fit <- g2sbart(Y, Graphs, 200, 15, 200, "poisson", graphs_weight, hyperpar, modelpar, offset_train = 2*exp(offset_train), 
                        offset_test = 2*exp(offset_test), nthreads = 1, verbose = T, seed = 1234)
USCancer_Time = difftime(Sys.time(), USCancer_Time, units = "secs")
Y_pred_ho <- colMeans(t(t(exp(USCancer_Fit$phi.test))* 2 * exp(offset_test)))
USCancer_RMSPE = sqrt(mean((Y_ho - Y_pred_ho)^2))

real_count_res = data.frame(RMSPE = c(USCancer_RMSPE))
rownames(real_count_res) <- c('US Cancer')

save(real_count_res, file = 'data/real_count_res.Rdata', compress = 'xz')

