rm(list=ls())
library(G2SBart)

# US Cancer
load("data/USFlood.Rdata")
USFlood_Time = Sys.time()
USFlood_Fit <- G2SBart::g2sbart(Y, Graphs, 200, 15, 200, family = "binary", graphs_weight,
                                hyperpar = hyperpar, modelpar = modelpar,
                                nthreads = 9, verbose = F, seed = 1234)
USFlood_Time = difftime(Sys.time(), USFlood_Time, units = "secs")
Y_ho_pred <- as.numeric(USFlood_Fit$prob.test.mean > 0.5)
USFlood_ACC = mean(Y_ho==Y_ho_pred)



