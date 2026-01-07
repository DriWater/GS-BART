rm(list=ls())
library(G2SBart)

# US Cancer
load("data/USFlood.Rdata")
USFlood_Time = Sys.time()
USFlood_Fit <- G2SBart::g2sbart(Y, Graphs, 200, 15, 200, family = "binary", graphs_weight,
                                hyperpar = hyperpar, modelpar = modelpar,
                                nthreads = 1, verbose = F, seed = 1234)
USFlood_Time = difftime(Sys.time(), USFlood_Time, units = "secs")
Y_ho_pred <- as.numeric(USFlood_Fit$prob.test.mean > 0.5)
USFlood_ACC = mean(Y_ho==Y_ho_pred)


load("data/USAir.Rdata")
USAir_Time = Sys.time()
USAir_Fit <- G2SBart::g2sbart(Y_bin, Graphs, 200, 15, 200, family = "binary", graphs_weight,
                              nthreads = 1, verbose = F, seed = 1234)
USAir_Time = difftime(Sys.time(), USAir_Time, units = "secs")
Y_ho_pred <- as.numeric(USAir_Fit$prob.test.mean > 0.5)
Y_bin_ho_tmp <- na.omit(Y_bin_ho)
remove.id = which(is.na(Y_bin_ho))
Y_ho_pred <- Y_ho_pred[-remove.id]
USAir_ACC = mean(Y_bin_ho_tmp==Y_ho_pred)

real_classification_res = data.frame(ACC = c(USFlood_ACC, USAir_ACC))
rownames(real_classification_res) <- c('Flood', 'Air Pollution')

save(real_classification_res, file = 'data/real_classification_res.Rdata', compress = 'xz')
