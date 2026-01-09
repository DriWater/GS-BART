rm(list=ls())
library(GSBart)

# NYC Education
load("data/NYEducation.Rdata")
NYEdu_Time=Sys.time()
NYEdu_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, 
                   verbose = F, nthreads = 6, seed = 1234) 
NYEdu_Time=difftime(Sys.time(), NYEdu_Time, units = "secs")
NYEdu_MSPE = mean((NYEdu_GSB$yhat.test.mean - y.test.unstandardized)^2)
NYEdu_MAPE = mean(abs(NYEdu_GSB$yhat.test.mean - y.test.unstandardized))

# King County House
load("data/KingHouse.Rdata")
KingHouse_Time=Sys.time()
KingHouse_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, verbose = F, 
                       nthreads = 8, seed = 1234)
KingHouse_Time=difftime(Sys.time(), KingHouse_Time, units = "secs")
KingHouse_MSPE = mean((KingHouse_GSB$yhat.test.mean - y.test.unstandardized)^2)
KingHouse_MAPE = mean(abs(KingHouse_GSB$yhat.test.mean - y.test.unstandardized))

# US Election
load("data/USelection.Rdata")
USElection_Time=Sys.time()
USElection_GSB = gsbart(y.train.unstandardized, Graphs, 200, 15, 200, graphs_weight = graphs_weight, 
                        nthreads = 8, seed = 1239)
USElection_Time=difftime(Sys.time(), USElection_Time, units = "secs")
USElection_MSPE = mean((USElection_GSB$yhat.test.mean - y.test.unstandardized)^2)
USElection_MAPE = mean(abs(USElection_GSB$yhat.test.mean - y.test.unstandardized))

real_continuous_res = data.frame(MSPE = c(NYEdu_MSPE, KingHouse_MSPE, USElection_MSPE),
                                 MAPE = c(NYEdu_MAPE, KingHouse_MAPE, USElection_MAPE))
rownames(real_continuous_res) <- c('NYC Education', 'KingHouse', 'US Election')

save(real_continuous_res, file = 'data/real_continuous_res.Rdata', compress = 'xz')
