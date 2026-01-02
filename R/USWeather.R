rm(list=ls())
library(BART)
library(data.tree)
library(igraph)
library(dplyr)
library(ggplot2)
library(sf)
library(network)
library(ComplexDomain)
library(caret)
library(matrixStats)
library(caret)
library(meshed)
library(glmmfields)
# library(tidycensus)
set.seed(1234)
setwd("/Users/shurenhe/Downloads/GSBart_rebuttal/Airpollution/")
load("AsthmaCountyProcessed.Rdata")
source("/Users/shurenhe/Downloads/code/GSBart/GSBart_graph_fun.R")

Y <- UScountydata_processed$Unhealthy.Days
non.na.id <- which(!is.na(Y))
in_train = createDataPartition(non.na.id, p = .8, list = F)
in_train = non.na.id[in_train]

X <- UScountydata_processed[,c(6:34)]
X <- X %>% st_drop_geometry()

sf_bnd=st_geometry(UScountydata_processed)
sf_coords_all=st_coordinates(st_centroid(sf_bnd))
X[,c(1,2)] <- sf_coords_all
nn.list= st_relate(st_geometry(UScountydata_processed$geometry), pattern = "****1****")
g0 = graph.adjlist(nn.list) %>% as.undirected() %>% igraph::simplify()

Y_bin <- Y
Y_bin[which(Y>0)] <- 1

Y_bin_ho <- Y_bin[-in_train]
Y_bin <- Y_bin[in_train]
X_ho <- X[-in_train, ]
X <- X[in_train, ]
X[,3:29] <- X[,3:29]/UScountydata_processed$population[in_train]
X_ho[,3:29] <- X_ho[,3:29]/UScountydata_processed$population[in_train]

n_ref = 100
n_bin = 100
n_rounds = 100

Graphs = vector(mode ='list',length =  n_rounds)

for( i in 1:n_rounds){

  g0_binned=StructurizeGraph (g0, in_train, bins = n_ref,

                              method = "greedy",

                              connect_method = "random")



  output=Generate_GraphCandidates(g0_binned, X, X_ho, num.ost = 5, n_bin)

  Graphs[[i]] <- output

}

# save(sf_bnd, g0, X, X_ho, Y_bin, Y_bin_ho, Graphs, non.na.id, in_train, file = "USWeather.RData")
load("USWeather.RData")

n_rounds = 50
sigmasq_mu=(3/(2*sqrt(n_rounds)))^2
# sigmasq_y=var(Y)
alpha=0.95
beta=2
nu = 3; q = 0.9
rb=.9
split_eps=1
quant = qchisq(1-q, nu)
n_ost = 7
n_chain = 27
# lambda_s = (var(Y)*quant)/nu
a = 3;
b = 0.03
# a = 3;
# b = 0.3*3/n_rounds
p = mean(Y_bin==1)
mu_mu = qlogis(p)/n_rounds
max_depth = 6
iter_pre_tree = 15
prob_split_by_x = min(n_chain / (n_chain + 2), 0.85)

hyperpar=c(
  n_rounds=n_rounds,
  n_ost=n_ost,
  alpha=alpha,
  beta=beta,
  nu = nu, q = q,
  rb=rb,
  split_eps=split_eps,
  n_chain = n_chain,
  max_depth = max_depth,
  iter_pre_tree = iter_pre_tree,
  prob_split_by_x = prob_split_by_x,
  a = a,
  b = b,
  mu_mu = mu_mu
)


modelpar=c(sigmasq_mu = sigmasq_mu)


Gs_fn = function(y,lambda){
  tmp <- exp(lambda)/(1+exp(lambda))
  tmp[is.na(tmp)] <- 1
  return(y - tmp)
}
Hs_fn = function(y,lambda){
  tmp <- -exp(lambda)/(1+exp(lambda))^2
  tmp[is.na(tmp)] <- 0
  return(tmp)
}

source("/Users/shurenhe/Downloads/code/GSBart/Gen_GSBart_new.R")

GSBart_fit = Gen_GSBart(Gs_fn, Hs_fn, Graphs, Y_bin,
                        Y_bin_ho, modelpar, hyperpar,
                        40, 15, verbose = T, seed = 1234)

# GSBart_fit = GSBart::lgsbart(Graphs, Y_bin, modelpar, hyperpar,
#                         5, 3, verbose = T, seed = 1234)

save(GSBart_fit, file = 'Weather_GS_Bart_res.RData')

y_hat <- exp(GSBart_fit$lambda.hat.test)/(1+exp(GSBart_fit$lambda.hat.test))
y_hat <- as.numeric(y_hat > 0.5)

Y_bin_ho_tmp <- na.omit(Y_bin_ho)
remove.id = which(is.na(Y_bin_ho))
y_hat_tmp <- y_hat[-remove.id]
print(mean(y_hat_tmp == Y_bin_ho_tmp))

MCMC = 8000
BURNIN = 2000
THIN = 10
ndpost = (MCMC - BURNIN)/THIN
nskip = BURNIN
n_test = length(Y_bin_ho)
n_rounds = 100

BART_Fit = pbart(X, Y_bin, X_ho, ntree = n_rounds, ndpost = ndpost, nskip = nskip, keepevery = THIN)

tmp <- BART_Fit$prob.test.mean
BART_y_test = as.numeric(tmp>0.5)
BART_ACC = mean(BART_y_test[-remove.id] == Y_bin_ho_tmp)

xgb_grid = expand.grid(
  # nrounds = n_rounds,
  max_depth = 20,
  eta = 5 / seq(10, 100, length = 10),
  gamma = seq(0, 1, length = 10),
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)


library(parallel)
library(doParallel)
library(xgboost)
nworkers <- detectCores()
cl <- makeCluster(nworkers)
registerDoParallel(cl)

XGB_Time = Sys.time()

sim_list <- foreach(i = 1:nrow(xgb_grid), .packages= c('caret', 'xgboost')) %dopar% {
  cv <- xgb.cv(data = as.matrix(X), label = Y_bin, nrounds = 100, params = as.list(xgb_grid[i,]),
               nfold = 5, metrics = "error", objective='binary:logistic', verbose = F)
  RMSE <- mean(cv$evaluation_log$test_error_mean[20:50])
}
para_select = as.list(xgb_grid[which.min(unlist(sim_list)),])
XGB_Fit<- xgboost(data = as.matrix(X), label = Y_bin, nrounds = 50, params = para_select,
                  metrics = "error", objective='binary:logistic', verbose = F)
XGB_Time = Sys.time() - XGB_Time
b = predict(XGB_Fit, as.matrix(X_ho))
XGBoost_y_test <- as.numeric(b>0.5)
XGB_ACC = mean(XGBoost_y_test[-remove.id] == Y_bin_ho_tmp)

#
# XGBoost_CM_1 = confusionMatrix(factor(XGBoost_y_test_1, levels = (0:(K - 1))),
#                                factor(multinom_resp_test_1, levels = (0:(K - 1))))


meshout<-spmeshed(as.matrix(Y_bin, nrow = length(Y_bin)), x= X[,-c(1,2)], coords = X[,c(1,2)], family="binomial", axis_partition=c(4,4), n_samples = MCMC, n_burn=BURNIN, n_thin=THIN, prior=list(phi=c(1,15)), verbose= 0, n_threads = 1)

#posteriormeans best_post_mean<-meshout$beta_mcmc%>%apply(1:2,mean) #processmeans wmesh<-data.frame(w_mgp=meshout$w_mcmc%>%summary_list_mean()) #predictions ymesh<-data.frame(y_mgp=meshout$yhat_mcmc%>%summary_list_mean())

X_tmp_all <- rbind(X, X_ho)
in_train_X = 1:nrow(X)

X_tmp_all <- X_tmp_all %>% mutate_all(~(scale(.) %>% as.vector))
X_tmp_all <- X_tmp_all %>% select(-c(E_OTHERRACE))
X_tmp_new <- X_tmp_all[in_train_X, ]
X_ho_tmp_new <- X_tmp_all[-in_train_X,]

m_glm <- glm(Y_bin ~ as.matrix(X_tmp_new), family = binomial(link = "logit"))
confint(m_glm)
# tmp <- lm(Y ~ as.matrix(X_tmp_new))
# summary(tmp)
# b <- predict(tmp, X_ho_tmp_new)
# sqrt(mean((Y_ho - b)^2))

m_glm_residuals <- residuals(m_glm)
ggplot(X_tmp_new, aes(X_tmp_new$lon, X_tmp_new$lat, colour = m_glm_residuals)) +
  scale_color_gradient2() +
  geom_point(size = 3)

d <- cbind(Y_bin, X_tmp_new)
options(mc.cores=5)
m_spatial <- glmmfields(Y_bin ~ .,family = binomial(link = "logit"),
                        data = d, lat = "lat", lon = "lon", nknots = 20, iter = 2000, chains = 5,
                        prior_intercept = student_t(3, 0, 10),
                        prior_beta = student_t(3, 0, 11),
                        prior_sigma = half_t(3, 0, 10),
                        prior_gp_theta = half_t(3, 0, 5),
                        prior_gp_sigma = half_t(3, 0, 10),
                        seed = 1234 # passed to rstan::sampling()
)

# m_spatial <- glmmfields(Y ~ .,family = poisson(link = "log"),
#                         data = d, lat = "lat", lon = "lon", nknots = 20, iter = 2000, chains = 5,
#                         prior_intercept = student_t(3, 0, 10),
#                         prior_beta = student_t(3, 0, 10),
#                         prior_sigma = half_t(3, 0, 10),
#                         prior_gp_theta = half_t(3, 0, 5),
#                         prior_gp_sigma = half_t(3, 0, 10),
#                         seed = 1234 # passed to rstan::sampling()
#)
Y_pred <- predict(
  m_spatial,
  newdata = X_ho_tmp_new,
  type = "response"
)$estimate
Y_pred <- as.numeric(Y_pred>0.5)
remove.id = which(is.na(Y_bin_ho))

mean(Y_bin_ho[-remove.id] == Y_pred[-remove.id])
