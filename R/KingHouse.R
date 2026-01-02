rm(list=ls())
library(ggplot2)
library(ggpubr)
library(GSBart)
load('~/Downloads/code/test_code/KingHouse.Rdata')
setwd("/Users/shurenhe/Downloads/code/GSBart_test/KingHouse/")


load('KingHouseBGGBoost.Rdata')
load('KingHouse.Rdata')
load('/Users/shurenhe/Downloads/code/GSBart_test/KingHouse/traindata.RData')

remove(Graphs)
remove(InputGraphs)
InputGraphs = suppressWarnings(
lapply(Graphs_tmp, function(x) {
  lapply(x, function(y) {
    return(list(root = y$root,
                children = y$children,
                mesh2tr = igraph::vertex.attributes(y$g)$mesh2tr,
                mesh2ho = igraph::vertex.attributes(y$g)$mesh2ho,
                tr2mesh = igraph::graph.attributes(y$g)$tr2mesh))
  }) 
})
)

InputGraphs = lapply(InputGraphs, function(x) {
  lapply(x, function(y) {
    list(root = y$root - 1,
         children = unname(lapply(y$children, function(z){unname(z) - 1})),
         mesh2tr = unname(lapply(y$mesh2tr, function(z){unname(z) - 1})),
         mesh2ho = unname(lapply(y$mesh2ho, function(z){unname(z) - 1})),
         tr2mesh = unname(y$tr2mesh - 1))
  }) 
})


Unstandardize = function(x, std_par){ return(x * std_par['scale'] + std_par['mean']) }
y.train.unstandardized <- Unstandardize(y.train, stdY$std_par)

GSBart_Time = Sys.time()
KingGSB <- gsbart(y.train.unstandardized, Graphs = InputGraphs, 200, 15, 200,  graphs_weight = graphs_weight, nthreads = 9, verbose = T, sparse = F)
GSBart_Time = difftime(Sys.time(), GSBart_Time, units = "secs")
KingGSB$time <- GSBart_Time

plot(KingGSB$sigma[150:350], type = 'l')

plot(KingGSB$tau[150:350], type = 'l')

geweke_test=function(sigma_out){
  geweke_results=(coda::geweke.diag(sigma_out))
  # Extract z-scores
  zvals <- geweke_results$z
  
  # Compute two-sided p-values
  pvals <- 2 * (1 - pnorm(abs(zvals)))
  return(pvals)
}

geweke_test(KingGSB$sigma[150:350])

# save(KingGSB, file = "/Users/shurenhe/Downloads/GSBart_rebuttal/KingHouse/KingGSBres350.RData")

load("/Users/shurenhe/Downloads/GSBart_rebuttal/KingHouse/KingGSBres350.RData")

set.seed(1234)
subInputGraphs = vector(mode = "list", length = length(InputGraphs))
Graph_len = length(InputGraphs[[1]])
for( i in 1:length(InputGraphs)){
  graphs_id = sample(1:Graph_len, 18, prob = graphs_weight)
  graphs_id = graphs_id[order(graphs_id)]
  subInputGraphs[[i]] = InputGraphs[[i]][graphs_id]
}

hyperpar['mu0'] = NULL
graphs_weight = rep(1/18, 18)

save(stdY, Y, Y_ho, hyperpar, modelpar, graphs_weight, subInputGraphs, y.test.unstandardized, file = "~/Downloads/code/test_code/KingHouse18sub.Rdata")



# Flex-BART
in_train = which(y.train %in% stdY$Y)
X = X_all[in_train,]
X_ho = X_all[-in_train, ]
coords = coords_all[in_train, ]
coords_ho = coords_all[-in_train, ]
X_tmp = cbind(coords, X)
X_ho_tmp = cbind(coords_ho,X_ho)
Y = y.train
Y_ho = y.test
y.test.unstandardized = Unstandardize(Y_ho, stdY$std_par)

coords_std <- Standardize(as.matrix(coords))


n_train = nrow(X_tmp); n_test = nrow(X_ho_tmp)
coords_all = rbind(coords, coords_ho)
df_coords <- data.frame(lat = coords_all[,1], long = coords_all[,2])
sf_coords <- st_as_sf(df_coords, coords = c(1:2))

nn.list = st_nn(sf_coords, sf_coords, k = 4)
g0 = graph.adjlist(nn.list) %>% as.undirected() %>% igraph::simplify()
E(g0)$weight = runif(ecount(g0), 0, 1)
components(g0)$membership
adj_matrix <- as_adjacency_matrix(g0, sparse = FALSE)

X_all <- rbind(X_tmp, X_ho_tmp)
vertex_id_train = c(1:n_train)
vertex_id_test = c((n_train+1):(n_train+n_test))
X_all_std <- apply(X_all, 2, function(x){(x - min(x)) / (max(x) - min(x))})
X_train_std = X_all_std[vertex_id_train,]
X_test_std = X_all_std[vertex_id_test,]


Flex_BART_Time=Sys.time()
Flex_BART_Fit = flexBART::network_BART(Y, vertex_id_train, X_train_std, vertex_id_test,
                                       X_test_std, A = adj_matrix, M = 200, nd = 400, burn = 2000, thin = 10)
Flex_BART_Time = difftime(Sys.time(), Flex_BART_Time, units = "secs")
Flex_BART_Y_out = t(apply(Flex_BART_Fit$yhat.test, 1, function(X){Unstandardize(X, stdY$std_par)}))
Flex_BART_pred_unstandardized = colMeans(Flex_BART_Y_out)
Flex_BART_MSPE = mean((Flex_BART_pred_unstandardized - y.test.unstandardized)^2)
Flex_BART_MAPE = mean(abs(Flex_BART_pred_unstandardized - y.test.unstandardized))
print(Flex_BART_MSPE, Flex_BART_MAPE, Flex_BART_Time)


library(stochtree)
nburnin = 15; ndpost = 200
XBART_Time=Sys.time()
XBART_Fit <- stochtree::bart(X_tmp, Y, X_test = X_ho_tmp, num_gfr = (nburnin+ndpost), num_burnin = 0, num_mcmc = 0, 
                             mean_forest_params = list(num_trees = 200), general_params = list(keep_gfr = T))
XBART_Time = difftime(Sys.time(), XBART_Time, units = "secs")
XBART_Fit$time = XBART_Time
XBART_Y_out = t(apply(XBART_Fit$y_hat_test, 2, function(X){Unstandardize(X, stdY$std_par)}))
XBART_pred_unstandardized = colMeans(XBART_Y_out[(nburnin+1):(ndpost+nburnin),])
XBART_MSPE = mean((XBART_pred_unstandardized - y.test.unstandardized)^2)
XBART_MAPE = mean(abs(XBART_pred_unstandardized - y.test.unstandardized))
print(c(XBART_MSPE, XBART_MAPE, XBART_Time))


# load("/Users/shurenhe/Downloads/GSBart_rebuttal/KingHouse/KingHouse6subGSB.Rdata")
# 
Unstandardize = function(x, std_par){ return(x * std_par['scale'] + std_par['mean']) }
# GSB.pred.unstandardized = Unstandardize(colMeans(KingHouseres$phi.test), stdY$std_par)
# GSBart_MSPE = mean((GSB.pred.unstandardized - y.test.unstandardized)^2)
# GSBart_MAPE = mean(abs(GSB.pred.unstandardized - y.test.unstandardized))


#load("~/Downloads/GSBart_rebuttal/KingHouse/KingHouseGSBres.Rdata")
# GSBart_Y_out = t(apply(KingHouseres$phi.test, 1, function(X){Unstandardize(X, stdY$std_par)}))
# GSBart_pred_unstandardized = colMeans(GSBart_Y_out)
GSB_trace <- colMeans((y.test.unstandardized - t(KingGSB$yhat.test))^2)
GS_iter = c(1:length(GSB_trace))

# GSB_trace <- colMeans((y.test.unstandardized - t(GSBart_Y_out))^2)
# GS_iter = c(1:length(GSB_trace))

df_mse <- data.frame(
  Iteration = c(GS_iter),
  MSPE = c(GSB_trace),
  Method = c(rep("bold(GS-BART)",length(GSB_trace)))
)

text_data <- data.frame(
  Method = c("bold(GS-BART)"),
  Iteration = c(30),
  RMSE = c(0.32),
  label = c("MSPE:0.094; ESS:23.77")
)

library(coda)

mean(coda::effectiveSize(KingGSB$yhat.test))
mean(pnorm(coda::geweke.diag(KingGSB$yhat.test)$z) > 0.05)

pdf(file = "/Users/shurenhe/Downloads/GSBart_rebuttal/KingHouse/KingHouse_traceplot.pdf", width  = 6, height = 4)

# 
# install.packages("ggbreak") 
# library(ggbreak)

p1 <- ggplot(df_mse, aes(x = Iteration, y = MSPE, color = Method)) +
  geom_line(col = 'black', alpha = 1) + # xlab('Tree Updates') +
  # geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
  #                              xintercept = c(14)),  aes(xintercept = xintercept),
  #            linetype = "dashed", color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  # geom_text(data = text_data,
  #           aes(x = 35, y = 0.3205, label = label),
  #           color = "black", size = 5,
  #           inherit.aes = FALSE) +
  geom_label(data = text_data,
             aes(x = 35, y = 0.1023, label = label),
             color = "black", size = 4,
             fill = "white", # 背景色
             label.size = 0.4, # 边框粗细，设为0可以去掉边框
             inherit.aes = FALSE)+
  # coord_cartesian(ylim = c(0.309, 0.321)) + 
  # scale_y_break(c(0.32, 0.5)) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15)  # Adjust size as needed
  ) +
  # ggtitle("KingHouse: GS-BART (ESS = 23.77/per 100 iter.)")
  ggtitle("KingHouse: GS-BART")
p1

# p1 <- ggplot(df_mse, aes(x = Iteration, y = RMSE, color = Method)) +
#   geom_line(alpha = 0.8) +
#   geom_vline(data = data.frame(Method = c("bold(GS-BART)"),
#                                xintercept = c(800)),  aes(xintercept = xintercept),
#              linetype = "dashed", color = "red", linewidth = 0.7) +
#   geom_text(data = text_data,
#             aes(x = Iteration, y = RMSE, label = label),
#             color = "blue", size = 6,
#             inherit.aes = FALSE) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
# p1
dev.off()

load("~/Downloads/GSBart_rebuttal/KingHouse/KingHouseGSBres.Rdata")
sum((colSums(KingHouseres$varcount)/sum(KingHouseres$varcount))[1:7])
