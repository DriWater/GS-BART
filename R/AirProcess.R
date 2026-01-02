rm(list=ls())

set.seed(1234)
setwd("/Users/shurenhe/Downloads/GSBart_rebuttal/Airpollution/")
load("AsthmaCountyProcessed.Rdata")
load("USWeather.RData")

InputGraphs = lapply(Graphs, function(x) {
  lapply(x, function(y) {
    #   stopifnot(length(igraph::graph.attributes(y$g)$tr2mesh) == n_train,
    #             length(igraph::graph.attributes(y$g)$ho2mesh) == n_test,
    #             length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2tr))) == n_train,
    #             length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2ho))) == n_test)
    return(list(root = y$root,
                children = y$children,
                mesh2tr= igraph::vertex.attributes(y$g)$mesh2tr,
                mesh2ho = igraph::vertex.attributes(y$g)$mesh2ho,
                tr2mesh = igraph::graph.attributes(y$g)$tr2mesh))
  }) 
})
# Reduce every index by 1 to make everything 0 indexed
Graphs = lapply(InputGraphs, function(x) {
  lapply(x, function(y) {
    list(root = y$root - 1,
         children = unname(lapply(y$children, function(z){unname(z) - 1})),
         mesh2tr = unname(lapply(y$mesh2tr, function(z){unname(z) - 1})),
         mesh2ho = unname(lapply(y$mesh2ho, function(z){unname(z) - 1})),
         tr2mesh = unname(y$tr2mesh - 1))
  }) 
})

n_rounds = M = 50
n_ost=7
n_train = length(Y_bin)
n_test = length(Y_bin_ho)
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
  nu = nu, rb=rb,
  split_eps=split_eps,
  max_depth = 6,
  mu0 = qlogis(mean(Y_bin==1))/n_rounds,
  tree_iter = 15,
  a = 3,
  b = 0.5 * 6 / n_rounds
)

graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/n_ost,n_ost), rep(min(p / (p + 2), 0.85)/p,p))

modelpar=c(sigmasq = (3/(2*sqrt(n_rounds)))^2, tausq = (6/(2*sqrt(n_rounds)))^2)

save(Y_bin, Y_bin_ho, in_train, hyperpar, modelpar, graphs_weight,
     Graphs, file = "/Users/shurenhe/Documents/GitHub/GS-BART/data/USAir.Rdata", compress = "bzip2")

load("~/Downloads/GSBart_rebuttal/Airpollution/USAir.Rdata")

# y_hat <- exp(colMeans(USAir$phi.test))/(1+exp(colMeans(USAir$phi.test)))
y_hat <- as.numeric(USAir$prob.test.mean > 0.5)

Y_bin_ho_tmp <- na.omit(Y_bin_ho)
remove.id = which(is.na(Y_bin_ho))
y_hat_tmp <- y_hat[-remove.id]
print(mean(y_hat_tmp == Y_bin_ho_tmp))

library(G2SBart)
load('~/Downloads/code/test_code/AirPollution2.Rdata')
hyperpar['tree_iter'] = 30; hyperpar['max_depth'] = 15
n_rounds = length(InputGraphs)
modelpar['tausq'] = (3/(2*sqrt(length(InputGraphs))))^2
hyperpar['b'] = 0.5 * 6 / n_rounds
hyperpar <- as.list(hyperpar); modelpar <- as.list(modelpar)
USAir <- g2sbart(Y_bin, InputGraphs, 200, 15, 200, family = "binary", graphs_weight,
                hyperpar, modelpar, nthreads = 8, verbose = T, seed = 1234)
