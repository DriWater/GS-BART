rm(list=ls())

set.seed(1234)
setwd("/Users/shurenhe/Downloads/GSBart_rebuttal/USFlood/")
load("AsthmaCountyProcessed.Rdata")
load("AsthmaCounty.Rdata")
load('/Users/shurenhe/Downloads/code/GSBart_test/UScounty/USFlood.RData')

InputGraphs = lapply(Graphs, function(x) {
  lapply(x, function(y) {
    #   stopifnot(length(igraph::graph.attributes(y$g)$tr2mesh) == n_train,
    #             length(igraph::graph.attributes(y$g)$ho2mesh) == n_test,
    #             length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2tr))) == n_train,
    #             length(unique(do.call(c, igraph::vertex.attributes(y$g)$mesh2ho))) == n_test)
    return(list(root = y$root,
                children = y$children,
                mesh2trList = igraph::vertex.attributes(y$g)$mesh2tr,
                mesh2hoList = igraph::vertex.attributes(y$g)$mesh2ho,
                tr2mesh = igraph::graph.attributes(y$g)$tr2mesh,
                ho2mesh = igraph::graph.attributes(y$g)$ho2mesh))
  }) 
})
# Reduce every index by 1 to make everything 0 indexed
InputGraphs = lapply(InputGraphs, function(x) {
  lapply(x, function(y) {
    list(root = y$root - 1,
         children = unname(lapply(y$children, function(z){unname(z) - 1})),
         mesh2trList = unname(lapply(y$mesh2trList, function(z){unname(z) - 1})),
         mesh2hoList = unname(lapply(y$mesh2hoList, function(z){unname(z) - 1})),
         tr2mesh = unname(y$tr2mesh - 1),
         ho2mesh = unname(y$ho2mesh - 1))
  }) 
})

n_rounds = M = 100 
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
  nu = nu, rb=rb,
  split_eps=split_eps,
  max_depth = 6,
  mu0 = qlogis(mean(Y==1))/n_rounds,
  tree_iter = 15,
  a = 3,
  b = 0.5 * 3 / n_rounds
)

graphs_weight = c(rep((1-min(p / (p + 2), 0.85))/n_ost,n_ost), rep(min(p / (p + 2), 0.85)/p,p))

modelpar=c(sigmasq = (3/(2*sqrt(n_rounds)))^2, tausq = (3/(2*sqrt(n_rounds)))^2)

save(Y, Y_ho, hyperpar, modelpar, graphs_weight, InputGraphs, file = "~/Downloads/code/test_code/USFlood.Rdata")

