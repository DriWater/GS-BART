source('R/Ushape_GenData_fun.R')

library(caret)

n_bin = 100
n_rounds = 50
nref = 100
num_ost = 5

## Generate U-shape dataset
hyperpar = c("n_rounds" = n_rounds,
             "nref" = nref,
             "num_ost" = num_ost,
             "n_bin"=n_bin)

n_train = 800; n_test = 200; p = 5
sim_Ushape = genUshape_fG(n_train, n_test, p = p, hyperpar, seed=1234)
sim_Ushape[['graphs_weight']] = c(rep((1-min(5/7, 0.85))/7, 7), rep(min(5/7, 0.85)/5, 5))
sim_Ushape$train_X = NULL; sim_Ushape$test_X = NULL

## Generate Friedman Example
set.seed (1234)
sigma <- 1
n <- 500
p = 5
f <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", 1:p)
Ey <- f(x)
in_train <- createDataPartition(Ey, p = .8, list = FALSE)
X = x[in_train,]; X_ho = x[-in_train,]
n_train = nrow(X); n_test = nrow(X_ho)

Graphs = lapply(1:n_rounds, function(x){Generate_GraphCandidates_graphbin(X = X, X_ho = X_ho, n_bin = 100, chain.type = 'regular')})
sim_Friedman = list(X = X, X_ho = X_ho, Graphs = Graphs, f_true = Ey[in_train], f_ho_true = Ey[-in_train])
sim_Friedman[['graphs_weight']] = rep(0.2, 5)


## Generate Torus Example
load("data/TorusSim.Rdata")
Graphs= lapply(InputGraphs, function(x) {
  lapply(x, function(y){
    list(root = y$root,
         children = y$children, 
         mesh2tr = y$mesh2tr,
         mesh2ho = y$mesh2ho,
         tr2mesh = y$tr2mesh)
  }) 
})
sim_Torus = sim
sim_Torus$Graphs = Graphs
sim_Torus[['graphs_weight']] = c(rep((1-min(6/8, 0.85))/7, 7), rep(min(6/8, 0.85)/6, 6))
sim_Torus$train_X = NULL; sim_Torus$test_X = NULL



save(sim_Ushape, sim_Friedman, sim_Torus, file = "data/sim_input.RData", compress = "bzip2")


