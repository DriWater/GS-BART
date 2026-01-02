source('/Users/shurenhe/Documents/GitHub/GS-BART/R/Ushape_GenData_fun.R')
n_bin = 100
n_rounds = 50
nref = 100
num_ost = 5

## Generate U-shape dataset
Standardize = function(Y, std_par = NULL) {
  if(is.null(std_par)) {
    ymean = mean(Y)
    yscale = 2 * max(abs(Y))
  } else {
    ymean = std_par['mean']
    yscale = std_par['scale']
  }
  Y = (Y - ymean) / yscale
  std_par = c('mean' = ymean, 'scale' = yscale)
  return(list(Y = Y, std_par = std_par))
}

Unstandardize = function(x, std_par) {
  return(x * std_par['scale'] + std_par['mean']) 
} 

hyperpar = c("n_rounds" = n_rounds,
             "nref" = nref,
             "num_ost" = num_ost,
             "n_bin"=n_bin)
n_train = 800; n_test = 200; p = 5
sim_Ushape = genUshape_fG(n_train, n_test, p = p, hyperpar, seed=1234)

save(coords, coords_ho, X, X_ho, Y, Y_ho, f, f_ho, cluster_true, cluster_ho_true, n, n_ho,
     phi_uni, sigmasq_uni, beta0_uni,
     file = "data/sim_input.RData", compress = "bzip2")

## Generate Friedman Example
p = 5
f <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

in_train <- createDataPartition(Ey, p = .8, list = FALSE)
X = x[in_train,]; X_ho = x[-in_train,]
n_train = nrow(X); n_test = nrow(X_ho)

Graphs = lapply(1:n_rounds, function(x){Generate_GraphCandidates_graphbin(X = X, X_ho = X_ho, n_bin = 100, chain.type = 'regular')})


## Generate Torus Example

