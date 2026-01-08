rm(list=ls())

library(sf)
library(igraph)
library(network)
library(caret)

source('R/graph_gen_fun.R')

## Section 3.2 

n_rounds = 50 # number of trees
nbin = 100 # number of bins for chain graphs
nref = 100  # number of bins for arborescences
nost = 5 # number of arborescences for each graph set

### Example 1 (only unstructured numerical features)

n = 500 # number of observations
p = 5 # number of covariates

x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", 1:p)

in_train <- createDataPartition(1:nrow(x), p = .8, list = FALSE)

X = x[in_train,] # training data 
X_ho = x[-in_train,] # testing data 

Graphs = lapply(1:n_rounds, function(x){Generate_GraphCandidates_graphbin(X = X, X_ho = X_ho, 
                                        n_bin = 100, chain.type = 'stratified')})

### Example 2 (network features)

load("data/NYEduData.Rdata")

nn.list= st_relate(st_geometry(sf_df2), pattern = "****1****")
g0 = graph.adjlist(nn.list) %>% as.undirected() %>% igraph::simplify()
E(g0)$weight = runif(ecount(g0), 0, 1)
Edu_dat = sf_df2[components(g0)$membership==1,] # save only connected components

sf_bnd=st_geometry(Edu_dat)
Edu_dat = Edu_dat%>%st_drop_geometry()
plot(sf_bnd) # plot mesh
g0 = graph.adjlist(st_relate(sf_bnd,pattern="****1****"))%>%as.undirected()%>%simplify() # build network based on mesh

X = Edu_dat[,c(2:6)]
in_train <- createDataPartition(1:nrow(X), p = .8, list = FALSE) # the index of testing and training data 

Graphs = Generate_GraphList_parallel(g0 = g0, X = X, in_train = in_train, nrounds = n_rounds, nost = nost, 
                                     nbin = nbin, nref = nref, nthreads = 5, chain.type = "stratified")

### Example 3 (spatial features on manifold)

# generate X from rectangle
n = 25
x<-seq(-0.99,0.99,length=n)
X=expand.grid(x=x,y=x)
colnames(X) <- c("x1", "x2")

# transfer coordinates to sfc form
sf_coords <- st_as_sf(as.data.frame(X), coords = c("x1", "x2"))
sf_coords=st_geometry(sf_coords)

# boundary is a 2D rectangle 
sf_bnd=sf::st_as_sfc(sf::st_bbox(c(xmin=-1, ymin=-1, xmax=1, ymax=1),
                                 crs=sf::st_crs(sf_coords)))
plot(sf_bnd)

in_train <- createDataPartition(1:nrow(X), p = .8, list = FALSE) # the index of testing and training data 
Graphs = Generate_GraphList_parallel(sf_coords = sf_coords, sf_bnd = sf_bnd, X = X, in_train = in_train, nrounds = n_rounds, 
                                     nost = nost, nbin = nbin, nref = nref, nthreads = 5, mesh.type = "stratified", chain.type = "stratified")

