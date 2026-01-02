### Simulation data on U-shape domain ###
library(ComplexDomain)
library(ggnet)
library(igraph)
library(dplyr)
library(network)
library(parallel)
library(doRNG)
library(sf)
source('/Users/shurenhe/Documents/GitHub/GS-BART/R/graph_gen_fun.R')

genVMesh <- function(sf_bnd, coords_ref = NULL, n_ref = 100, graph = TRUE, max_retry = 10, ...) {
  
  attempt <- 1
  sf_mesh <- NULL
  seeds=sample(1:10000,max_retry)
  while (attempt <= length(seeds)) {
    tryCatch({
      # If no coords_ref, sample points inside boundary
      if (is.null(coords_ref)) {
        set.seed(seeds[attempt])
        coords_ref <- st_sample(sf_bnd, size = n_ref, ...)
      }
      
      # Voronoi + intersection
      voronoi <- st_voronoi(st_union(coords_ref), sf_bnd)
      sf_mesh <- st_intersection(st_cast(voronoi), sf_bnd)
      
      # If we made it here, success → exit loop
      break
    },
    error = function(e) {
      #message(sprintf("Attempt %d failed: %s", attempt, e$message))
      coords_ref <<- NULL  # reset for next attempt
      attempt <<- attempt + 1
    })
  }
  
  if (is.null(sf_mesh)) {
    stop("Mesh generation failed after ", max_retry, " attempts.")
  }
  
  if (!graph) {
    return(sf_mesh)
  } else {
    neighboring <- st_relate(sf_mesh, pattern = "****1****")
    g <- graph.adjlist(neighboring) %>%
      as.undirected() %>%
      igraph::simplify()
    return(list(mesh = sf_mesh, g = g))
  }
}

genUshape_fG_vis = function(n, n_ho, p,hyperpar,seed=1234) {
  # rate is the expansion rate
  #Ushape=gensfUbnd(rate=1,rot_angle=45)
  set.seed(seed)
  n_bin = hyperpar['n_bin']
  n_rounds=hyperpar['n_rounds']
  nref=hyperpar['nref']
  num_ost=hyperpar['num_ost']
  
  rot_angle = 45
  ## Generate a ushape domain.
  # if cluster.value is provided, generate 3 clusters on u-shape
  Ushape = gensfUbnd(
    rate = 1,
    rot_angle = rot_angle,
    cluster.value = c(-1, 0.5, 1)
  )
  sf_bnd = Ushape$bnd
  sf_cluster = Ushape$cluster
  bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2])
  colnames(bnd) = c('x', 'y')
  ## Generate n uniform training locations
  message("Sampling locations")
  sf_coords = st_sample(sf_bnd, n)
  coords = st_coordinates(sf_coords)
  colnames(coords) = c('lon', 'lat')
  
  ## generate n_ho uniform hold-out locations
  mesh_ho = genVMesh(sf_bnd,
                     n_ref = n_ho,
                     graph = TRUE,
                     type = 'regular')
  
  sf_coords_ho= st_centroid(mesh_ho$mesh)
  coords_ho=st_coordinates(sf_coords_ho)
  colnames(coords_ho) = c('lon', 'lat')
  cluster_true = unlist(st_within(sf_coords, sf_cluster))
  cluster_ho_true = unlist(st_within(sf_coords_ho, sf_cluster))
  coords_all = rbind(coords, coords_ho)
  
  # generate X from GP
  message("Simulating X")
  X_all = simGpU(p, coords_all, angle = rot_angle)
  colnames(X_all) = paste("X", 1:p, sep = "")
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  # generate true function
  message("Evaluating f at X")
  f_true = evalFunU(coords, X, angle = rot_angle)
  f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)
  
  # modify true function to create discontinuity
  f_true[cluster_true == 1] = f_true[cluster_true == 1] - 4
  f_true[cluster_true == 2] =  f_true[cluster_true == 2] + 4
  f_true[cluster_true == 3] = -0.5 * f_true[cluster_true == 3]
  
  f_ho_true[cluster_ho_true == 1] = f_ho_true[cluster_ho_true == 1] - 4
  f_ho_true[cluster_ho_true == 2] =  f_ho_true[cluster_ho_true == 2] + 4
  f_ho_true[cluster_ho_true == 3] = -0.5 * f_ho_true[cluster_ho_true == 3]
  
  train_X = cbind(coords, X)
  test_X = cbind(coords_ho, X_ho)
  
  num.ost = num_ost
  nrounds = n_rounds
  Graphs = list()
  M = nrounds
  message("Generating Graphs")
  Graphs = foreach::foreach(m = 1:M,
                            # .packages = c("sf", "igraph", "ComplexDomain"),
                            .packages = c("sf", "igraph"),
                            .export = c("Generate_GraphCandidates", "Generate_OST",
                                        "Generate_Chain", "x2sfpoint", "GetPredecessors")) %dorng% {
                                          # mesh = ComplexDomain::genVMesh(sf_bnd,
                                          #                                n_ref = nref,
                                          #                                graph = TRUE,
                                          #                                type = 'random')
                                          
                                          mesh = genVMesh(sf_bnd,
                                                          n_ref = nref,
                                                          graph = TRUE,
                                                          type = 'random')
                                          
                                          g0 = mesh$g
                                          sf_mesh = mesh$mesh
                                          n_ho_tmp = length(sf_coords_ho)
                                          
                                          igraph::graph_attr(g0, 'tr2mesh') = unlist(st_within(sf_coords, sf_mesh))
                                          igraph::graph_attr(g0, 'ho2mesh') = unlist(st_within(sf_coords_ho, sf_mesh))
                                          
                                          
                                          vertex_attr(g0)=list(mesh2tr=split(1:n,factor(graph_attr(g0,'tr2mesh'),1:vcount(g0))),
                                                               mesh2ho=split(1:n_ho_tmp,factor(graph_attr(g0,'ho2mesh'),1:vcount(g0))))
                                          Generate_GraphCandidates(g0,
                                                                   X = train_X,
                                                                   X_ho = test_X,
                                                                   num.ost = num.ost,
                                                                   n_bin = n_bin, chain.type = 'random')
                                        }
  
  return(list(
    X = X,
    X_ho = X_ho,
    train_X = train_X,
    test_X = test_X,
    coords = coords,
    coords_ho = coords_ho,
    mesh_ho = mesh_ho,
    Graphs = Graphs,
    f_true = f_true,
    f_ho_true = f_ho_true,
    sf_bnd = sf_bnd,
    sf_cluster = sf_cluster
  ))
  
}

### Generate full data on U shape
genUshape_fG = function(n, n_ho, p, hyperpar, seed=1234, nthreads = 1) {
  # rate is the expansion rate
  #Ushape=gensfUbnd(rate=1,rot_angle=45)
  set.seed(seed)
  n_bin = hyperpar['n_bin']
  n_rounds=hyperpar['n_rounds']
  nref=hyperpar['nref']
  num_ost=hyperpar['num_ost']
  
  rot_angle = 45
  ## Generate a ushape domain.
  # if cluster.value is provided, generate 3 clusters on u-shape
  Ushape = gensfUbnd(
    rate = 1,
    rot_angle = rot_angle,
    cluster.value = c(-1, 0.5, 1)
  )
  sf_bnd = Ushape$bnd
  sf_cluster = Ushape$cluster
  bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2])
  colnames(bnd) = c('x', 'y')
  ## Generate n uniform training locations
  message("Sampling locations")
  sf_coords = st_sample(sf_bnd, n)
  coords = st_coordinates(sf_coords)
  colnames(coords) = c('lon', 'lat')
  ## generate n_ho uniform hold-out locations
  sf_coords_ho = st_sample(sf_bnd, n_ho)
  coords_ho = st_coordinates(sf_coords_ho)
  colnames(coords_ho) = c('lon', 'lat')
  cluster_true = unlist(st_within(sf_coords, sf_cluster))
  cluster_ho_true = unlist(st_within(sf_coords_ho, sf_cluster))
  
  sf_coords_all = c(sf_coords, sf_coords_ho)
  coords_all = rbind(coords, coords_ho)
  
  # generate X from GP
  message("Simulating X")
  X_all = simGpU(p, coords_all, angle = rot_angle)
  colnames(X_all) = paste("X", 1:p, sep = "")
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  
  # generate true function
  message("Evaluating f at X")
  f_true = evalFunU(coords, X, angle = rot_angle)
  f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)
  
  # modify true function to create discontinuity
  f_true[cluster_true == 1] = f_true[cluster_true == 1] - 4
  f_true[cluster_true == 2] =  f_true[cluster_true == 2] + 4
  f_true[cluster_true == 3] = -0.5 * f_true[cluster_true == 3]
  
  f_ho_true[cluster_ho_true == 1] = f_ho_true[cluster_ho_true == 1] - 4
  f_ho_true[cluster_ho_true == 2] =  f_ho_true[cluster_ho_true == 2] + 4
  f_ho_true[cluster_ho_true == 3] = -0.5 * f_ho_true[cluster_ho_true == 3]
  
  train_X = as.matrix(cbind(coords, X))
  test_X = as.matrix(cbind(coords_ho, X_ho))
  
  num.ost = num_ost
  nrounds = n_rounds
  Graphs = list()
  M = nrounds
  message("Generating Graphs")
  
  if(nthreads > 1){
    Graphs = Generate_GraphList_parallel(sf_coords = sf_coords_all, sf_bnd = sf_bnd, g0 = NULL, 
                                         rbind(train_X, test_X), in_train = c(1:n), nrounds, num_ost, 
                                         n_bin, nref, nthreads = 5, mesh.type = "stratified", chain.type = "stratified")
  }else{
    Graphs = Generate_GraphList(sf_coords = sf_coords_all, sf_bnd = sf_bnd, g0 = NULL, 
                                rbind(train_X, test_X), in_train = c(1:n), nrounds, num_ost, 
                                n_bin, nref, mesh.type = "stratified", chain.type = "stratified")
  }
  
  return(list(
    X = X,
    X_ho = X_ho,
    coords = coords,
    coords_ho = coords_ho,
    Graphs = Graphs,
    f_true = f_true,
    f_ho_true = f_ho_true,
    sf_bnd = sf_bnd
  ))
  
}
####
genUshape_YG = function(n, n_ho, p, hyperpar,
                        noise = 0.1, standardize = TRUE,testloc='random',seed=1234) {
  # rate is the expansion rate
  #Ushape=gensfUbnd(rate=1,rot_angle=45)
  set.seed(seed)
  n_bin = hyperpar['n_bin']
  n_rounds=hyperpar['n_rounds']
  nref=hyperpar['nref']
  num_ost=hyperpar['num_ost']
  rot_angle = 45
  ## Generate a ushape domain.
  # if cluster.value is provided, generate 3 clusters on u-shape
  Ushape = gensfUbnd(
    rate = 1,
    rot_angle = rot_angle,
    cluster.value = c(-1, 0.5, 1)
  )
  sf_bnd = Ushape$bnd
  sf_cluster = Ushape$cluster
  #ggplot(Ushape$cluster) + geom_sf(aes(fill = as.factor(beta)))
  bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2])
  colnames(bnd) = c('x', 'y')
  ## Generate n uniform training locations
  message("Sampling locations")
  sf_coords = st_sample(sf_bnd, n)
  coords = st_coordinates(sf_coords)
  colnames(coords) = c('lon', 'lat')
  ## generate n_ho uniform hold-out locations
  
  ## generate n_ho uniform hold-out locations
  if(testloc=='grid'){
    mesh_grid=genVMesh(sf_bnd,n_ref=n_ho,type='regular')$mesh
    sf_coords_ho = st_centroid(mesh_grid)
  }else{
    sf_coords_ho = st_sample(sf_bnd, n_ho)
  }
  coords_ho = st_coordinates(sf_coords_ho)
  colnames(coords_ho) = c('lon', 'lat')
  n_ho=nrow(coords_ho)
  cluster_true = unlist(st_within(sf_coords, sf_cluster))
  cluster_ho_true = unlist(st_within(sf_coords_ho, sf_cluster))
  
  coords_all = rbind(coords, coords_ho)
  message("Simulating X")
  X_all = simGpU(p, coords_all, angle = rot_angle)
  colnames(X_all) = paste("X", 1:p, sep = "")
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  # generate true function
  message("Evaluating f at X")
  f_true = evalFunU(coords, X, angle = rot_angle)
  f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)
  
  # modify true function to create discontinuity
  f_true[cluster_true == 1] = f_true[cluster_true == 1] - 4
  f_true[cluster_true == 2] =  f_true[cluster_true == 2] + 4
  f_true[cluster_true == 3] = -0.5 * f_true[cluster_true == 3]
  
  f_ho_true[cluster_ho_true == 1] = f_ho_true[cluster_ho_true == 1] - 4
  f_ho_true[cluster_ho_true == 2] =  f_ho_true[cluster_ho_true == 2] + 4
  f_ho_true[cluster_ho_true == 3] = -0.5 * f_ho_true[cluster_ho_true == 3]
  
  message("Sampling Y")
  sigmasq_y = noise
  Y0 = f_true + rnorm(n, 0, sigmasq_y)
  Y0_ho = f_ho_true + rnorm(n_ho, 0, sigmasq_y)
  
  ## standardize data
  if(standardize) {
    stdY = Standardize(Y0)
    Y = stdY$Y
    
    Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  } else {
    Y = Y0
    Y_ho = Y0_ho
  }
  
  
  num.ost = num_ost
  nrounds = n_rounds
  Graphs = list()
  M = nrounds
  message("Generating Graphs")
  Graphs = foreach::foreach(m = 1:M,
                            # .packages = c("sf", "igraph", "ComplexDomain"),
                            .packages = c("sf", "igraph"),
                            .export = c("Generate_GraphCandidates", "Generate_OST",
                                        "Generate_Chain", "x2sfpoint", "GetPredecessors")) %dorng% {
                                          mesh = genVMesh(sf_bnd,
                                                          n_ref = nref,
                                                          graph = TRUE,
                                                          type = 'random')
                                          g0 = mesh$g
                                          sf_mesh = mesh$mesh
                                          
                                          
                                          igraph::graph_attr(g0, 'tr2mesh') = unlist(st_within(sf_coords, sf_mesh))
                                          igraph::graph_attr(g0, 'ho2mesh') = unlist(st_within(sf_coords_ho, sf_mesh))
                                          
                                          
                                          igraph::vertex_attr(g0)=list(mesh2tr=split(1:n,factor(graph_attr(g0,'tr2mesh'),1:vcount(g0))),
                                                                       mesh2ho=split(1:n_ho,factor(graph_attr(g0,'ho2mesh'),1:vcount(g0))))
                                          
                                          Generate_GraphCandidates(g0,
                                                                   X = X,
                                                                   X_ho = X_ho,
                                                                   num.ost = num.ost,
                                                                   n_bin = n_bin)
                                        }
  
  train_X = as.matrix(cbind(coords, X))
  test_X = as.matrix(cbind(coords_ho, X_ho))
  
  return(list(
    Y = Y,
    Y_ho = Y_ho,
    X = X,
    X_ho = X_ho,
    coords = coords,
    coords_ho = coords_ho,
    train_X = train_X,
    test_X = test_X,
    Graphs = Graphs,
    f_true = f_true,
    f_ho_true = f_ho_true,
    sf_bnd = sf_bnd,
    sf_cluster = sf_cluster,
    stdY=stdY
  ))
  
}



##### Generate U-shape data for classification data
#### sf functions and sf point/polygon types are used extensively
genUshapeExample_C_YG = function(n, n_ho, p, n_bin = 20, n_rounds, nref, num_ost,
                                 noise = 0.1, standardize = TRUE,testloc='random',seed=1234) {
  # rate is the expansion rate
  #Ushape=gensfUbnd(rate=1,rot_angle=45)
  set.seed(seed)
  rot_angle = 45
  ## Generate a ushape domain.
  # if cluster.value is provided, generate 3 clusters on u-shape
  ubnd = genURegion(angle = rot_angle)
  outer = as.matrix(cbind(ubnd$x, ubnd$y))
  sfc_ubnd = st_polygon(list(as.matrix(rbind(outer, outer[1,
  ])))) %>% st_geometry()
  
  centroid = st_point(c(0, 0))
  sfc_line=st_linestring(matrix(c(10,-9.6,-10,10.4),nrow=2,byrow=TRUE))
  sfc_circle = centroid %>% st_buffer( 0.91) %>%
    st_geometry()
  sfc_hcircle=st_collection_extract(st_split(sfc_circle,sfc_line))[2]
  sfc_clust12 = st_cast(st_difference(sfc_ubnd, sfc_hcircle),
                        "POLYGON")
  sfc_clust3 = st_intersection(sfc_ubnd, sfc_hcircle)
  beta_true_uniq = c(1,2,3)
  sf_uclust = st_sf(beta = beta_true_uniq, geometry = c(sfc_clust12,
                                                        sfc_clust3))
  Ushape = list(bnd=sfc_ubnd,cluster=sf_uclust)
  sf_bnd = Ushape$bnd
  sf_cluster = Ushape$cluster
  #ggplot(Ushape$cluster) + geom_sf(aes(fill = as.factor(beta)))
  bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2])
  colnames(bnd) = c('x', 'y')
  ## Generate n uniform training locations
  message("Sampling locations")
  sf_coords = st_sample(sf_bnd, n)
  coords = st_coordinates(sf_coords)
  colnames(coords) = c('lon', 'lat')
  ## generate n_ho uniform hold-out locations
  if(testloc=='grid'){
    mesh_grid=genVMesh(sf_bnd,n_ref=n_ho,type='regular')$mesh
    sf_coords_ho = st_centroid(mesh_grid)
  }else{
    sf_coords_ho = st_sample(sf_bnd, n_ho)
  }
  coords_ho = st_coordinates(sf_coords_ho)
  colnames(coords_ho) = c('lon', 'lat')
  n_ho=nrow(coords_ho)
  
  cluster_true = unlist(st_within(sf_coords, sf_cluster))
  cluster_ho_true = unlist(st_within(sf_coords_ho, sf_cluster))
  
  ## Generate piecewise constant functions
  
  # estimate geodesic distance between observed/testing locations
  # NOTE: this may take a while
  coords_all = rbind(coords, coords_ho)
  #dist_res = gdist(coords_all, bnd)# this can potentially be solved using Nan Wu and Dunson's JRSSB paper
  # #dist_mat_all = dist_res[1:(n + n_ho), 1:(n + n_ho)]
  # #rm(dist_res)
  #
  # generate X from GP
  message("Simulating X")
  X_all = simGpU(p, coords_all, angle = rot_angle)
  colnames(X_all) = paste("X", 1:p, sep = "")
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  # check what X[, 1] looks like
  # ggplot() +
  #   geom_boundary(bnd) +
  #   geom_point(aes(x = lon, y = lat, color = X[, 1]), data = as.data.frame(coords)) +
  #   scale_color_gradientn(colours = rainbow(5))  +
  #   labs(x = 's_1', y = 's_2', title = 'X[, 1]')
  
  # generate true function
  message("Evaluating f at X")
  f_true = evalFunU(coords, X, angle = rot_angle)
  f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)
  
  # modify true function to create discontinuity
  f_true[cluster_true == 1] = f_true[cluster_true == 1] - 4
  f_true[cluster_true == 2] =  f_true[cluster_true == 2] + 4
  f_true[cluster_true == 3] = -7#*f_true[cluster_true == 3]
  
  f_ho_true[cluster_ho_true == 1] = f_ho_true[cluster_ho_true == 1] - 4
  f_ho_true[cluster_ho_true == 2] =  f_ho_true[cluster_ho_true == 2] + 4
  f_ho_true[cluster_ho_true == 3] = -7# * f_ho_true[cluster_ho_true == 3]
  
  message("Sampling Y")
  sigmasq_y = noise
  Y0 = f_true + rnorm(n, 0, sigmasq_y)
  Y0_ho = f_ho_true + rnorm(n_ho, 0, sigmasq_y)
  
  ## standardize data
  if(standardize) {
    stdY = Standardize(Y0)
    Y = stdY$Y
    
    Y_ho = (Y0_ho - stdY$std_par['mean']) / stdY$std_par['scale']
  } else {
    Y = Y0
    Y_ho = Y0_ho
  }
  
  
  ggplot() +
    geom_boundary(bnd) +
    geom_point(aes(x = lon, y = lat, color = f_true), data = as.data.frame(coords)) +
    scale_color_gradientn(colours = rainbow(5))  +
    labs(x = 's_1', y = 's_2', title = 'f_true')
  
  ############
  
  
  num.ost = num_ost
  nrounds = n_rounds
  Graphs = list()
  M = nrounds
  message("Generating Graphs")
  Graphs = foreach::foreach(m = 1:M,
                            # .packages = c("sf", "igraph", "ComplexDomain"),
                            .packages = c("sf", "igraph"),
                            .export = c("Generate_GraphCandidates", "Generate_OST",
                                        "Generate_Chain", "x2sfpoint", "GetPredecessors")) %dorng% {
                                          # mesh = ComplexDomain::genVMesh(sf_bnd,
                                          #                                n_ref = nref,
                                          #                                graph = TRUE,
                                          #                                type = 'random')
                                          mesh = genVMesh(sf_bnd,
                                                          n_ref = nref,
                                                          graph = TRUE,
                                                          type = 'random')
                                          g0 = mesh$g
                                          sf_mesh = mesh$mesh
                                          
                                          
                                          igraph::graph_attr(g0, 'tr2mesh') = unlist(st_within(sf_coords, sf_mesh))
                                          igraph::graph_attr(g0, 'ho2mesh') = unlist(st_within(sf_coords_ho, sf_mesh))
                                          
                                          
                                          igraph::vertex_attr(g0)=list(mesh2tr=split(1:n,factor(graph_attr(g0,'tr2mesh'),1:vcount(g0))),
                                                                       mesh2ho=split(1:n_ho,factor(graph_attr(g0,'ho2mesh'),1:vcount(g0))))
                                          Generate_GraphCandidates(g0,
                                                                   X = X,
                                                                   X_ho = X_ho,
                                                                   num.ost = num.ost,
                                                                   n_bin = n_bin)
                                        }
  
  train_X = as.matrix(cbind(coords, X))
  test_X = as.matrix(cbind(coords_ho, X_ho))
  
  return(list(
    Y = Y,
    Y_ho = Y_ho,
    X = X,
    X_ho = X_ho,
    coords = coords,
    coords_ho = coords_ho,
    train_X = train_X,
    test_X = test_X,
    Graphs = Graphs,
    f_true = f_true,
    f_ho_true = f_ho_true,
    sf_bnd = sf_bnd
  ))
  
}

#######
genUData_fG_BAMDT=function(n, n_ho, p, M, n_ref,seed=1234){
  set.seed(seed)
  rot_angle = 45
  ## Generate a ushape domain.
  # if cluster.value is provided, generate 3 clusters on u-shape
  Ushape = gensfUbnd(
    rate = 1,
    rot_angle = rot_angle,
    cluster.value = c(-1, 0.5, 1)
  )
  sf_bnd = Ushape$bnd
  sf_cluster = Ushape$cluster
  #ggplot(Ushape$cluster) + geom_sf(aes(fill = as.factor(beta)))
  ubnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2])
  colnames(ubnd) = c('x', 'y')
  ## Generate n uniform training locations
  message("Sampling locations")
  sf_coords = st_sample(sf_bnd, n)
  coords = st_coordinates(sf_coords)
  colnames(coords) = c('lon', 'lat')
  ## generate n_ho uniform hold-out locations
  sf_coords_ho = st_sample(sf_bnd, n_ho)
  coords_ho = st_coordinates(sf_coords_ho)
  colnames(coords_ho) = c('lon', 'lat')
  cluster_true = unlist(st_within(sf_coords, sf_cluster))
  cluster_ho_true = unlist(st_within(sf_coords_ho, sf_cluster))
  
  ## Generate piecewise constant functions
  
  # estimate geodesic distance between observed/testing locations
  # NOTE: this may take a while
  coords_all = rbind(coords, coords_ho)
  #dist_res = gdist(coords_all, bnd)# this can potentially be solved using Nan Wu and Dunson's JRSSB paper
  # #dist_mat_all = dist_res[1:(n + n_ho), 1:(n + n_ho)]
  # #rm(dist_res)
  #
  # generate X from GP
  message("Simulating X")
  X_all = simGpU(p, coords_all, angle = rot_angle)
  colnames(X_all) = paste("X", 1:p, sep = "")
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  # check what X[, 1] looks like
  # ggplot() +
  #   geom_boundary(bnd) +
  #   geom_point(aes(x = lon, y = lat, color = X[, 1]), data = as.data.frame(coords)) +
  #   scale_color_gradientn(colours = rainbow(5))  +
  #   labs(x = 's_1', y = 's_2', title = 'X[, 1]')
  
  # generate true function
  message("Evaluating f at X")
  f_true = evalFunU(coords, X, angle = rot_angle)
  f_ho_true = evalFunU(coords_ho, X_ho, angle = rot_angle)
  
  # modify true function to create discontinuity
  f_true[cluster_true == 1] = f_true[cluster_true == 1] - 4
  f_true[cluster_true == 2] =  f_true[cluster_true == 2] + 4
  f_true[cluster_true == 3] = -0.5 * f_true[cluster_true == 3]
  
  f_ho_true[cluster_ho_true == 1] = f_ho_true[cluster_ho_true == 1] - 4
  f_ho_true[cluster_ho_true == 2] =  f_ho_true[cluster_ho_true == 2] + 4
  f_ho_true[cluster_ho_true == 3] = -0.5 * f_ho_true[cluster_ho_true == 3]
  
  
  coords_all = rbind(coords, coords_ho)
  dist_res = gdist(coords_all, ubnd)
  dist_mat_all = dist_res[1:(n + n_ho), 1:(n + n_ho)]
  rm(dist_res)
  
  
  coords_knots = list() # knot coordinates for each tree
  graphs = list()       # spatial graph for each tree
  knot_idx = list()
  for (m in 1:M) {
    # subsample training locations as knots
    knot_idx[[m]] = sample.int(n, n_ref)
    coords_knots[[m]] = coords[knot_idx[[m]], ]
    dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
    
    # get CDT graph on knots
    mesh = gen2dMesh(coords_knots[[m]], ubnd)
    graph0 = constrainedDentri(n_ref, mesh,
                               gaurantee_connected = T, dist_mat = dist_mat_knot)
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    E(graph0)$weight=runif(ecount(graph0))
    graphs[[m]] = graph0
  }
  
  # ensuring spatial graphs are connected
  for (m in 1:M) {
    if (components(graphs[[m]])$no != 1)
      stop(paste("Disconnected graph:", m))
  }
  
  # assign observations to their nearest knots
  # projections[i, j] is the index of the nearest knot of obs. i in weak learner j
  # similirly, projections_ho is for hold-out locations
  projections = array(0, dim = c(n, M))
  projections_ho = array(0, dim = c(n_ho, M))
  for (m in 1:M) {
    # get distance between training locations and knots
    cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
    projections[, m] = apply(cdist_mat_knots, 1, which.min)
    
    # get distance between hold-out locations and knots
    cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
    projections_ho[, m] = apply(cdist_mat_ho_knots, 1, which.min)
  }
  
  
  return(list(X=X, coords=X, n=n, p=p, ubnd=ubnd,
              X_ho=X_ho, coords_ho=coords_ho, n_ho=n_ho, dist_mat_all=dist_mat_all,
              f_true=f_true,
              f_ho_true=f_ho_true,
              projections=projections,
              graphs=graphs,
              projections_ho=projections_ho))
}

GenUshape_G_BAMDT=function(coords,coords_ho,hyperpar,dist_mat_all,coords_ref=NULL){
  
  n=nrow(coords);n_ho=nrow(coords_ho)
  #rm(dist_res)
  M=hyperpar['n_rounds']
  if(is.null(coords_ref)){
    coords_knots = list()
  }else{
    coords_knots=coords_ref
  }
  # knot coordinates for each tree
  graphs = list()       # spatial graph for each tree
  knot_idx = list()
  for (m in 1:M) {
    # subsample training locations as knots
    if(is.null(coords_ref)){
      knot_idx[[m]] = sample.int(n, hyperpar['nref'])
      coords_knots[[m]] = coords[knot_idx[[m]], ]
    }
    
    dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
    
    # get CDT graph on knots
    mesh = gen2dMesh(coords_knots[[m]], bnd)
    graph0 = constrainedDentri(hyperpar['nref'], mesh,
                               gaurantee_connected = T, dist_mat = dist_mat_knot)
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    E(graph0)$weight=runif(ecount(graph0))
    graphs[[m]] = graph0
  }
  
  # ensuring spatial graphs are connected
  for (m in 1:M) {
    if (components(graphs[[m]])$no != 1)
      stop(paste("Disconnected graph:", m))
  }
  
  # assign observations to their nearest knots
  # projections[i, j] is the index of the nearest knot of obs. i in weak learner j
  # similirly, projections_ho is for hold-out locations
  projections = array(0, dim = c(n, M))
  projections_ho = array(0, dim = c(n_ho, M))
  for (m in 1:M) {
    # get distance between training locations and knots
    cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
    projections[, m] = apply(cdist_mat_knots, 1, which.min)
    
    # get distance between hold-out locations and knots
    cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
    projections_ho[, m] = apply(cdist_mat_ho_knots, 1, which.min)
  }
  
  return(list(graphs=graphs,projections=projections,projections_ho=projections_ho))
}

GenUshape_G_BAMDT_Sacramento=function(sf_bnd,coords,coords_ho,hyperpar,coords_ref=NULL){
  
  coords_all = rbind(coords, coords_ho)
  n=nrow(coords);n_ho=nrow(coords_ho)
  M=hyperpar['n_rounds']
  n_ref=hyperpar['nref']
  graphs = list()       # spatial graph for each tree
  projections = array(0, dim = c(n, M))
  projections_ho = array(0, dim = c(n_test, M))
  
  for (m in 1:M) {
    
    if(is.null(coords_ref)){
      Vmesh = genVMesh(sf_bnd,
                       n_ref=hyperpar['nref'],
                       graph = TRUE,
                       type = 'random')
    }else{
      Vmesh = genVMesh(sf_bnd,
                       coords_ref=coords_ref[[m]][1:hyperpar['nref']],
                       graph = TRUE,
                       type = 'random')
    }
    
    graph0 = Vmesh$g
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    E(graph0)$weight=runif(ecount(graph0))
    graphs[[m]] = graph0
    projections[, m] = as.numeric(st_within(sf_coords,Vmesh$mesh))
    projections_ho[, m] = as.numeric(st_within(sf_coords_ho,Vmesh$mesh))
  }
  
  return(list(graphs=graphs,projections=projections,projections_ho=projections_ho))
}




### Generate Candidate Graphs for BAMDT given dist_all

GenUshape_G_BAMDT=function(coords,coords_ho,hyperpar,dist_mat_all,coords_ref=NULL){
  
  n=nrow(coords);n_ho=nrow(coords_ho)
  #rm(dist_res)
  M=hyperpar['n_rounds']
  if(is.null(coords_ref)){
    coords_knots = list()
  }else{
    coords_knots=coords_ref
  }
  # knot coordinates for each tree
  graphs = list()       # spatial graph for each tree
  knot_idx = list()
  for (m in 1:M) {
    # subsample training locations as knots
    if(is.null(coords_ref)){
      knot_idx[[m]] = sample.int(n, hyperpar['nref'])
      coords_knots[[m]] = coords[knot_idx[[m]], ]
    }
    
    dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
    
    # get CDT graph on knots
    mesh = gen2dMesh(coords_knots[[m]], bnd)
    graph0 = constrainedDentri(hyperpar['nref'], mesh,
                               gaurantee_connected = T, dist_mat = dist_mat_knot)
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    E(graph0)$weight=runif(ecount(graph0))
    graphs[[m]] = graph0
  }
  
  # ensuring spatial graphs are connected
  for (m in 1:M) {
    if (components(graphs[[m]])$no != 1)
      stop(paste("Disconnected graph:", m))
  }
  
  # assign observations to their nearest knots
  # projections[i, j] is the index of the nearest knot of obs. i in weak learner j
  # similirly, projections_ho is for hold-out locations
  projections = array(0, dim = c(n, M))
  projections_ho = array(0, dim = c(n_ho, M))
  for (m in 1:M) {
    # get distance between training locations and knots
    cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
    projections[, m] = apply(cdist_mat_knots, 1, which.min)
    
    # get distance between hold-out locations and knots
    cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
    projections_ho[, m] = apply(cdist_mat_ho_knots, 1, which.min)
  }
  
  return(list(graphs=graphs,projections=projections,projections_ho=projections_ho))
}

### Generate Candidate Graphs for BAMDT given locations

GenUshape_G_BAMDT_V2=function(sf_bnd,coords,coords_ho,hyperpar,coords_ref=NULL){
  
  coords_all = rbind(coords, coords_ho)
  n=nrow(coords);n_ho=nrow(coords_ho)
  M=hyperpar['n_rounds']
  n_ref=hyperpar['nref']
  graphs = list()       # spatial graph for each tree
  projections = array(0, dim = c(n, M))
  projections_ho = array(0, dim = c(n_test, M))
  
  for (m in 1:M) {
    
    if(is.null(coords_ref)){
      Vmesh = genVMesh(sf_bnd,
                       n_ref=hyperpar['nref'],
                       graph = TRUE,
                       type = 'random')
    }else{
      Vmesh = genVMesh(sf_bnd,
                       coords_ref=coords_ref[[m]][1:hyperpar['nref']],
                       graph = TRUE,
                       type = 'random')
    }
    
    graph0 = Vmesh$g
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    E(graph0)$weight=runif(ecount(graph0))
    graphs[[m]] = graph0
    projections[, m] = as.numeric(st_within(sf_coords,Vmesh$mesh))
    projections_ho[, m] = as.numeric(st_within(sf_coords_ho,Vmesh$mesh))
  }
  
  return(list(graphs=graphs,projections=projections,projections_ho=projections_ho))
}



### Generate Torus Example

evalFunTorusBAMDT=function(coords_polar, X) {
  theta = coords_polar[, 3]; phi = coords_polar[, 4]
  # f = phi + X[, 1] * sin(theta)
  f = X[, 1] * phi + 0.3 * sin(theta)
  return(f)
}

genTorus_fG = function(p, hyperpar,seed=1234) {
  load('torus_dist_mat.RData')
  ## or run the following code to generate torus_dist_mat (caution: slow)
  ##
  ## define a torus
  # R = 6; r = 4
  # phi_min = -pi / 6; phi_max = (2 - 1/3) * pi
  #
  # # generate n uniform training locations
  # n = 500
  # coords_p = genLocationsTorus(n, R, r, phi_min, phi_max)$coords_p
  # coords_c = polar2Cartesian(coords_p)
  #
  # # generate n_ho uniform hold-out locations
  # n_ho = 200
  # coords_p_ho = genLocationsTorus(n_ho, R, r, phi_min, phi_max)$coords_p
  # coords_c_ho = polar2Cartesian(coords_p_ho)
  #
  # # get geodesic distance matrices
  # coords_p_all = rbind(coords_p, coords_p_ho)
  # dist_mat_all = gdistTorus(coords_p_all, R, r, phi_min, phi_max,
  #                           n_theta = 200, n_phi = 200, k_nn = 8,
  #                           return_both = F)
  # save(n, n_ho, coords_p,coords_c,coords_p_ho, coords_c_ho,dist_mat_all, file = 'torus_dist_mat.RData')
  # ##
  ##
  ##
  set.seed(seed)
  X_all = matrix(0, nrow = n + n_ho, ncol = p)
  X_all = simGpTorus(p, rbind(coords_p,coords_p_ho))
  colnames(X_all) = paste("X", 1:p, sep = "")
  
  X = X_all[1:n, , drop = F]
  X_ho = X_all[-(1:n), , drop = F]
  
  colnames(X) = paste("X", 1:p, sep = "")
  colnames(X_ho) = paste("X", 1:p, sep = "")
  
  # true partition
  cluster_true = ifelse(coords_p[, 4] < 2/3*pi, 1, 0)
  idx = cluster_true == 0
  cluster_true[idx] = ifelse(coords_p[idx, 4] < 7/6*pi, 2, 3)
  
  cluster_ho_true = ifelse(coords_p_ho[, 4] < 2/3*pi, 1, 0)
  idx = cluster_ho_true == 0
  cluster_ho_true[idx] = ifelse(coords_p_ho[idx, 4] < 7/6*pi, 2, 3)
  
  # cluster_true = ifelse(coords_p[, 3] < coords_p[, 4], 1, 2)
  # cluster_ho_true = ifelse(coords_p_ho[, 3] < coords_p_ho[, 4], 1, 2)
  
  # true function
  f_true = evalFunTorusBAMDT(coords_p, X)
  f_ho_true = evalFunTorusBAMDT(coords_p_ho, X_ho)
  
  f_true[cluster_true == 2] = -f_true[cluster_true == 2]
  f_ho_true[cluster_ho_true == 2] = -f_ho_true[cluster_ho_true == 2]
  
  
  Graphs = list()
  message("Generating Graphs")
  
  coords_knots = list()
  knot_idx = list()
  for (m in 1:hyperpar['n_rounds']) {
    # subsample training locations as knots
    knot_idx[[m]] = sample.int(n, hyperpar['nref'])
    coords_knots[[m]] = coords_p[knot_idx[[m]], ]
    dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
    
    graph0 = KNNGraph(dist_mat_knot, k_nn = 5)
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    g0 = graph0
    
    cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
    tr2mesh = apply(cdist_mat_knots, 1, which.min)
    
    # get distance between test locations and knots
    cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
    ho2mesh= apply(cdist_mat_ho_knots, 1, which.min)
    graph_attr(g0,'tr2mesh')=tr2mesh
    graph_attr(g0,'ho2mesh')=ho2mesh
    vertex_attr(g0)=list(mesh2tr=split(1:n,factor(graph_attr(g0,'tr2mesh'),1:hyperpar['nref'])),
                         mesh2ho=split(1:n_ho,factor(graph_attr(g0,'ho2mesh'),1:hyperpar['nref'])))
    
    graphs = Generate_GraphCandidates(g0=g0,X=X,X_ho=X_ho,
                                      num.ost=hyperpar['num_ost'],
                                      n_bin=20)
    Graphs[[m]]=graphs
  }
  
  
  train_X = as.matrix(cbind(coords_c, X))
  test_X = as.matrix(cbind(coords_c_ho, X_ho))
  
  return(list(
    # Y = Y,
    #  Y_ho = Y_ho,
    X = X,
    X_ho = X_ho,
    coords = coords_c,
    coords_ho = coords_c_ho,
    coords_p=coords_p,
    coords_p_ho=coords_p_ho,
    train_X = train_X,
    test_X = test_X,
    Graphs = Graphs,
    dist_mat_all=dist_mat_all,
    f_true = f_true,
    f_ho_true = f_ho_true
  ))
  
}

GenTorus_G_BAMDT=function(coords_p,hyperpar,dist_mat_all){
  
  n=nrow(coords_p);n_ho=nrow(dist_mat_all)-n
  coords_knots = list()
  graphs = list()
  knot_idx = list()
  M=hyperpar['n_rounds']
  n_knot=hyperpar['nref']
  for (m in 1:M) {
    # subsample training locations as knots
    knot_idx[[m]] = sample.int(n, n_knot)
    coords_knots[[m]] = coords_p[knot_idx[[m]], ]
    dist_mat_knot = dist_mat_all[knot_idx[[m]], knot_idx[[m]]]
    
    graph0 = KNNGraph(dist_mat_knot, k_nn = 5)
    E(graph0)$eid = as.integer(1:ecount(graph0))  # edge id
    V(graph0)$vid = as.integer(1:vcount(graph0))  # vertex id
    graphs[[m]] = graph0
  }
  
  # check connectivity
  for (m in 1:M) {
    # print( plotGraph(coords_knots[[m]], graphs[[m]], title = m) )
    if (components(graphs[[m]])$no != 1)
      stop(paste("Disconnected graph:", m))
  }
  
  # make mapping from observations to knots
  projections = matrix(0, nrow = n, ncol = M)
  projections_ho = array(0, dim = c(n_ho,  M))
  # cdist_mat_ho_knots = list()
  for (m in 1:M) {
    cdist_mat_knots = dist_mat_all[ 1:n, knot_idx[[m]] ]
    projections[, m] = apply(cdist_mat_knots, 1, which.min)
    
    # get distance between test locations and knots
    cdist_mat_ho_knots = dist_mat_all[ (n + 1):(n + n_ho), knot_idx[[m]] ]
    projections_ho[,  m] = apply(cdist_mat_ho_knots, 1, which.min)
  }
  
  
  return(list(graphs=graphs,projections=projections,projections_ho=projections_ho))
}


logistic_func = function(x) {
  1/(1+exp(-x))
}
logit_func = function(x) {
  log(x/(1-x))
}
softmax_func = function(x) {
  exp(x)/rowSums(exp(x))
}







