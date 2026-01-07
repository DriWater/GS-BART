# Compute RMSE at each sweep

library(deldir)
library(igraph)
library(foreach)
library(doParallel)
library(doRNG) 

x2sfpoint=function(x){
  sf::st_as_sf(data.frame(x=unname(as.data.frame(x)),y=rep(0,length(x))),coords=c('x','y'))
}

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

RMSETrace <- function(Y_out, Y_ho) {
  rmse_trace <- apply(Y_out, 1, function(x) {
    sqrt(mean((x - Y_ho)^2))
  }) 
}

vcount_sub_ost = function(graph,root){
  n=vcount(graph$g)
  vcount = rep(0,n)
  vcount = vcount_pass(graph,root,n)
  #vcount = vcount + branch
  return(vcount)
}

############## orient the tree by specifying a root
as.ost = function (g.st, root=NULL){
  n.v=vcount(g.st)
  # randomly pick a root vertex if root vertex is not specified
  if(is.null(root)){
    root=sample(n.v,1)
  }
  ## Find parent of each vertex
  to=GetPredecessors(g.st,root)
  to[to==0]=NA
  from = 1:n.v
  # get directed edge list
  e=data.frame(from=from,to=to) 
  # generate an igraph directed oriented tree
  g = graph_from_edgelist(as.matrix(e)[-root,],directed=TRUE)  
  vertex_attr(g)=vertex_attr(g.st) 
  graph_attr(g)=graph_attr(g.st) 
  ## terminal leaves
  leaves=which(igraph::degree(g,mode='in')==0)
  # get children vertices for each vertex
  children=split(e$from,factor(e$to,levels=1:n.v))
  children[leaves] <- list(rep(NULL,length(leaves)))
  
  return(list(root=root, leaves=leaves, children=children, e=data.frame(from=from,to=to), g=g)) 
} 

compute_coverage=function(posterior_samples, y_true,...){

upper.lvl = NULL
lower.lvl = NULL
for(i in 1:ncol(posterior_samples)){
  tmp <- ci(posterior_samples[,i], ...)
  upper.lvl <- append(upper.lvl, tmp$CI_high)
  lower.lvl <- append(lower.lvl, tmp$CI_low)
}

coverage <- mean((y_true < upper.lvl) & (y_true > lower.lvl))
ci_len <- mean(upper.lvl - lower.lvl)
return(list(coverage=coverage,ci_len=ci_len))
}
 

geweke_test=function(sigma_out){
  geweke_results=(geweke.diag(sigma_out))
  # Extract z-scores
  zvals <- geweke_results$z
  
  # Compute two-sided p-values
  pvals <- 2 * (1 - pnorm(abs(zvals)))
  return(pvals)
}
plot_predCI=function(Y_out,Y_ho){
  posterior_samples=t(Y_out)
  posterior_mean <- apply(posterior_samples, 1, mean)
  lower_95 <- apply(posterior_samples, 1, quantile, probs = 0.025)
  upper_95 <- apply(posterior_samples, 1, quantile, probs = 0.975)
  df <- data.frame(
    y = Y_ho,
    posterior_mean = posterior_mean,
    lower_95 = lower_95,
    upper_95 = upper_95
  )
  
  g=ggplot(df, aes(x = y, y = posterior_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.1, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")+
    theme_minimal()
  return(g)
}

######################################
GetPredecessors <- function(g.st, root){
  pathdata <- igraph::get.shortest.paths(graph = g.st, from = root)
  pathlist <- prelist <- lapply(pathdata$vpath, function(x){
    y = as.vector(x)
    y[(length(y)-1)]
  })
  prelist[[root]] <- 0
  return(unlist(prelist))
}

############## 'mapping' is now called 'children'
Generate_OST = function (g0, root=NULL){
  n.v=igraph::vcount(g0)
  ## generate a spanning tree of graph0 by finding mst with random edge weights
  g.st=igraph::mst(g0,weight=runif(ecount(g0)))
  # randomly pick a root vertex if root vertex is not specified
  if(is.null(root)){
    root=sample(n.v,1)
  }
  ## Find parent of each vertex
  to=GetPredecessors(g.st,root)
  to[to==0]=NA
  from = 1:n.v
  # get directed edge list
  e=data.frame(from=from,to=to)
  # generate an igraph directed oriented tree
  g = igraph::graph_from_edgelist(as.matrix(e)[-root,],directed=TRUE)
  # vertex_attr(g)=igraph::vertex_attr(g.st)
  # graph_attr(g)=igraph::graph_attr(g.st)
  ## terminal leaves
  leaves=which(igraph::degree(g,mode='in')==0)
  # get children vertices for each vertex
  children=split((e$from-1),factor((e$to-1),levels=0:(n.v-1)))
  children[leaves] <- list(rep(numeric(0),length(leaves)))
  children = unname(children)
  mesh2tr = unname(igraph::vertex.attributes(g.st)$mesh2tr)
  mesh2ho = unname(igraph::vertex.attributes(g.st)$mesh2ho)
  tr2mesh = unname(igraph::graph.attributes(g.st)$tr2mesh)
  # ho2mesh = unname(igraph::graph.attributes(g.st)$ho2mesh)
  if(is.null(mesh2ho)){
    return(list(root = (root - 1), children = children,
                mesh2tr = mesh2tr, tr2mesh = tr2mesh))
  }
  return(list(root = (root - 1), children = children,
              mesh2tr = mesh2tr, mesh2ho = mesh2ho,
              tr2mesh = tr2mesh))
}

Generate_Chain = function(x, x_ho, n_bin=20, type = 'regular'){
  pred_flag = !is.null(x_ho)
  rangex=range(c(x))
  if(is.null(ncol(x))) {
    n_train = length(x)
  } else {
    n_train = nrow(x)
  }
  if(!is.null(x_ho)){
    if(is.null(ncol(x_ho))) {
      n_test = length(x_ho)
    } else {
      n_test = nrow(x_ho)
    }
  }
  if(rangex[1] == rangex[2]) {
    warning("An input feature has no variation; generating a singleton chain\n")
    root = 1
    leaves = 1
    children = list(1)
    children[1] = list(numeric(0))
    mesh2tr= list(seq(1, n_train))
    tr2mesh = rep(1, n_train)
    if(pred_flag){
      mesh2ho = list(seq(1, n_test))
      ho2mesh = rep(1, n_test)
    }else{
      mesh2ho = vector(mode = 'list', length = 1)
    }
  }
  xline=sf::st_linestring(matrix(c(rangex[1],rangex[2],0,0),2,2))
  uniq_x_all=unique(c(x)) # still awaiting to be solved ?
  is.cat=length(uniq_x_all)<n_bin
  if(is.cat){
    n_bin=length(uniq_x_all)
    sort_uniq = sort(uniq_x_all)
    midpoints <- (head(sort_uniq, -1) + tail(sort_uniq, -1)) / 2
    multipoints=c(-Inf, midpoints, Inf)
  }else{
    
    ## change type to regular if equal space
    if(type == 'stratified'){
      # Breakpoints of bins (length n_bin + 1)
      breaks <- seq(rangex[1], rangex[2], length.out = n_bin)
      
      # Sample one point uniformly from each bin
      cutpoints <- sapply(1:(n_bin-1), function(i) {
        runif(1, min = breaks[i], max = breaks[i + 1])
      }) 
      multipoints=c(-Inf,cutpoints,Inf)
      }else{
    multipoints=c(-Inf, sort(sf::st_coordinates(sf::st_line_sample(xline, n=n_bin-1,type=type))[,1]), Inf)
    }
  }
  # lines <- lapply(1:(length(multipoints)-1), function(i) {
  #   sf::st_linestring(matrix(c(multipoints[i], multipoints[i+1],0,0),2,2)) %>% sf::st_sfc(.)
  # })
  # sf_lmesh <- do.call('c', lines)
  root = n_bin-1
  children = as.list(c(NA,0:(n_bin-2)))
  leaves=1
  children[leaves]= list(numeric(0))
  children = unname(children)
  tr2mesh = findInterval(x,multipoints) - 1
  mesh2tr = unname(split(0:(n_train-1),factor(tr2mesh, 0:(n_bin-1))))
  # root = 0
  # children = as.list(c(1:(n_bin-1),NA))
  # leaves=n_bin; children[leaves]= list(numeric(0))
  # children = unname(children)
  # tr2mesh = findInterval(x,multipoints) - 1
  # mesh2tr = unname(split(0:(n_train-1),factor(tr2mesh, 0:(n_bin-1))))
  if(pred_flag){
    ho2mesh = findInterval(x_ho, multipoints) - 1
    mesh2ho =  unname(split(0:(n_test-1),factor(ho2mesh, 0:(n_bin-1))))
    return(list(root = root, children = children, mesh2tr = mesh2tr,
                mesh2ho = mesh2ho,tr2mesh = tr2mesh))
  }else{
    return(list(root = root, children = children, mesh2tr = mesh2tr, tr2mesh = tr2mesh))
  }
}

x2sfpoint = function(x){
  sf::st_as_sf(data.frame(x=unname(as.data.frame(x)),y=rep(0,length(x))),coords=c('x','y'))
}

### Generate candidate graphs (num.ost trees + p chains) to split
Generate_GraphCandidates = function(g0=NULL, X, X_ho=NULL, num.ost=NULL,n_bin=20, chain.type = 'regular'){
  ## X is the unstructured features
  ## g0 is the structural graph
  ## num.ost is the number of
  if(is.null(X)){
    output = vector(mode='list',length=num.ost)
    for(i in 1:num.ost){
      output[[i]] = Generate_OST(g0)
    }
  }else if(is.null(g0)){
    p = ncol(X)
    output = vector(mode='list',length=p)
    for(i in 1:p){
      output[[i]]=Generate_Chain(X[,i],X_ho[,i],n_bin, type = chain.type);
    }
  }else{
    p = ncol(X)
    output = vector(mode='list',length=p+num.ost)
    for(i in 1:num.ost){
      output[[i]] = Generate_OST(g0)
    }
    for(i in 1:p){
      output[[num.ost+i]]=Generate_Chain(X[,i],X_ho[,i],n_bin, type = chain.type);
    }
  }
  return(output)
}

## Recursively generate balanced partition given a general graph
constructBalanceCluster=function(graph, nclust) {
  #graph=g0
  # First, get the number of vertices
  N = igraph::vcount(graph)
  # Check that nclust is valid
  if(nclust > N) {
    stop("nclust cannot be larger than the number of vertices in graph")
  }
  if(nclust < 1) {
    stop("nclust must be positive and at least 1")
  }
  if(!igraph::is.connected(graph)) {
    warning("Input graph is not connected; clusters may not be contiguous\n")
  }
  # Allocate the membership vector
  membership = rep(0, N)
  # Construct minimum spanning tree on the graph
  minspantree = igraph::mst(graph,weights=runif(ecount(graph))) 
  ## convert to OST
  minspantree.ost=as.ost(minspantree)
  
  # target cluster size
  avg_size=floor(N/nclust)
  g.forest=minspantree.ost 
  edgelist=matrix(NA,nrow=nclust-1,ncol=2)
  ## recursively bipartition the graph
  for (i in 1:(nclust-1)){ 
    ## find which edge to remove such that it results in a desired cluster size
    vi=which.min(abs(vcount_sub_ost(g.forest,minspantree.ost$root)-avg_size))[1]
    vj=minspantree.ost$e$to[vi]
    g.forest=Split_OST(g.forest,vi)
    edgelist[i,]=c(vi,vj)
  }    
  ## membership of nodes
  membership = igraph::components(g.forest$g)$membership
  edgelist=cbind(membership[edgelist[,1]],membership[edgelist[,2]])
  return(list(g=g.forest,membership=membership,binedgelist=edgelist))
}

### Generate candidate graphs (num.ost trees + p chains) to split
Generate_GraphCandidates_graphbin = function(g0=NULL, X=NULL,X_ho=NULL, trainid, testid, n_ost=NULL, n_ref = 20, n_bin=20, chain.type = 'regular'){
  ## X is the unstructured features
  ## g0 is the structural graph
  ## num.ost is the number of     
  if(is.null(X)){ 
    output = vector(mode='list',length=num.ost)
    for(i in 1:num.ost){
      output[[i]] = Generate_OST_graphbin(g0,n_ref,trainid,testid)
    }
  }else{
    p = ncol(X)
    if(is.null(g0)){
      output = vector(mode='list',length=p)
      for(i in 1:p){
        output[[i]]=Generate_Chain(X[,i], X_ho[,i], n_bin, type = chain.type);
      }
    }else{
      output = vector(mode='list',length=p+n_ost)
      for(i in 1:n_ost){
        output[[i]] = Generate_OST_graphbin(g0,n_ref,trainid,testid)
      }
      p = ncol(X)
      for(i in 1:p){
        output[[n_ost+i]]=Generate_Chain(X[,i],X_ho[,i], n_bin, type = chain.type);
      }
    }
  }
  return(output)
}


### Creating an appropriate structural graph when starting with a general graph
# Output is connected, undirected, simple, with vertex attributes mesh2tr and mesh2ho;
# graph attributes tr2mesh and ho2mesh training/testing respectively
# One way to connect is to consider each disconnected component to be its own structured feature
# The default (easier) way is to randomly add edges to force it to start as connected
# Constructing the bins can be done randomly or in a better way like considering balancedness
StructurizeGraph = function(inp_graph, training_ids, bins = 0.05,
                            method = "random", connect_method = "random") {
  n_orig = vcount(inp_graph)
  testing_ids = (1:n_orig)[-training_ids]
  pred_flag = (length(testing_ids)>0)
  if(bins >= 1) {
    if(bins >= n_orig) {
      bins = n_orig
    }
  } else {
    bins = floor(bins*n_orig)
  }
  if(pred_flag){
    simplified_graph = inp_graph %>%
      igraph::as.undirected(mode = "collapse") %>%
      igraph::simplify() %>%
      igraph::set_vertex_attr("mesh2tr", index = training_ids, value = (1:length(training_ids)-1)) %>%
      igraph::set_vertex_attr("mesh2tr", index = testing_ids, value = c()) %>%
      igraph::set_vertex_attr("mesh2ho", index = testing_ids, value = (1:length(testing_ids)-1)) %>%
      igraph::set_vertex_attr("mesh2ho", index = training_ids, value = c())
  }else{
    simplified_graph = inp_graph %>%
      igraph::as.undirected(mode = "collapse") %>%
      igraph::simplify() %>%
      igraph::set_vertex_attr("mesh2tr", index = training_ids, value = (1:length(
        training_ids)-1)) %>%
      igraph::set_vertex_attr("mesh2tr", index = testing_ids, value = c())
  }
  simplified_graph_components = igraph::components(simplified_graph)
  if(simplified_graph_components$no > 1) {
    clust_sizes = simplified_graph_components$csize
    n_components = simplified_graph_components$no
    membership_vec = simplified_graph_components$membership
    vtx_seq = c()
    if(connect_method == "random") {
      biggest_cluster = which.max(clust_sizes)
      biggest_cluster_vertices = which(membership_vec == biggest_cluster)
      clust_ids_to_merge = (1:n_components)[-biggest_cluster]
      for(clust_id in clust_ids_to_merge) {
        merge_vertices = which(membership_vec == clust_id)
        merge_vertex_small = sample(merge_vertices, 1)
        merge_vertex_big = sample(biggest_cluster_vertices, 1)
        vtx_seq = c(vtx_seq, merge_vertex_small, merge_vertex_big)
      }
    }
    if(connect_method == "breakoff") {
      stop("Feature breakout connection not yet implemented")
    }
    
    simplified_graph = igraph::add_edges(simplified_graph, vtx_seq)
  }
  
  if(method == "betweenness") {
    mapping = igraph::cluster_edge_betweenness(simplified_graph,
                                               directed = FALSE,
                                               modularity = FALSE,
                                               membership = FALSE) %>%
      igraph::cut_at(no = bins)
    
  } else if(method == "random") {
    igraph::E(simplified_graph)$weight = runif(ecount(simplified_graph))
    random_mst = igraph::mst(simplified_graph)
    edges_to_remove = sample(E(random_mst), bins - 1, replace = FALSE)
    mapping = random_mst %>%
      igraph::delete.edges(edges_to_remove) %>%
      igraph::components() %>%
      `$`("membership")
  } else if(method == "greedy") {
    mapping = igraph::cluster_fast_greedy(simplified_graph) %>%
      igraph::cut_at(no = bins)
  } else if(method == "balanced"){
    mapping = constructBalanceCluster(simplified_graph, bins) %>%
      `$`("membership")
  } else {
    stop("Unknown method provided")
  }
  if(pred_flag){
    structured_binned_graph = igraph::contract(simplified_graph, mapping,
                                               vertex.attr.comb = function(x){as.numeric(na.omit(c(x)))}) %>%
      igraph::simplify() %>%
      igraph::set_graph_attr("tr2mesh", mapping[training_ids] - 1) %>%
      igraph::set_graph_attr("ho2mesh", mapping[testing_ids] - 1)
  }else{
    structured_binned_graph = igraph::contract(simplified_graph, mapping,
                                               vertex.attr.comb = function(x){as.numeric(na.omit(c(x)))}) %>%
      igraph::simplify() %>%
      igraph::set_graph_attr("tr2mesh", mapping[training_ids] - 1)
  }
  return(structured_binned_graph)
}
 
Generate_GraphList_parallel = function(sf_coords = NULL, sf_bnd = NULL, g0 = NULL, X, in_train, nrounds=50,
                              nost=5, nbin=100, nref=100, nthreads = 1, mesh.type = "random", chain.type = "random"){
  nthreads = min(nthreads, nrounds, parallel::detectCores(logical = FALSE))
  cl <- makeCluster(nthreads)
  registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  n = length(in_train); n_ho = 0
  X_ho = NULL
  if(!is.null(X)){
    n_ho = nrow(X) -n;
    if(max(in_train) > nrow(X)){
      stop("training data index out of boundary")
    }
    if(length(in_train) < nrow(X)){ X_ho = X[-in_train, ,drop = FALSE] }
  }
  if(!is.null(sf_bnd)){
    if(is.null(sf_coords)){
      stop("Coordinates of data can not be null once boundary is provided!")
    }
    
    n_ho = length(sf_coords) - n;
    if (mesh.type %in% c('regular', 'hexagonal')) {
  
  mesh <- genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = mesh.type)
  g0 <- mesh$g
  sf_mesh <- mesh$mesh
  vcnt <- vcount(g0)

  igraph::graph_attr(g0, 'tr2mesh') <- unlist(st_within(sf_coords[in_train], sf_mesh)) - 1
  igraph::vertex_attr(g0, 'mesh2tr') <- split(0:(n-1), factor(graph_attr(g0, 'tr2mesh'), 0:(vcnt - 1)))

  if (n_ho > 0) {
    igraph::graph_attr(g0, 'ho2mesh') <- unlist(st_within(sf_coords[-in_train], sf_mesh)) - 1
    igraph::vertex_attr(g0, 'mesh2ho') <- split(0:(n_ho - 1), factor(graph_attr(g0, 'ho2mesh'), 0:(vcnt - 1)))
  }

  Graphs <- replicate(nrounds, Generate_GraphCandidates(
    g0, X = X[in_train, ,drop = FALSE], X_ho = X_ho,
    num.ost = nost, n_bin = nbin, chain.type = chain.type
  ), simplify = FALSE)

} else if (mesh.type %in% c('stratified', 'random')) {
  
  mesh0 <- genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = 'regular')
  Graphs <- foreach::foreach(m = 1:nrounds,
    .packages = c("sf", "igraph", "ComplexDomain", "purrr"),
    .export = c("Generate_GraphCandidates", "Generate_OST", "Generate_Chain", "x2sfpoint", "GetPredecessors")
  ) %dorng% {

    mesh <- switch(mesh.type,
      stratified = {
        points_list <- st_geometry(mesh0$mesh) %>%
          map(~ st_sample(.x, size = 1))

        coords_ref <- st_as_sf(data.frame(
          geometry = st_sfc(unlist(points_list, recursive = FALSE))
        ))
        
        tmp = try(genVMesh(sf_bnd, coords_ref, n_ref = nref, graph = TRUE))
        while("try-error" %in% class(tmp)){
          tmp = try(genVMesh(sf_bnd, coords_ref, n_ref = nref, graph = TRUE))
        }
        tmp
      },
      random = genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = "random"),
      stop("Unknown mesh.type: ", mesh.type)
    )

    g0 <- mesh$g
    sf_mesh <- mesh$mesh
    vcnt <- vcount(g0)

    igraph::graph_attr(g0, 'tr2mesh') <- unlist(st_within(sf_coords[in_train], sf_mesh)) - 1
    igraph::vertex_attr(g0, 'mesh2tr') <- split(0:(n-1), factor(graph_attr(g0, 'tr2mesh'), 0:(vcnt - 1)))

    if (n_ho > 0) {
      igraph::graph_attr(g0, 'ho2mesh') <- unlist(st_within(sf_coords[-in_train], sf_mesh)) - 1
      igraph::vertex_attr(g0, 'mesh2ho') <- split(0:(n_ho - 1), factor(graph_attr(g0, 'ho2mesh'), 0:(vcnt - 1)))
    }

    Generate_GraphCandidates(g0, X = X[in_train, ,drop = FALSE], X_ho = X_ho, num.ost = nost,
                             n_bin = nbin, chain.type = chain.type)
  }
}
} else if (!is.null(sf_coords)) {

  coords <- st_coordinates(sf_coords)
  tri <- deldir::deldir(coords[,1], coords[,2])
  edges_df <- tri$delsgs[, c("ind1", "ind2")]
  g0 <- igraph::graph_from_edgelist(as.matrix(edges_df), directed = FALSE)

  Graphs <- foreach::foreach(m = 1:nrounds,
    .packages = c("sf", "igraph", "ComplexDomain"),
    .export = c("Generate_GraphCandidates", "Generate_OST", "StructurizeGraph", "Generate_Chain", "x2sfpoint", "GetPredecessors")
  ) %dorng% {
    g0_binned <- StructurizeGraph(g0, in_train, bins = nref, method = "greedy", connect_method = "random")
    Generate_GraphCandidates(g0_binned, X = X[in_train, ,drop =FALSE], X_ho = X_ho, num.ost = nost,
                             n_bin = nbin, chain.type = chain.type)
  }

} else if (!is.null(g0)) {

  Graphs <- foreach::foreach(m = 1:nrounds,
    .packages = c("sf", "igraph", "ComplexDomain"),
    .export = c("Generate_GraphCandidates", "Generate_OST", "StructurizeGraph", "Generate_Chain", "x2sfpoint", "GetPredecessors",
                "constructBalanceCluster", "as.ost", "vcount_sub_ost")
  ) %dorng% {
    g0_binned <- StructurizeGraph(g0, in_train, bins = nref, method = "greedy", connect_method = "random")
    Generate_GraphCandidates(g0_binned, X = X[in_train, ,drop =FALSE], X_ho = X_ho, num.ost = nost,
                             n_bin = nbin, chain.type = chain.type)
  }

} else {
  # No mesh, no sf_coords, no g0 — fallback
  Graphs <- foreach::foreach(m = 1:nrounds,
    .packages = c("sf", "igraph", "ComplexDomain"),
    .export = c("Generate_GraphCandidates", "Generate_OST", "StructurizeGraph", "Generate_Chain", "x2sfpoint", "GetPredecessors")
  ) %dorng% {
    Generate_GraphCandidates(g0 = NULL, X = X[in_train, ,drop =FALSE], X_ho = X_ho,
                             num.ost = 0, n_bin = nbin, chain.type = chain.type)
  }
}
}

Generate_GraphList = function(sf_coords = NULL, sf_bnd = NULL, g0 = NULL, X, in_train,
                              nrounds = 50, nost = 5, nbin = 100, nref = 100,
                              mesh.type = "random", chain.type = "random") {
  n  = length(in_train)
  n_ho = 0
  X_ho = NULL
  
  if (!is.null(X)) {
    n_ho = nrow(X) - n
    if (max(in_train) > nrow(X)) stop("training data index out of boundary")
    if (length(in_train) < nrow(X)) X_ho = X[-in_train, ,drop = FALSE]
  }
  
  # --- Case 1: boundary provided ---
  if (!is.null(sf_bnd)) {
    if (is.null(sf_coords)) stop("Coordinates of data cannot be NULL once boundary is provided!")
    
    n_ho = length(sf_coords) - n
    
    if (mesh.type %in% c("regular", "hexagonal")) {
      # Build mesh once, reuse
      mesh <- genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = mesh.type)
      g0   <- mesh$g
      sf_mesh <- mesh$mesh
      vcnt <- igraph::vcount(g0)
      
      igraph::graph_attr(g0, "tr2mesh") <- unlist(sf::st_within(sf_coords[in_train], sf_mesh)) - 1
      igraph::vertex_attr(g0, "mesh2tr") <- split(0:(n-1),
                                                  factor(igraph::graph_attr(g0, "tr2mesh"), 0:(vcnt - 1)))
      
      if (n_ho > 0) {
        igraph::graph_attr(g0, "ho2mesh") <- unlist(sf::st_within(sf_coords[-in_train], sf_mesh)) - 1
        igraph::vertex_attr(g0, "mesh2ho") <- split(0:(n_ho-1),
                                                    factor(igraph::graph_attr(g0, "ho2mesh"), 0:(vcnt - 1)))
      }
      
      Graphs <- replicate(
        nrounds,
        Generate_GraphCandidates(
          g0, X = X[in_train,  ,drop =FALSE], X_ho = X_ho,
          num.ost = nost, n_bin = nbin, chain.type = chain.type
        ),
        simplify = FALSE
      )
      
    } else if (mesh.type %in% c("stratified", "random")) {
      # Build per-iteration mesh
      mesh0 <- genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = "regular")
      
      Graphs <- lapply(seq_len(nrounds), function(m) {
        mesh <- switch(
          mesh.type,
          stratified = {
            points_list <- sf::st_geometry(mesh0$mesh) |>
              purrr::map(~ sf::st_sample(.x, size = 1))
            coords_ref <- sf::st_as_sf(data.frame(
              geometry = sf::st_sfc(unlist(points_list, recursive = FALSE))
            ))
            tmp = try(genVMesh(sf_bnd, coords_ref, n_ref = nref, graph = TRUE))
            while("try-error" %in% class(tmp)){
              tmp = try(genVMesh(sf_bnd, coords_ref, n_ref = nref, graph = TRUE))
            }
            tmp
          },
          random = genVMesh(sf_bnd, n_ref = nref, graph = TRUE, type = "random"),
          stop("Unknown mesh.type: ", mesh.type)
        )
        
        g0   <- mesh$g
        sf_mesh <- mesh$mesh
        vcnt <- igraph::vcount(g0)
        
        igraph::graph_attr(g0, "tr2mesh") <- unlist(sf::st_within(sf_coords[in_train], sf_mesh)) - 1
        igraph::vertex_attr(g0, "mesh2tr") <- split(0:(n-1),
                                                    factor(igraph::graph_attr(g0, "tr2mesh"), 0:(vcnt - 1)))
        
        if (n_ho > 0) {
          igraph::graph_attr(g0, "ho2mesh") <- unlist(sf::st_within(sf_coords[-in_train], sf_mesh)) - 1
          igraph::vertex_attr(g0, "mesh2ho") <- split(0:(n_ho-1),
                                                      factor(igraph::graph_attr(g0, "ho2mesh"), 0:(vcnt - 1)))
        }
        
        Generate_GraphCandidates(
          g0, X = X[in_train, ,drop =FALSE], X_ho = X_ho,
          num.ost = nost, n_bin = nbin, chain.type = chain.type
        )
      })
    } else {
      stop("Unknown mesh.type: ", mesh.type)
    }
    
    # --- Case 2: no boundary, but have coordinates ---
  } else if (!is.null(sf_coords)) {
    
    coords <- sf::st_coordinates(sf_coords)
    tri <- deldir::deldir(coords[,1], coords[,2])
    edges_df <- tri$delsgs[, c("ind1", "ind2")]
    g0 <- igraph::graph_from_edgelist(as.matrix(edges_df), directed = FALSE)
    
    Graphs <- lapply(seq_len(nrounds), function(m) {
      g0_binned <- StructurizeGraph(g0, in_train, bins = nref, method = "greedy", connect_method = "random")
      Generate_GraphCandidates(
        g0_binned, X = X[in_train, ,drop =FALSE], X_ho = X_ho,
        num.ost = nost, n_bin = nbin, chain.type = chain.type
      )
    })
    
    # --- Case 3: no boundary/coords, but have g0 ---
  } else if (!is.null(g0)) {
    
    Graphs <- lapply(seq_len(nrounds), function(m) {
      g0_binned <- StructurizeGraph(g0, in_train, bins = nref, method = "greedy", connect_method = "random")
      Generate_GraphCandidates(
        g0_binned, X = X[in_train, ,drop =FALSE], X_ho = X_ho,
        num.ost = nost, n_bin = nbin, chain.type = chain.type
      )
    })
    
    # --- Case 4: full fallback ---
  } else {
    Graphs <- lapply(seq_len(nrounds), function(m) {
      Generate_GraphCandidates(
        g0 = NULL, X = X[in_train, ,drop = FALSE], X_ho = X_ho,
        num.ost = 0, n_bin = nbin, chain.type = chain.type
      )
    })
  }
  
  return(Graphs)
}

