#' Calculates functional evenness
#' 
#' Calculates functional evenness from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @param approxEvenness a logical. Is evenness calculated precisely or estimated?
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @section Details:
#' When using 'approxEvenness=T', the minimum spanning tree is constructed on a KNN. KNN is computed iteratively, with increasing K, until a minimum spanning tree can be constructed.
#' This is not how minimum spanning tree is supposed to be computed efficiently, but I couldn't find any implementation in any R-packages, that do not use the distance matrix to find the minimum spanning tree.
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @importFrom spdep knearneigh
#' @importFrom spdep knn2nb
#' @importFrom spdep nb2listw
#' @importFrom ape mst
#' @export

functional_evenness <- function(x, a = rep(1, nrow(x)), approxEvenness=F) {
  x <- x[n>0,] # Subset present species
  n <- n[n>0]  # Subset present species
  n <- n/sum(n)# Relative abundances
  
  S <- length(n)# Number of species
  Ed <- 1/(S-1) # Expected edge length in minimum spanning tree
  
  if (isFALSE(approxEvenness)) {
    dx <- dist(x) %>% as.matrix 
    l <- unclass(ape::mst(dx))==1 # Logical encoding of minimum spanning tree
    li <- which(l, arr.ind = T) # Indices of mst-vertices
    li.lower <- li[diff(t(li))<0,] # Subset lower triangle
    
    EW <- apply(li.lower, 1, function(y) dx[y[1],y[2]]/sum(n[y]))
  }
  else {
    ## Slight change to spdep::mstree, since it apparently cannot handle listwdist(?)
    mod_spdep_mstree <- function (nbw, ini = NULL) {
      np <- length(nbw$neighbours) # Only change is nbw[[2]] --> nbw$neighbours...
      nodes <- cbind(FALSE, 0, rep(Inf, np)) # (and "n"-->"np" to avoid name collision)
      if (is.null(ini)) 
        ini <- sample(1:np, 1)
      nodes[ini, 1] <- TRUE
      nodes[nbw$neighbours[[ini]], 2] <- ini
      nodes[nbw$neighbours[[ini]], 3] <- nbw$weights[[ini]]
      mst <- matrix(0, np - 1, 3)
      for (i in 1:(np - 1)) {
        id.min <- which.min(nodes[, 3])
        if (!is.finite(nodes[id.min, 3])) 
          stop("Graph is not connected!")
        nodes[id.min, 1] <- TRUE
        mst[i, ] <- c(nodes[id.min, 2], 
                      id.min, 
                      nodes[id.min,3])
        id.out <- !nodes[nbw$neighbours[[id.min]], 1]
        node.can <- nbw$neighbours[[id.min]][id.out]
        node.cost <- nbw$weights[[id.min]][id.out]
        id.best <- node.cost < nodes[node.can, 3]
        nodes[node.can[id.best], 2] <- id.min
        nodes[node.can[id.best], 3] <- node.cost[id.best]
        nodes[id.min, 3] <- Inf
      }
      attr(mst, "class") <- c("mst", "matrix")
      mst
    }
    
    # The approximation will almost always overestimate the evenness, 
    # since it is based on KNN
    dnn <- seq(5*ncol(x),nrow(x),length.out = 5) %>% ceiling # Numbers of neighbours to try
    
    li.lower <- NULL
    for (i in dnn) { # Select a number of neighbours
      if (!is.null(li.lower)) break # (If graph is disconnected, try new neighbour number.)
      try({ # Try to find a minimum spanning trree
        li.lower <- x %>% 
          spdep::knearneigh(i) %>% 
          spdep::knn2nb() %>% 
          spdep::nb2listw() %>% # Generate a knn-list
          {
            .$weights <- lapply(1:length(.$weights), function(y) { # Set the weights equal to the dissimilarity
              sapply(.$neighbours[[y]], function(z) sum((x[y,]-x[z,])^2)) %>% sqrt 
            })
            . 
          } %>%  
          mod_spdep_mstree # Find minimum spanning tree
      },T)
    }
    if (is.null(li.lower)) stop("No connected graph found.")
    EW <- apply(li.lower, 1, function(y) sqrt(sum((x[y[1],]-x[y[2],])^2))/sum(n[1:2])) # Weighted length of mst-edge
  }
  
  PEW <- EW/sum(EW) # Relative weighted length of mst-edge
  
  out <- (sum(pmin(PEW,Ed))-Ed)/(1-Ed)
  
  return(out)
}