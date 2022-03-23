#' Calculates Rao's quadratic entropy
#' 
#' Calculates Rao's quadratic entropy from a species-trait matrix and possibly an species-abundance matrix/vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @param raoScale a logical. Is quadratic entropy divided by maximum value (ade4::divcmax)?
#' @param approxRao an integer. If NOT 0, number of species pairs to estimate quadratic entropy on.
#' @param gower a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.
#' @return a number.
#' 
#' @section Details:
#' This function implements Rao's quadratic entropy, which is defined as the expected dissimilarity between a random pair of individuals from a population.
#' 
#' As such it is quite similar to ade4::divc, however it has 3 major differences; 
#' 
#' 1) It allows the user to approximate the quadratic entropy by sampling a specified number of pairs. (Which unfortunately also makes calculating the maximum possible value impractical.)
#' 2) It handles two dissimilarity measures; Euclidean distance and Gower dissimilarity, which allows the flexibility of using the index on non-ordinated traits.
#' 
#' It is also considerably faster on large numbers of communities. 
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#'  
#'  
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @export
#' 
#' @examples 
#' # An example using the data also supplied in the package. Note that the first part of the example, just shows how to create species-abundance and species-trait tables from the data, as well as subsetting species in the intersection of both data sources.
#' library(hash)
#' FIA_dict <- hash(PLANTS_meta$plants_code, PLANTS_meta$fia_code)
#' PLANTS_dict <- invert(FIA_dict)
#' 
#' PLANTS_tG <- PLANTS_traits %>%
#'   gower_traits(T, NA.tolerance = .1)
#' 
#' tSP <- FIA %>%
#'   filter(INVYR>2000) %>%
#'   group_by(ID,SPCD) %>%
#'   summarise(
#'     dens = mean(individuals/samples, na.rm = T),
#'     .groups = "drop"
#'   ) %>%
#'   filter(SPCD %in% values(FIA_dict,PLANTS_tG$traits$plants.code)) %>%
#'   pivot_wider(ID,SPCD,values_from=dens,names_prefix = "SP",values_fill = 0)
#' 
#' PLANTS_tG <- PLANTS_traits %>%
#'   filter(plants_code %in% values(PLANTS_dict, str_remove(names(tSP)[-1],"^SP"))) %>% 
#'   gower_traits(T)
#' 
#' PCoA_traits <- ape::pcoa(gowerDissimilarity(PLANTS_tG$traits[,-1],PLANTS_tG$weights))
#' 
#' nPC <- PCoA_traits$values[,3:4] %>% 
#'   apply(1,diff) %>% 
#'   sign %>% 
#'   diff %>% 
#'   {
#'     which(.!=0)[1]
#'   } %>% 
#'   names %>% 
#'   as.numeric()
#' 
#' ordTraits <- PCoA_traits$vectors[,1:nPC]
#' 
#' par(mfrow = c(2,3))
#' quadratic_entropy(ordTraits,,tSP[,-1], F, gower = F) %>% hist(25,main = "Absolute Euclidean",xlim=c(0,0.05))
#' quadratic_entropy(ordTraits,,tSP[,-1], F, 10, gower = F) %>% hist(25,main = "Approximate Euclidean",xlim=c(0,0.05))
#' quadratic_entropy(ordTraits,,tSP[,-1], T, gower = F) %>% hist(25,main = "Relative Euclidean",xlim=0:1)
#' quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], F, gower = T) %>% hist(25,main = "Absolute Gower",xlim=c(0,0.05))
#' quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], F,10, gower = T) %>% hist(25,main = "Approximate Gower",xlim=c(0,0.05))
#' quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], T, gower = T) %>% hist(25,main = "Relative Gower",xlim=0:1)
#' 
#' 
#' # Computation time comparison between algorithms.
#' ordDist <- dist(ordTraits)
#' 
#' compRes <- lapply((1:5)^4, function(n) {
#'   wC <- sample(nrow(tSP),n)
#'   
#'   tibble(
#'     communities = n,
#'     mine = system.time(quadratic_entropy(ordTraits,,tSP[wC,-1],T,gower = F))[3],
#'     ade4 = system.time(ade4::divc(as.data.frame(t(tSP[wC,-1])),ordDist,T))[3]
#'   )
#' }) 
#' 
#' compRes %>% 
#'   tibble %>% 
#'   unnest(1) %>%
#'   pivot_longer(c(mine,ade4),names_to="implementation",values_to="time") %>% 
#'   ggplot(aes(communities,time,color=implementation)) +
#'   geom_path() +
#'   geom_path()

quadratic_entropy <- function(x, w = rep(1,ncol(x)), a = rep(1, nrow(x)), raoScale=T, approxRao=0, gower=T) {
  # Adjusted ade4::divcmax to forego the check for euclidean nature, 
  # as well as simply returning the abundance vector which maximizes the entropy.
  divcmax_spec <- function (dis, w=1, epsilon = 1e-08, comment = FALSE, gower=F) {
    if (!inherits(dis, "dist")) 
      stop("Distance matrix expected")
    if (epsilon <= 0) 
      stop("epsilon must be positive")
    D2 <- if (gower) as.matrix(dis)*sum(w)/2 else D2 <- as.matrix(dis)^2/2
    n <- dim(D2)[1]
    result <- data.frame(matrix(0, n, 4))
    names(result) <- c("sim", "pro", "met", 
                       "num")
    relax <- 0
    x0 <- apply(D2, 1, sum)/sum(D2)
    result$sim <- x0
    objective0 <- t(x0) %*% D2 %*% x0
    if (comment == TRUE) 
      print("evolution of the objective function:")
    xk <- x0
    repeat {
      repeat {
        maxi.temp <- t(xk) %*% D2 %*% xk
        if (comment == TRUE) 
          print(as.character(maxi.temp))
        deltaf <- (-2 * D2 %*% xk)
        sature <- (abs(xk) < epsilon)
        if (relax != 0) {
          sature[relax] <- FALSE
          relax <- 0
        }
        yk <- (-deltaf)
        yk[sature] <- 0
        yk[!(sature)] <- yk[!(sature)] - mean(yk[!(sature)])
        if (max(abs(yk)) < epsilon) {
          break
        }
        alpha.max <- as.vector(min(-xk[yk < 0]/yk[yk < 0]))
        alpha.opt <- as.vector(-(t(xk) %*% D2 %*% yk)/(t(yk) %*% 
                                                         D2 %*% yk))
        if ((alpha.opt > alpha.max) | (alpha.opt < 0)) {
          alpha <- alpha.max
        }
        else {
          alpha <- alpha.opt
        }
        if (abs(maxi.temp - t(xk + alpha * yk) %*% D2 %*% 
                (xk + alpha * yk)) < epsilon) {
          break
        }
        xk <- xk + alpha * yk
      }
      if (prod(!sature) == 1) {
        if (comment == TRUE) 
          print("KT")
        break
      }
      vectD2 <- D2 %*% xk
      u <- 2 * (mean(vectD2[!sature]) - vectD2[sature])
      if (min(u) >= 0) {
        if (comment == TRUE) 
          print("KT")
        break
      }
      else {
        if (comment == TRUE) 
          print("relaxation")
        satu <- (1:n)[sature]
        relax <- satu[u == min(u)]
        relax <- relax[1]
      }
    }
    return(xk)
  }
  
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  if (inherits(x,"data.frame")) x <- as.matrix(x)
  if (inherits(a,"data.frame")) a <- as.matrix(a)
  
  if (!is.matrix(a)) stop("Unable to coerce 'a' to matrix.")
  if (!is.matrix(x)) stop("Unable to coerce 'x' to matrix.")
  
  if (missing(w)) w <- rep(1,ncol(x))
  if (!is.vector(w) | !is.numeric(w)) stop("'w' must be a numeric vector.")
  
  a <- a %>% 
    inset(is.na(.),0) %>% 
    divide_by(rowSums(.))

  if (approxRao!=0) {
    if (raoScale) stop("Unable to use divcmax on less than full dissimilarity matrix")
    out <- apply(a, 1, function(c) {
      wSP <- which(!is.na(c))
      cX <- x[wSP,]
      c <- c[wSP]
      
      sapply(1:approxRao, function(y) {
        # Sample a indice-pair
        ind <- sample.int(nrow(cX),2,T,prob=c) 
        
        # Calculate trait-wise differences between the indices
        delta <- abs(cX[ind[1],]-cX[ind[2],]) 
        
        # Calculate Gower dissimilarity
        if (gower) {
          delta %>% 
            inset(is.na(.),0) %>% 
            multiply_by_matrix(w) %>% 
            divide_by((!is.na(delta)) %*% w) %>% 
            raise_to_power(2)
        } else {
          delta %>% 
            raise_to_power(2) %>% 
            sum 
        }
      }) %>% 
        # Mean pair dissimilarity (expected dissimilarity)
        mean / 2
    })
    
    return(out)
  }
  else {
  
  D <- if (gower) gowerDissimilarity(x,w) else as.matrix(dist(x,"euclidean"))
  D[is.na(D)] <- 0
  # D[upper.tri(D,T)] <- NA
  
  scaleD <- if (raoScale) {
    MEB <- as.vector(divcmax_spec(as.dist(D), gower = gower))
    
    as.vector(t(MEB) %*% D^2 %*% MEB) / 2
  } else 1
  
  # For each community
  out <- apply(a, 1, function(c) {
    # Abundance weighted sum of lower triangle of dissimilarity matrix
    (t(c) %*% D^2 %*% c) / 2
    }) 
  return(out / scaleD)
  }
}
