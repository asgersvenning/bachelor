#' Calculates the weighted empirical cumulative density function
#' 
#' Calculates the weighted empirical cumulative density of a vector, given a vector or matrix of weights.
#' The ECDF is evaluated at the observed values of 'x'.
#' 
#' @param x numeric vector. A vector of observations to calculate the ECDF on.
#' @param w a numeric vector, matrix or data frame. A vector, matrix or data frame of weights to calculate the ECDF on. If a matrix or data frame, rows are considered vectors of weights.
#' @param scale a logical. If true 'dis' is range-scaled, such that 'dis' is in [0,1].
#' @return a table containing the weighted ECDF of 'dis' weighted by 'SP' evaluated at intervals of 'quantileSpacing'.
#' 
#' @section Details:
#' This function calculates the weighted ECDF efficiently by constructing a matrix, with columns corresponding sorted values in 'x' and row correspond to values of 'x', where the values indicate whether the given value of 'x' is larger than the column value.
#' 
#' m_{i,j} = 1_{\[sort(x)_j >= x_i\]}
#' 
#' By computing:
#' 
#' ecdf = m \cdot w^T 
#' 
#' where w is matrix or vector of weights for the ECDF.
#' 
#' @importFrom magrittr %>% 
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_by 
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom Rdpack reprompt
#' @export
#' 
#' @examples
#' Unif <- runif(100)
#' Norm <- rnorm(100)
#' tw <- matrix(rpois(100*50,sample(5,100,T)^2),ncol=100,byrow=T)
#' 
#' tibble(
#'   ecdf = list(w.ECDF(Unif,tw,scale=T),
#'               w.ECDF(Norm,tw,scale=T)),
#'   distribution = c("Uniform", "Normal (1-D)")
#' ) %>%
#'   unnest(ecdf) %>%
#'   ggplot(aes(x,ecdf,fill=distribution,color=distribution)) +
#'   # geom_point(size = .5) +
#'   stat_summary_bin(geom = "ribbon", alpha = .5,
#'                    fun.data = mean_sdl,bins = 50) +
#'   coord_cartesian(xlim = 0:1, ylim = 0:1, clip = "off", expand = F)
#'   
#' # This example shows the complexity of computation as a function of the length of x and the rows of w.
#' # This can be understood as how fast the computation is carried out given a number of species and a number of sites.
#' res <- lapply((2:15)^3, function(x) {
#'   sapply((2:15)^3, function(y) {
#'     Norm <- rnorm(x)
#'     tw <- matrix(rpois(x*y,sample(5,x,T)^2),ncol=x,byrow=T)
#'     
#'     c(species = x, communities = y, time = system.time(w.ECDF(Norm,tw,scale=T))[3])
#'   })
#' })
#' 
#' res %>% 
#'   lapply(t) %>% 
#'   lapply(as_tibble) %>% 
#'   tibble %>% 
#'   unnest(1) %>% 
#'   ggplot(aes(species^(1/3),communities^(1/3),fill=time.elapsed)) +
#'   geom_raster() +
#'   scale_fill_viridis_c(option = "A",
#'                        trans = "log10",
#'                        na.value = "black",
#'                        limits = c(0.01,50),
#'                        labels = function(x) scales::number(x,accuracy = .01),
#'                        n.breaks = 7) +
#'   scale_x_continuous(labels = function(x) sapply(x, function(y) as.numeric(y)^3),
#'                      n.breaks = 15) +
#'   scale_y_continuous(labels = function(x) sapply(x, function(y) as.numeric(y)^3),
#'                      n.breaks = 15) +
#'   coord_cartesian(expand = F) +
#'   labs(x = "species", y = "communities")


w.ECDF <- function(x,
                   w=rep(1,length(x)),
                   scale = F) {
  if (is.vector(w)) w <- matrix(w,nrow=1)
  if (inherits(w,"data.frame")) w <- as.matrix(w)
  if (!is.vector(x) && dim(x)[1] == 1) x <- as.vector(unlist(x))
  if (!is.matrix(w) | !is.numeric(w)) stop("'w' must be a numeric vector, matrix or data frame.")
  if (!is.vector(x) | !is.numeric(x)) stop("'x' must be a numeric vector.")
  w <- w[,!is.na(x)]
  w <- w / rowSums(w)
  x <- x[!is.na(x)]
  
  
  if (scale) x <- (x - min(x))/diff(range(x))
  sobs <- sort(x)
  
  ECDF <- outer(sobs,x,">=") %*% t(w)
  
  ECDF %>% 
    as.data.frame %>% 
    as_tibble %>% 
    pivot_longer(everything(),names_to="dummy",values_to="ecdf") %>% 
    group_by(dummy) %>% 
    mutate(x = sobs) %>% 
    ungroup %>% 
    select(!dummy) 
}


