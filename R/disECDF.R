#' Calculates the weighted empirical cumulative density function
#' 
#' Calculates the weighted empirical cumulative density of a vector, given a vector or matrix of weights.
#' The ECDF is evaluated at equal intervals, with a user specified spacing.
#' 
#' @param dis numeric vector. A vector of observations to calculate the ECDF on.
#' @param SP a numeric vector or matrix. A vector or matrix of weights to calculate the ECDF on.
#' @param spacing a double in [0,max(dis)]. A double indicating the spacing of intervals to evaluate the ECDF on.
#' @param scale a logical. If true 'dis' is range-scaled, such that 'dis' is in [0,1].
#' @return a table containing the weighted ECDF of 'dis' weighted by 'SP' evaluated at intervals of 'quantileSpacing'.
#' 
#' @section Details:
#' This function calculates the weighted ECDF efficiently by constructing a matrix, with columns corresponding increasing values in the range of 'dis' and row correspond to entries in 'dis', where the values indicate whether the given value of 'dis' is larger than the column value.
#' 
#' m_{i,j} = 1_{\[v_i,dis_j\]}
#' 
#' where "1_{f(x)}" is a indicator function. 'v' is then a vector of values:
#' 
#' v_i = min(dis) + spacing * i
#' 
#' By computing:
#' 
#' ecdf = w \cdot m 
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
#' disUnif <- runif(100)
#' disNorm <- rnorm(100)
#' tsp <- matrix(rpois(100*50,sample(5,100,T)^2),ncol=100,byrow=T)
#' 
#' tibble(
#'   ecdf = list(disECDF(disUnif,tsp,scale=T),
#'               disECDF(disNorm,tsp,scale=T)),
#'   distribution = c("Uniform", "Normal (1-D)")
#' ) %>%
#'   unnest(ecdf) %>%
#'   ggplot(aes(distance,ecdf,fill=distribution,color=distribution)) +
#'   geom_point(size = .5) +
#'   stat_summary_bin(geom = "ribbon", alpha = .5,
#'                    fun.data = mean_sdl,bins = 50) +
#'   coord_cartesian(xlim = 0:1, ylim = 0:1, clip = "off", expand = F)


disECDF <- function(dis,
                    SP,
                    spacing = .01,
                    scale = F) {
  if (scale) dis <- (dis - min(dis,na.rm=T))/diff(range(dis,na.rm=T))
  
  breaks <- seq(min(dis,na.rm=T),max(dis,na.rm=T),spacing)
  
  disECDF <- sapply(breaks, function(x) dis<x)
  
  spECDF <- (SP / rowSums(SP)) %*% disECDF
  
  spECDF %>% 
    as.data.frame %>% 
    rownames_to_column("ID") %>% 
    as_tibble %>% 
    pivot_longer(!ID,names_to="dummy",values_to="ecdf") %>% 
    group_by(ID) %>% 
    mutate(distance = breaks) %>% 
    select(!dummy) %>%
    ungroup 
}




