#' Create a regular grid from a bounding box. 
#' 
#' @param bbox a vector giving the bounding box of the grid.
#' @param delta numeric. the width of the grid cells.
#' @return a tibble 
#' 
#' @importFrom tibble tibble
#' @export

bboxToRegularGrid <- function(bbox, delta) {
  if (!is.numeric(bbox) | length(bbox) != 4) stop("Bounding box must be a numeric of length 4.")
  if (!is.numeric(delta) | length(delta) != 1) stop("Delta must be a numeric of length 1.")
  range <- (bbox/delta) %>% {
      ifelse(sign(.)*c(1,1,-1,-1)==1,floor(.),ceiling(.))
    }
  
  xB <- delta*(range[1]:range[3])
  yB <- delta*(range[2]:range[4])
  
  tibble(
    x = rep(xB, length(yB)),
    y = rep(yB, each = length(xB))
  )
}