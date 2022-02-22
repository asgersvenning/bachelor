#' Create a regular aggregation grid 
#' 
#' \code{createAggregationGrid} Creates a regular stars grid, which can be used for aggregating an sf polygon sfc. 
#' 
#' The created grid is centered on (0,0)
#' 
#' @param object a sf object with a bounding box which can be obtained via sf::st_bbox.
#' @param delta numeric. Width of cells in the grid (in the units of the supplied crs).
#' @param crs an object of class crs. Desired coordinate reference system (CRS) of the output grid.
#' 
#' @return An object of class stars with 2 dimensions and 1 attribute, ID, which is unique for each grid cell.
#' 
#' @section Warnings:
#' It is important that the supplied CRS is uses a equal area projection, to avoid spatial biases in the size of the aggregation cells.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export

createAggregationGrid <- function(object = asgerbachelor::ecoregions_L0,
                                  delta = 50000,
                                  crs = st_crs(2163)) {
  if (class(crs)!="crs") stop("crs must be of class 'crs'.")
  
  object %>% 
    sf::st_bbox(crs=crs) %>% 
    asgerbachelor::bboxToRegularGrid(delta) %>%
    mutate(ID = dplyr::row_number()) %>% 
    stars::st_as_stars() %>% 
    sf::st_set_crs(crs)
}
