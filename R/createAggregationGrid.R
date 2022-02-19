#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export

createAggregationGrid <- function(object = asgerbachelor::ecoregions_L0,
                                  delta = 50000,
                                  crs = st_crs(2163),
                                  dir) {
  if (class(crs)!="crs") stop("crs must be of class 'crs'.")
  
  object %>% 
    sf::st_bbox(crs=crs) %>% 
    asgerbachelor::bboxToRegularGrid(delta) %>%
    mutate(ID = dplyr::row_number()) %>% 
    stars::st_as_stars() %>% 
    sf::st_set_crs(crs)
}
