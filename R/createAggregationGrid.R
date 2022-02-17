createAggregationGrid <- function(delta = 50000,
                                  crs = st_crs(2163),
                                  dir) {
  if (class(crs)!="crs") stop("crs must be of class 'crs'.")
  
  boundary <- if (!missing(dir)) st_read(dir) else asgerbachelor::ecoregions_L0
  
  boundary %>% 
    st_bbox(crs=crs) %>% 
    asgerbachelor::bboxToRegularGrid(delta) %>%
    mutate(ID = row_number()) %>% 
    st_as_stars() %>% 
    st_set_crs(crs)
}
