#' Converts the output of \code{ecoregionsToGrid} to a RasterStack.
#' 
#' @param ecoregionsGrid a tibble created by \code{ecoregionsToGrid}. 
#' @param grid the same grid which was supplied to \code{ecoregionsToGrid} to create ecoregionsGrid.
#' 
#' @return a RasterStack with a layer for each class in ecoregionsGrid. Values correspond to coverage proportion.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom tidyr complete
#' @importFrom tidyr nest
#' @importFrom tidyselect contains
#' @importFrom purrr map
#' @importFrom methods as
#' @export

ecoregionsGridToRasterstack <- function(ecoregionsGrid, grid) {
  out <- ecoregionsGrid %>% 
    select(c(ID,tidyselect::contains("_KEY"),proportion)) %>% 
    group_by(dplyr::across(tidyselect::contains("_KEY"))) %>% 
    complete(
      ID = grid %>% 
        dim %>% 
        prod %>% 
        seq
    ) %>% 
    ungroup %>% 
    nest(rasterData = !contains("_KEY")) %>% 
    mutate(raster = map(rasterData, function(x) {
      grid %>%
        mutate(proportion = x$proportion) %>% 
        select(!ID) %>% 
        as("Raster")
    })) %>% 
    summarize(
      rasterStack = raster::stack(raster) %>%
        magrittr::set_names(unlist(across(contains("_KEY"), ~paste0("C_", .x)))) %>% 
        list
    ) 
  
  return(out$rasterStack[[1]])
}
