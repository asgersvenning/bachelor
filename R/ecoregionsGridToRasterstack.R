ecoregionsGridToRasterstack <- function(ecoregionsGrid, grid) {
  out <- ecoregionsGrid %>% 
    select(c(ID,contains("_KEY"),proportion)) %>% 
    group_by(across(contains("_KEY"))) %>% 
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
