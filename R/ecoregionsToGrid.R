ecoregionsToGrid <- function(grid, 
                             level, 
                             simplify, 
                             complete = F) {
  if (!(level %in% 1:4)) stop("Ecoregion level must be between 1-4.")
  lcol <- paste0("L", level, "_KEY") %>% 
    rlang::parse_expr()
  
  regions <- if (level == 1) {
    asgerbachelor::ecoregions_L1 
    } else if (level == 2) {
      asgerbachelor::ecoregions_L2
    } else if (level == 3) {
      asgerbachelor::ecoregions_L
    } else if (level == 4) {
      asgerbachelor::ecoregions_L4} 
  
  regions <- regions %>% 
    st_transform(crs=st_crs(grid))
  
  if (!missing(simplify)) {
    if (!is.numeric(simplify) | length(simplify) != 1 | simplify < 0) stop("Simplify must be a numeric of length 1 and positive.")
    regions <- st_simplify(regions, dTolerance = simplify)
  }
  
  st_intersection(grid %>% 
                    st_as_sf,
                  regions) %>% 
    as_tibble %>% 
    mutate(
      area = units::drop_units(st_area(geometry)),
      geometry = arrange(as_tibble(grid), ID)[ID,1:2],
      gridArea = units::drop_units(st_area(grid)$area[ID])
    ) %>% 
    group_by(ID, {{lcol}}, geometry) %>% 
    summarize(proportion = sum(area)/gridArea,
              .groups = "drop") %>% 
    mutate(across({{lcol}}, factor)) %>% 
    {
      if (complete) group_by(., ID, geometry) %>% 
        complete({{lcol}}, fill = list(proportion = 0)) %>% 
        ungroup else .
    }
}