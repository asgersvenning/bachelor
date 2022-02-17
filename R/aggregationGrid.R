#' Aggregation grid used for aggregating both FIA records and EPA ecoregions.
#'
#' Species are encoded by FIA species codes, and coordinates are in EPSG:2163 (for details = "sf::st_crs(2163)$proj4string"). Density is species per FIA plot. 
#'
#' @format Aggregation grid with two dimensions (x = 92, y = 58) and 1 attribute:
#' \describe{
#'   \item{ID}{dbl A unique grid cell ID.}
#' }
"aggregationGrid"