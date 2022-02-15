#' FIA records aggregated to a equal area grid. 
#'
#' Species are encoded by FIA species codes, and geographical bounding boxes in latitude/longititude degrees. Density is species per FIA plot. 
#'
#' @format A tibble with 243 rows and 10 variables:
#' \describe{
#'   \item{INVYR}{dbl Year of inventory}
#'   \item{LON_int}{chr Longitudinal bounding box in decimal degrees.} 
#'   \item{LAT_int}{chr Latitudinal bounding box in decimal degrees.}
#'   \item{SPCD}{dbl FIA species code}
#'   \item{density}{dbl Number of species per plot in the bounding box. Calculated simply as number of individuals divided by number of plots.}
#' }
"FIA"