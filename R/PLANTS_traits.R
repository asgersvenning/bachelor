#' Species traits for species in the FIA from PLANTS 
#'
#' Table containing the traits of all species present in the FIA database, sourced from the PLANTS database. Is mainly used to create species trait matrices, which is very simple; first subset the species codes of interest then subset the trait columns of interest.
#'
#' @format A tibble with 2257 rows and 83 variables:
#' \describe{
#'   \item{plants_code}{chr species code from the PLANTS database - USDA National Plant Data Team.}
#'   \item{...}{species traits.} 
#' }
"PLANTS_traits"