#' Data for connections between species codes from FIA and PLANTS
#'
#' A table that contains the FIA and PLANTS species codes, as well as common species names, for the species present in the FIA dataset. Mostly used for connecting FIA species codes with species traits from the PLANTS database. 
#'
#' @format A tibble with 2257 rows and 3 columns:
#' \describe{
#'   \item{fia_code}{dbl species code used in the FIA database - USDA Forest service}
#'   \item{common_name}{chr common species name, US-English} 
#'   \item{plants_code}{chr species code used in the PLANTS database - USDA National Plant Data Team}
#' }
"PLANTS_meta"