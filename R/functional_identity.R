#' Calculates functional identity
#' 
#' Calculates functional identity from a species-trait matrix and possibly an abundance vector.
#' Functional identity is simply a fancy way of saying (community-weighted) mean traits.
#' 
#' @param x numeric matrix. Species-trait matrix. 
#' @param a optional numeric vector. Species-abundances (matrix), can be missing.
#' @return A table of community functional identity (community weighted mean traits). In case 'a' is missing, simply the mean a table with one row, containing the mean traits.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#'  
#' @section Details:
#' Names of functional traits are derived from x if present, otherwise they are formatted as "X"+"column number" to be aligned with FD::dbFD. 
#'  
#' In case of some of or all of the traits are factors, the trait levels must be converted to dummy columns. This can be done using the 'gower_traits' function in this package.
#'  
#' @importFrom magrittr %>% 
#' @importFrom tibble as_tibble
#' @importFrom Rdpack reprompt
#' @export
#' 
#' @examples 
#' \dontrun{
#' # An example using the data also supplied in the package. 
#' # Note that the first part of the example, 
#' # just shows how to create species-abundance and species-trait tables from the data, 
#' # as well as subsetting species in the intersection of both data sources.
#' library(hash)
#' FIA_dict <- hash(PLANTS_meta$plants_code, PLANTS_meta$fia_code)
#' PLANTS_dict <- invert(FIA_dict)
#' 
#' tSP <- FIA %>% 
#'   filter(INVYR>2000) %>% 
#'   group_by(ID,SPCD) %>% 
#'   summarise(
#'     dens = mean(individuals/samples, na.rm = T),
#'     .groups = "drop"
#'   ) %>% 
#'   filter(SPCD %in% values(FIA_dict,PLANTS_tG$traits$plants.code)) %>% 
#'   pivot_wider(ID,SPCD,values_from=dens,names_prefix = "SP",values_fill = 0)
#' 
#' PLANTS_tG <- PLANTS_traits %>% 
#'   filter(plants_code %in% values(PLANTS_dict, str_remove(names(tSP)[-1],"^SP"))) %>% 
#'   gower_traits(T)
#' 
#' tFI <- functional_identity(PLANTS_tG$traits[,-1],rep(1,372))
#' }

functional_identity <- function(x, a = rep(1, nrow(x))) {
  # Attempt to coerce x and a to matrices.
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  if (inherits(x,"data.frame")) x <- as.matrix(x)
  if (inherits(a,"data.frame")) a <- as.matrix(a)
  
  if (!is.matrix(a)) stop("Unable to coerce 'a' to matrix.")
  if (!is.matrix(x)) stop("Unable to coerce 'x' to matrix.")
  
  # Check if the number of species x and a match.
  if (ncol(a) != nrow(x)) stop("Incompatible number of species do not match! Number of columns in 'a' must match number of rows in 'x'.")
  
  # Set missing abundances to 0, amounts to removing them from the means.
  a[is.na(a)] <- 0
  
  # Abundance weighted column means
  out <- ((t(replace(x,is.na(x),0)) %*% t(a)) / (t(!is.na(x)) %*% t(a))) %>% t 
  
  # Set infinite means, which arise when the sum of abundances for non-missing species for a trait is 0, to NA. 
  out[is.infinite(out)] <- NA
  
  if (is.null(colnames(x))) { # Set the output names based on the input trait names
    colnames(out) <- paste0("X",1:ncol(out)) # Or standardize with FD::dbFD output
  } 
  
  return(as_tibble(out))
}



