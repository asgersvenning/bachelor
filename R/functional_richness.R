#' Calculates functional richness
#' 
#' Calculates functional richness from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @param relative a logical. Calculate functional richness as a proportion of maximum value.
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @section Details:
#' When using 'relative=T', the maximum volume is computed as the product of the functional trait axis ranges, and standardization is carried out by dividing the convex hull volume by this hypercube volume.
#'
#' OBS: The computation of the convex hull, is extremely slow for large numbers of species and traits. In that case, consider ordinating the traits to a lower dimensionality.
#'
#' @importFrom Rdpack reprompt
#' @importFrom geometry convhulln
#' @export
#' 
#' @examples 
#' # An example using the data also supplied in the package. Note that the first part of the example, just shows how to create species-abundance and species-trait tables from the data, as well as subsetting species in the intersection of both data sources.
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
#' PCoA_traits <- ape::pcoa(FD::gowdis(PLANTS_tG$traits[,-1],PLANTS_tG$weights))
#' 
#' nPC <- PCoA_traits$values[,3:4] %>% 
#'   apply(1,diff) %>% 
#'   sign %>% 
#'   diff %>% 
#'   {
#'     which(.!=0)[1]
#'   } %>% 
#'   names %>% 
#'   as.numeric()
#' 
#' ordTraits <- PCoA_traits$vectors[,1:nPC]
#' 
#' par(mfrow = c(1,1))
#' functional_richness(ordTraits,tSP[,-1]) %>% hist(main = "Distribution of functional richnesses\nat a 50x50 km scale in the FIA dataset",
#'                                                  xlab = "Functional richness")


functional_richness <- function(x, a = rep(1, nrow(x)), relative = F) {
  # Attempt to coerce x and a to matrices.
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  if (inherits(a,"data.frame")) a <- as.matrix(a)
  if (inherits(x,"data.frame")) x <- as.matrix(x)
  
  # Return resonable error messages.
  if (!is.matrix(x)) stop("Cannot coerce 'x' to matrix.")
  if (!is.matrix(a)) stop("Cannot coerce 'a' to matrix.")
  
  out <- apply(a, 1, function(c) {
    # If a community doesn't have any species, set functional richness to 0.
    if (all(c == 0)) return(0)

    # The number of species must be greater than the number of traits, to compute the convex hull.
    if (!(sum(c > 0) > ncol(x))) return(NA)
    
    # Subset present species
    xC <- x[c>0,]
    
    # Volume of convex hull
    geometry::convhulln(xC,output.options = "FA")$vol /
      # Possibly standardized by the volume of the minimum encapsulating (hyper-)cube
      ifelse(relative, prod(abs(diff(apply(xC,2,range)))),1) 
  })
  
  return(out)
}