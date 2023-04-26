#' Calculates functional richness
#' 
#' Calculates functional richness from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param w an optional numeric vector. Trait weights for use in gower dissimilarity computation (if 'gower' = T).
#' @param a optional numeric vector. Species-abundances.
#' @param relative a logical. Calculate functional richness as a proportion of maximum value.
#' @param ndim an optional integer. If supplied the trait space will be ordinated using principal coordinate analysis (PCoA) to a lower dimensional space equal to the value of this parameter. This can speed up the computational time significantly, and is almost mandatory if the number of traits is > 10.
#' @param gower a logical. If true the trait dissimilarity metric will be gower dissimilarity, otherwise euclidean distance will be used.
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @section Details:
#' When using 'relative=T', the maximum volume is computed as the product of the functional trait axis ranges, and standardization is carried out by dividing the convex hull volume by this hypercube volume.
#'
#' OBS: The computation of the convex hull, is extremely slow for large numbers of species and traits. In that case, consider ordinating the traits to a lower dimensionality using the 'ndim' parameter.
#'
#' @importFrom Rdpack reprompt
#' @importFrom geometry convhulln
#' @importFrom stats dist
#' @importFrom stats cmdscale
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
#' functional_richness(ordTraits,tSP[,-1]) %>% 
#'   hist(main = "Distribution of functional richnesses\nat a 50x50 km scale in the FIA dataset",
#'        xlab = "Functional richness")
#' }


functional_richness <- function (x, w, a = rep(1, nrow(x)), relative = F, ndim = NULL, gower = T) {
  if (is.vector(a)) 
    a <- matrix(a, nrow = 1)
  if (inherits(a, "data.frame")) 
    a <- as.matrix(a)
  if (inherits(x, "data.frame")) 
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("Cannot coerce 'x' to matrix.")
  if (!is.matrix(a)) 
    stop("Cannot coerce 'a' to matrix.")
  
  if (ncol(x)>=nrow(x)) stop(paste0("Number of traits must be larger than the number of species, use ndim < ", nrow(x),"!"))
  
  if (!is.null(ndim)) {
    if (!is.numeric(ndim) | length(ndim) != 1 | (ndim %% 1) != 0) stop("'ndim' must be null or a whole number.")
    
    x <- cmdscale({if (gower) gowerDissimilarity(x,w) else dist(x)},ndim)
  }
  
  
  out <- apply(a, 1, function(c) {
    if (all(c == 0)) return(0)
    if (sum(c > 0) <= ncol(x)) return(NA)
    xC <- x[c > 0, ]
    geometry::convhulln(xC, output.options = "FA")$vol/ifelse(relative, 
                                                              prod(abs(diff(apply(xC, 2, range)))), 1)
  })
  return(out)
}