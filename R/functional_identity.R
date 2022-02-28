#' Calculates functional identity
#' 
#' Calculates functional identity from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#'  
#' @section Details:
#' Names of functional traits are derived from x if present, otherwise they are formatted as "X"+"column number" to be aligned with FD::dbFD. 
#'  
#' @importFrom magrittr %>% 
#' @importFrom tibble as_tibble
#' @importFrom Rdpack reprompt
#' @export

functional_identity <- function(x, a = rep(1, nrow(x))) {
  a <- a/sum(a) # Relative abundances
  
  out <- (t(x) %*% a) %>% as.vector %>% t %>% as_tibble # Abundance weighted column means
  
  if (!is.null(colnames(x))) { # Set the output names based on the input trait names
    colnames(out) <- colnames(x)
  } else {
    colnames(out) <- paste0("X",1:ncol(out)) # Or standardize with FD::dbFD output
  } 
  
  return(out)
}