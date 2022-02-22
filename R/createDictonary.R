#' Create a dictionary 
#' 
#' \code{createDictionary} creates a dictionary for a set of values given a set of unique integer codes.
#' 
#' @param meta a data frame or named list. An object containing the codes and values.
#' @param code character or integer. Name or index of the location of codes in meta.
#' @param value character or integer. Name or index of the location of values in meta.
#' 
#' @return A vector where the value at the index position equals to code, matches the value corresponding to code.
#' 
#' @section Warning:
#' Given that the function creates the dictionary by hashing the values at the literal position of code, the length of the output vector will be equal to the maximum code value.
#'
#' @export

createDictionary <- function(meta, code, value) {
  codes <- meta[[code]]
  values <- meta[[value]]
  if (length(codes) != length(values)) stop("code must match the length of value.")
  if (length(codes) != length(unique(codes))) stop("code must contain unique values.")
  if (!is.numeric(codes) | any(codes<1)) stop("Species codes must be numeric and positive to create a dictionary.")
  
  dict <- rep(NA, max(codes))
  dict[codes] <- values
  
  return(dict)
}
