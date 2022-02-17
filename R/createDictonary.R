createDictionary <- function(meta, code, value) {
  codes <- meta[[code]]
  values <- meta[[value]]
  
  if (!is.numeric(codes) | any(codes<1)) stop("Species codes must be numeric and positive to create a dictionary.")
  
  dict <- rep(NA, max(codes))
  dict[codes] <- values
  
  return(dict)
}
