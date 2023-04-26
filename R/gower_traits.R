#' Transforms a traits matrix to conform to Gower dissimilarity calculations
#' 
#' Creates dummy variables for multi-level factor/character traits, binary variables for bi-level factor/character traits and scales all variables by their range.
#' Furthermore this function supplies weights, which are inversely proportional to the number of dummy variables created for a given trait.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param hasNames a logical. Does the first column contain species names/codes?
#' @param NA.tolerance a number between 0 and 1. Subset species which have a lower proportion of missing traits than 'NA.tolerance'.
#' @return a list which includes the transformed trait table and the weights for each column in the table.
#' 
#' @section Details:
#' This functions performs automatic transformations of a trait-matrix/-table, to "Gower variables". 
#' This means that the variables are appropriate for usage in Gower Dissimilarity calculations, without further transformations.
#' Combined with the weights, which are also automatically calculated, and in conjunction with the function "gowerDissimilarity" (also found in this package), this makes the calculation of diverse Gower Dissimilarities for all other purposes than pairwise dissimilarities, much easier than what was previously available in R. 
#' 
#' @importFrom magrittr %>% 
#' @importFrom magrittr set_colnames
#' @importFrom magrittr set_names
#' @importFrom magrittr inset
#' @importFrom magrittr is_weakly_greater_than
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @importFrom dplyr where
#' @importFrom dplyr bind_cols
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @importFrom tidyr replace_na
#' @importFrom Rdpack reprompt
#' @export


gower_traits <- function(x, hasNames = F, NA.tolerance = .25) {
  # Attempt to coerce the traits to a table (tibble).
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  if (is.null(colnames(x))) colnames(x) <- paste0("var",1:ncol(x))
  if (inherits(x,"data.frame") | inherits(x,"matrix")) x <- as_tibble(x) else stop("'x' must be a matrix or dataframe.")
  
  # Replace "_" in variable names with ".", to ensure that the function is able to separate the variable name from the suffix.
  colnames(x) <- str_replace(colnames(x),"_",".")
  # Temporarily remove species names from the data table, if present.
  xV <- if (hasNames) x[,-1] else x
  
  mSp <- xV %>% 
    is.na %>% 
    as.matrix %>% 
    rowMeans %>%
    is_weakly_greater_than(NA.tolerance) %>% 
    not %>% 
    which

  xT <- xV %>% 
    # Coerce columns to numeric, if the operation is succesful on more than 75% of non-missing values.
    mutate(across(where(~mean(!is.na(as.numeric(.x[!is.na(.x)])))>.25),as.numeric)) %>% 
    # Subset variables which contains more than one non-missing unique value.
    select(where(~!(all(is.na(.x)) | all(.x[!is.na(.x)] == .x[!is.na(.x)][1])))) %>% 
    summarise(across(everything(), ~tibble(value = as.vector(.x)) %>% list)) %>% 
    mutate(across(everything(), ~map(.x, function(x) {
      if (is.logical(x) | all(x[!is.na(x)] %in% c(0:1))) {
        return(
          x %>% 
            mutate(
              "Binary" = as.numeric(value)
            ) %>% 
            select(!value)
        )
      } 
      else if (is.numeric(x$value)) {
        return(x)
      } 
      else {
        uniq <- unique(x$value[!is.na(x$value)])
        nUniq <- length(uniq)
        if (nUniq < 2) {
          return(tibble())
        }
        else if (nUniq == 2) {
          return(
            x %>% 
              mutate(
                "{paste0(uniq, collapse='_or_')}" := as.numeric(factor(value)) - 1
              ) %>% 
              select(!value)
          )
        }
        else {
          return(
            x %>%
              mutate(out = outer(value,uniq,"==") %>%
                       inset(,as.numeric(.)) %>% 
                       set_colnames(uniq) %>% 
                       as_tibble) %>% 
              unnest(out) %>% 
              select(!value)
          )
        }
      }
    })
    )
    ) %>% 
    unnest(everything(),names_sep = "_") %>% 
    set_colnames(colnames(.) %>% str_remove("_value")) %>% 
    mutate(across(everything(), ~(.x-min(.x,na.rm=T))/diff(range(.x,na.rm=T)))) 
  
  # Calculate variable weights
  xW <- names(xT) %>% 
    str_remove("_(?<=_).+") %>% 
    table %>% 
    rep(1/.,.)
  
  # Add species names back, if they are present
  xT <- if (hasNames) bind_cols(x[,1],xT) else xT
  
  return(list(traits = xT[mSp,], weights = xW))
}


