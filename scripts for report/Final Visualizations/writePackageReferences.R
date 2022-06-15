library(rcrossref)
library(rrapply)
library(magrittr)
library(future.apply)
c("base", "tidyverse","ggplot2","patchwork","magrittr","hash","gamm4","kernlab","sf","stars","terra","mgcv","FD","geometry") %>% sapply(require,character.only=T)

query_crossref <- function(citation) {
  out <- citation %>%
    toBibtex()
  
  newAt <- attributes(out)
  
  new_url <- NULL
  if (!any(newAt$names == "url") & tolower(citation$bibtype) == "article") {
    query <- c(
      "query.bibliographic" = paste0(unlist(str_extract_all(citation$title,"[[:alnum:][ ]]+")),collapse=""),
      "query.author" = paste0(unlist(citation$author),collapse=" "),
      "query.container-title" = citation$journal
    )  
    
    res <- cr_works(limit = 1, flq = query)$data
    
    new_url <- c(res$url,"url")
  }
  
  newAt$names <- c(newAt$names[1:(length(newAt$names)-1)],new_url[2],"key","")
  
  key <- str_extract(out[1],"(?<=[{]).+(?=,$)")
  
  c(paste0("@",out[1] %>% 
             str_extract("(?<=^@).+(?=\\{,$)") %>% 
             tolower %>% 
             ifelse(. %in% c(
               "article",
               "book","booklet",
               "proceedings","inproceedings",
               "manual","techreport",
               "mastersthesis","phdthesis",
               "misc","unpublished"
             ), ., "misc"),"{",key,","),
    out[2:(length(out)-1)],
    if (!is.null(new_url)) paste0("  url = {",new_url[1],"},") else NULL,
    paste0("  key = {",key,"}"),
    out[length(out)]) %>%
    set_attributes(newAt)
}

citation_with_keys <- function(package,key,prefix="R-") {
  if (!missing(key) && length(package) != length(key)) stop("'package' and 'key' must have same length.")
  if (missing(key)) key <- vector("list",length=length(package))
  # plan("multisession")
  out <- map2(package,key, function(p,k) {
    p %>% 
      citation %>% 
      unclass %>% 
      set_attr("bibtype",tolower(attr(.,"bibtype"))) %>% 
      rrapply(function(x) !is.list(x), function(node) sapply(node, function(x) paste0(str_extract_all(str_remove_all(x,"[ ](?=[ ])"),"[[:alnum:]():-[ ]\\./'=]+") %>% unlist,collapse=""))) %>% 
      map2(paste0(prefix,
                  ifelse(is.null(k),p,k),
                  if (length(.) > 1) paste0("-",1:length(.)) else NULL), 
           function(c,fk) {
             c %>% 
               inset2("bibtype",ifelse(tolower(attr(.,"bibtype")) %in% c(
                 "article",
                 "book","booklet",
                 "proceedings","inproceedings",
                 "manual","techreport",
                 "mastersthesis","phdthesis",
                 "misc","unpublished"
               ),tolower(attr(.,"bibtype")),"misc")) %>%
               inset2("key",fk) %>% 
               do.call(bibentry,.)
           }
      ) 
  }) %>% 
    Reduce(c,.) %>% 
    lapply(query_crossref) %>% 
    Reduce(c,.) %>%
    set_class("Bibtex")
  # plan("sequential")
  return(out)
}


c("base", "tidyverse","ggplot2","patchwork","magrittr","hash","gamm4","kernlab","sf","stars","terra","mgcv","FD","geometry") %>%
  citation_with_keys %>%
  c(readLines("../references_from_mendeley.bib")) %>%
  writeLines("../references.bib",useBytes = T)
