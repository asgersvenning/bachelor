% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLANTS_meta.R
\docType{data}
\name{PLANTS_meta}
\alias{PLANTS_meta}
\title{Data for connections between species codes from FIA and PLANTS}
\format{
A tibble with 2257 rows and 3 columns:
\describe{
\item{fia_code}{dbl species code used in the FIA database - USDA Forest service}
\item{common_name}{chr common species name, US-English}
\item{plants_code}{chr species code used in the PLANTS database - USDA National Plant Data Team}
}
}
\usage{
PLANTS_meta
}
\description{
A table that contains the FIA and PLANTS species codes, as well as common species names, for the species present in the FIA dataset. Mostly used for connecting FIA species codes with species traits from the PLANTS database.
}
\keyword{datasets}
