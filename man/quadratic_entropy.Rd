% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quadratic_entropy.R
\name{quadratic_entropy}
\alias{quadratic_entropy}
\title{Calculates Rao's quadratic entropy}
\usage{
quadratic_entropy(
  x,
  w = rep(1, ncol(x)),
  a = rep(1, nrow(x)),
  raoScale = T,
  approxRao = 0,
  gower = T
)
}
\arguments{
\item{x}{numeric matrix. Species-trait matrix.}

\item{a}{optional numeric vector. Species-abundances.}

\item{raoScale}{a logical. Is quadratic entropy divided by maximum value (ade4::divcmax)?}

\item{approxRao}{an integer. If NOT 0, number of species pairs to estimate quadratic entropy on.}

\item{gower}{a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.}
}
\value{
a number.
}
\description{
Calculates Rao's quadratic entropy from a species-trait matrix and possibly an species-abundance matrix/vector.
}
\section{Details}{

This function implements Rao's quadratic entropy, which is defined as the expected dissimilarity between a random pair of individuals from a population.

As such it is quite similar to ade4::divc, however it has 3 major differences;
\enumerate{
\item It allows the user to approximate the quadratic entropy by sampling a specified number of pairs. (Which unfortunately also makes calculating the maximum possible value impractical.)
\item It handles two dissimilarity measures; Euclidean distance and Gower dissimilarity, which allows the flexibility of using the index on non-ordinated traits.
}

It is also considerably faster on large numbers of communities.
}

\examples{
# An example using the data also supplied in the package. Note that the first part of the example, just shows how to create species-abundance and species-trait tables from the data, as well as subsetting species in the intersection of both data sources.
library(hash)
FIA_dict <- hash(PLANTS_meta$plants_code, PLANTS_meta$fia_code)
PLANTS_dict <- invert(FIA_dict)

PLANTS_tG <- PLANTS_traits \%>\%
  gower_traits(T, NA.tolerance = .1)

tSP <- FIA \%>\%
  filter(INVYR>2000) \%>\%
  group_by(ID,SPCD) \%>\%
  summarise(
    dens = mean(individuals/samples, na.rm = T),
    .groups = "drop"
  ) \%>\%
  filter(SPCD \%in\% values(FIA_dict,PLANTS_tG$traits$plants.code)) \%>\%
  pivot_wider(ID,SPCD,values_from=dens,names_prefix = "SP",values_fill = 0)

PLANTS_tG <- PLANTS_traits \%>\%
  filter(plants_code \%in\% values(PLANTS_dict, str_remove(names(tSP)[-1],"^SP"))) \%>\% 
  gower_traits(T)

PCoA_traits <- ape::pcoa(gowerDissimilarity(PLANTS_tG$traits[,-1],PLANTS_tG$weights))

nPC <- PCoA_traits$values[,3:4] \%>\% 
  apply(1,diff) \%>\% 
  sign \%>\% 
  diff \%>\% 
  {
    which(.!=0)[1]
  } \%>\% 
  names \%>\% 
  as.numeric()

ordTraits <- PCoA_traits$vectors[,1:nPC]

par(mfrow = c(2,3))
quadratic_entropy(ordTraits,,tSP[,-1], F, gower = F) \%>\% hist(25,main = "Absolute Euclidean",xlim=c(0,0.05))
quadratic_entropy(ordTraits,,tSP[,-1], F, 10, gower = F) \%>\% hist(25,main = "Approximate Euclidean",xlim=c(0,0.05))
quadratic_entropy(ordTraits,,tSP[,-1], T, gower = F) \%>\% hist(25,main = "Relative Euclidean",xlim=0:1)
quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], F, gower = T) \%>\% hist(25,main = "Absolute Gower",xlim=c(0,0.05))
quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], F,10, gower = T) \%>\% hist(25,main = "Approximate Gower",xlim=c(0,0.05))
quadratic_entropy(PLANTS_tG$traits[,-1],PLANTS_tG$weights,tSP[,-1], T, gower = T) \%>\% hist(25,main = "Relative Gower",xlim=0:1)


# Computation time comparison between algorithms.
ordDist <- dist(ordTraits)

compRes <- lapply((1:5)^4, function(n) {
  wC <- sample(nrow(tSP),n)
  
  tibble(
    communities = n,
    mine = system.time(quadratic_entropy(ordTraits,,tSP[wC,-1],T,gower = F))[3],
    ade4 = system.time(ade4::divc(as.data.frame(t(tSP[wC,-1])),ordDist,T))[3]
  )
}) 

compRes \%>\% 
  tibble \%>\% 
  unnest(1) \%>\%
  pivot_longer(c(mine,ade4),names_to="implementation",values_to="time") \%>\% 
  ggplot(aes(communities,time,color=implementation)) +
  geom_path() +
  geom_path()
}
\references{
\insertRef{Villeger2008}{asgerbachelor}
}
