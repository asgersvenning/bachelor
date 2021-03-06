% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createAggregationGrid.R
\name{createAggregationGrid}
\alias{createAggregationGrid}
\title{Create a regular aggregation grid}
\usage{
createAggregationGrid(
  object = asgerbachelor::ecoregions_L0,
  delta = 50000,
  crs = st_crs(2163)
)
}
\arguments{
\item{object}{a sf object with a bounding box which can be obtained via sf::st_bbox.}

\item{delta}{numeric. Width of cells in the grid (in the units of the supplied crs).}

\item{crs}{an object of class crs. Desired coordinate reference system (CRS) of the output grid.}
}
\value{
An object of class stars with 2 dimensions and 1 attribute, ID, which is unique for each grid cell.
}
\description{
\code{createAggregationGrid} Creates a regular stars grid, which can be used for aggregating an sf polygon sfc.
}
\details{
The created grid is centered on (0,0)
}
\section{Warnings}{

It is important that the supplied CRS is uses a equal area projection, to avoid spatial biases in the size of the aggregation cells.
}

