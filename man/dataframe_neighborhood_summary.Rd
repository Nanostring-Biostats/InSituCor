% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeighborhoodCalculations.R
\name{dataframe_neighborhood_summary}
\alias{dataframe_neighborhood_summary}
\title{Summarize a data frame's values over a neighborhood}
\usage{
dataframe_neighborhood_summary(df, neighbors)
}
\arguments{
\item{df}{A data frame, with rows aligned to the rows of neighbors. Can contain numeric and/or character columns.}

\item{neighbors}{A neighbors adjacency matrix}
}
\value{
A list with elements corresponding to data frame variables. Entries are either numeric vectors
summarizing numeric variablies in conditionon or count matrices containing the tabulated values from character vectors.
}
\description{
Numeric and character columns are respectively summed or tabulated across each cell's neighbors.
}
