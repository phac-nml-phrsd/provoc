% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuse.R
\name{lineage_similarity}
\alias{lineage_similarity}
\title{Finds and prints all similarities among lineage}
\usage{
lineage_similarity(data, simplify = FALSE, almost = 1)
}
\arguments{
\item{data}{A lineage matrix}
}
\value{
A list of length 4 containing information on which lineages differ by one,
the Jaccard similarity between lineages, which lineages are subsets and almost subsets
of each other. in is_subset and is_almost_subset a value is true if the lineages of the
column name is a subset/almost a subset of the lineages of the row name.
}
\description{
Finds and prints all similarities among lineage
}
