% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/provoc.R
\name{remove_identical_lineages}
\alias{remove_identical_lineages}
\title{Remove Identical Lineages}
\usage{
remove_identical_lineages(fused_df, annihilate)
}
\arguments{
\item{fused_df}{The fused data frame from the \code{fuse()} function}

\item{annihilate}{if TRUE will remove all duplicate lineages, if FALSE will warn user if there is duplicate lineages}
}
\value{
A data frame with no mutations that have a duplicate combination of lineages
}
\description{
To ensure the predictor matrix is singular, if annihilate is TRUE the function returns all the mutations with a unique combination of lineages. If FALSE it warns user of duplicate lineages
}
\examples{
# This function is internally used and not typically called by the user.
}
