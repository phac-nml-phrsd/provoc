% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/provoc.R
\name{constant_with_by}
\alias{constant_with_by}
\title{Finds all columns of the data that are constant with the specified by group}
\usage{
constant_with_by(data, by)
}
\arguments{
\item{data}{Data frame containing count, coverage, and lineage columns.}

\item{by}{Column name to group and process data.}
}
\value{
A dataframe with the columns by and all that are constant with by
}
\description{
Finds all columns of the data that are constant with the specified by group
}
