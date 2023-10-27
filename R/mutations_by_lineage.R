#' Lists of mutations by lineage
#' 
#' Top 700 most common mutations for each sequence with lineage information in NextStrain GenBank database as of May 05 2022. Useful for creating a variant matrix (for publications, it is highly recommended that you fully understand how and why you are constructing a variant matrix).
#' 
#' The mutation lists are imperfect, and it is a good idea to further subset the mutation lists prior to analysis. For instance, by using only lineages with more than 1000 sequences in the database (the largest count for a lineage is a proxy measure for this).
#' 
#' A python script (written by Dr. Art Poon and his lab) was used to apply \href{https://github.com/lh3/minimap2}{minimap2} to all data in the \href{https://nextstrain.org/sars-cov-2/#datasets}{NextStrain GenBank} collection of sequences. This python script is in the \code{data-raw} folder in the source version of this repo on GitHub. 
#' TODO: document the usage of retrieve-nsgb.py.
#' 
#' The sequence database can be manually downloaded by clicking this link: \url{https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz}
#' 
#' The metadata (including assigned lineage) can be downloaded by clicking this link: \url{https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz}
#' 
#' @format A data frame with four columns.
#' \describe{
#'      \item{lineage}{The designated lineage, as labelled in GenBank. Note that these are imperfect and sequence lineage assignment has changed over time.}
#'      \item{mutation}{The name of the mutation (relative to \code{NC_045512}), in the format "type|position|alt".}
#'      \item{count}{As of March 2022, the number of GenBank records labelled with this lineage that contained this mutation.}
#' }
#' 
#' @source \url{https://nextstrain.org}
"mutations_by_lineage"


#' Facts about the lineages
#' 
#' Especially useful for checking the earliest/latest sequence date
#' 
#' @format A data frame
#' 
#' @source \url{https://nextstrain.org}
"lineage_facts"

#' Lists of mutations per lineage, as per Constellations
#' 
#' The names of mutations according to Pangolin's "Constellations" files. These are current as of October 2024. The constellation files for pre-Omicron variants are not expected to change.
#' 
#' @format A data frame
#' 
#' @source \url{https://github.com/cov-lineages/constellations}
"constellation_lists"
