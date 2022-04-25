#' Gather information about lineages
#' 
#' Downloads and parses the table at \url{https://cov-lineages.com/lineages.org} using the \code{rvest} package.
#' 
#' @param path Path to save (or load) the lineage information. Default \code{data/}.
#' 
#' @return A data frame, with the side effect of creating a file \code{path/lineage-list.csv} if it does not exist.
#' @export
get_lineage_info <- function(path = "data"){
    if(!dir.exists(path)) {
        stop(paste0("I cannot find the directory ''", path, "''"))
    }
    if (!file.exists(paste0(path, "/lineage_list.csv"))) {
        library(rvest)
        covlin <- rvest::read_html("https://cov-lineages.org/lineage_list.html")
        covlin_tab <- rvest::html_table(covlin)[[1]]
        names(covlin_tab) <- c("Lineage", "Common_Countries", "Earliest", 
            "Designated", "Assigned", "Description", "WHO_name")
        write.csv(covlin_tab, paste0(path, "/lineage_list.csv"))
        return(covlin_tab)
    } else {
        return(read.csv(paste0(path, "/lineage_list.csv")))
    }
}

