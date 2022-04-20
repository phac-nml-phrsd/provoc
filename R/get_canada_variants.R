#' Download weekly breakdown of SARS-CoV-2 variants in Canada
#' 
#' @param path Where to read the file (or store it if not already downloaded). Relative paths are fine. I suggest using the \code{here} package to get full path names so this function can do better error checking.
#' 
#' @return A data frame with weekly variant numbers
#' @export
get_canada_variants <- function(path = "data/canada_variants.csv"){
    directories <- strsplit(path, "/")[[1]]
    file_dir <- paste(directories[1:(length(directories) - 1)])
    if(!dir.exists(file_dir)) {
        stop(paste0("I cannot find the directory ''", file_dir, "''"))
    }
    if(file.exists(path)) {
        canada_variants <- read.csv(path)
    } else {
        canada_variants <- read.csv("https://health-infobase.canada.ca/src/data/covidLive/covid19-epiSummary-variants.csv")
        names(canada_variants) <- c("Variant_Grouping", "Identifier", "Lineage", "Percent_CT_of_Sample", "Collection_Week")
        name_rms <- "( \\(Gamma\\))|( \\(Beta\\))|( \\(Eta\\))|( \\(Delta\\))|( \\(Alpha\\))"
        canada_variants$Lineage <- gsub(name_rms, "", canada_variants$Lineage)
        write.csv(canada_variants, file = path, row.names = FALSE)
    }
    return(canada_variants)
}
