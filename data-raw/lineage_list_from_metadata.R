handle <- file("metadata.tsv", "r")
header <- readLines(handle, n = 1)

# For loops in R are not slow!
# (Unless, of course, you're growing a list, which is exactly what I'm doing.)
lineage_list <- list()
while(TRUE){
    mtpi <- readLines(handle, n = 1)
    if(length(mtpi) == 0) {
        break
    }
    mtpi <- strsplit(mtpi, "\t")[[1]]

    lineage <- mtpi[20]
    mutations <- paste(unlist(strsplit(mtpi[46:48], ",")), collapse = ",")
    date <- mtpi[28]
    country <- mtpi[8]


    if(!lineage %in% names(lineage_list)) {
        lineage_list[[lineage]] <- list(
            mutations = mutations,
            dates = date,
            countries = country
        )
    } else {
        lineage_list[[lineage]]$mutations <- c(lineage_list[[lineage]]$mutations, mutations)
        lineage_list[[lineage]]$dates <- c(lineage_list[[lineage]]$dates, date)
        lineage_list[[lineage]]$countries <- c(lineage_list[[lineage]]$countries, country)
    }
}
close(handle)

#deprecated
#saveRDS(lineage_list, file = here("data-raw", "lineage_list.rds"))
