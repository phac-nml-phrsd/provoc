handle <- file("metadata.tsv", "r")

# For loops in R are not slow!
# (Unless you're growing a list, which is what I'm doing.)
lineage_list <- list()
while(TRUE){
    mtpi <- readLines(handle, n = 1)
    if(length(mtpi) == 0) {
        break
    }
    mtpi <- strsplit(mtpi, "\t")[[1]]

    lineage <- mtpi[20]
    mutations <- unlist(strsplit(mtpi[46:48], ","))
    date <- mtpi[28]
    country <- mtpi[8]


    if(!lineage %in% names(lineage_list)) {
        lineage_list[[lineage]] <- list(
            mutations = mutations,
            date = date,
            countries = data.frame(country = country, count = 1)
        )
    } else {
        old_list <- lineage_list[[lineage]]
        old_list$mutations <- union(old_list$mutations, mutations)
        old_list$date <- c(old_list$date, date)
        if(country %in% old_list$countries$country) {
            old_list$countries$count[which(country == old_list$countries$country)] <- old_list$countries$count[which(country == old_list$countries$country)] + 1
        } else {
            old_list$countries <- rbind(old_list$countries, data.frame(country = country, count = 1))
        }
        lineage_list[[lineage]] <- old_list
    }
}

sapply(lapply(lineage_list, function(x) x$date)[order(names(lineage_list))], length)

close(handle)

length(lineage_list)

saveRDS(lineage_list, file = "lineage_list.rds")