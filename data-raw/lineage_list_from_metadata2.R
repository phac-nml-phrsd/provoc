library(lubridate)
library(dplyr)

# wget https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
# gunzip -k metadata.tsv.gz
handle <- file("metadata.tsv", "r")
header <- readLines(handle, n = 1)

# For loops in R are not slow!
# (Unless, of course, you're growing a list, which is exactly what I'm doing.)
ll <- list()
counter <- 0
while(TRUE){
    counter <- counter + 1
    if (!counter %% 5000) print(counter)
    mtpi <- readLines(handle, n = 1)
    if(length(mtpi) == 0) {
        break
    }
    mtpi <- strsplit(mtpi, "\t")[[1]]

    lineage <- mtpi[20]
    mutations <- unlist(strsplit(mtpi[46:48], ","))
    date <- as.character(mtpi[28])
    country <- mtpi[8]


    if(!lineage %in% names(ll)) {
        ll[[lineage]] <- vector(mode = "list", length = length(mutations))
        names(ll[[lineage]]) <- mutations
        for(mut in mutations) {
            ll[[lineage]][[mut]] <- list(
                count = 1, 
                date1 = date, date2 = date, daten2 = date, daten1 = date)
        }
        ll[[lineage]]$seqs <- 1
    } else {
        for(mut in mutations) {
            if(mut %in% names(ll[[lineage]])) {
                old <- ll[[lineage]][[mut]]
                ll[[lineage]][[mut]] <- list(
                    count = old$count + 1,
                    date1 = ifelse(date < old$date1, date, old$date1),
                    date2 = ifelse(date < old$date1, old$date1,
                        ifelse(date < old$date2, date, old$date2)
                    ),
                    daten2 = ifelse(date > old$daten1, old$daten1,
                        ifelse(date > old$daten2, date, old$daten2)
                    ),
                    daten1 = ifelse(date > old$daten1, date, old$daten1))
            } else {
            ll[[lineage]][[mut]] <- list(
                count = 1, 
                date1 = date, date2 = date, daten2 = date, daten1 = date)
            }
        }
        ll[[lineage]]$seqs <- ll[[lineage]]$seqs + 1
    }
}
close(handle)

saveRDS(ll, here("data-raw", "ll.rds"))

ll_coolj <- lapply(1:length(ll), function(i) {
    if(names(ll)[i] != "?"){
        x <- ll[[i]]
        x2 <- bind_rows(sapply(x, as.data.frame))
        x2$lineage <- as.character(names(ll)[i])
        x2$mut <- names(x)
        seq_col <- which(names(x2) == "X[[i]]")
        x2$seqs <- as.numeric(x2[x2$mut == "seqs", seq_col])
        x2 <- x2[x2$mut != "seqs", -seq_col]
        enough <- which(x2$count > quantile(x2$count, 0.1))
        if(length(enough) > 0.33*nrow(x2)) x2 <- x2[enough,]
        x2 <- x2[rev(order(x2$count)),]
        if(nrow(x2) > 500) x2 <- x2[1:500,]
        x2
    }
})

sum(unlist(lapply(ll_coolj, nrow)))

ll2 <- bind_rows(ll_coolj)

mutations_by_lineage <- ll2 

usethis::use_data(mutations_by_lineage, overwrite = TRUE)
