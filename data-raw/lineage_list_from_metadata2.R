library(lubridate)
library(dplyr)
library(here)

# GET the latest metadata
# wget --show-progress -O "data-raw/metadata_$(date +"%Y-%m-%d").tsv.gz" https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
# USE the latest metadata
mfiles <- list.files("data-raw", pattern = "^metadata")
handle <- gzfile(here("data-raw", rev(sort(mfiles))[1]), "r")

header <- readLines(handle, n = 1) 
header <- strsplit(header, split = "\t")[[1]]
pango_lineage <- which(header == "pango_lineage")
mutations <- which(header %in% c("insertions", "deletions", "substitutions"))
date_submitted <- which(header == "date_submitted")
country_col <- which(header == "country")
sra_accession <- which(header == "sra_accession")

# For loops in R are not slow!
# (Unless, of course, you're growing a list, which is exactly what I'm doing.)
ll <- list()
seq_ids <- c()
if(file.exists(here("data-raw", "ll.rds"))) {
    ll <- readRDS(here("data-raw", "ll.rds"))
    seq_ids <- readRDS(here("data-raw", "seq_ids.rds"))
}
counter <- 0
t0 <- Sys.time()
while(TRUE){
    mtp <- readLines(handle, n = 1000)
    for(i in 1:1000){
        counter <- counter + 1
        mtpi <- mtp[i]
        if(length(mtpi) == 0) {
            break
        }
        mtpi <- strsplit(mtpi, "\t")[[1]]

        this_seq <- mtpi[sra_accession]
        if(this_seq %in% seq_ids) {
            next
        }
        if (!counter %% 5000) print(counter)

        seq_ids <- c(seq_ids, this_seq)

        lineage <- mtpi[pango_lineage]
        mutations <- unlist(strsplit(mtpi[mutations], ","))
        date <- as.character(mtpi[date_submitted])
        country <- mtpi[country_col]

        if(!lineage %in% names(ll)) {
            ll[[lineage]] <- vector(mode = "list", length = length(mutations))
            names(ll[[lineage]]) <- mutations
            for(mut in mutations) {
                ll[[lineage]][[mut]] <- list(
                    count = 1, 
                    date1 = date, daten1 = date)
            }
            ll[[lineage]]$seqs <- 1
        } else {
            for(mut in mutations) {
                if(mut %in% names(ll[[lineage]])) {
                    old <- ll[[lineage]][[mut]]
                    ll[[lineage]][[mut]] <- list(
                        count = old$count + 1,
                        date1 = ifelse(date < old$date1, date, old$date1),
                        daten1 = ifelse(date > old$daten1, date, old$daten1))
                } else {
                ll[[lineage]][[mut]] <- list(
                    count = 1, 
                    date1 = date, daten1 = date)
                }
            }
            ll[[lineage]]$seqs <- ll[[lineage]]$seqs + 1
        }
    }
}
close(handle)
Sys.time() - t0

saveRDS(ll, here("data-raw", "ll.rds"))
saveRDS(seq_ids, here("data-raw", "seq_ids.rds"))

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
