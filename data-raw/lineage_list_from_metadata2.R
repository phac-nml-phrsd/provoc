library(lubridate)
library(dplyr)
library(here)
library(parallel)

ncores <- detectCores() - 2
cl <- makeCluster(ncores)

# GET the latest metadata
# wget --show-progress -O "data-raw/metadata_$(date +"%Y-%m-%d").tsv.gz" https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
# USE the latest metadata
mfiles <- list.files("data-raw", pattern = "^metadata")
handle <- gzfile(here("data-raw", rev(sort(mfiles))[1]), "r")

header <- readLines(handle, n = 1) 
header <- strsplit(header, split = "\t")[[1]]
pango_lineage <- which(header == "pango_lineage")
mutation_cols <- which(header %in% c("insertions", "deletions", "substitutions"))
date_submitted <- which(header == "date_submitted")
country_col <- which(header == "country")
sra_accession <- which(header == "sra_accession")
clusterExport(cl, c("pango_lineage", "mutation_cols", "date_submitted"))

kill_switch <- FALSE
ll <- data.frame()
counter <- 0
N <- 500000
t1 <- Sys.time()
while(!kill_switch) {
    t0 <- Sys.time()
    mtp <- strsplit(readLines(handle, n = N), split = "\t")
    if(length(mtp) == 0) kill_switch <- TRUE
    ml <- bind_rows(parLapply(cl, mtp, function(x) {
        if(length(x) > 0) {
            mutations <- unlist(strsplit(x[mutation_cols], ","))
            if(length(mutations) > 0) {
                data.frame(
                    lineage = x[pango_lineage],
                    mutation = mutations,
                    date1 = as.character(x[date_submitted]),
                    daten = as.character(x[date_submitted]),
                    count = 1,
                    stringsAsFactors = FALSE
                )
            }
        }}))

    ll <- summarise(group_by(bind_rows(ll, ml), lineage, mutation),
        count = sum(count), date1 = min(date1), daten = max(date1),
        .groups = "drop")

    saveRDS(ll, here("data-raw", "ll-part.rds"))
    counter <- counter + N
    print(paste0(counter, ", ", difftime(Sys.time(), t0, units = "secs")/N, " seconds per line."))
}
print(difftime(Sys.time(), t1, units = "mins"))

close(handle)
saveRDS(ll, here("data-raw", "ll.rds"))

mutations_by_lineage <- ll %>%
    filter(!lineage %in% c("?", "Unassigned", "unclassifiable")) %>%
    group_by(lineage, mutation) %>%
    mutate(total_count = max(count)) %>%
    filter(count > 0.2*total_count) %>%
    ungroup %>%
    select(-total_count) %>%
    group_by(lineage) %>%
    arrange(count) %>%
    slice(1:700) %>%
    ungroup()

#usethis::use_data(mutations_by_lineage, overwrite = TRUE)
