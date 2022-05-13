library(jsonlite)
library(dplyr)
library(lubridate)

# Requires the output of retrieve-nsgb.py, which downloads and parses nextstrain genbank data.
# retrieve-nsgb.py requires a Unix-like CLI and requires the correct command line arguments to work.
# This will be better documented in a future version of this package (possibly in a bonus vignette; it will NOT be added as package functionality because it requires downloading and processing a file that is currently dozens of gigabytes and growing rapidly). 

mutations <- read_json("../CoCoVoC/data/nsgb-2022-03-11.json")

for (i in seq_along(mutations)) { 
    d1 <- data.frame(
        mutation = names(mutations[[i]]$mutation), 
        count = as.numeric(mutations[[i]]$mutation))
    d1$frequency <- d1$count / mutations[[i]]$count
    d1 <- d1[order(d1$count),]
    d1 <- tail(d1, 700)
    d1$lineage <- names(mutations)[i]

    if (i == 1) {
        mutations_by_lineage <- d1
    } else {
        mutations_by_lineage <- bind_rows(mutations_by_lineage, d1)
    }
}

low_count <- sapply(unique(mutations_by_lineage$lineage),
    function(x) max(mutations_by_lineage$count[mutations_by_lineage$lineage == x]))
low_count <- low_count[low_count <= 100]

mutations_by_lineage <- mutations_by_lineage[!mutations_by_lineage$lineage %in% names(low_count),]

format(object.size(mutations_by_lineage), units = "MB")
head(mutations_by_lineage)

mutations_by_lineage <- mutations_by_lineage[, c("lineage", "mutation", "count")]
usethis::use_data(mutations_by_lineage, overwrite = TRUE)



# Gather metadata ---------------------------------------------------
# wget https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
# gunzip -k metadata.tsv.gz

# WARNING: 4.7GB file. Not for the faint of memory
mt <- read.csv("metadata.tsv", sep = "\t")
names(mt)

facts <- mt %>% 
    filter(pango_lineage != "?")) %>%
    group_by(pango_lineage) %>%
    summarise(count = n(),
        earliest_date = min(ymd(date), na.rm = TRUE),
        latest_date = max(ymd(date), na.rm = TRUE),
        )

lins <- unique(mt$pango_lineage)
lins <- lins[!lins %in% c("?", "", "unclassifiable", "None")]

lineage_facts <- sapply(lins, function(lin) {
    sublins <- mt[mt$pango_lineage == lin,]

    country_tab <- sort(table(sublins$country[!is.na(sublins$country)]), decreasing = TRUE)
    division_tab <- sort(table(sublins$division[!is.na(sublins$division)]), decreasing = TRUE)
    dates <- ymd(sublins$date)


    c(
        pango_lineage = lin, 
        count = nrow(sublins),
        Canada_count = as.numeric(country_tab["Canada"]),
        earliest_date = as.character(min(dates, na.rm = TRUE)),
        median_date = as.character(median(dates, na.rm = TRUE)),
        latest_date = as.character(max(dates, na.rm = TRUE)),
        top_10_countries = paste(names(country_tab)[1:10], collapse = ";"),
        top_10_divisions = paste(names(division_tab)[1:10], collapse = ";")
    )
})

lineage_facts <- as.data.frame(t(lineage_facts))
rownames(lineage_facts) <- NULL

head(lineage_facts)
format(object.size(lineage_facts), units = "MB")

usethis::use_data(lineage_facts, overwrite = TRUE)
