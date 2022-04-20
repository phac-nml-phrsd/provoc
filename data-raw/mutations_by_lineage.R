library(jsonlite)
library(dplyr)

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
