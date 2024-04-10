#' Create a lineage definition matrix from a list of lineages with their mutations
#'
#' The input should be a named list where the names represent lineages and the values are vectors of mutations. The output will be a valid lineage matrix.
#'
#' @param lineage_list A named list of vectors of mutations. Mutation names should be in the same format as those in your data, otherwise post-processing will be required.
#'
#' @return A lineage matrix (rownames are lineages, colnames are mutations, entry i,j is 1 if lineage i contains mutation j, 0 otherwise).
#' @export
lineage_defs_from_list <- function(lineage_list) {
    if (!is.list(lineage_list)) {
        stop("Input must be a list")
    } else if (is.null(names(lineage_list))) {
        stop("Input must be a named list.")
    } else if (any(sapply(lineage_list, function(x) class(x) != "character"))) {
        stop("Input must be a named list of character vectors")
    }

    all_mutations <- unique(unlist(lineage_list))
    lineage_defs <- as.data.frame(t(1 * sapply(lineage_list,
                function(x) all_mutations %in% x)))
    names(lineage_defs) <- all_mutations
    lineage_defs
}


#' Given a vector of lineage names, uses GenBank data to determine the lineage matrix.
#'
#' Requires one of a few potentially arbitrary choices: Use the \code{max_n} most common mutations for a given lineage or use the mutations in the top \code{top_quantile} quantile.
#'
#' See \code{?mutations_by_lineage} for more information about the GenBank data used and see \code{unique(mutations_by_lineage$lineage)} to see available lineage names.
#'
#' @param lineages A character vector of lineage names in the same format as the
#' @param max_n The maximum number of mutations from a given lineage to retain.
#' @param top_quantile Only take mutations that are in the top \code{top_quantile} quantile. E.g. 0.05 gives the mutations in more than 95% of the sequences of that lineage.
#' @param mutation_format "tpa" (default) for \code{type|pos|alt} (e.g. \code{~|2832|G}; note that it includes the pipes "|") or "aa" for amino acid (e.g. \code{aa:orf1a:K856R}).
#'
#' @return A lineage matrix (rownames are lineage, colnames are mutations, entry i,j is 1 if lineage i contains mutation j, 0 otherwise).
#' @export
lineage_defs_from_lineage <- function(lineages,
    max_n = NULL, top_quantile = NULL, mutation_format = c("tpa", "aa"),
    mutations = provoc::mutations_by_lineage) {
    if(is.null(max_n) & is.null(top_quantile)) {
        warning("Parameters not set, using all mutations (lineage_defs will be massive and estimation will be slow and fraught with correlations)")
        top_quantile <- 0
        max_n <- 700
    } else if (is.null(top_quantile)) {
        top_quantile <- 0
    } else if (is.null(max_n)) {
        max_n <- 700
    }

    lineage_list <- lapply(lineages, function(x) {
        m <- mutations[mutations$lineage == x,]
        m <- m[order(-m$count),]
        m$top_n <- FALSE
        m$top_n[1:min(nrow(m), max_n)] <- TRUE
        m$top_quantile <- m$count < quantile(m$count, 1 - top_quantile)
        m$mutation[m$top_n & m$top_quantile]
    })
    names(lineage_list) <- lineages

    lineage_def <- lineage_defs_from_list(lineage_list)
    if(mutation_format[1] == "aa"){
        for(i in 1:ncol(lineage_defs)){
            cn <- strsplit(colnames(lineage_defs
)[i], split = "\\|")[[1]]
            colnames(lineage_defs
)[i] <- parse_mutation(cn[1], as.numeric(cn[2]) + 1, cn[3])
        }
    }
    lineage_def
}


#' Use mutations found in the data to determine the lineage matrix
#'
#' Uses the included GenBank mutation list (\code{?mutations_by_lineage}) to determine a lineage matrix that is guaranteed to have the mutations from your data (assuming you provide the correct mutation name format).
#'
#' @param type A vector containing either +, -, or ~, representing insertions, deletions, and mutations, respectively. Can also be i, d, or m.
#' @param pos A vector containing the 1-indexed positions.
#' @param alt A vector containing the alternative (nucleotide(s) for mutations and insertions or the number of deletions).
#' @param aa The vector of mutations in amino acid format (like those create by \code{parse_mutation}). TODO: Not yet implemented.
#' @param data A data frame containing columns labelled "type", "pos", and "alt", each with entries in the expected format. TODO: Not yet implemented.
#' @param mutation_format "tpa" (default) for \code{type|pos|alt} (e.g. \code{~|2832|G}; note that it includes the pipes "|") or "aa" for amino acid (e.g. \code{aa:orf1a:K856R}).
#' @param max_n The maximum number of mutations from a given lineage to retain.
#' @param top_quantile Only take mutations that are in the top \code{top_quantile} quantile. E.g. 0.95 gives the mutations that are observed in more than 95% of the sequences of that lineage.
#' @param matches The minimum number of matching mutations in the NextStrain GenBank data.
#' @param lineage_names A character vector of lineage names. Must match the names in \code{mutations_by_lineage}.
#' @param start_date Earliest date in the study. Must be in ISO-8601 format (as should all dates with no exceptions).
#' @param check_after Check if there's an observed sequence after the start of your study? Set FALSE if the start date is recent.
#' @param check_canada Checks if the lineage was ever observed in Canada. Default FALSE.
#'
#' @return A lineage matrix (rownames are lineage, colnames are mutations, entry i,j is 1 if lineage i contains mutation j, 0 otherwise).
#' @export
lineage_defs_from_data <- function(type = NULL, pos = NULL, alt = NULL,
    aa = NULL, data = NULL,
    mutation_format = c("tpa", "aa"),
    max_n = 80, top_quantile = 0.05,
    matches = 3,
    start_date = NULL, check_after = TRUE, check_canada = FALSE,
    mutations = provoc::mutations_by_lineage) {

    # TODO: Error checking (same size, not all null, what if type and aa are both specified, etc.)
    if(is.null(max_n) & is.null(top_quantile)) {
        warning("Parameters not set, using all mutations")
        top_quantile <- 1
        max_n <- 700
    } else if (is.null(top_quantile)) {
        top_quantile <- 1
    } else if (is.null(max_n)) {
        max_n <- 700
    }

    mutations2 <- mutations
    if("m" %in% type) {
        mutations2$mutation <- gsub("~", "m", mutations2$mutation)
        mutations2$mutation <- gsub("\\+", "i", mutations2$mutation)
        mutations2$mutation <- gsub("-", "d", mutations2$mutation)
    }

    mutations <- unique(paste(type, pos, alt, sep = "|"))

    all_lineages <- unique(mutations$lineage)
    lineages <- lapply(all_lineages, function(x) {
        sum(mutations %in% mutations$mutation[mutations$lineage == x])
    })

    lineages <- all_lineages[lineages >= matches]

    if (!is.null(start_date)) {
        lineages <- extant_lineages(lineage_names = lineages,
            start_date = start_date,
            check_after = check_after, check_canada = check_canada)
    }

    lineage_defs_from_lineage(lineages = lineages,
        max_n = max_n, top_quantile = top_quantile,
        mutation_format = mutation_format)
}

#' Filter lineages active on a given date.
#'
#' Using NGSB data (see \code{?mutations_by_lineage}), checks that the earliest sequence in a lineage was observed by \code{start_date}. Optionally checks if the latest sequence in the lineage was observed after start_date. Optionally checks if the lineage was ever observed in Canada.
#'
#' @param lineage_names A character vector of lineage names. Must match the names in \code{mutations_by_lineage}.
#' @param start_date Earliest date in the study. Must be in ISO-8601 format (as should all dates with no exceptions).
#' @param check_after Check if there's an observed sequence after the start of your study? Set FALSE if the start date is recent.
#' @param check_canada Checks if the lineage was ever observed in Canada. Default FALSE.
#'
#' @details The code also ignores any + symbols and anything after them, so lineages such as B.1.617.2+K417N (Delta+) will be treated as B.1.617.2 (Delta).
#'
#' @return A character vector.
extant_lineages <- function(lineage_names, start_date, check_after = TRUE, check_canada = FALSE) {
    start_date <- lubridate::ymd(start_date)
    lineage_names2 <- sapply(strsplit(lineage_names, "\\+"), `[`, 1)
    include_lineage <- logical(length(lineage_names2))

    for(i in seq_along(include_lineage)) {
        if(lineage_names2[i] %in% lineage_facts$pango_lineage) {
            lin_date <- lubridate::ymd(lineage_facts$earliest_date[lineage_facts$pango_lineage == lineage_names2[i]])
            include_lineage[i] <- lin_date <= start_date
        }

        if(check_after & include_lineage[i]) {
            late_date <- lubridate::ymd(lineage_facts$latest_date[lineage_facts$pango_lineage == lineage_names2[i]])
            include_lineage[i] <- include_lineage[i] & (late_date >= start_date)
        }
        if(check_canada & include_lineage[i]) {
            include_lineage[i] <- include_lineage[i] &
                !is.na(lineage_facts$Canada_count[lineage_facts$pango_lineage == lineage_names2[i]])
        }
    }

    lineage_names[include_lineage]
}
