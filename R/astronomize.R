#' re_findall
#'
#' Emulate behaviour of Python's re.findall() function. Lovingly stolen from \url{https://github.com/PoonLab/gromstole}.
#'
#' @param pat:  regex pattern
#' @param s:  character, a single string
#' @return character, vector of all matching substrings
re_findall <- function(pat, s) {
    if (!is.character(s)) {
        stop("re.findall() requires a character object for input 's'")
    }
    matches <- gregexpr(pat, s)
    index <- as.integer(matches[[1]])
    match_length <- attr(matches[[1]], "match.length")

    sapply(seq_along(index), function(i) {
        start <- index[i]
        stop <- start + match_length[i] - 1
        substr(s, start, stop)
    })
}

#' Extract mutation list from a directory of constellation files.
#'
#' "Constellations" are files produced from \url{https://github.com/cov-lineages/constellations}, which represent the cov-lineage team's best knowledge about which mutations define a lineage. These are updated frequently, please \code{git pull} in your cloned copy accordingly. See details.
#'
#' Code is adapted from \code{scripts/estimate-freqs.R} in \url{https://github.com/PoonLab/gromstole}.
#'
#' @param path Path to the constellations folder in the cov-lineages/constellations repository. The default assumes that the current project and the constellations repo are in the same directory.
#'
#' @return A lineage matrix for use with provoc.
#' @export
#'
#' @details From the repo, a constellation "a collection of mutations which are functionally meaningful, but which may arise independently a number of times".
#'
#' This function requires a local clone of the constellations repository. By default, I assume that your current project is in the same parent directory as the constellation file (e.g. your working directory is \code{.../parent/current_project} and the constellations repo is cloned as \code{.../parent/constellations}).
#'
#' If this is not the case, a path to the root folder of the constellations repo is required. The path is to the root folder, not the folder with the constellation files.
#'
#' There are no options to specify which lineages to include. I am working on a \code{annihilate(coco, lineage_defs)} function to remove mutations that aren't shared between coco and lineage_defs and squash lineages that have too few observed mutations or are too similar to other lineages.
#'
#' @examples
#' if(dir.exists("../constellations")) {
#'     lineage_defs <- astronomize()
#' }
#'
astronomize <- function(path = NULL) {

    if (is.null(path)) {
        return(provoc::lineage_defs_from_list(provoc::constellation_lists))
    }

    if (!dir.exists(path)) {
        return(provoc::lineage_defs_from_list(provoc::constellation_lists))
    }

    orfs <- list(
        "orf1a" = c(265, 13468),
        "orf1b" = c(13467, 21555),
        "S" = c(21562, 25384),
        "orf3a" = c(25392, 26220),
        "E" = c(26244, 26472),
        "M" = c(26522, 27191),
        "orf6" = c(27201, 27387),
        "orf7a" = c(27393, 27759),
        "orf7b" = c(27755, 27887),
        "orf8" = c(27893, 28259),
        "N" = c(28273, 29533),
        "orf10" = c(29557, 29674),
        "nsp2" = c(806, 2719),
        "nsp3" = c(2720, 8554),
        "nsp4" = c(8555, 10054),
        "nsp5" = c(10055, 10972),
        "nsp6" = c(10973, 11842),
        "nsp7" = c(11263, 11509),
        "nsp12" = c(13442, 16236),
        "nsp13" = c(16237, 18039),
        "nsp15" = c(19621, 20658)
    )

    # Removes posssible trailing slash for consistency, then adds it back
    path <- gsub("/$", "", path)
    path <- paste0(path, "/constellations/definitions/")
    sitelist <- lapply(list.files(path, full.names = TRUE), function(stelfile) {
        constellation <- jsonlite::read_json(stelfile, simplifyVector = TRUE)
        constellation$sites <- unique(constellation$sites)

        len_1a <- (orfs[["orf1a"]][2] - orfs[["orf1a"]][1]) / 3 + 1

        # convert constellation to label notation in the mapped files
        sites <- lapply(unique(constellation$sites), function(d) {
            toks <- toupper(strsplit(d, ":")[[1]])

            if (toks[1] != "DEL" && toks[1] != "NUC")
                toks <- c("aa", toks)

            if (toks[2] == "S" || toks[2] == "SPIKE") {
                toks[[2]] <- "S"
            } else if (toks[1] == "DEL") {
                toks[[1]] <- "del"
            } else if (toks[1] == "NUC") {
                toks <- toks[-1]
            } else if (toks[2] == "8") {
                toks[[2]] <- "orf8"
            } else if (toks[2] == "ORF1AB" || toks[2] == "1AB") {
                num <- as.numeric(re_findall("\\d+", toks[3]))

                if (num <= len_1a) {
                    toks[[2]] <- "orf1a"
                } else {
                    # Determine nucleotide position relative to start of orf1b
                    new_pos <- (((num - 1) * 3 + orfs[["orf1a"]][1]) -
                            orfs[["orf1b"]][1]) / 3
                    toks[[3]] <- gsub(num, floor(new_pos) + 1,
                        toks[[3]])
                    toks[[2]] <- "orf1b"
                }

            } else if (nchar(toks[2]) >= 3 &&
                    substring(toks[2], 1, 3) == "ORF") {
                toks[[2]] <- tolower(toks[2])
            } else if (substring(toks[2], 1, 3) == "NSP") {
                start_pos <- orfs[[tolower(toks[2])]][1]
                codon <- as.integer(re_findall("\\d+", toks[3]))
                nuc_pos <- start_pos + (codon - 1) * 3
                if (nuc_pos >= orfs[["orf1a"]][1] &&
                        nuc_pos <= orfs[["orf1a"]][2]) {
                    toks[[2]] <- "orf1a"
                } else if (nuc_pos >= orfs[["orf1b"]][1] &&
                        nuc_pos <= orfs[["orf1b"]][2]) {
                    toks[[2]] <- "orf1b"
                } else {
                    stop("Could not convert nsp to orf1a/b")
                }
                new_pos <- ((nuc_pos - orfs[[toks[2]]][1]) / 3) + 1
                toks[[3]] <- gsub(codon, floor(new_pos), toks[[3]])
            }

            if (grepl("+", toks[1], fixed = TRUE)) {
                ins <- strsplit(toks[[1]], split = "[+]")[[1]]
                toks[[1]] <- gsub(" ", "",
                    paste("+", ins[1], ".", ins[2]))
            }
            toks <- paste(toks, collapse = ":")
        })

        sites <- unlist(sites, recursive = FALSE)
    })

    lineage_defs <- as.matrix(dplyr::bind_rows(sapply(sitelist, function(sites) {
        x <- rep(1, length(sites))
        names(x) <- sites
        x
    })))

    lineage_defs[is.na(lineage_defs)] <- 0
    rownames(lineage_defs) <- gsub(".json", "", list.files(path, full.names = FALSE))
    rownames(lineage_defs) <- gsub("_constellation", "", rownames(lineage_defs),
        fixed = TRUE)
    rownames(lineage_defs) <- gsub("c", "", rownames(lineage_defs))

    lineage_defs <- lineage_defs[, colSums(lineage_defs) > 0]

    # Manual fixes based on known naming anomalies
    colnames(lineage_defs)[which(colnames(lineage_defs) == "+22205.GAGCCAGAA")] <-
        "ins:22205:9"
    colnames(lineage_defs)[which(colnames(lineage_defs) == "+28262.AACA")] <- "ins:28262:4"
    colnames(lineage_defs)[which(colnames(lineage_defs) == "28271-")] <- "del:28271:1"
    colnames(lineage_defs)[which(colnames(lineage_defs) == "A28271-")] <- "del:28271:1"
    # Spike protein starts at position 21562
    # Position reported as index of amino acids, hence 3*246
    colnames(lineage_defs)[which(colnames(lineage_defs) == "aa:S:RSYLTPG246-")] <-
        paste0("del:", 21562 + 3 * 246, ":21")
    colnames(lineage_defs)[which(colnames(lineage_defs) == "aa:S:Y144-")] <-
        paste0("del:", 21562 + 3 * 144, ":1")
    colnames(lineage_defs)[which(colnames(lineage_defs) == "aa:S:HV69-")] <-
        paste0("del:", 21562 + 3 * 69, ":2")
    # ORF 1a starts at 265
    colnames(lineage_defs)[which(colnames(lineage_defs) == "aa:orf1a:SGF3675-")] <-
        paste0("del:", 265 + 3 * 3675 - 1, ":9")
    as.matrix(lineage_defs)
}

#' Filter lineage_defs for specific lineages, keeping mutations that are present in at least one lineage
#'
#' @param lineage_defs The result of \code{astronomize()}. If NULL, tries to run \code{astronoimize}.
#' @param lineages Vector of lineage names (must be in \code{rownmaes(lineage_defs)}). Defaults to lineages circulating in 2021-2022.
#' @param return_df Should the function return a data frame? Note that returned df is transposed compared to lineage_defs. Default FALSE.
#' @param path Passed on to \code{astronomize} if \code{lineage_defs} is NULL.
#' @param shared_order Put shared mutations first? Default TRUE.
#'
#' @return A lineage matrix with fewer rows and columns than \code{lineage_defs}. If \code{return_df}, the columns represent lineage names and a \code{mutations} column is added.
#' @export
#'
#' @details After removing some lineage, the remaining mutations might not be present in any of the remaining lineage. This function will remove mutations that no longer belong to any lineage.
#'
#' Note that return_df will
filter_lineages <- function(
    lineage_defs = NULL,
    lineages = c("B.1.526", "B.1.1.7", "B.1.351", "B.1.617.2",
        "B.1.427", "B.1.429", "P.1"),
    return_df = FALSE,
    path = NULL,
    shared_order = TRUE) {

    if (is.null(lineage_defs)) {
        lineage_defs <- astronomize(path = path)
    }

    lineage_defs <- lineage_defs[lineages, ]
    lineage_defs <- lineage_defs[, apply(lineage_defs, 2, sum) > 0]

    if (shared_order) {
        lineage_defs <- lineage_defs[rev(order(apply(lineage_defs, 1, sum))),
            rev(order(apply(lineage_defs, 2, sum)))]
    }

    if (return_df) {
        mutnames <- colnames(lineage_defs)
        lineage_defs <- as.data.frame(t(lineage_defs))
        lineage_defs$mutation <- mutnames
        rownames(lineage_defs) <- NULL
    }

    return(lineage_defs)
}

#' Obtain and clean barcodes file from usher_barcores in Freyja (or elsewhere)
#' 
#' @param path The location to store the barcodes file. Tries data/clean, data/, then the current directory.
#' @param write Should the file be written to disk to avoid downloading? If TRUE, uses the first path that exists.
#' @param url The URL (or file path) for the barcodes file. Defaults 
#' @param update If TRUE, overwrite the existing barcodes file
#' 
#' @export
usher_barcodes <- function(
    path = c("data/clean/", "data/", "./")[1],
    write = TRUE,
    url = "https://raw.githubusercontent.com/andersen-lab/Freyja/main/freyja/data/usher_barcodes.csv",
    update = FALSE) {

    if (substr(path, nchar(path), nchar(path)) != "/") {
        path <- paste0(path, "/")
    }

    barcodes_exists <- FALSE
    for (pathname in c(path, "data/clean/", "data/", "./")) {
        if (file.exists(paste0(pathname, "usher_barcodes.csv")) && !update) {
            barcodes <- read.csv(paste0(pathname, "usher_barcodes.csv"))
            barcodes_exists <- TRUE
            break
        }
    }

    if (!barcodes_exists) {
        cat(ifelse(update, "", "usher_barcodes.csv not found."),
            "Downloading from Freyja repository.\n")
        barcodes <- read.csv(url)
        # First column is lineage names
        rownames(barcodes) <- barcodes[, 1]
        barcodes <- barcodes[, -1]
        for (pathname in c(path, "data/clean/", "data/", "./")) {
            if (dir.exists(pathname)) {
                cat("Writing to", paste0(pathname, "usher_barcodes.csv"), "\n")
                write.csv(barcodes,
                    paste0(pathname, "usher_barcodes.csv"),
                    row.names = TRUE)
                break
            }
        }
    }

    if ("X" %in% colnames(barcodes)) {
        rownames(barcodes) <- barcodes[, "X"]
        barcodes[, "X"] <- NULL
        barcodes <- as.matrix(barcodes)
    }

    colnames(barcodes) <- paste0("~",
        sapply(colnames(barcodes), function(x) {
            substr(x, 2, nchar(x))
        })) |> provoc:::parse_mutations()
    barcodes
}
