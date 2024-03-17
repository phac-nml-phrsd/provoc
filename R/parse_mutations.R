
#' Parse mutations from (type, pos, alt) to amino acid
#' 
#' As with many things in this repo, lovingly stolen from the \code{seq_utils.py} script in \url{https://github.com/PoonLab/gromstole} (as well as PoonLab/covizu).
#' 
#' @param type Either "~", "+", or "-" for mutation, insertion, or deletion, respectively.
#' @param pos The one-indexed position relative to the reference.
#' @param alt The alternate nucleotide, nucleotides, or number of deletions, respectively.
#' @param reffile The path to the reference file that was used to define mutations. Usually NC_045512.fa except in very particular circumstances.
#' 
#' @return e.g. aa:orf1a:K856R
parse_mutation <- function(type, pos, alt, 
    reffile = system.file("extdata/NC_045512.fa", package = "provoc")) {
    gcode <- list(
        'TTT'= 'F', 'TTC'= 'F', 'TTA'= 'L', 'TTG'= 'L',
        'TCT'= 'S', 'TCC'= 'S', 'TCA'= 'S', 'TCG'= 'S',
        'TAT'= 'Y', 'TAC'= 'Y', 'TAA'= '*', 'TAG'= '*',
        'TGT'= 'C', 'TGC'= 'C', 'TGA'= '*', 'TGG'= 'W',
        'CTT'= 'L', 'CTC'= 'L', 'CTA'= 'L', 'CTG'= 'L',
        'CCT'= 'P', 'CCC'= 'P', 'CCA'= 'P', 'CCG'= 'P',
        'CAT'= 'H', 'CAC'= 'H', 'CAA'= 'Q', 'CAG'= 'Q',
        'CGT'= 'R', 'CGC'= 'R', 'CGA'= 'R', 'CGG'= 'R',
        'ATT'= 'I', 'ATC'= 'I', 'ATA'= 'I', 'ATG'= 'M',
        'ACT'= 'T', 'ACC'= 'T', 'ACA'= 'T', 'ACG'= 'T',
        'AAT'= 'N', 'AAC'= 'N', 'AAA'= 'K', 'AAG'= 'K',
        'AGT'= 'S', 'AGC'= 'S', 'AGA'= 'R', 'AGG'= 'R',
        'GTT'= 'V', 'GTC'= 'V', 'GTA'= 'V', 'GTG'= 'V',
        'GCT'= 'A', 'GCC'= 'A', 'GCA'= 'A', 'GCG'= 'A',
        'GAT'= 'D', 'GAC'= 'D', 'GAA'= 'E', 'GAG'= 'E',
        'GGT'= 'G', 'GGC'= 'G', 'GGA'= 'G', 'GGG'= 'G',
        '---'= '-', 'XXX'= '?'
    )
    orfs = list(
        'orf1a'= c(265, 13468),
        'orf1b'= c(13467, 21555),
        'S'= c(21562, 25384),
        'orf3a'= c(25392, 26220),
        'E'= c(26244, 26472),
        'M'= c(26522, 27191),
        'orf6'= c(27201, 27387),
        'orf7a'= c(27393, 27759),
        'orf7b'= c(27755, 27887),
        'orf8'= c(27893, 28259),
        'N'= c(28273, 29533),
        'orf10'= c(29557, 29674)
    )

    pos <- as.numeric(pos) - 1 # convert to 0-indexed

    refseq <- unlist(strsplit(readLines(reffile)[-1], split = ""))

    if (type == "~") {
        this_orf <- "None"
        this_left <- this_right <- "None"

        orfl <- sapply(orfs, function(x) x[1] < pos & x[2] > pos)
        if(any(orfl)) {
            this_orf <- names(orfl)[orfl][1]
            this_left <- orfs[[this_orf]][1]
            this_right <- orfs[[this_orf]][2]
        }

        if(this_orf != "None") {
            # 3 nucleotides per codon
            codon_left <- 3 * ((pos-this_left) %/% 3)
            codon_pos <- (pos-this_left) %% 3

            rcodon <- refseq[this_left:this_right + 1][codon_left:(codon_left+2) + 1]
            ramino <- gcode[[paste0(rcodon, collapse = "")]]

            qcodon <- rcodon
            qcodon[codon_pos + 1] <- alt
            qcodon <- paste0(qcodon, collapse = "")
            qamino <-  gcode[[paste0(qcodon, collapse = "")]]

            if (ramino != qamino)
                return(paste0("aa:", this_orf, ":", ramino, round(1+codon_left/3), qamino))
        }
    } else if (type == "+") {
        # Revert to 1-indexing
        return(paste0("ins:", pos + 1, ":", nchar(alt)))
    } else if (type == "-") {
        return(paste0("del:", pos + 1, ":", alt))
    }

    # Revert to 1-indexing
    return(paste0(refseq[pos + 1], pos + 1, alt))
}


#' Parse all unique mutations in a vector
#'
#' @param muts A vector of mutations in the format "+50535C", "-43234.2", and "+342234.AC".
#'
#' @return A data frame with columns `label` and `mutation`.
parse_unique_mutations <- function(muts) {
    unique_muts <- unique(muts)
    new_muts <- vector("character", length(unique_muts))  # new_muts vector

    for (i in seq_along(unique_muts)) {
        thismut <- unique_muts[i]
        first_char <- substr(thismut, 1, 1)

        if (first_char == "~") {
            # Extract position and alternate nucleotide
            pos <- as.numeric(substr(thismut, 2, nchar(thismut) - 1))
            alt <- substr(thismut, nchar(thismut), nchar(thismut))
            # Call parse_mutation with "~"
            new_muts[i] <- provoc:::parse_mutation("~", pos, alt)
        } else if (first_char %in% c("-", "+")) {
            # Split mutation string
            splits <- strsplit(thismut, "[+-]", perl=TRUE)[[1]]
            type <- first_char
            if (type == "+") {
              pos <- as.numeric(gsub(".*?([0-9]+).*", "\\1", splits[2]))
              alt <- gsub("[^a-zA-Z]", "", splits[2])
            }
            else {
              split_numbers <- strsplit(as.character(splits[2]), "[.]")[[1]]
              pos <- as.numeric(split_numbers[1])
              alt <- as.numeric(split_numbers[2])
            }
            # Call parse_mutation with type "-" or "+"
            new_muts[i] <- provoc:::parse_mutation(type, pos, alt)
        } else {
            new_muts[i] <- thismut
        }
    }

    # Turn new_muts into a dataframe
    data.frame(label = unique_muts, mutation = new_muts, stringsAsFactors = FALSE)
}


#' Parse output of the Gromstole pipeline, from SNVs to AAs
#' 
#' This function can be slow.
#' 
#' @param labels The `labels` column from GromStole output.
#' 
#' @return A vector of the same length of `labels`.
#' @export
parse_mutations <- function(labels) {
    unique_aa <- provoc:::parse_unique_mutations(unique(labels))
    return(unique_aa$mutation[match(labels, unique_aa$label)])
}
