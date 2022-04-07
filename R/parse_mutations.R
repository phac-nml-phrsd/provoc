#' Amino acid one-letter codes
#' 
#' As with many things in this repo, lovingly stolen from \url{https://github.com/PoonLab/gromstole}.
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

#' Open reading frame positions
#' 
#' As with many things in this repo, lovingly stolen from \url{https://github.com/PoonLab/gromstole}.
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

#' Parse mutations from (type, pos, alt) to amino acid
#' 
#' As with many things in this repo, lovingly stolen from \url{https://github.com/PoonLab/gromstole}.
#' 
#' @param type Either "~", "+", or "-" for mutation, insertion, or deletion, respectively.
#' @param pos The one-indexed position relative to the reference.
#' @param alt The alternate nucleotide, nucleotides, or number of deletions, respectively.
#' @param reffile The path to the reference file that was used to define mutations. Usually NC_045512.fa except in very particular circumstances.
#' 
#' @return e.g. aa:orf1a:K856R
#' @export
parse_mutation <- function(type, pos, alt, 
    reffile = system.file("inst/extdata/NC_045512.fa", package = "provoc")) {

    pos <- as.numeric(pos) - 1 # convert to 0-indexed

    refseq <- unlist(strsplit(readLines(reffile)[-1], split = ""))

    if (type == "~") {
        this_orf <- "None"
        this_left <- this_right <- "None"

        i <- 1
        while(pos >= orfs[[i]][1]) {
            i <- i + 1
        }
        i <- i - 1
        if (pos < orfs[[i]][2]) {
            this_orf <- names(orfs)[i]
            this_left <- orfs[[i]][1] 
            this_right <- orfs[[i]][2] 
        }

        if(this_orf != "None") {
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
        return(paste0("ins:", pos + 1, ":", nchar(alt)))
    } else if (type == "-") {
        return(paste0("del:", pos + 2, ":", alt))
    }

    return(paste0(refseq[pos + 1], pos + 1, alt))
}
