# ProVoC

PROportions of Variants of Concern using counts, coverage, and a variant matrix.

- Counts: The number of times a given mutation was observed.
- Coverage: The number of times the position of a given variant was read.
- Variant Matrix: A matrix where the row names are variants, the column names are mutations, and the entries are 1 if the row variant has the column mutation and 0 otherwise.

# Usage

There are two steps to using this software: create the variant matrix and run the model(s).

The variant matrix can be created several ways:

1. User-specified. The columns must be labelled with mutations (in the same format as the ones in your data), and the row names must be the names of the variant.
2. Specify variants of concern, and use Nextstrain Genbank data to determine the mutations and "nuisance" lineages.
    - A nuisance lineage is a lineage that shares a pre-specified percent of mutations with your VoCs.

Data must be specified as follows.
There must be a column labelled `count` and a column labelled `coverage`.
There can be an optional column labelled `sample`, which will be interpreted as a grouping variable and the analysis will be run once for each unique value of `sample`.
Any other columns will be ignored, so you are free to leave them in your data for future purposes (e.g. location and date will still be tied to unique sample ids).

The package currently can only be installed by time travellers who travel forward to a time when I have completed the package.
For all non Time Lords, you can source the `provoc.R` script and use the individual functions in their current state.
