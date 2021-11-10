# biostat_significance
R package to calculate the statistical significance of the alignment score of two sequences 

## Installation
Run the following command on the RStudio console
`devtools::install_github('arr15334/biostat_significance')`

## Usage
`library(statSignificance)`
`x <- stat_significance('file1.fasta', 'file2.fasta', '[DNA/RNA/Protein]'], substitutionMatrix, gapOpen, gapExtension, shuffles, shuffledSequence)`

The function will return:

+ an histogram pointing the original score compared with the n permutation scores -> x$h
+ A numerical summary of the scores -> x$scores.summary
+ Original sequences alignment score -> x$S
+ Gumbel k constant -> x$k
+ Standardized score -> x$s.st
+ Probability that the score can be achieved randomly -> x$P
