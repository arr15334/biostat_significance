% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_significance.R
\name{stat_sign}
\alias{stat_sign}
\title{Statistical significance}
\usage{
stat_sign(
  seq1,
  seq2,
  seq.type,
  alignment,
  sub.matrix,
  gap.open = 3,
  gap.ext = 1,
  shuffles = 100,
  shuffle.seq = 1
)
}
\arguments{
\item{seq1}{First sequence}

\item{seq2}{Second sequence}

\item{seq.type}{Sequence type (DNA/RNA/Protein)}

\item{alignment}{Alignment type (global or local)}

\item{sub.matrix}{Substitution matrix to use for score calculation (BLOSUM62, custom mat,...)}

\item{gap.open}{Penalty for opening a new gap}

\item{gap.ext}{Penalty for extending current gap}

\item{shuffles}{Number of shuffles to calculate}

\item{shuffle.seq}{Sequence to be shuffled (1 or 2)}
}
\value{
Histogram of scores, original score, probability, standardized score and k constant of gumbell distribution
}
\description{
Calculates probability of an alignment score to have been obtained randomly
}
