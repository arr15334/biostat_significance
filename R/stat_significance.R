#library(Biostrings)
#library(seqinr)
#library(ggplot2)
#library(VGAM)


#' Statistical significance
#'
#' Calculates probability of an alignment score to have been obtained randomly
#' @param seq1 First sequence
#' @param seq2 Second sequence
#' @param seq.type Sequence type (DNA/RNA/Protein)
#' @param alignment Alignment type (global or local)
#' @param sub.matrix Substitution matrix to use for score calculation (BLOSUM62, custom mat,...)
#' @param gap.open Penalty for opening a new gap
#' @param gap.ext Penalty for extending current gap
#' @param shuffles Number of shuffles to calculate
#' @param shuffle.seq Sequence to be shuffled (1 or 2)
#' @return Histogram of scores, original score, probability, standardized score and k constant of gumbell distribution
#' @example
#' scores <- stat_sign('../sequence_B.fasta', '../NC_000001.11.fasta', 'DNA','global', 'BLOSUM62', shuffles=50)
#' @export
stat_sign <- function(seq1, seq2, seq.type, alignment, sub.matrix,
                      gap.open = 3, gap.ext = 1, shuffles = 100,
                      shuffle.seq = 1) {
  randallscore <- double()
  x <- read_sequence(seq1, seq2, seq.type)$seq1.read
  y <- read_sequence(seq1, seq2, seq.type)$seq2.read
  # score original sequences
  randallscore[1] <- Biostrings::pairwiseAlignment(x, y,
                                       substitutionMatrix = sub.matrix,
                                       gapOpening = -gap.open,
                                       gapExtension = -gap.ext,
                                       type=alignment,
                                       scoreOnly = TRUE)
  n1 <- length(x)
  n2 <- length(y)
  for (i in 2:shuffles) {
    if (shuffle.seq == 1) {
      x <- seqinr::c2s(sample(x, n1, replace=FALSE))
      x <- Biostrings::DNAString(x)
    } else {
      y <- seqinr::c2s(sample(y,n2, replace=FALSE))
      y <- Biostrings::DNAString(y)
    }
    randallscore[i] <- Biostrings::pairwiseAlignment(x, y,
                                         substitutionMatrix = sub.matrix,
                                         gapOpening = -gap.open,
                                         gapExtension = -gap.ext,
                                         type=alignment,
                                         scoreOnly = TRUE)
  }
  xmean <- mean(randallscore)
  xs <- sd(randallscore)
  lambda <- 1.2825 / xs
  u <- xmean - 0.45 * xs
  k <- exp(lambda * u) / (n1*n2)
  s.st <- lambda * randallscore[1] - log(k*n1*n2)
  p.s <- 1 - exp(-exp(-s.st))
  g <- 0.57721566
  s <- (xmean-u)/g  #no estoy segura sobre el uso de los parametros en el histograma
  #s <- 0.78*s.st
  #Adjusting the gumbel distribution
  #dens <- dgumbel(randallscore, location = u, scale = s)
  #coef <- sum(dens)/length(randallscore)
  #dens2 <- dens/ coef
  randallscore2 <- data.frame(randallscore)
  histograma <- ggplot2::ggplot(data=randallscore2, ggplot2::aes(x=randallscore)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..) , binwidth=5 ) +
    ggplot2::stat_function(fun= dgumbel, args= list(loc= u, scale= s)) +
    ggplot2::geom_vline(xintercept = randallscore[1], color = "red") +
    ggplot2::labs(caption = paste("Alignment", alignment, "\nSequence1: ", seq1,
                         "\nSequence2: ", seq2, "\nType:", seq.type,
                         "\nMatrix:", sub.matrix, "\nGap open:", gap.open,
                         "\nGap extension:", gap.ext,
                         "\n# shuffles:", shuffles,
                         "\nShuffled sequence:", shuffle.seq,
                         "\nMode:", u,
                         "\nScale parameter:", s,
                         "\nOriginal Score:", randallscore[1],
                         "\nProbability:", p.s))

  return(list(scores.summary=summary(randallscore),
              h=histograma,
              S=randallscore[1],
              k=k,
              sprima=s.st,
              P=p.s
              ))
}
#mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)

#' Statistical significance
#'
#' Calculates probability of an alignment score to have been obtained randomly
#' @param seq1 First sequence as .fasta file
#' @param seq2 Second sequence as .fasta file
#' @param seq.type Type of sequences (DNA/RNA/Protein)
#' @return A list of the two sequences converted to given type
#' @export
read_sequence <- function(seq1, seq2, seq.type) {
  if (seq.type == 'RNA') {
    # read first fasta then convert to string object
    seq1.read <- Biostrings::readRNAStringSet(seq1)
    seq2.read <- Biostrings::readRNAStringSet(seq2)
    x <- Biostrings::RNAString(seq1.read)
    y <- Biostrings::RNAString(seq2.read)
  } else if (seq.type == 'PROTEIN') {
    seq1.read <- Biostrings::readAAStringSet(seq1)
    seq2.read <- Biostrings::readAAStringSet(seq2)
    x <- Biostrings::AAString(seq1.read)
    y <- Biostrings::AAString(seq2.read)
  } else {
    seq1.read <- Biostrings::readDNAStringSet(seq1)
    seq2.read <- Biostrings::readDNAStringSet(seq2)
    x <- Biostrings::DNAString(toString(seq1.read))
    y <- Biostrings::DNAString(toString(seq2.read))
  }
  return(list(seq1.read=x, seq2.read=y))
}


#scores <- stat_sign('./sequence_B.fasta', '../task3/fasta/NC_000001.11.fasta', 'DNA','global', mat, shuffles=50)






