% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_scoring_big_matrices.r
\name{exp_environ_missingness}
\alias{exp_environ_missingness}
\title{Expression matrix (in environ object) missingness evaluation}
\usage{
exp_environ_missingness(
  obj,
  col_max_miss = 0.1,
  row_max_miss = 0.1,
  elim = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{obj}{Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).}

\item{col_max_miss}{Column maximum missingness (default = 0.1)}

\item{row_max_miss}{Row maximum missingness (default = 0.1)}

\item{elim}{Whether to eliminate rows & cols with too many NAs (default = TRUE), alternative is to just throw a warning.}

\item{verbose}{Print what function has done.}
}
\value{
List (of two vectors: bad row & column indexes)
}
\description{
Warning: this function modifies any object it is given. This function does two things at once:
1. automatically removes completely-NA rows and columns.
2. finds out if some rows or columns have more missingness than a threshold (default = 0.1 (otherwise written, 10\%)), and eliminates them (if requested).
}
\details{
It then returns a list of two vectors, containing the bad row & column indexes respectively.
}
