% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_scoring_big_matrices.r
\name{compute_t_stat_matrix_environ}
\alias{compute_t_stat_matrix_environ}
\title{Compute t-stat matrix - from environment-stored expression data}
\usage{
compute_t_stat_matrix_environ(
  obj,
  number_of_threads = 1,
  col_max_miss = 0.1,
  row_max_miss = 0.1
)
}
\arguments{
\item{obj}{Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).}

\item{number_of_threads}{Number of threads to parallelize calculation over. Default = 1.}

\item{col_max_miss}{Maximum missingness (is.na) for the column (samples), between 0 and 1. Default is 10\% (col_max_miss=0.1). If a column contains more than 10\% NAs, it is eliminated from further analyses.}

\item{row_max_miss}{Maximum missingness (is.na) for the row (genes), between 0 and 1. Default is 10\% (row_max_miss=0.1). If a row contains more than 10\% NAs, it is eliminated from further analyses.}
}
\value{
Dataframe
}
\description{
Compute gene-by-celltype specificity matrix with the t-statistic. The input must be an environment containing an exp & annot.
}
