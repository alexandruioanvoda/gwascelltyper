% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_scoring_big_matrices.r
\name{read_exp_and_annot_environ}
\alias{read_exp_and_annot_environ}
\title{Read exp and annot}
\usage{
read_exp_and_annot_environ(obj, exp_rds, annot_rds)
}
\arguments{
\item{obj}{Environment to read exp and annotLevels to.}

\item{exp_rds}{Path to RDS file containing the expression matrix.}

\item{annot_rds}{Path to RDS file containing the annotation dataframe.}
}
\value{
None
}
\description{
Read expression and annotation data into the environment, then check the sample/cell tags for correspondence between exp and annot.
}
