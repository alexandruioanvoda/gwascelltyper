% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_testing.r
\name{finalStep_bjk_joint_test}
\alias{finalStep_bjk_joint_test}
\title{Compute the block jackknife joint test discrete}
\usage{
finalStep_bjk_joint_test(mat)
}
\arguments{
\item{mat}{Dataframe containing blocked enrichment estimates (block == N) of enrichment and from whole data (block == 0). Needs to have 5 columns: Name (aka gene-set), P_value, BLOCK (1:200 as index of deleted block, 0 for whole data), method (in this case, LDSC or MAGMA) and Z (Z-score obtained from P-value).}
}
\value{
Dataframe
}
\description{
Compute GWAS enrichment in discrete lists of genes with the block jackknife joint test.
}
