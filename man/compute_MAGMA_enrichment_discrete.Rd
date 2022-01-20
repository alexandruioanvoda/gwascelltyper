% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_testing.r
\name{compute_MAGMA_enrichment_discrete}
\alias{compute_MAGMA_enrichment_discrete}
\title{Compute MAGMA enrichment}
\usage{
compute_MAGMA_enrichment_discrete(
  gene_sets,
  sumstats_file,
  output_dir = NULL,
  upstream_kb = 10,
  downstream_kb = 1.5,
  gene_nomenclature = "hgnc",
  population = "eur",
  gwas_sample_size = NULL,
  genome_ref_path = paste0(system.file(package = "gwascelltyper"), "/extdata/g1000_",
    population),
  magma_path = paste0(system.file(package = "gwascelltyper"), "/extdata/magma-linux64")
)
}
\arguments{
\item{gene_sets}{List of genesets. Each geneset must be a character vector, containing HGNC gene names.}

\item{sumstats_file}{Address of properly formatted MAGMA GWAS file.}

\item{output_dir}{Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.}

\item{upstream_kb}{Number of kb upstream of gene, SNP-to-gene mapping window. Default = 10.}

\item{downstream_kb}{Number of kb downstream of gene, SNP-to-gene mapping window. Default = 1.5.}

\item{gene_nomenclature}{Type of gene IDs fed in from the expression data (default = hgnc; alternatives = entrez).}

\item{population}{Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)}

\item{gwas_sample_size}{Sample size of GWAS.}

\item{genome_ref_path}{Path to genome references for MAGMA (default = packageLocation/extdata/magma/g1000_eur).}

\item{magma_path}{Path to MAGMA executable (default is in package extdata folder, where download_dependencies() should store it).}
}
\value{
Dataframe
}
\description{
Compute GWAS enrichment in discrete lists of genes with MAGMA.
}