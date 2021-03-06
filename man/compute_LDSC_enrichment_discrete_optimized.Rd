% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_testing.r
\name{compute_LDSC_enrichment_discrete_optimized}
\alias{compute_LDSC_enrichment_discrete_optimized}
\title{Compute LDSC enrichment - optimized}
\usage{
compute_LDSC_enrichment_discrete_optimized(
  gene_sets,
  sumstats_file,
  output_dir = NULL,
  number_of_threads = 1,
  window_size = 1e+05,
  population = "eur",
  ld_wind_cm = 1,
  ldsc_path = paste0(system.file(package = "gwascelltyper"), "/extdata/"),
  gene_coord = paste0(system.file(package = "gwascelltyper"),
    "/extdata/refGene_coord.txt"),
  print_snps = paste0(system.file(package = "gwascelltyper"),
    "/extdata/hapmap3_snps/hm."),
  bim_prefix = paste0(system.file(package = "gwascelltyper"), "/extdata/1000G_",
    toupper(population), "_Phase3_plink/1000G.", toupper(population), ".QC."),
  ref_ld_chr = paste0(system.file(package = "gwascelltyper"), "/extdata/1000G_",
    toupper(population), "_Phase3_baseline/baseline."),
  w_ld_chr = paste0(system.file(package = "gwascelltyper"), "/extdata/1000G_",
    toupper(population), "_Phase3_weights_hm3_no_hla/weights.")
)
}
\arguments{
\item{gene_sets}{List of genesets. Each geneset must be a character vector, containing HGNC gene names.}

\item{sumstats_file}{Address of properly formatted LDSC GWAS file.}

\item{output_dir}{Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.}

\item{number_of_threads}{Number of threads to parallelize calculation over. Default = 1.}

\item{window_size}{Size of window surrounding gene, for SNP-to-gene mapping purposes. Default = 100000 (100kb).}

\item{population}{Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)}

\item{ld_wind_cm}{Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM). Default = 1.}

\item{ldsc_path}{Path to folder containing all LDSC scripts.}

\item{gene_coord}{Adress of UCSC bed file containing all the gene coordinates.}

\item{print_snps}{Path of file containing HapMap3 SNPs.}

\item{bim_prefix}{Path + prefix of plink BIM files required by LDSC to compute LD scores.}

\item{ref_ld_chr}{Path of file containing LD reference panel. Tell LDSC which LD Scores to use as the predictors in the LD Score regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. LDSC will automatically concatenate .l2.ldscore files split across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. If the filename prefix contains the symbol (at), LDSC will replace the (at) symbol with chromosome numbers. Otherwise, LDSC will append chromosome numbers to the end of the filename prefix. Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz Example 2: --ref-ld-chr ld/(at)_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz}

\item{w_ld_chr}{Filename prefix for file with LD Scores with sum r^2 taken over SNPs included in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz. LDSC will read files split into 22 chromosomes in the same manner as --ref-ld-chr.}
}
\value{
Dataframe
}
\description{
Compute GWAS enrichment in discrete lists of genes with LDSC. Optimized to load all dependency files into Python-allocated RAM only once.
}
