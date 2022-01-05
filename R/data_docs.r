#' GTEX EWCE-computed tissue-gene specificity scores
#'
#' Data from The Genotype-Tissue Expression (GTEX) project (ref. 1);
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by GTEX,
#' ranked by EWCE (ref. 2)
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_gtex_ewce)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets gtex ewce
#'
#' @references 1. The Genotype-Tissue Expression (GTEX) project. Nat. Genet. 2013.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010069/}{PubMed})\cr
#' 2. Nathan Skene et al. Front. Neurosci. 2016.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730103/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_gtex_ewce)
#' head(gene_score_matrix_gtex_ewce)
#'
"gene_score_matrix_gtex_ewce"


#' GTEX t-stat-computed tissue-gene specificity scores
#'
#' Data from The Genotype-Tissue Expression (GTEX) project (ref. 1);
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by GTEX,
#' ranked by t-statistic (same method as that used in ref. 2)
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_gtex_t_stat)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets gtex t-stat t statistic
#'
#' @references 1. The Genotype-Tissue Expression (GTEX) project. Nat. Genet. 2013.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010069/}{PubMed})\cr
#' 2. H. Finucane et al. Nat. Gen. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5896795/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_gtex_t_stat)
#' head(gene_score_matrix_gtex_t_stat)
#'
"gene_score_matrix_gtex_t_stat"


#' HPCA EWCE-computed tissue-gene specificity scores
#'
#' Data from Mabbott et al.;
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by HPCA (Human Primary Cell Atlas - HPCA, ref. 1),
#' ranked by EWCE (ref. 2); processed expression data obtained from the SingleR package by Dr. Dvir Aran.
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_hpca_ewce)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets hpca ewce
#'
#' @references 1. Neil A. Mabbott et al. BMC Genomics. 2013.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24053356}{PubMed})\cr
#' 2. Nathan Skene et al. Front. Neurosci. 2016.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730103/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_hpca_ewce)
#' head(gene_score_matrix_hpca_ewce)
#'
"gene_score_matrix_hpca_ewce"


#' HPCA t-stat-computed tissue-gene specificity scores
#'
#' Data from Mabbott et al.;
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by HPCA (Human Primary Cell Atlas - HPCA, ref. 1),
#' ranked by t-statistic (ref. 2); processed expression data obtained from the SingleR package by Dr. Dvir Aran.
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_hpca_t_stat)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets hpca t statistic
#'
#' @references 1. Neil A. Mabbott et al. BMC Genomics. 2013.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24053356}{PubMed})\cr
#' 2. H. Finucane et al. Nat. Gen. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5896795/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_hpca_t_stat)
#' head(gene_score_matrix_hpca_t_stat)
#'
"gene_score_matrix_hpca_t_stat"


#' BPE EWCE-computed tissue-gene specificity scores
#'
#' Data from Mabbott et al.;
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by the Blueprint & ENCODE projects (ref. 1 & 2, combined & preprocessed data taken from SingleR - ref. 3),
#' ranked by EWCE (ref. 4); processed expression data obtained from the SingleR package by Dr. Dvir Aran.
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_bpe_ewce)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets hpca ewce
#'
#' @references 1. Stunnenberg & Hirst 2016\cr
#' 2. ENCODE Project Consortium 2012\cr
#' 3. Aran et al. 2019\cr
#' 4. Nathan Skene et al. Front. Neurosci. 2016.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730103/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_bpe_ewce)
#' head(gene_score_matrix_bpe_ewce)
#'
"gene_score_matrix_bpe_ewce"


#' BPE t-stat-computed tissue-gene specificity scores
#'
#' Data from Mabbott et al.;
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by the Blueprint & ENCODE projects (ref. 1 & 2, combined & preprocessed data taken from SingleR - ref. 3),
#' ranked by t-statistic (ref. 4); processed expression data obtained from the SingleR package by Dr. Dvir Aran.
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_bpe_t_stat)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets hpca t statistic
#'
#' @references 1. Stunnenberg & Hirst 2016\cr
#' 2. ENCODE Project Consortium 2012\cr
#' 3. Aran et al. 2019\cr
#' 4. H. Finucane et al. Nat. Gen. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5896795/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_bpe_t_stat)
#' head(gene_score_matrix_bpe_t_stat)
#'
"gene_score_matrix_bpe_t_stat"


#' ImmGen EWCE-computed tissue-gene specificity scores
#'
#' Data from Shay & Kang;
#' Contains matrix of tissue (columns) specificity scores for each gene (row) profiled by the ImmGen project (ref. 1),
#' ranked by EWCE (ref. 2).
#'
#' @docType data
#'
#' @usage data(gene_score_matrix_immgen_ewce)
#'
#' @format An object of class \code{"data matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords datasets immgen ewce
#'
#' @references 1. Shay & Kang 2013\cr
#' 2. Nathan Skene et al. Front. Neurosci. 2016.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730103/}{PubMed})
#'
#' @examples
#' data(gene_score_matrix_immgen_ewce)
#' head(gene_score_matrix_immgen_ewce)
#'
"gene_score_matrix_immgen_ewce"


#' Standard summstats column headers
#'
#' Equivalency table written by Nathan Skene et al. (2018) for their MAGMA.Celltyping package. Used to standardize GWAS summary statistics headers.
#'
#' @docType data
#'
#' @usage data(standard_sumstats_column_headers)
#'
#' @format An object of class \code{"data frame"}; see \code{\link[data.frame]{data.frame}}.
#'
#' @keywords standard sumstats column headers gwas
#'
#' @references 1. Nathan Skene et al. Nat. Genet. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/29785013}{PubMed})
#'
#' @examples
#' data(standard_sumstats_column_headers)
#' head(standard_sumstats_column_headers)
#'
"standard_sumstats_column_headers"


#' SNP location data
#'
#' SNP location written by Nathan Skene et al. (2018) for their MAGMA.Celltyping package. Used to find rsIDs for SNPs when they're not present.
#'
#' @docType data
#'
#' @usage data(snp_loc)
#'
#' @format An object of class \code{"data frame"}; see \code{\link[data.frame]{data.frame}}.
#'
#' @keywords snp location gwas
#'
#' @references 1. Nathan Skene et al. Nat. Genet. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/29785013}{PubMed})
#'
#' @examples
#' data(standard_sumstats_column_headers)
#' head(standard_sumstats_column_headers)
#'
"snp_loc"


#' Annotations for HPCA
#'
#' Cluster / cell-type annotation for the columns of the expression matrix of the Human Primary Cell Atlas (HPCA) data. Taken and processed from the SingleR package by Dvir Aran.
#'
#' @docType data
#'
#' @usage data(annot_hpca)
#'
#' @format An object of class \code{"list"}; see \code{\link[list]{list}}.
#'
#' @keywords annotation HPCA
#'
#' @references 1. Nathan Skene et al. Nat. Genet. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/29785013}{PubMed})
#'
#' @examples
#' data(annot_hpca)
#' head(annot_hpca)
#'
"annot_hpca"


#' The HPCA expression matrix
#'
#' Expression matrix of the Human Primary Cell Atlas (HPCA) data. Columns are sample IDs, rows are HGNC gene IDs. Processed data taken from the SingleR package by Dvir Aran.
#'
#' @docType data
#'
#' @usage data(exp_hpca)
#'
#' @format An object of class \code{"matrix"}; see \code{\link[matrix]{matrix}}.
#'
#' @keywords expression HPCA
#'
#' @references 1. Nathan Skene et al. Nat. Genet. 2018.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/29785013}{PubMed})
#'
#' @examples
#' data(exp_hpca)
#' head(exp_hpca)
#'
"exp_hpca"
