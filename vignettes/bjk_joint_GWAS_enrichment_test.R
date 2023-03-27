## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("gwascelltyper", quietly = TRUE))
#      devtools::install_github("alexandruioanvoda/gwascelltyper")
#  
#  # once installed:
#  library("gwascelltyper")
#  download_dependencies(pop = "eur")

## ----eval=FALSE---------------------------------------------------------------
#  download.file(url = "https://s3.tebi.io/bjk-vignette/bjk_vignette.tgz",
#                destfile = "bjk_vignette.tgz")
#  utils::untar("bjk_vignette.tgz") # This unzips into a folder called "To_share".
#  
#  library("gwascelltyper")
#  mat <- readRDS("./To_share/bjk_last_step_data.rds")
#  head(mat)
#  finalStep_bjk_joint_test(mat)

## ----eval=FALSE---------------------------------------------------------------
#  library("gwascelltyper")
#  gene_sets <- readRDS("./To_share/geneset.rds")[c("DC", "All_genes_control")]
#  # We're subsetting to just one gene-set (DC -> dendritic cells)
#  # and the control (background) gene-set.
#  
#  results <- compute_bjk_enrichment_discrete(gene_sets = gene_sets,
#                                  sumstats_files = c("./To_share/magma_CD.gwas", "./To_share/ldsc_CD.gwas.gz"),
#                                  output_dir = "./To_share/interm", # where we'll store our intermediary files.
#                                  pop = "eur",
#                                  n_blocks = 200,
#                                  number_of_threads = 1) # On macbooks, this can be increased to 2 or 3 to increase the computation speed
#  
#  results[[2]] # The matrix with delete-n values
#  results[[1]] # The final bjk, magma and ldsc values

