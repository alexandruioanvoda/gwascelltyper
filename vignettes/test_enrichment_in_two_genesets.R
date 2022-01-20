## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  print("hi")
#  

## ----eval=TRUE----------------------------------------------------------------
library(gwascelltyper)
munge_sumstats_for_MAGMA(preview=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(gwascelltyper)
#  
#  path = paste0(getwd(), "/cd_build37_40266.txt") # it could also be "/path/to/your/GWASsumstats.txt"
#  
#  munge_sumstats_for_MAGMA(sumstats_raw = path)
#  
#  file.exists(paste0(path, "_processed.sumstats"))

