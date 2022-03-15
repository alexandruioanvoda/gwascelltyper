#' Detect OS
#'
#' Function that outputs "windows", "unix" or "osx" depending on what system it is run on. Created by Will from R-bloggers here: https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
#'
#' @return character vector of size 1
#'
#' @export
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


#' Download dependencies
#'
#' Downloads data dependencies (e.g. linkage disequilibrium, gene annotations, etc.) pre-formatted for each of the 3 packages (LDSC, MAGMA & SNPsea)
#'
#' @param address Path where dependencies are downloaded. Do not modify, parameter meant for developers.
#'
#' @importFrom utils askYesNo
#' @importFrom utils download.file
#'
#' @return Path
#'
#' @export
download_dependencies <- function(address = paste0(system.file(package="gwascelltyper"), "/extdata/")) {
  # Because some of these downloads may take quite a while to finish,
  # and R timeouts downloads after 60 seconds by default, we're recalibrating that parameter to 1 hour for users with slow internet.
  options(timeout=60*60)
  possible_pops = c("eur", "eas", "afr", "amr", "sas", "subpop") # For which there would be readily-available, pre-processed MAGMA LD reference datasets.
  
  # MAGMA data --------
  if (utils::askYesNo("Do you want to download the MAGMA dependency data? (recommended yes)")) {
    utils::download.file(url="https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip", destfile = paste0(address, "NCBI38.zip"))
    utils::unzip(paste0(address,"NCBI38.zip"), exdir = address)
    utils::download.file(url="https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip", destfile = paste0(address, "NCBI37.3.zip"))
    utils::unzip(paste0(address,"NCBI37.3.zip"), exdir = address)
    # reformatting the MAGMA gene-loc files a bit:
    df <- utils::read.table(paste0(address, "NCBI38.gene.loc")); df <- df[,c(6, 2:5, 1)]
    utils::write.table(df, file = paste0(address, "NCBI38.hgnc.gene.loc"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    file.rename(from = paste0(address, "NCBI38.gene.loc"), to = paste0(address, "NCBI38.entrez.gene.loc"))
    df <- utils::read.table(paste0(address, "NCBI37.3.gene.loc")); df <- df[,c(6, 2:5, 1)]
    utils::write.table(df, file = paste0(address, "NCBI37.3.hgnc.gene.loc"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    file.rename(from = paste0(address, "NCBI37.3.gene.loc"), to = paste0(address, "NCBI37.3.entrez.gene.loc"))
    # Now, let's download LD data:
    population <- tolower(readline(prompt="Enter which population (from: eur, eas, afr, amr, sas, subpop): "))
    if (!population %in% possible_pops) {stop("Wrong population name. Choose from eur, eas, afr, amr, sas or subpop.")}
    utils::download.file(paste0("https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_", population, ".zip"), destfile = paste0(address,"g1000_", population,".zip"))
    utils::unzip(paste0(address,"g1000_", population,".zip"), exdir = address)
    #file.remove(paste0(address,"g1000_", population, ".zip"))
  }
  
  # MAGMA executable --------
  if (askYesNo("Do you want to download the SNPsea executable?")) {
    if (get_os() == "windows") {
      utils::download.file(url = "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_win.zip", destfile = paste0(address, "magma_v1.10_win.zip"))
      utils::unzip(paste0(address,"magma_v1.10_win.zip"), exdir = address)
    } else {
      if (get_os() == "unix") {
        utils::download.file(url = "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_static.zip", destfile = paste0(address, "magma_v1.10_static.zip"))
        utils::unzip(paste0(address,"magma_v1.10_static.zip"), exdir = address)
      } else {
        if (get_os() == "osx") {
          utils::download.file(url = "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_mac.zip", destfile = paste0(address, "magma_v1.10_mac.zip"))
          utils::unzip(paste0(address,"magma_v1.10_mac.zip"), exdir = address)
        } else {
          stop("gwascelltyper cannot detect operating system properly. Please raise an issue on github mentioning this and sessionInfo() outputs.")
        }
      }
    }
  }
  
  # SNPsea data --------
  required_files_for_snpsea <- c("Lango2010.txt.gz", "TGP2011.bed.gz", "NCBIgenes2013.bed.gz")
  if (sum(required_files_for_snpsea %in% list.files(address)) != 3) {
    if (askYesNo("Do you want to download the SNPsea dependencies for the eur population? (recommended yes)")) {
      cat("Downloading the SNPsea data.\n")
      utils::download.file(url = "http://files.figshare.com/1504037/SNPsea_data_20140520.zip",
                           destfile = paste0(address, "SNPsea_data_20140520.zip"))
      utils::unzip(zipfile = paste0(address, "SNPsea_data_20140520_bak.zip"),
                   exdir = address)
      # Reformatting some dependency data, won't take long.
      x <- read.table(file = paste0(address, "NCBIgenes2013.bed.gz"))[,c(1,2,3,5,4)]
      xfile <- gzfile(paste0(address, "NCBIgenes2013_hgnc.bed.gz"), "w")
      write.table(x = x, file = xfile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      close(xfile)
      rm(xfile)
      file.copy(from = paste0(address, "NCBIgenes2013.bed.gz"), to = paste0(address, "NCBIgenes2013_entrez.bed.gz"))
    }
  } else {
    cat("All SNPsea dependency files are downloaded for the eur population, in the following path:", address,"\n")
  }
  # SNPsea executable --------
  if (askYesNo("Do you want to download the SNPsea executable?")) {
    if (get_os() == "windows") {
      stop("SNPsea does not have a publicly available executable for Windows, but you can compile the C++ code from https://github.com/slowkow/snpsea\nIt needs to be placed in gwascelltyper/extdata/snpsea-executable")
    } else {
      if (get_os() == "unix") {
        download.file(url = "https://github.com/slowkow/snpsea/blob/master/bin/snpsea-linux64?raw=true",
                      destfile = paste0(address, "snpsea-executable"))
      } else {
        if (get_os() == "osx") {
          stop("SNPsea does not have a publicly available executable for OSX, but you can compile the C++ code from https://github.com/slowkow/snpsea\nIt needs to be placed in gwascelltyper/extdata/snpsea-executable")
        } else {
          stop("gwascelltyper cannot detect operating system properly.")
        }
      }
    }
  }
  
  # LDSC data --------
  if (utils::askYesNo("Do you want to download the LDSC dependency data?")) {
    population <- tolower(readline(prompt="Enter which population (from: eur, eas, afr, amr, sas, subpop): "))
    if (!population %in% possible_pops) {stop("Wrong population name. Choose from eur, eas, afr, amr, sas or subpop.")}
    if (!population %in% c("eur", "eas")) {stop("LDSC does not have pre-processed downloadable LD scores on the public server for afr, amr, sas and subpop, but instructions are provided in the wiki (https://github.com/bulik/ldsc/wiki/) for how to make such files. Please download eur or eas to the data folder (to see how they're formatted and placed), and repeat this for your own LD data.")}
    
    # Population-agnostic files that you should download and process:
    utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz", destfile = paste0(address, "hapmap3_snps.tgz"))
    utils::untar(paste0(address, "hapmap3_snps.tgz"), exdir = address)
    utils::download.file(url = "https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2", destfile = paste0(address, "w_hm3.snplist.bz2"))
    utils::write.table(x = utils::read.table(file = paste0(address, "w_hm3.snplist.bz2"), header = TRUE), quote = FALSE, row.names = FALSE, file = paste0(address, "w_hm3.snplist"))
    utils::download.file(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz", destfile = paste0(address, "refGene.txt.gz"))
    refGene <- read.table(file = paste0(address, "refGene.txt.gz"))
    refGene <- refGene[,c(13,3,5,6)] # These are the only columns we need.
    refGene <- refGene[-which(duplicated(refGene[,1])),] # There are duplicate gene names due to isoform annotations. To keep it simple, we'll just remove duplicated gene IDs.
    refGene <- refGene[order(refGene[,2]),]
    colnames(refGene) <- c("GENE", "CHR", "START", "END")
    write.table(x = refGene, file = paste0(address, "refGene_coord.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
    
    if (population == "eur") {
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz", destfile = paste0(address,"1000G_Phase3_baseline_ldscores.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz", destfile = paste0(address," weights_hm3_no_hla.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz", destfile = paste0(address,"1000G_Phase3_plinkfiles.tgz"))
      temp <- lapply(paste0(address, c("1000G_Phase3_baseline_ldscores.tgz", "weights_hm3_no_hla.tgz", "1000G_Phase3_plinkfiles.tgz")),
                     function(x) {utils::untar(x, exdir = address); return(1)})
      # Renaming some of the files for proper standard input and deleting temp archives
      file.rename(from = paste0(address,"weights_hm3_no_hla"), to = paste0(address,"1000G_EUR_Phase3_weights_hm3_no_hla"))
    } else if (population == "eas") {
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz",destfile = paste0(address,"1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz",destfile = paste0(address,"1000G_Phase3_EAS_weights_hm3_no_hla.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_plinkfiles.tgz",destfile = paste0(address,"1000G_Phase3_EAS_plinkfiles.tgz"))
      temp <- lapply(paste0(address, c("1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz", "1000G_Phase3_EAS_weights_hm3_no_hla.tgz", "1000G_Phase3_EAS_plinkfiles.tgz")),
                     function(x) {utils::untar(x, exdir = address); return(1)})
      file.rename(paste0(address, "1000G_Phase3_EAS_plinkfiles"), paste0(address, "1000G_Phase3_EAS_plink"))
      # Renaming some of the files for proper standard input
      for (i in list.files(paste0(address, "1000G_Phase3_eas_weights_hm3_no_hla"), full.names = TRUE)) {
        file.rename(from = i, to = paste0(dirname(i), "/", to = gsub(pattern = "weights.EAS.hm3_noMHC", replacement = "weights", basename(i))))
      }
    }
  }
  
  # LDSC executable scripts & environments --------
  if (utils::askYesNo("Do you want to download the LDSC scripts?")) {
    system(paste0("cd ", address, " && git clone https://github.com/bulik/ldsc.git && mv ./ldsc/* ./ && conda env create --file environment.yml && source activate ldsc && conda deactivate"))
  }

  return(address)
}

#' Find version of MAGMA binary
#'
#' Returns the version number of a MAGMA binary
#'
#' @param path Path to MAGMA binary
#'
#' @return version (character)
#'
#' @export
find_magma_version <- function(path) {
  version = gsub("[[:space:]]", "", x = substr(x = gsub(x = suppressWarnings(system(paste0(path, " --h"), intern = TRUE))[1],
                                                        pattern = "Welcome to MAGMA v",
                                                        replacement = "", fixed = TRUE),
                                               start = 1, stop = 5))
  return(version)
}
