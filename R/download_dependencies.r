#' Detect OS
#'
#' Function that outputs "windows", "linux" or "osx" depending on what system it is run on. Created by Will from R-bloggers here: https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
#'
#' @return character vector of size 1
#'
#' @export
get_os <- function() {
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

#' Dependency data paths
#'
#' Checks package-specific paths for where the dependency-data is usually stored. Dependency data is mostly data on LD regions, but other required stats too, and are usually population- (e.g. eas, afr, eur, etc.) and package- (LDSC, MAGMA & SNPsea) specific.
#'
#' @param pop Reference population shortcode (available right now: eur; Soon to be added: afr, amr, eas, eur, sas).
#' @param for_which Which softwares to check (available right now: all, magma, snpsea or ldsc).
#'
#' @importFrom utils askYesNo
#'
#' @return Named list (of paths)
#'
#' @export
dep_path <- function(pop = c("eur", "afr", "amr", "eas", "sas")[1],
                     for_which = c("all", "magma", "snpsea", "ldsc")[1]) {
  
  popu <- toupper(pop)
  
  # First, defining the relative file names
  deps <- list()
  if (for_which %in% c("all", "ldsc")) {
    # LDSC executables
    deps[["ldsc_scripts"]] <- c("ldsc.py", "make_annot.py")
    
    # LDSC data
    deps[["1000g_bim_ldsc"]] <- paste0("1000G_",popu,"_Phase3_plink/1000G.",popu,".QC.",1:22,".bim")
    deps[["1000g_fam_ldsc"]] <- paste0("1000G_",popu,"_Phase3_plink/1000G.",popu,".QC.",1:22,".fam")
    deps[["1000g_bed_ldsc"]] <- paste0("1000G_",popu,"_Phase3_plink/1000G.",popu,".QC.",1:22,".bed")
    deps[["1000g_frq_ldsc"]] <- paste0("1000G_",popu,"_Phase3_frq/1000G.",popu,".QC.",1:22,".frq")
    deps[["1000g_ldscore_universal_ldsc"]] <- paste0("1000G_",popu,"_Phase3_weights_hm3_no_hla/weights.",1:22,".l2.ldscore.gz")
    deps[["1000g_ldscore_category_ldsc"]] <- paste0("1000G_",popu,"_Phase3_baseline/baseline.", 1:22, ".l2.ldscore.gz")
    deps[["1000g_m_m550_ldsc"]] <- c(paste0("1000G_",popu,"_Phase3_baseline/baseline.", 1:22, ".l2.M"),
                                               paste0("1000G_",popu,"_Phase3_baseline/baseline.", 1:22, ".l2.M_5_50"))
    deps[["1000g_annot_ldsc"]] <- paste0("1000G_",popu,"_Phase3_baseline/baseline.", 1:22, ".annot.gz")
    # I think this may actually not be needed, need to check later ^
    deps[["hm3snps_ldsc"]] <- paste0("hapmap3_snps/hm", 1:22, ".snp")
    deps[["geneloc_ldsc"]] <- "refGene_coord.txt"
  }
  if (for_which %in% c("all", "magma")) {
    if (get_os()=="linux") {
      deps[["magma_binary"]] <- "magma"
    } else if (get_os()=="osx") {
      deps[["magma_binary"]] <- "magma"
    } else if (get_os()=="windows") {
      deps[["magma_binary"]] <- "magma.exe"
    }
    
    # MAGMA data
    deps[["1000g_bim_magma"]] <- paste0("g1000_",pop,".bim")
    deps[["1000g_fam_magma"]] <- paste0("g1000_",pop,".fam")
    deps[["1000g_bed_magma"]] <- paste0("g1000_",pop,".bed")
    deps[["geneloc_magma"]] <- c("NCBI37.3.hgnc.gene.loc", "NCBI38.hgnc.gene.loc",
                                           "NCBI37.3.entrez.gene.loc", "NCBI38.entrez.gene.loc")
  }
  if (for_which %in% c("all", "snpsea")) {
    if (get_os()=="linux") {
      deps[["snpsea_binary"]] <- "snpsea-linux64"
    } else if (get_os()=="osx") {
      deps[["snpsea_binary"]] <- "snpsea-osx"
    } else if (get_os()=="windows") {
      deps[["snpsea_binary"]] <- "snpsea.exe"
    }
    
    # SNPsea data
    deps[["1000g_6colbed_snpsea"]] <- "TGP2011.bed.gz" # 1000G EUR linkage intervals # No other population provided yet, but am interested in figuring how to best do this
    deps[["txt_snpsea"]] <- "Lango2010.txt.gz" # EUR LD-pruned SNPs
    deps[["geneloc_snpsea"]] <- c("NCBIgenes2013_hgnc.bed.gz", "NCBIgenes2013_entrez.bed.gz")
  }
  
  # Now, checking that all the requested (all/magma/ldsc/snpsea) files exist in an absolute path:
  i=1
  while (i <= length(deps)) {
    if (!all(file.exists(paste0(system.file(package="gwascelltyper"), "/extdata/", deps[[i]])))) {
      answer <- utils::askYesNo(paste0("The dependency files in category: ",
                                       names(deps)[i], ", are not downloaded or missing,",
                                       "\nand thus you either need to download the dependency data (Yes),",
                                       "\nor cancel this job for the moment (No/Cancel):"),
                                default = FALSE)
      if (isTRUE(answer)) {
        download_dependencies()
        i=0 # This will be increased to one at the end of the loop structure, so the loop resets (all files are checked again)
      } else {
        stop(paste0("ERROR: The dependency files in category: ", names(deps)[i], ", are missing for population ",pop,"."))
      }
    }
    i = i+1
  }
  
  # Processing some paths a bit before sending to the wrapper function:
  if (for_which %in% c("all", "ldsc")) {
    # deps$gene_coord is now deps$geneloc_ldsc
    deps$bim_prefix <- substr(deps[["1000g_bim_ldsc"]][1], 1,
                              nchar(deps[["1000g_bim_ldsc"]][1])-5)
    deps$print_snps <- substr(deps[["1000g_bim_ldsc"]][1], 1,
                              nchar(deps[["hm3snps_ldsc"]][1])-5)
    deps$ref_ld_chr <- substr(deps[["1000g_ldscore_category_ldsc"]][1], 1,
                              nchar(deps[["1000g_ldscore_category_ldsc"]][1])-15)
    deps$w_ld_chr <- substr(deps[["1000g_ldscore_universal_ldsc"]][1], 1,
                            nchar(deps[["1000g_ldscore_universal_ldsc"]][1])-15)
  }
  if (for_which %in% c("all", "magma")) {
    deps$genome_ref_path <- paste0("g1000_", pop)
    #deps$magma_path <- 
  }
  if (for_which %in% c("all", "snpsea")) {
  }
  
  
  
  # Now, return the named list of the dependency data
  deps <- lapply(deps, function(x) {return(paste0(system.file(package="gwascelltyper"), "/extdata/", x))})
  return(deps)
}

#' Download dependencies
#'
#' Downloads data dependencies (e.g. linkage disequilibrium, gene annotations, etc.) pre-formatted for each of the 3 packages (LDSC, MAGMA & SNPsea)
#'
#' @param address Path where dependencies are downloaded. Do not modify, parameter meant for developers.
#' @param pop Population to download LD data for (available right now: eur; available soon: eas, afr, amr, sas, subpop).
#'
#' @importFrom utils askYesNo
#' @importFrom utils download.file
#'
#' @return Path
#'
#' @export
download_dependencies <- function(address = paste0(system.file(package="gwascelltyper"), "/extdata/"),
                                  pop = c("eur", "eas", "afr", "amr", "sas", "subpop")[1]) {
  # Because some of these downloads may take quite a while to finish,
  # and R timeouts downloads after 60 seconds by default, we're recalibrating that parameter to 1 hour for users with slow internet.
  if (!dir.exists(address)) {dir.create(address)}
  options(timeout=60*60)
  #possible_pops = c("eur", "eas", "afr", "amr", "sas", "subpop") # For which there would be readily-available, pre-processed MAGMA LD reference datasets.
  #population <- tolower(readline(prompt="Enter which population (from: eur, eas, afr, amr, sas, subpop): "))
  #if (!population %in% possible_pops) {stop("Wrong population name. Choose from eur, eas, afr, amr, sas or subpop.")}
  
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
    utils::download.file(paste0("https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_", pop, ".zip"), destfile = paste0(address,"g1000_", pop,".zip"))
    utils::unzip(paste0(address,"g1000_", pop,".zip"), exdir = address)
    #file.remove(paste0(address,"g1000_", pop, ".zip"))
  }
  
  # MAGMA executable --------
  if (askYesNo("Do you want to download the MAGMA executable?")) {
    if (get_os() == "windows") {
      utils::download.file(url = "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_win.zip", destfile = paste0(address, "magma_v1.10_win.zip"))
      utils::unzip(paste0(address,"magma_v1.10_win.zip"), exdir = address)
    } else {
      if (get_os() == "linux") {
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
  required_files_for_snpsea <- c("Lango2010.txt.gz", "TGP2011.bed.gz", "NCBIgenes2013_hgnc.bed.gz", "NCBIgenes2013_entrez.bed.gz")
  if (sum(required_files_for_snpsea %in% list.files(address)) != 4) {
    if (askYesNo("Do you want to download the SNPsea dependencies for the eur population? (recommended yes)")) {
      cat("Downloading the SNPsea data.\n")
      utils::download.file(url = "http://files.figshare.com/1504037/SNPsea_data_20140520.zip",
                           destfile = paste0(address, "SNPsea_data_20140520.zip"))
      utils::unzip(zipfile = paste0(address, "SNPsea_data_20140520.zip"),
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
    
    
    # Download the file
    cat("Please wait while this large file downloads. It takes about 5 minutes on a fast internet connection, but can take much more on a weaker connection.\n")
    utils::download.file("https://s3.tebi.io/ldscdatazip/EUR_complete.zip",
                         destfile = paste0(address, "EUR_complete.zip"))
    utils::unzip(paste0(address, "EUR_complete.zip"),
                 exdir = address)
    
    
    # Population-agnostic files that you should download and process:
    utils::download.file(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz", destfile = paste0(address, "refGene.txt.gz"))
    refGene <- read.table(file = paste0(address, "refGene.txt.gz"))
    refGene <- refGene[,c(13,3,5,6)] # These are the only columns we need.
    refGene <- refGene[-which(duplicated(refGene[,1])),] # There are duplicate gene names due to isoform annotations. To keep it simple, we'll just remove duplicated gene IDs.
    refGene <- refGene[order(refGene[,2]),]
    colnames(refGene) <- c("GENE", "CHR", "START", "END")
    write.table(x = refGene, file = paste0(address, "refGene_coord.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
  
  # LDSC executable scripts & environments --------
  if (utils::askYesNo("Do you want to download the LDSC scripts?")) {
    system(paste0("cd ", address, " && git clone https://github.com/bulik/ldsc.git && mv ./ldsc/*",
                  "./ && conda env create --file environment.yml && source activate ldsc",
                  "&& conda deactivate"))
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
