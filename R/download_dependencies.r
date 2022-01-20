#' Download baseline data
#'
#' Downloads baseline data (e.g. LD regions & others) pre-formatted for each of the 3 packages (LDSC, MAGMA & SNPsea)
#'
#' @param population Population to download reference data for (fast downloads available for afr, amr, eas, eur, sas).
#' @param address Path where dependencies are downloaded. Do not modify, parameter meant for developers.
#'
#' @importFrom utils askYesNo
#' @importFrom utils download.file
#'
#' @return Path
#'
#' @export
download_dependencies <- function(population = "eur", address = paste0(system.file(package="gwascelltyper"), "/extdata/")) {
  require(data.table)
  require(R.utils)
  if (!population %in% c("afr", "amr", "eas", "eur", "sas", "subpop")) {
    stop("Population options: afr amr eas eur sas subpop.", population, "is not a publicly available choice for said packages yet. Please consult the MAGMA, LDSC & SNPsea websites to create your own LD reference data.")
  } else {
    cat("Population options: afr amr eas eur sas subpop. Downloading the", population, "population.\n")
  }


  # LDSC
  if (utils::askYesNo("Do you want to download the LDSC scripts?")) {
    system(paste0("cd ", address, " && git clone https://github.com/bulik/ldsc.git && mv ./ldsc/* ./ && conda env create --file environment.yml && source activate ldsc && conda deactivate"))
  }

  if (askYesNo(paste0("Do you want to download the LDSC dependencies for the ", population, " population?"))) {
    if (!file.exists(paste0(address, "refGene_coord.txt"))) {
      cat("Downloading refGene annotations from the UCSC server.\n")
      utils::download.file(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz", destfile = paste0(address, "refGene.txt.gz"))
      refGene <- read.table(file = paste0(address, "refGene.txt.gz"))
      refGene <- refGene[,c(13,3,5,6)] # These are the only columns we need.
      refGene <- refGene[-which(duplicated(refGene[,1])),] # There are duplicate gene names due to isoform annotations. To keep it simple, we'll just remove duplicated gene IDs.
      refGene <- refGene[order(refGene[,2]),]
      colnames(refGene) <- c("GENE", "CHR", "START", "END")
      write.table(x = refGene, file = paste0(address, "refGene_coord.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
      file.remove(paste0(address, "refGene.txt.gz"))
    }
    if (!dir.exists(paste0(address, "hapmap3_snps"))) {
      cat("Downloading LDSC's HapMap3 SNP data.\n")
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz",
                    destfile = paste0(address,"hapmap3_snps.tgz"))
      system(paste0("cd ", address," && tar -xzvf hapmap3_snps.tgz && rm hapmap3_snps.tgz"))
    }
    if (!file.exists(paste0(address, "w_hm3.snplist"))) {
      cat("Downloading another HapMap3 dependency for formatting sumstats.\n")
      utils::download.file(url = "https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2", destfile = paste0(address, "w_hm3.snplist.bz2"))
      cat("Decompressing the file.\n")
      write.table(x = read.table(file = paste0(address, "w_hm3.snplist.bz2"), header = TRUE), quote = FALSE, row.names = FALSE, file = paste0(address, "w_hm3.snplist"))
      file.remove(paste0(address, "w_hm3.snplist.bz2"))
    }
    if (population == "eur") {
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz",
                    destfile = paste0(address,"1000G_Phase3_baseline_ldscores.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz",
                    destfile = paste0(address,"weights_hm3_no_hla.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz",
                    destfile = paste0(address,"1000G_Phase3_plinkfiles.tgz"))
      system(paste0("cd ", address , " && tar -xzvf ./1000G_Phase3_baseline_ldscores.tgz && tar -xzvf ./weights_hm3_no_hla.tgz && tar -xzvf ./1000G_Phase3_plinkfiles.tgz"))
      # Renaming some of the files for proper standard input and deleting temp archives
      file.rename(from = paste0(address,"weights_hm3_no_hla"), to = paste0(address,"1000G_EUR_Phase3_weights_hm3_no_hla"))
      file.remove(paste0(address,"1000G_Phase3_baseline_ldscores.tgz"))
      file.remove(paste0(address,"1000G_Phase3_plinkfiles.tgz"))
      file.remove(paste0(address,"weights_hm3_no_hla.tgz"))
    } else if (population == "eas") {
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz",
                    destfile = paste0(address,"1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz",
                    destfile = paste0(address,"1000G_Phase3_EAS_weights_hm3_no_MHC.tgz"))
      utils::download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_plinkfiles.tgz",
                    destfile = paste0(address,"1000G_Phase3_EAS_plinkfiles.tgz"))
      system(paste0("cd ", address, " && tar -xvzf 1000G_Phase3_EAS_plinkfiles.tgz && rm 1000G_Phase3_EAS_plinkfiles.tgz && mv 1000G_Phase3_EAS_plinkfiles 1000G_EAS_Phase3_plink && mkdir ./1000G_EAS_Phase3_baseline && tar -C ./1000G_EAS_Phase3_baseline/ -xvzf 1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz && rm 1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz && tar -xvzf 1000G_Phase3_EAS_weights_hm3_no_MHC.tgz && rm 1000G_Phase3_EAS_weights_hm3_no_MHC.tgz && mv 1000G_Phase3_EAS_weights_hm3_no_MHC 1000G_EAS_Phase3_weights_hm3_no_hla"))
      # Renaming some of the files for proper standard input
      for (i in list.files(paste0(address, "1000G_Phase3_eas_weights_hm3_no_hla"), full.names = TRUE)) {
        file.rename(from = i, to = paste0(dirname(i), "/", to = gsub(pattern = "weights.EAS.hm3_noMHC", replacement = "weights", basename(i))))
      }
    } else {
      cat("LDSC does not have pre-processed downloadable LD scores on the public server for afr, amr, sas and subpop, but instructions are provided in the wiki (https://github.com/bulik/ldsc/wiki/) for how to make such files. Please download eur or eas to the data folder (to see how they're formatted and placed), and repeat this for your own LD data.\n")
    }
  }

  # MAGMA
  if (askYesNo(paste0("Do you want to download the MAGMA dependencies for the ", population, " population?"))) {
    utils::download.file(paste0("https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_", population, ".zip"), destfile = paste0(address,"g1000_", population,".zip"))
    utils::unzip(paste0(address,"g1000_", population,".zip"), exdir = address)
    file.remove(paste0(address,"g1000_", population, ".zip"))
  }


  cat("Checking if SNPsea has all of its dependencies downloaded.\n")
  required_files_for_snpsea <- c("Lango2010.txt.gz", "TGP2011.bed.gz", "NCBIgenes2013.bed.gz")
  if (sum(required_files_for_snpsea %in% list.files(address)) != 3) {
    if (askYesNo("Do you want to download the SNPsea dependencies for the eur population?")) {
      cat("Downloading the SNPsea data.\n")
      system(paste0("cd ", address, " && curl -LOk http://files.figshare.com/1504037/SNPsea_data_20140520.zip"))
      system(paste0("cd ", address, " && unzip SNPsea_data_20140520.zip && rm SNPsea_data_20140520.zip"))
      x <- read.table(file = paste0(address, "NCBIgenes2013.bed.gz"))[,c(1,2,3,5,4)]
      xfile <- gzfile(paste0(address, "NCBIgenes2013_hgnc.bed.gz"), "w")
      write.table(x = x, file = xfile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      close(xfile)
      rm(xfile)
      file.rename(from = paste0(address, "NCBIgenes2013.bed.gz"), to = paste0(address, "NCBIgenes2013_entrez.bed.gz"))
    }
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
