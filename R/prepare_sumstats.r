#' Process summary statistics for MAGMA
#'
#' Process genome-wide association study summary statistics file into the required MAGMA format for analyses.
#'
#' @param sumstats_raw File containing the raw summary statistics file.
#' @param sumstats_munged File to output into. Default = paste0(sumstats_raw, "_processed.sumstats").
#' @param preview Show ideal MAGMA sumstats file then exit. If false (default), it proceeds to the processing.
#'
#' @importFrom utils data
#' @import data.table
#'
#' @return Path
#'
#' @export
munge_sumstats_for_MAGMA <- function(sumstats_raw, sumstats_munged = paste0(sumstats_raw,"_processed.sumstats"), preview = FALSE) {
  # Things left to do: use another CHR:BP-to-SNP matching rsID database than the one from MAGMA.Celltyping (which only has 2mil SNPs)
  if (preview) {
    cat("The perfect sumstats file for MAGMA should look something like this:

    SNP	CHR	BP	A1	A2	BETA	P
    rs3094315	1	752566	G	A	-0.0558	0.02774
    rs3131972	1	752721	A	G	0.0549	0.03149

    Where the columns are standardized.
    If the SNPs have variable sample size (e.g. sumstats come from meta-analysis), then that can be stored as a column called N, with integers.\n")
    return(NULL)
  }
  # This is a modified version (ugly) of the MAGMA.Celltyping package sumstats-munging functions (ugly as well),
  # that works no matter which version of sed & awk is used. It overwrites the input path, beware.
  if (is.null(sumstats_raw)) {
    stop("No raw sumstats file supplied in arguments.")
  }


  cat("This function can take some time if the summary statistics file is large (e.g. 10-15 mins if sumstats have 11+ million lines). Please have a bit of patience.\n")


  cat("Loading reference header names.\n")
  standard_sumstats_column_headers <- gwascelltyper::standard_sumstats_column_headers


  # Helper functions
  ask_genome_build <- function() {
    genomebuild <- as.numeric(readline("Which genome build is the data from? 1 for GRCh37, 2 for GRCh38: "))
    if (!genomebuild %in% c(1, 2)) stop("Genome build must be entered as either 1 (for GRCh37) or 2 (for GRCh38)")
    if (genomebuild == 1) genomebuild <- "GRCh37"
    else genomebuild <- "GRCh38"
    return(genomebuild)
  }


  cat("Reading summary statistics file.\n")
  sumstats <- fread(sumstats_raw)
  cat("Standardizing header.\n")
  colnames(sumstats) <- toupper(colnames(sumstats))
  for (i in 1:ncol(sumstats)) {
    if (colnames(sumstats)[i] %in% standard_sumstats_column_headers$Uncorrected) {
      colnames(sumstats)[i] <- standard_sumstats_column_headers$Corrected[which(standard_sumstats_column_headers$Uncorrected == colnames(sumstats)[i])]
    }
  }


  # Check if key headers are present
  if (!sum(c("A1", "A2") %in% colnames(sumstats)) == 2) stop("Either or both A1 & A2 columns missing from the sumstats header.")
  if (!sum(c("CHR", "BP") %in% colnames(sumstats)) == 2) stop("Either or both CHR & BP columns missing from the sumstats header.")
  if (sum(c("Z", "OR", "BETA", "LOG_ODDS", "SIGNED_SUMSTAT") %in% colnames(sumstats)) < 1) stop(paste0("Error: cannot find a column name representing signed statistic in GWAS sumstats file (Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT). Header of file:\n", colnames(sumstats), "\n"))


  # Make sure the CHR column contents are properly formatted
  sumstats$CHR = as.factor(sumstats$CHR)
  if (length(grep("chr", sumstats$CHR[1])) != 0) sumstats$CHR = gsub("chr", "", sumstats$CHR)


  # Make sure that there's an rsID column (SNP) and make one if there isn't any.
  if (sum(c("CHR", "BP") %in% colnames(sumstats)) == 2 & sum("SNP" %in% colnames(sumstats)) == 0) {
    cat("There is no SNP column found within the data. It must be inferred from CHR and BP information. Note: this process drops any SNPs which are not from Hapmap.\nInitial number of rows:", nrow(sumstats), "\n")

    # Load genome builds and ask which one of these the SNPs are from
    snp_loc <- gwascelltyper::snp_loc
    genomebuild <- ask_genome_build()

    # SNP matching
    sumstats = merge(sumstats, snp_loc[snp_loc$Build == genomebuild,][, -4], by = c("CHR", "BP"))
    cat("Number of rows remaining:", nrow(sumstats), ".\n")
  } else {
    cat("Checking if all of the entries in the SNP column have rs IDs.\n")
    no_rsids <- union(which(substr(sumstats$SNP, 1, 2) != "rs"), which(grepl("\\D", substr(sumstats$SNP, 3, nchar(sumstats$SNP)))))

    sumstats_norsid <- sumstats[no_rsids,]
    if (nrow(sumstats_norsid)>0) {
      cat("Number of rows that don't have an rsID in the SNP column:", nrow(sumstats_norsid), "out of", nrow(sumstats),".\nAttempting to get SNP rsIDs for those.\n")
      sumstats <- sumstats[-no_rsids,]

      # Load genome builds and ask which one of these the SNPs are from
      snp_loc <- gwascelltyper::snp_loc
      genomebuild <- ask_genome_build()

      sumstats_norsid <- sumstats_norsid[, !"SNP"]
      sumstats_norsid = merge(sumstats_norsid, snp_loc[snp_loc$Build == genomebuild,][, -4], by = c("CHR", "BP")) # SNP matching
      cat("Number of rows remaining (for which we managed to get rsIDs out of the ones that didn't have one):", nrow(sumstats_norsid), ".\n")

      # Make sure that the colnames are matched when column-binding sumstats and sumstats_norsid
      setcolorder(sumstats_norsid, colnames(sumstats))

      # Merge them
      sumstats <- rbind(sumstats, sumstats_norsid)
    } else {
      cat("All SNPs have rsIDs.\n")
    }
  }

  if (!sum(colnames(sumstats)[1:3] == c("SNP", "CHR", "BP")) == 3) {
    cat("Reordering columns so that SNP, CHR & BP are the first 3, in this order.\n")
    setcolorder(sumstats, c("SNP", "CHR", "BP", setdiff(colnames(sumstats), c("SNP", "CHR", "BP"))))
  }

  cat("Writing processed summary statistics.\n")
  fwrite(sumstats, file = sumstats_munged, quote = FALSE, sep = "\t", showProgress = TRUE)

}

#' Process summary statistics for LDSC
#'
#' Process genome-wide association study summary statistics file into the required LDSC format (given MAGMA-formatted sumstats) for analyses.
#'
#' @param sumstats_magma File containing the raw summary statistics file.
#' @param sumstats_munged File to output into. Default = paste0(sumstats_raw, "_processed.sumstats").
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param gwas_sample_size If GWAS sample size is not a constant integer for all SNPs, then leave NULL. Otherwise, must be an integer.
#' @param merge_alleles Leave NULL if you do not want to merge the SNPs with those in HapMap3 (as recommended in LDSC guidelines). Otherwise, should be the path of the w_hm3.snplist file.
#' @param ldsc_path Path to LDSC python script (default is extdata folder in package path).
#' @param preview Show ideal LDSC sumstats file then exit. If false (default), it proceeds to the processing.
#'
#' @return Path
#'
#' @importFrom R.utils compressFile
#' @importFrom utils read.table
#' @importFrom utils write.table
#'
#' @export
munge_sumstats_for_LDSC_from_MAGMA <- function(sumstats_magma = NULL, sumstats_munged = paste0(sumstats_magma,"_processed.sumstats.gz"),
                                               number_of_threads = 1, gwas_sample_size = NULL,
                                               merge_alleles = paste0(system.file(package="gwascelltyper"), "/extdata/w_hm3.snplist"),
                                               ldsc_path = paste0(system.file(package="gwascelltyper"), "/extdata/"),
                                               preview = FALSE) {
  if (preview) {
    cat("The perfect sumstats file for LDSC should look something like this:

    SNP	A1	A2	Z	N
    rs3094315	G	A	-2.201	13974
    rs3131972	A	G	2.151	13974

    Where the columns are standardized.
    The N is the SNP sample size column. The Z column contains the Z score of the SNP. A1 and A2 refer to allele 1 and allele 2.\n")
    return(NULL)
  }
  if ((!is.null(merge_alleles)) && (!file.exists(merge_alleles))) stop("Error: Supplied merge_alleles param is has to be null or it has to point to a file that exists.")
  if (is.null(sumstats_magma)) stop("No magma sumstats file supplied in arguments.")
  if (is.null(ldsc_path)) stop("No path towards LDSC's scripts provided.")


  intermed_file <- tempfile(tmpdir = dirname(sumstats_magma), pattern = "intermed_", fileext = "_file.gwas")
  file.copy(from = sumstats_magma, to = intermed_file, overwrite = TRUE)


  command <- ""
  if (is.null(merge_alleles)) { # See why here: https://github.com/bulik/ldsc/issues/146
    if (is.null(gwas_sample_size)) {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads,
                        " python ", ldsc_path, "munge_sumstats.py --sumstats ", intermed_file, " --merge-alleles ", merge_alleles,
                        " --out ", intermed_file)
    } else {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads,
                        " python ", ldsc_path, "munge_sumstats.py --sumstats ", intermed_file,
                        " --out ", intermed_file, " --N ", gwas_sample_size)
    }
  } else {
    if (is.null(gwas_sample_size)) {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads,
                        " python ", ldsc_path, "munge_sumstats.py --sumstats ", intermed_file, " --merge-alleles ", merge_alleles,
                        " --out ", intermed_file)
    } else {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads,
                        " python ", ldsc_path, "munge_sumstats.py --sumstats ", intermed_file, " --merge-alleles ", merge_alleles,
                        " --out ", intermed_file, " --N ", gwas_sample_size)
    }
  }
  print(command)
  system(command = command)


  if (!is.null(merge_alleles)) {
    cat("Often, the --merge_alleles parameter of the munge_sumstats.py script of LDSC can create missingness in
    non-rsID-columns (see https://github.com/bulik/ldsc/issues/170#issue-507105766), like this:
    SNP	A1	A2	Z	N\nrs3094315\nrs3131972\nrs3131969
    This following operation will eliminate this minor problem by taking out such missing-beta rows.\n")
    sumstats <- read.table(file = paste0(intermed_file, ".sumstats.gz"), header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
    sumstats <- sumstats[which(sumstats[,2]!=""),]
    sumstats[,which(colnames(sumstats)=="N")] <- as.integer(sumstats[,which(colnames(sumstats)=="N")])
    write.table(file = paste0(intermed_file, ".sumstats"),
                x=sumstats, quote = FALSE,row.names = FALSE, sep = "\t")
    file.remove(paste0(intermed_file, ".sumstats.gz"))
    compressFile(filename = paste0(intermed_file, ".sumstats"),
                 destname = paste0(intermed_file, ".sumstats.gz"),
                 ext = "gz", FUN=gzfile)
  }


  file.rename(from = paste0(intermed_file, ".sumstats.gz"), to = sumstats_munged)
  file.remove(intermed_file) # Deleting the intermediary MAGMA-munged sumstats
  return(sumstats_munged)
}

#' Process summary statistics for SNPsea
#'
#' Process genome-wide association study summary statistics file into the required SNPsea format (given MAGMA-formatted sumstats) for analyses. Leaves only the genome-wide significant SNPs (because that is the SNPsea approach).
#'
#' @param sumstats_magma File containing GWAS summary statistics file properly formatted for MAGMA.
#' @param sumstats_munged File to output SNPsea-format into. Default = paste0(sumstats_magma, "_processed.sumstats").
#' @param preview Show ideal SNPsea sumstats file then exit. If false (default), it proceeds to the processing.
#'
#' @return Path
#' 
#' @import data.table
#'
#' @export
munge_sumstats_for_SNPSEA_from_MAGMA <- function(sumstats_magma = NULL, sumstats_munged = paste0(sumstats_magma,"_processed.sumstats"), preview = FALSE) {
  if (preview) {
    cat("The perfect summary statistics for SNPsea should look like below:

    CHR     POS     SNP     P
    chr1    40069939        rs3916164       3e-10
    chr1    158575729       rs857684        4e-16

    Where the P-value column needs to be in the Ne-NN format, where N is a digit.\n")
    return(NULL)
  }
  file.copy(from = sumstats_magma, to = sumstats_munged)
  cat("Munging sumstats function for SNPSEA is still in development. Guideline for manual formatting below.
  # head of Red_blood_cell_count-Harst2012-45_SNPs.gwas
  # Harst et al. 2012
  # doi:10.1038/nature11677
  # PMID: 23222517
  # 45 SNPs associated with red blood cell count (RBC) taken from Table 1.
  # Positions are on hg19. SNPs are included if P <= 5e-8.
  CHR     POS     SNP     P
  chr1    40069939        rs3916164       3e-10
  chr1    158575729       rs857684        4e-16
  chr1    199007208       rs7529925       8e-09
  chr1    248039451       rs3811444       5e-10\n\n")
  # The part where you read the file
  df <- read.table(file = sumstats_munged, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[,c("CHR", "BP", "SNP", "P")]
  df$CHR <- paste0("chr", df$CHR)
  colnames(df) <- c("CHR","POS","SNP", "P")
  cat("Eliminating", nrow(df[which(as.numeric(df$P) > 5e-8),]), "SNPs from", nrow(df),"total (based on filter: Pval <= 5e-8), making the final sumstat", nrow(df[which(as.numeric(df$P) <= 5e-8),]), "SNPs long.\n")
  df <- df[which(as.numeric(df$P) <= 5e-8),]
  df$P <- formatC(as.numeric(df$P), format = "e", digits = 0)
  write.table(df, file = sumstats_munged, row.names=FALSE, sep="\t", quote=FALSE)
  return(sumstats_munged)
}
