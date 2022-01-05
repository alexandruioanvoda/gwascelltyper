#' Compute LDSC enrichment - optimized
#'
#' Compute GWAS enrichment in discrete lists of genes with LDSC. Optimized to load all dependency files into Python-allocated RAM only once.
#'
#' @param gene_sets List of genesets. Each geneset must be a character vector, containing HGNC gene names.
#' @param sumstats_file Address of properly formatted LDSC GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param population Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param ldsc_path Path to folder containing all LDSC scripts.
#' @param gene_coord Adress of UCSC bed file containing all the gene coordinates.
#' @param window_size Size of window surrounding gene, for SNP-to-gene mapping purposes. Default = 100000 (100kb).
#' @param bim_prefix Path + prefix of plink BIM files required by LDSC to compute LD scores.
#' @param print_snps Path of file containing HapMap3 SNPs.
#' @param ref_ld_chr Path of file containing LD reference panel. Tell LDSC which LD Scores to use as the predictors in the LD Score regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. LDSC will automatically concatenate .l2.ldscore files split across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. If the filename prefix contains the symbol (at), LDSC will replace the (at) symbol with chromosome numbers. Otherwise, LDSC will append chromosome numbers to the end of the filename prefix. Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz Example 2: --ref-ld-chr ld/(at)_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz
#' @param w_ld_chr Filename prefix for file with LD Scores with sum r^2 taken over SNPs included in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz. LDSC will read files split into 22 chromosomes in the same manner as --ref-ld-chr.
#' @param ld_wind_cm Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM). Default = 1.
#'
#' @return Dataframe
#'
#' @export
compute_LDSC_enrichment_discrete_optimized <- function(gene_sets, sumstats_file, output_dir = NULL, number_of_threads = 1, window_size = 100000, population = "eur", ld_wind_cm = 1,
                                                       ldsc_path = paste0(system.file(package="gwascelltyper"), "/extdata/"),
                                                       gene_coord = paste0(system.file(package="gwascelltyper"), "/extdata/refGene_coord.txt"),
                                                       print_snps = paste0(system.file(package="gwascelltyper"), "/extdata/hapmap3_snps/hm."),
                                                       bim_prefix = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_plink/1000G.",toupper(population),".QC."),
                                                       ref_ld_chr = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_baseline/baseline."),
                                                       w_ld_chr = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_weights_hm3_no_hla/weights.")) {
  # Libraries
  require(doParallel)
  require(data.table)
  require(R.utils)
  # Checks
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("Sumstats file does not exist at provided path.")}
  if (is.null(ldsc_path)) {stop("What LDSC path?")}
  if (is.null(gene_coord)) {stop("What gene_coord param?")}
  if (is.null(bim_prefix)) {stop("What bim_prefix param?")}
  if (is.null(print_snps)) {stop("What print_snps param?")}
  if (is.null(ref_ld_chr)) {stop("What ref_ld_chr param?")}
  if (is.null(w_ld_chr)) {stop("What w_ld_chr param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  file_dependencies <- c(paste0(ldsc_path, c("ldsc.py", "make_annot.py")),
                         gene_coord,
                         paste0(print_snps, 1:22, ".snp"),
                         paste0(bim_prefix, do.call(paste0, expand.grid(paste0(1:22, "."), c("bed", "fim", "fam")))),
                         c(paste0(ref_ld_chr, do.call(paste0, expand.grid(paste0(1:22, "."), c("annot.gz", "l2.ldscore.gz", "l2.M", "l2.M_5_50")))),
                           paste0(gsub(pattern = "/baseline.", replacement = "/", x = ref_ld_chr, fixed = TRUE), "print_snps.txt")),
                         paste0(w_ld_chr, 1:22, ".l2.ldscore.gz"))
  cat("Checking LDSC file dependencies.\n")
  if (any(!file.exists(file_dependencies))) {
    cat("\n\n", file_dependencies[!file.exists(file_dependencies)], "\n")
    stop("The files above were passed as arguments but do not exist.")
  } else {
    rm(file_dependencies)
  }
  # A few other boilerplate checks:
  window_size <- format(window_size, scientific = FALSE)
  if (names(gene_sets)[length(gene_sets)] != "All_genes_control") stop("No <All_genes_control> contained within the last geneset from the gene_sets object.")
  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = names(gene_sets)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    names(gene_sets) <- format_celltype_names(names(gene_sets))
  }
  names_of_gene_sets <- names(gene_sets) # Retaining this in memory as the gene_sets object will get modified into a df.


  cat("LDSC enrichment analysis started at:", as.character(Sys.time()), "\n")


  # Boiler plate code in case user wants to store intermediary and output files somewhere permanent instead of the default path (default = a temporary directory, that is erased after the run).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for LDSC files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }


  cat("Creating gene-set file to make the thin-annot out of.\n")
  geneset_file <- paste0(output_dir, "/geneset_file.tsv")
  n <- max(lengths(gene_sets))
  for (i in 1:length(gene_sets)) {
    length(gene_sets[[i]]) <- n
  }
  gene_sets <- do.call(cbind, gene_sets)
  write.table(gene_sets, file = geneset_file, quote = FALSE, sep = "\t", na = "", row.names = FALSE, col.names = FALSE)


  cat("Making a thin-annot of all gene-sets (one per column), per each chromosome.\n")
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads," MKL_NUM_THREADS=", number_of_threads,
                    " python ", ldsc_path, "make_annot_optimized.py --gene-set-file ", geneset_file,
                    " --gene-coord-file ", gene_coord,
                    " --windowsize ", window_size,
                    " --bimfile ", bim_prefix,
                    " --annot-file ", output_dir, "/out.")
  system(command = command)
  cat("Thin annot making finished at:", as.character(Sys.time()), "\n")


  cat("LDSC consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.\n")
  gc()
  cat("Initializing R multithreading.\n")
  registerDoParallel(cores=number_of_threads)


  foreach(i=1:22, .options.snow=list(preschedule=TRUE)) %do% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", ldsc_path, "ldsc.py --l2 --bfile ",
                      bim_prefix, as.character(i), " --ld-wind-cm ", ld_wind_cm,
                      " --annot ", output_dir, "/out.chr", as.character(i), ".annot.gz",
                      " --thin-annot --out ", output_dir, "/merged_ld_scores.", as.character(i),
                      " --print-snps ", print_snps, as.character(i), ".snp")
    system(command = command)

    cat("Reading in the LD score, M & M_5_50 files.\n")
    merged_ld <- read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.ldscore.gz"), header = TRUE)
    merged_m <- read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.M"), header = FALSE)
    merged_m_5_50 <- read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.M_5_50"), header = FALSE)

    cat("Separating the LD score, M & M_5_50 files.\n")
    for (j in 1:(ncol(merged_ld)-3)) {
      de_merged <- merged_ld[,c(1,2,3,j+3)]
      write.table(x = de_merged, file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
      gzip(paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore"), destname=paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore.gz"), overwrite=TRUE, remove=TRUE)
      write.table(x = merged_m[,j], file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.M"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(x = merged_m_5_50[,j], file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.M_5_50"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }

  cat("LD scoring finished at:", as.character(Sys.time()), "\n")

  registerDoSEQ() # This finishes the parallel backend
  print("Calling garbage collection now, to free memory after the LD/geneset parallel operations done by DoParallel.")
  gc()


  print("Writing LDCTS file, line for each cell type")
  suppressWarnings(file.remove(paste0(output_dir, "/main.ldcts")))
  for (i in 1:(length(names_of_gene_sets)-1)) {
    ldcts_line <- paste(gsub(" ", "_", names_of_gene_sets[i]), # This is the celltype name
                        paste0(output_dir, "/geneset.", as.character(i), ".,", output_dir, "/geneset.", length(names_of_gene_sets), "."), sep = "\t") # This is the files for the cell-type and the control gene-set
    write(ldcts_line, file=paste0(output_dir, "/main.ldcts"), append=TRUE)
  }


  # This command actually returns the P-values for each celltype against the GWAS trait
  print("This command actually returns the P-values for each cell type against the GWAS trait")
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads, " python ", ldsc_path, "ldsc.py --h2-cts ", sumstats_file,
                    " --ref-ld-chr ", ref_ld_chr, " ",
                    "--out ", output_dir, "/ldsc_final_output ",
                    "--ref-ld-chr-cts ", output_dir, "/main.ldcts ",
                    "--w-ld-chr ", w_ld_chr)
  print(command)
  system(command)
  print("Done. Returning results dataframe.")


  # Here we build the results dataframe for plotting.
  results <- read.table(paste0(output_dir, "/ldsc_final_output.cell_type_results.txt"), header = T, stringsAsFactors = FALSE)
  colnames(results)[4] <- "P_value"


  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)

  cat("LDSC enrichment analysis finished at:", as.character(Sys.time()), "\n")

  return(results)
}


#' Compute LDSC enrichment
#'
#' Compute GWAS enrichment in discrete lists of genes with LDSC.
#'
#' @param gene_sets List of genesets. Each geneset must be a character vector, containing HGNC gene names.
#' @param sumstats_file Address of properly formatted LDSC GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param population Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param ldsc_path Path to folder containing all LDSC scripts.
#' @param gene_coord Adress of UCSC bed file containing all the gene coordinates.
#' @param window_size Size of window surrounding gene, for SNP-to-gene mapping purposes. Default = 100000 (100kb).
#' @param bim_prefix Path + prefix of plink BIM files required by LDSC to compute LD scores.
#' @param print_snps Path of file containing HapMap3 SNPs.
#' @param ref_ld_chr Path of file containing LD reference panel. Tell LDSC which LD Scores to use as the predictors in the LD Score regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. LDSC will automatically concatenate .l2.ldscore files split across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. If the filename prefix contains the symbol (at), LDSC will replace the (at) symbol with chromosome numbers. Otherwise, LDSC will append chromosome numbers to the end of the filename prefix.Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz Example 2: --ref-ld-chr ld/(at)_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz
#' @param w_ld_chr Filename prefix for file with LD Scores with sum r^2 taken over SNPs included in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz. LDSC will read files split into 22 chromosomes in the same manner as --ref-ld-chr.
#' @param ld_wind_cm Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM). Default = 1.
#'
#' @return Dataframe
#'
#' @export
compute_LDSC_enrichment_discrete <- function(gene_sets, sumstats_file, output_dir = NULL, number_of_threads = 1, window_size = 100000, population = "eur", ld_wind_cm = 1,
                                             ldsc_path = paste0(system.file(package="gwascelltyper"), "/extdata/"),
                                             gene_coord = paste0(system.file(package="gwascelltyper"), "/extdata/refGene_coord.txt"),
                                             print_snps = paste0(system.file(package="gwascelltyper"), "/extdata/hapmap3_snps/hm."),
                                             bim_prefix = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_plink/1000G.",toupper(population),".QC."),
                                             ref_ld_chr = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_baseline/baseline."),
                                             w_ld_chr = paste0(system.file(package="gwascelltyper"), "/extdata/1000G_",toupper(population),"_Phase3_weights_hm3_no_hla/weights.")) {
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("Sumstats file does not exist at provided path.")}
  if (is.null(ldsc_path)) {stop("What LDSC path?")}
  if (is.null(gene_coord)) {stop("What gene_coord param?")}
  if (is.null(bim_prefix)) {stop("What bim_prefix param?")}
  if (is.null(print_snps)) {stop("What print_snps param?")}
  if (is.null(ref_ld_chr)) {stop("What ref_ld_chr param?")}
  if (is.null(w_ld_chr)) {stop("What w_ld_chr param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  file_dependencies <- c(paste0(ldsc_path, c("ldsc.py", "make_annot.py")),
                         gene_coord,
                         paste0(print_snps, 1:22, ".snp"),
                         paste0(bim_prefix, do.call(paste0, expand.grid(paste0(1:22, "."), c("bed", "fim", "fam")))),
                         c(paste0(ref_ld_chr, do.call(paste0, expand.grid(paste0(1:22, "."), c("annot.gz", "l2.ldscore.gz", "l2.M", "l2.M_5_50")))),
                           paste0(gsub(pattern = "/baseline.", replacement = "/", x = ref_ld_chr, fixed = TRUE), "print_snps.txt")),
                         paste0(w_ld_chr, 1:22, ".l2.ldscore.gz"))
  cat("Checking LDSC file dependencies.\n")
  if (any(!file.exists(file_dependencies))) {
    cat("\n\n", file_dependencies[!file.exists(file_dependencies)], "\n")
    stop("The files above were passed as arguments but do not exist.")
  } else {
    rm(file_dependencies)
  }
  # A few other boilerplate checks:
  window_size <- format(window_size, scientific = FALSE)
  if (names(gene_sets)[length(gene_sets)] != "All_genes_control") stop("No <All_genes_control> contained within the last geneset from the gene_sets object.")
  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = names(gene_sets)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    names(gene_sets) <- format_celltype_names(names(gene_sets))
  }


  cat("LDSC enrichment analysis started at:", as.character(Sys.time()), "\n")


  # Boiler plate code in case user wants to store intermediary and output files somewhere permanent instead of the default path (default = a temporary directory, that is erased after the run).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for LDSC files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }


  print("LDSC consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.")
  gc()
  require(doParallel)
  registerDoParallel(cores=number_of_threads)


  print("Writing lists of genes to files and making thin-annots for each gene set...")
  foreach(i=1:(length(gene_sets)-1), .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    geneset_file <- paste0(output_dir, "/genes_", gsub(" ", "_", names(gene_sets)[i]), "_list.txt") # The gsub command is the cell name
    con <- file(geneset_file)
    writeLines(gene_sets[[i]], con)
    close(con)
    for (j in seq(from = 1, to = 22, by = 1)) { # Compute annot file / each chromosome
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", ldsc_path, "make_annot.py --gene-set-file ",
                        geneset_file, " --gene-coord-file ", gene_coord,
                        " --windowsize ", window_size,
                        " --bimfile ", bim_prefix, as.character(j), ".bim ",
                        "--annot-file ", output_dir, "/celltype.", as.character(i), ".", as.character(j), ".annot.gz")
      print(command)
      system(command = command)
    }
  }

  print("Calling garbage collection now, to free memory after the thin-annot parallel operations done by DoParallel.")
  gc()

  print("Write list of all genes (for control)")
  geneset_file <- paste0(output_dir, "/genes_controls_list.txt")
  con <- file(geneset_file)
  writeLines(gene_sets[[length(gene_sets)]], con) # gene_sets[[length(gene_sets)]] is the control set (vector with all genes)
  close(con)
  # Make thin-annot for control too
  foreach(j=1:22, .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", ldsc_path, "make_annot.py --gene-set-file ",
                      geneset_file, " --gene-coord-file ", gene_coord,
                      " --windowsize ", window_size,
                      " --bimfile ", bim_prefix, as.character(j), ".bim ",
                      "--annot-file ", output_dir, "/control.", as.character(j), ".annot.gz")
    print(command)
    system(command = command)
  }


  print("Calling garbage collection now, to free memory after the <control> thin-annot parallel operations done by DoParallel.")
  gc()

  cat("Thin annot making finished at:", as.character(Sys.time()), "\n")

  print("Writing LDCTS file, line for each cell type")
  for (i in 1:(length(gene_sets)-1)) {
    # Get cell type name
    cell <- names(gene_sets)[i]; cell <- gsub(" ", "_", cell); print(cell)
    # Write line to LDCTS
    ldcts_line <- paste(cell, paste0(output_dir, "/celltype.", as.character(i), ".,", output_dir, "/control."), sep = "\t")
    write(ldcts_line, file=paste0(output_dir, "/main.ldcts"), append=TRUE)
  }


  print("Computing LD for each gene-set / each chromosome.")
  foreach(i=1:(length(gene_sets)-1), .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    for (j in seq(from = 1, to = 22, by = 1)) {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", ldsc_path, "ldsc.py --l2 --bfile ",
                        bim_prefix, as.character(j), " --ld-wind-cm ", ld_wind_cm,
                        " --annot ", output_dir, "/celltype.", as.character(i), ".", as.character(j), ".annot.gz",
                        " --thin-annot --out ", output_dir, "/celltype.", as.character(i), ".", as.character(j),
                        " --print-snps ", print_snps, as.character(j), ".snp")
      print(command)
      system(command = command)
    }
  }


  print("Calling garbage collection now, to free memory after the LD/geneset parallel operations done by DoParallel.")
  gc()


  # Compute LD for the control too
  print("Compute LD for the control too")
  foreach(j=1:22, .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", ldsc_path, "ldsc.py --l2 --bfile ",
                      bim_prefix, as.character(j), " --ld-wind-cm ", ld_wind_cm,
                      " --annot ", output_dir, "/control.", as.character(j), ".annot.gz",
                      " --thin-annot --out ", output_dir, "/control.", as.character(j),
                      " --print-snps ", print_snps, as.character(j), ".snp")
    print(command)
    system(command = command)
  }


  print("Calling garbage collection now, to free memory after the LD/control parallel operations done by DoParallel.")
  gc()

  cat("LD scoring finished at:", as.character(Sys.time()), "\n")

  registerDoSEQ() # This finishes the parallel backend


  gc()


  # This command actually returns the P-values for each celltype against the GWAS trait
  print("This command actually returns the P-values for each cell type against the GWAS trait")
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads, " python ", ldsc_path, "ldsc.py --h2-cts ", sumstats_file,
                    " --ref-ld-chr ", ref_ld_chr, " ",
                    "--out ", output_dir, "/ldsc_final_output ",
                    "--ref-ld-chr-cts ", output_dir, "/main.ldcts ",
                    "--w-ld-chr ", w_ld_chr)
  print(command)
  system(command)
  print("Done. Returning results dataframe.")


  # Here we build the results dataframe for plotting.
  results <- read.table(paste0(output_dir, "/ldsc_final_output.cell_type_results.txt"), header = T, stringsAsFactors = FALSE)
  colnames(results)[4] <- "P_value"


  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)


  return(results)
}


#' Compute MAGMA enrichment
#'
#' Compute GWAS enrichment in discrete lists of genes with MAGMA.
#'
#' @param gene_sets List of genesets. Each geneset must be a character vector, containing HGNC gene names.
#' @param sumstats_file Address of properly formatted MAGMA GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param upstream_kb Number of kb upstream of gene, SNP-to-gene mapping window. Default = 10.
#' @param downstream_kb Number of kb downstream of gene, SNP-to-gene mapping window. Default = 1.5.
#' @param gene_nomenclature Type of gene IDs fed in from the expression data (default = hgnc; alternatives = entrez).
#' @param population Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param gwas_sample_size Sample size of GWAS.
#' @param genome_ref_path Path to genome references for MAGMA (default = packageLocation/extdata/magma/g1000_eur).
#' @param magma_path Path to MAGMA executable (default is in package extdata folder, where download_dependencies() should store it).
#'
#' @return Dataframe
#'
#' @export
compute_MAGMA_enrichment_discrete <- function(gene_sets, sumstats_file, output_dir = NULL, upstream_kb = 10, downstream_kb = 1.5,
                                              gene_nomenclature = "hgnc", population = "eur", gwas_sample_size = NULL,
                                              genome_ref_path = paste0(system.file(package = "gwascelltyper"),"/extdata/g1000_", population),
                                              magma_path = paste0(system.file(package = "gwascelltyper"),"/extdata/magma-linux64")) {
  if (sum(gene_nomenclature %in% c("hgnc", "entrez")) != 1) {stop("Gene nomenclature of your input expression data must be either of: hgnc, entrez.")}
  if (sum(population %in% c("eur", "eas")) != 1) {stop("Unsupported population param supplied. Possible entries: eur, eas.")} # You gotta enlarge this.
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(genome_ref_path)) {stop("What MAGMA input genome_ref_path param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}


  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = names(gene_sets)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    names(gene_sets) <- format_celltype_names(names(gene_sets))
  }


  # Checking if mapping of SNPs to genes exists.
  prefix <- get_magma_paths(sumstats_file, upstream_kb, downstream_kb, population, gene_nomenclature)
  if (!(file.exists(paste0(prefix, ".genes.annot")) & file.exists(paste0(prefix, ".genes.out")) & file.exists(paste0(prefix, ".genes.raw")))) {
    print(".genes.annot file containing the mapping of SNPs (from the supplied sumstats file) to genes not found. Generating it.")
    map_SNPs_to_genes_for_MAGMA(sumstats_file = sumstats_file,
                                upstream_kb = upstream_kb,
                                downstream_kb = downstream_kb,
                                N = gwas_sample_size,
                                population = population,
                                gene_nomenclature = gene_nomenclature,
                                genome_ref_path = genome_ref_path,
                                magma_path = magma_path)
  }
  else cat(prefix, ".genes.annot\n", prefix, ".genes.out\n", prefix, ".genes.raw\nAll files printed before exist. Proceeding with the analysis.\n")


  print("MAGMA consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.")
  gc()


  # Boiler plate code in case user wants to store intermediary and output files somewhere instead of the default (which is a temporary directory).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for MAGMA files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }

  analysis_name = basename(output_dir)
  sumstats_file = path.expand(sumstats_file)


  # Take out the controls (all genes) as this is a competitive analysis (not self-contained).
  if ("All_genes_control" %in% names(gene_sets)) {
    cat("Removing geneset called 'All_genes_control', as this is a competitive analysis (not self-contained).\n")
    gene_sets <- gene_sets[-which(names(gene_sets)=="All_genes_control")]
  }


  cat("Preparing name map.\n")
  name_map <- cbind(names(gene_sets), paste0("Name_", 1:length(gene_sets), "_geneset"))
  names(gene_sets) <- name_map[,2]


  # Preparation for the enrichment test
  cat("Generating gene_covars.\n")
  geneCovarFile_contents <- NULL
  for (set in 1:length(gene_sets)) {
    # Geneset length check
    if (length(gene_sets[[set]])==0) stop(paste0("Gene set #", set, " doesn't contain any gene."))
    cat("Geneset ", names(gene_sets)[set], " contains ", length(gene_sets[[set]]), " genes.\n")

    # Collapse with " ", while writing the gene-set name at the beginning
    geneCovarFile_contents[set] <- paste0(names(gene_sets)[set], " ", paste(gene_sets[[set]], collapse=" "))
  }
  geneCovarFile = tempfile(pattern = "geneCovarFile", tmpdir = output_dir)
  cat("Writing gene_covar file to", geneCovarFile, "\n")
  write.table(geneCovarFile_contents, file = geneCovarFile, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


  cat("Making sure MAGMA is chmod +x.\n")
  system(paste0("chmod +x ", magma_path))


  # Now you run the analysis.
  magma_cmd = sprintf("%s --gene-results '%s.genes.raw' --set-annot '%s' --out '%s/%s'",
                      magma_path, prefix, geneCovarFile, output_dir, analysis_name)
  print(magma_cmd)
  system(magma_cmd)

  # Read and properly format the results (e.g. missing name column)
  res = read.table(paste0(output_dir,"/",analysis_name, ".gsa.out"), stringsAsFactors = FALSE, header = TRUE)
  if ("FULL_NAME" %in% colnames(res)) {
    res <- res[,-which(colnames(res) == "VARIABLE")]
    colnames(res)[which(colnames(res)=="FULL_NAME")] <- "Name"
    res <- res[,7:1]
    colnames(res)[2] <- "P_value"
  } else {
    colnames(res)[which(colnames(res)=="VARIABLE")] <- "Name"
    colnames(res)[which(colnames(res)=="P")] <- "P_value"
    res <- res[,c(1,7,2:6)]
  }

  # This returns the true gene-set names from the map.
  res$Name <- name_map[match(res$Name, name_map[,2]),1]

  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)


  return(res)
}


#' Compute SNPsea enrichment linear
#'
#' Compute GWAS enrichment in gene specificity matrices (linear analysis) with SNPsea.
#'
#' @param gene_sets Matrix with HGNC genes as rows & cell-types as columns. Values need to be numeric.
#' @param sumstats_file Address of properly formatted SNPsea GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param number_of_threads Number of threads to parallelize calculation over (default = 1).
#' @param slop If a SNP overlaps no gene intervals, extend the SNP interval this many nucleotides further and try again (default: 10000).
#' @param gene_intervals BED file with gene intervals. The fourth column must contain the same gene identifiers as in --gene-matrix.
#' @param snp_intervals BED file with all known SNP intervals. The fourth column must contain the same SNP identifiers as in --snps and --null-snps.
#' @param null_snps Text file with SNP identifiers to sample when generating null matched or random SNP sets. These SNPs must be a subset of --snp-intervals.
#' @param null_snpsets Generate a distribution of scores with N null matched SNP sets to evaluate type 1 error (default: 0).
#' @param min_observations Stop testing a column in --gene-matrix after observing this many null SNP sets with specificity scores greater or equal to those obtained with the SNPs in --snps. Increase this value to obtain more accurate p-values (default: 25).
#' @param max_iterations Maximum number of null SNP sets tested against each column in --gene-matrix. Increase this value to resolve small p-values (default: 10000).
#' @param gene_nomenclature Type of gene names.
#'
#' @return Dataframe
#'
#' @export
compute_SNPSEA_enrichment_discrete <- function(gene_sets, sumstats_file, output_dir = NULL, slop = "10e3", number_of_threads = 1,
                                               null_snpsets = 0, min_observations = 100, max_iterations = "1e7", gene_nomenclature = "hgnc",
                                               snpsea_path = paste0(system.file(package = "gwascelltyper"),"/extdata/snpsea-linux64"),
                                               gene_intervals = paste0(system.file(package = "gwascelltyper"),"/extdata/NCBIgenes2013_", gene_nomenclature,".bed.gz"),
                                               snp_intervals = paste0(system.file(package = "gwascelltyper"),"/extdata/TGP2011.bed.gz"),
                                               null_snps = paste0(system.file(package = "gwascelltyper"),"/extdata/Lango2010.txt.gz")) {
  require(R.utils)
  if (is.null(gene_sets)) {stop("What gene score matrix?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(gene_intervals)) {stop("What SNPsea gene_intervals param?")}
  if (is.null(snp_intervals)) {stop("What SNPsea snp_intervals param?")}
  if (is.null(null_snps)) {stop("What SNPsea null_snps param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = names(gene_sets)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    names(gene_sets) <- format_celltype_names(names(gene_sets))
  }

  print("SNPsea consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.")
  gc()


  # Take out the controls (all genes) as this is a competitive analysis (not self-contained).
  if ("All_genes_control" %in% names(gene_sets)) {
    cat("Removing geneset called 'All_genes_control', as this is a competitive analysis (not self-contained).\n")
    gene_sets <- gene_sets[-which(names(gene_sets)=="All_genes_control")]
  }


  # Boiler plate code in case user wants to store intermediary and output files somewhere instead of the default (which is a temporary directory).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for LDSC files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }


  # Need to make .gct.gz file before running SNPsea
  print("Building 0/1 matrix from gene_sets list")
  gene_set_matrix <- matrix(0, nrow = length(unique(unlist(gene_sets))), ncol = length(gene_sets))
  colnames(gene_set_matrix) <- names(gene_sets)
  rownames(gene_set_matrix) <- unique(unlist(gene_sets))
  for (set in 1:length(gene_sets)) {
    present_genes <- which(rownames(gene_set_matrix) %in% gene_sets[[set]])
    if (length(present_genes) > 0) {
      gene_set_matrix[present_genes, set] <- 1
    }
  }


  cat("Creating the GCT-format matrix.\n")
  gene_set_matrix = suppressWarnings(data.frame(Name = rownames(gene_set_matrix), Description = rownames(gene_set_matrix), gene_set_matrix))
  rownames(gene_set_matrix) <- NULL


  print("Writing the SNPsea GCT file.")
  write.table(gene_set_matrix, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, file = paste0(output_dir, "/gene_set_matrix.gct"))
  # The line underneath adds a necessary header to our file
  system(paste0("( echo -e '#1.2\n",paste(nrow(gene_set_matrix), ncol(gene_set_matrix)-2, sep = "\t"),"'; cat ",output_dir,"/gene_set_matrix.gct) > ",output_dir,"/tmp && mv ",output_dir,"/tmp ",output_dir,"/gene_set_matrix.gct"))
  gzip(paste0(output_dir,"/gene_set_matrix.gct"), destname=paste0(output_dir,"/gene_set_matrix.gct.gz"), overwrite=TRUE, remove=FALSE)


  cat("Making sure SNPsea is chmod +x.\n")
  system(paste0("chmod +x ", snpsea_path))


  # Wrap the arguments for SNPSEA
  command <- paste0("options=(--snps ", sumstats_file,
                    " --gene-matrix ", output_dir, "/gene_set_matrix.gct.gz",
                    " --gene-intervals ", gene_intervals,
                    " --snp-intervals ", snp_intervals,
                    " --null-snps ", null_snps,
                    " --out ", output_dir, "/snpsea",
                    " --slop ", slop,
                    " --threads ", number_of_threads,
                    " --null-snpsets ", null_snpsets,
                    " --min-observations ", min_observations,
                    " --max-iterations ", max_iterations, ")")
  cat("Running SNPsea with the following parameters:\n", command, "\n")
  system(paste0(command, " && ", snpsea_path, " ${options[*]}"))


  # Read file located in out folder, then change headers from c("condition", "pvalue") to c("Name","P_value")
  results <- read.csv(file = paste0(output_dir, "/snpsea/condition_pvalues.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(results) <- c("Name", "P_value", "Nulls_observed", "Nulls_tested")


  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)


  return(results)
}


#' Compute MAGMA enrichment linear
#'
#' Compute GWAS enrichment in gene specificity matrices (linear analysis) with MAGMA.
#'
#' @param gene_score_matrix Matrix with HGNC genes as rows & cell-types as columns. Values need to be numeric.
#' @param sumstats_file Address of properly formatted MAGMA GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param upstream_kb Number of kb upstream of gene, SNP-to-gene mapping window. Default = 10.
#' @param downstream_kb Number of kb downstream of gene, SNP-to-gene mapping window. Default = 1.5.
#' @param gene_nomenclature (default = hgnc)
#' @param population Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param gwas_sample_size Sample size of GWAS.
#' @param genome_ref_path Path to genome references for MAGMA.
#' @param magma_path Path to MAGMA executable (default is in package extdata folder, where download_dependencies() should store it).
#'
#' @return Dataframe
#'
#' @export
compute_MAGMA_enrichment_linear <- function(gene_score_matrix, sumstats_file, output_dir = NULL, upstream_kb = 10, downstream_kb = 1.5,
                                            gene_nomenclature = "hgnc", population = "eur", gwas_sample_size,
                                            genome_ref_path = paste0(system.file(package = "gwascelltyper"),"/extdata/g1000_", population),
                                            magma_path = paste0(system.file(package = "gwascelltyper"),"/extdata/magma-linux64")) {
  if (sum(gene_nomenclature %in% c("hgnc", "entrez")) != 1) {stop("Gene nomenclature of your input expression data must be either of: hgnc, entrez.")}
  if (sum(population %in% c("eur", "eas")) != 1) {stop("Unsupported population param supplied. Possible entries: eur, eas.")} # You gotta enlarge this.
  if (is.null(gene_score_matrix)) {stop("What gene score matrix?")}
  if (any(is.na(gene_score_matrix))) {stop("Gene score matrix contains NAs. This input is not supported.")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(genome_ref_path)) {stop("What MAGMA input genome_ref_path param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = colnames(gene_score_matrix)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    colnames(gene_score_matrix) <- format_celltype_names(colnames(gene_score_matrix))
  }


  # Checking if mapping of SNPs to genes exists.
  prefix <- get_magma_paths(sumstats_file, upstream_kb, downstream_kb, population, gene_nomenclature)
  if (!(file.exists(paste0(prefix, ".genes.annot")) & file.exists(paste0(prefix, ".genes.out")) & file.exists(paste0(prefix, ".genes.raw")))) {
    print(".genes.annot file containing the mapping of SNPs (from the supplied sumstats file) to genes not found. Generating it.")
    map_SNPs_to_genes_for_MAGMA(sumstats_file = sumstats_file,
                                upstream_kb = upstream_kb,
                                downstream_kb = downstream_kb,
                                N = gwas_sample_size,
                                population = population,
                                gene_nomenclature = gene_nomenclature,
                                genome_ref_path = genome_ref_path,
                                magma_path = magma_path)
  }
  else cat(prefix, ".genes.annot\n", prefix, ".genes.out\n", prefix, ".genes.raw\nAll files printed before exist. Proceeding with the analysis.\n")


  print("MAGMA linear consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.")
  gc()


  # Boiler plate code in case user wants to store intermediary and output files somewhere instead of the default (which is a temporary directory).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for MAGMA files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }


  cat("Making quantiles.\n")
  if (nrow(gene_score_matrix) < 100) {
    stop("Less than 100 genes in the gene score matrix. Can't fit gene scores into percentiles (100 bins).")
  } else {
    gene_score_matrix <- bin_specificity_into_quantiles(gene_score_matrix, numberOfBins = 100)
  }


  # Initializing a few params
  analysis_name = basename(output_dir)
  sumstats_file = path.expand(sumstats_file)
  cat("Preparing name map.\n")
  name_map <- cbind(colnames(gene_score_matrix), paste0("Name_", 1:ncol(gene_score_matrix), "_geneset"))
  colnames(gene_score_matrix) <- name_map[,2]
  if (gene_nomenclature == "hgnc") {
    gene_score_matrix <- cbind(hgnc = rownames(gene_score_matrix), gene_score_matrix)
  } else {
    gene_score_matrix <- cbind(entrez = rownames(gene_score_matrix), gene_score_matrix)
  }


  cat("Generating gene_covars.\n")
  geneCovarFile = tempfile(pattern = "geneCovarFile", tmpdir = output_dir, fileext = ".txt")
  cat("Writing gene_covar file for MAGMA.\n")
  write.table(gene_score_matrix, file = geneCovarFile, quote = FALSE, row.names = FALSE, sep = "\t")


  cat("Making sure MAGMA is chmod +x.\n")
  system(paste0("chmod +x ", magma_path))


  magma_cmd = sprintf("%s --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos --out '%s/linear_%s'",
                      magma_path, prefix, geneCovarFile, output_dir, analysis_name)
  print(magma_cmd)
  system(magma_cmd)


  # Read and properly format the results (e.g. missing name column)
  res = read.table(paste0(output_dir,"/linear_",analysis_name, ".gsa.out"), stringsAsFactors = FALSE, header = TRUE)
  if ("FULL_NAME" %in% colnames(res)) {
    res <- res[,-which(colnames(res) == "VARIABLE")]
    colnames(res)[which(colnames(res)=="FULL_NAME")] <- "Name"
    res <- res[,7:1]
    colnames(res)[2] <- "P_value"
  } else {
    colnames(res)[which(colnames(res)=="VARIABLE")] <- "Name"
    colnames(res)[which(colnames(res)=="P")] <- "P_value"
    res <- res[,c(1,7,2:6)]
  }


  # This returns the true geneset names.
  res$Name <- name_map[match(res$Name, name_map[,2]),1]


  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)


  return(res)
}


#' Compute SNPsea enrichment linear
#'
#' Compute GWAS enrichment in gene specificity matrices (linear analysis) with SNPsea.
#'
#' @param gene_score_matrix Matrix with HGNC genes as rows & cell-types as columns. Values need to be numeric.
#' @param sumstats_file Address of properly formatted SNPsea GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param number_of_threads Number of threads to parallelize calculation over (default = 1).
#' @param slop If a SNP overlaps no gene intervals, extend the SNP interval this many nucleotides further and try again (default: 10000).
#' @param gene_intervals BED file with gene intervals. The fourth column must contain the same gene identifiers as in --gene-matrix.
#' @param snp_intervals BED file with all known SNP intervals. The fourth column must contain the same SNP identifiers as in --snps and --null-snps.
#' @param null_snps Text file with SNP identifiers to sample when generating null matched or random SNP sets. These SNPs must be a subset of --snp-intervals.
#' @param null_snpsets Generate a distribution of scores with N null matched SNP sets to evaluate type 1 error (default: 0).
#' @param min_observations Stop testing a column in --gene-matrix after observing this many null SNP sets with specificity scores greater or equal to those obtained with the SNPs in --snps. Increase this value to obtain more accurate p-values (default: 25).
#' @param max_iterations Maximum number of null SNP sets tested against each column in --gene-matrix. Increase this value to resolve small p-values (default: 10000).
#' @param gene_nomenclature Type of gene names.
#'
#' @return Dataframe
#'
#' @export
compute_SNPSEA_enrichment_linear <- function(gene_score_matrix, sumstats_file, output_dir = NULL, number_of_threads = 1, slop = "10e3",
                                             null_snpsets = 0, min_observations = 100, max_iterations = "1e7", gene_nomenclature = "hgnc",
                                             snpsea_path = paste0(system.file(package = "gwascelltyper"),"/extdata/snpsea-linux64"),
                                             gene_intervals = paste0(system.file(package = "gwascelltyper"),"/extdata/NCBIgenes2013_", gene_nomenclature,".bed.gz"),
                                             snp_intervals = paste0(system.file(package = "gwascelltyper"),"/extdata/TGP2011.bed.gz"),
                                             null_snps = paste0(system.file(package = "gwascelltyper"),"/extdata/Lango2010.txt.gz")) {
  require(R.utils)
  if (is.null(gene_score_matrix)) {stop("What gene score matrix?")}
  if (any(is.na(gene_score_matrix))) {stop("Gene score matrix contains NAs. This input is not supported.")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(gene_intervals)) {stop("What SNPsea gene_intervals param?")}
  if (is.null(snp_intervals)) {stop("What SNPsea snp_intervals param?")}
  if (is.null(null_snps)) {stop("What SNPsea null_snps param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  if (any(grepl(pattern = "[^A-Za-z0-9_]", x = colnames(gene_score_matrix)))) {
    cat("At least one of the cell type names contains at least a special character. Replacing with underscore.\n")
    colnames(gene_score_matrix) <- format_celltype_names(colnames(gene_score_matrix))
  }




  print("SNPsea linear consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.")
  gc()


  # Boiler plate code in case user wants to store intermediary and output files somewhere instead of the default (which is a temporary directory).
  store_results = FALSE
  if (is.null(output_dir)) {
    output_dir <- tempdir(check = TRUE)
    cat("Temporary directory for SNPsea files is", output_dir, "\nTemporary dir will be deleted after the analysis is done.\n")
  }
  else {
    store_results = TRUE

    # Make sure that the output dir ends without forward slash
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/") {
      output_dir <- gsub('.{1}$', '', output_dir)
    }

    if (!dir.exists(output_dir)) {
      cat("Output dir supplied does not exist. Creating ", output_dir, "\n")
      dir.create(output_dir)
    }
  }


  cat("Creating the GCT-format matrix.\n")
  gene_score_matrix = suppressWarnings(data.frame(Name = rownames(gene_score_matrix), Description = rownames(gene_score_matrix), gene_score_matrix))
  rownames(gene_score_matrix) <- NULL


  print("Writing the SNPsea GCT file.")
  write.table(gene_score_matrix, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, file = paste0(output_dir,"/gene_score_matrix.gct"))
  system(paste0("( echo -e '#1.2\n",paste(nrow(gene_score_matrix), ncol(gene_score_matrix)-2, sep = "\t"),"'; cat ",output_dir,"/gene_score_matrix.gct) > ",output_dir,"/tmp && mv ",output_dir,"/tmp ",output_dir,"/gene_score_matrix.gct"))
  gzip(paste0(output_dir,"/gene_score_matrix.gct"), destname=paste0(output_dir,"/gene_score_matrix.gct.gz"), overwrite=TRUE, remove=FALSE)


  cat("Making sure SNPsea is chmod +x.\n")
  system(paste0("chmod +x ", snpsea_path))


  # Wrap the arguments for SNPSEA
  command <- paste0("options=(--snps ", sumstats_file,
                    " --gene-matrix ", output_dir, "/gene_score_matrix.gct.gz",
                    " --gene-intervals ", gene_intervals,
                    " --snp-intervals ", snp_intervals,
                    " --null-snps ", null_snps,
                    " --out ", output_dir, "/snpsea",
                    " --slop ", slop,
                    " --threads ", number_of_threads,
                    " --null-snpsets ", null_snpsets,
                    " --min-observations ", min_observations,
                    " --max-iterations ", max_iterations, ")")
  cat("Running SNPsea with the following parameters:\n", command, "\n")
  system(paste0(command, " && ", snpsea_path, " ${options[*]}"))


  # Read file located in out folder, then change headers from c("condition", "pvalue") to c("Name","P_value")
  results <- read.csv(file = paste0(output_dir, "/snpsea/condition_pvalues.txt"), sep = "\t", header = TRUE)
  colnames(results) <- c("Name", "P_value", "Nulls_observed", "Nulls_tested")



  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)



  return(results)
}


#' Map SNPs to genes for MAGMA
#'
#' Compute GWAS enrichment in discrete lists of genes with MAGMA.
#'
#' @param gwas_sumstats_path Path to summary statistics file.
#' @param upstream_kb Default = 10
#' @param downstream_kb Default = 1.5
#' @param N Default = NULL
#' @param population Default = eur
#' @param gene_nomenclature Default = hgnc
#' @param genome_ref_path Path to g1000 dependency files. Depends on population used.
#' @param magma_path Path to MAGMA executable.
#'
#' @return Character
#'
#' @export
map_SNPs_to_genes_for_MAGMA <- function (sumstats_file, upstream_kb = 10, downstream_kb = 1.5, N = NULL,
                                         population = "eur", gene_nomenclature = "hgnc",
                                         genome_ref_path = paste0(system.file(package = "gwascelltyper"),"/extdata/g1000_", population),
                                         magma_path = paste0(system.file(package = "gwascelltyper"),"/extdata/magma-linux64")) {

  sumstats_file = path.expand(sumstats_file)
  prefix = get_magma_paths(sumstats_file, upstream_kb, downstream_kb, population, gene_nomenclature)


  cat("Figuring out the N.\n")
  if (is.null(N)) {
    if (!"N" %in% strsplit(readLines(con = sumstats_file, n = 1), "\t")[[1]]) {
      stop(paste0("Error: GWAS_sample_size param is NULL and there is no N column in the sumstats file: ", sumstats_file))
    } else {n_arg = "ncol=N"}
  }
  else {
    if (is.numeric(N)) {
      if (N < 1000) stop("Value of N provided is less than 1000. This seems unlikely?")
      if (N > 1e+08) stop("Value of N provided is over than 100000000. In 2019 this seems unlikely.")
      n_arg = paste0("N=", N) # This feeds in a unique sample size for all of the SNPs in the sumstats file.
    } else {stop("Error: The supplied GWAS_sample_size param is not numeric")}
  }


  cat("Figuring out the genome build used for these sumstats.\n")
  genome_build = get_genomebuild_for_sumstats(sumstats_file)
  gene_loc_dir = sprintf("%s/extdata/", system.file(package = "gwascelltyper"))
  if (genome_build == "GRCh37") genomeLocFile = paste0(gene_loc_dir, "/NCBI37.3.", gene_nomenclature,".gene.loc")
  if (genome_build == "GRCh38") genomeLocFile = paste0(gene_loc_dir, "/NCBI38.", gene_nomenclature,".gene.loc")

  cat("Making sure MAGMA is chmod +x.\n")
  system(paste0("chmod +x ", magma_path))

  cat("GWAS Sumstats appear to come from genome build:", genome_build, "\n")
  magma_cmd = sprintf("%s --annotate window=%s,%s --snp-loc '%s' --gene-loc '%s' --out '%s'",
                      magma_path, upstream_kb, downstream_kb, sumstats_file, genomeLocFile, prefix)
  system(magma_cmd)
  magma_cmd = sprintf("%s --bfile '%s' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'",
                      magma_path, path.expand(genome_ref_path), sumstats_file, n_arg, prefix, prefix)
  system(magma_cmd)

  return(sprintf("%s.genes.out", prefix))
}


#' Get genome build from sumstats
#'
#' Tries to estimate what genome build a sumstats file comes from.
#'
#' @param path Path of summary statistics file.
#'
#' @return Character
#'
#' @export
get_genomebuild_for_sumstats <- function (path) {
  topLines = read.table(path, nrows = 30000, header = TRUE, stringsAsFactors = FALSE)
  topSNPS = topLines$SNP
  topLOCs = sprintf("%s-%s-%s", topLines$SNP, topLines$CHR, topLines$BP)
  data("snp_loc")
  sub_SNP_LOC_DATA = snp_loc[sample(1:dim(snp_loc)[1], 1e+05),]
  topSNP_locs = sub_SNP_LOC_DATA[sub_SNP_LOC_DATA$SNP %in% topLines, ]
  topSNP_locs$locs = sprintf("%s-%s-%s", topSNP_locs$SNP, topSNP_locs$CHR, topSNP_locs$BP)
  return(names(sort(table(topSNP_locs[topSNP_locs$locs %in% topLOCs, ]$Build), decreasing = TRUE)[1]))
}


#' Get MAGMA paths
#'
#' Get or create MAGMA folders for the sumstats file provided. Outputs a path.
#'
#' @param sumstats_file Path for summary statistics file.
#' @param upstream_kb Kilobases upstream of gene to add.
#' @param downstream_kb Kilobases downstream of gene to add.
#' @param population Population (default = eur)
#' @param gene_nomenclature Gene nomenclature (default = hgnc)
#'
#' @return List
#'
#' @export
get_magma_paths <- function (sumstats_file, upstream_kb, downstream_kb, population, gene_nomenclature) {
  prefix <- paste0(dirname(sumstats_file), "/MAGMA_Files/", basename(sumstats_file), ".", upstream_kb, "UP.", downstream_kb, "DOWN_", population, "_", gene_nomenclature)
  if (!dir.exists(dirname(prefix))) dir.create(dirname(prefix), showWarnings = FALSE)
  return(prefix)
}


#' Bin specificity into quantiles
#'
#' Mostly an adjuvant function for MAGMA linear, which requires quantiles (natural numbers) as input.
#'
#' @param input_matrix Input matrix.
#' @param numberOfBins Number of bins.
#'
#' @return Dataframe
#'
#' @export
bin_specificity_into_quantiles <- function(input_matrix,numberOfBins){
  if (any(is.na(input_matrix))) {
    print("THIS MATRIX CONTAINS NAs.")
  }
  bin.columns.into.quantiles <- function(input_matrix,numberOfBins=100){
    quantileValues = rep(0,length(input_matrix))
    quantileValues[input_matrix>0] = as.numeric(cut(input_matrix[input_matrix>0],
                                                    breaks=unique(quantile(input_matrix[input_matrix>0], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)),
                                                    include.lowest=TRUE))
    return(quantileValues)
  }
  specificity_quantiles = apply(input_matrix,2,FUN=bin.columns.into.quantiles,numberOfBins=100)
  rownames(specificity_quantiles) = rownames(input_matrix)
  return(specificity_quantiles)
}


#' Seurat enrichment testing
#'
#' Takes a Seurat object with differentially expressed genes attached, and output results for enrichment testing.
#'
#' @param seurat_object
#'
#' @return Dataframe
#'
#' @export
Seurat_enrichment_testing <- function(seurat_object, magma_sumstat = NULL, ldsc_sumstat = NULL, snpsea_sumstat = NULL, gwas_sample_size = NULL) {

}
