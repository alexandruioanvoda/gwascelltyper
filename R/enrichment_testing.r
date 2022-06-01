#' Compute LDSC enrichment - optimized
#'
#' Compute GWAS enrichment in discrete lists of genes with LDSC. Optimized to load all dependency files into Python-allocated RAM only once.
#'
#' @param gene_sets List of genesets. Each geneset must be a character vector, containing HGNC gene names.
#' @param sumstats_file Address of properly formatted LDSC GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param window_size Size of window surrounding gene, for SNP-to-gene mapping purposes. Default = 100000 (100kb).
#' @param pop Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param ld_wind_cm Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM). Default = 1.
#' @param gene_nomenclature Type of gene IDs fed in from the expression data (default = hgnc; alternatives = entrez).
#' @param deps A named list, containing all of the dependencies for LDSC (e.g. Python script paths, SNP linkage disequilibrium data, etc.). Provided by default by another function, called gwascelltyper::dep_path().
#'
#' @import doParallel
#' @import data.table
#' @importFrom R.utils gzip
#' @importFrom R.utils compressFile
#' @importFrom utils write.table
#'
#' @return Dataframe
#'
#' @export
compute_LDSC_enrichment_discrete_optimized <- function(gene_sets, sumstats_file, output_dir = NULL, number_of_threads = 1, window_size = 100000, pop = "eur", ld_wind_cm = 1,
                                                       deps = dep_path(pop, "ldsc")) {
  # Checks
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("GWAS sumstats file does not exist at provided path.")}
  # dep_path() already checks the existence of the other required files!

  
  # A few other boilerplate checks: # This should be modified to a general gene-set renamer applicable to all enrichment-testing functions.
  window_size <- format(window_size, scientific = FALSE)
  if (names(gene_sets)[length(gene_sets)] != "All_genes_control") stop("No <All_genes_control> contained within the last geneset from the gene_sets object.")
  gobj <- gene_object_renamer(gene_sets)
  gene_sets <- gobj$gobj

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
  utils::write.table(gene_sets, file = geneset_file, quote = FALSE, sep = "\t", na = "", row.names = FALSE, col.names = FALSE)

  cat("Making a thin-annot of all gene-sets (one per column), per each chromosome.\n")
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads," MKL_NUM_THREADS=", number_of_threads,
                    " python ", deps$ldsc_scripts[2]," --gene-set-file ", geneset_file,
                    " --gene-coord-file ", deps$geneloc_ldsc,
                    " --windowsize ", window_size,
                    " --bimfile ", deps$bim_prefix,
                    " --annot-file ", output_dir, "/out.")
  system(command = command)
  cat("Thin annot making finished at:", as.character(Sys.time()), "\n")


  cat("LDSC consumes a lot of RAM memory, so we are calling garbage collection in R now, before proceeding to the analysis.\n")
  gc()
  cat("Initializing R multithreading.\n")
  registerDoParallel(cores=number_of_threads)

  foreach::foreach(i=1:22, .options.snow=list(preschedule=TRUE)) %do% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", deps$ldsc_scripts[1], " --l2 --bfile ",
                      deps$bim_prefix, as.character(i), " --ld-wind-cm ", ld_wind_cm,
                      " --annot ", output_dir, "/out.chr", as.character(i), ".annot.gz",
                      " --thin-annot --out ", output_dir, "/merged_ld_scores.", as.character(i),
                      " --print-snps ", deps$print_snps, as.character(i), ".snp")
    system(command = command)

    cat("Reading in the LD score, M & M_5_50 files.\n")
    merged_ld <- utils::read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.ldscore.gz"), header = TRUE)
    merged_m <- utils::read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.M"), header = FALSE)
    merged_m_5_50 <- utils::read.table(file = paste0(output_dir, "/merged_ld_scores.", i, ".l2.M_5_50"), header = FALSE)

    cat("Separating the LD score, M & M_5_50 files.\n")
    for (j in 1:(ncol(merged_ld)-3)) {
      de_merged <- merged_ld[,c(1,2,3,j+3)]
      utils::write.table(x = de_merged, file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
      R.utils::gzip(paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore"), destname=paste0(output_dir, "/geneset.", j, ".", i, ".l2.ldscore.gz"), overwrite=TRUE, remove=TRUE)
      utils::write.table(x = merged_m[,j], file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.M"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      utils::write.table(x = merged_m_5_50[,j], file = paste0(output_dir, "/geneset.", j, ".", i, ".l2.M_5_50"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
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
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads, " python ", deps$ldsc_scripts[1], " --h2-cts ", sumstats_file,
                    " --ref-ld-chr ", deps$ref_ld_chr, " ",
                    "--out ", output_dir, "/ldsc_final_output ",
                    "--ref-ld-chr-cts ", output_dir, "/main.ldcts ",
                    "--w-ld-chr ", deps$w_ld_chr)
  print(command)
  system(command)
  print("Done. Returning results dataframe.")


  # Here we build the results dataframe for plotting.
  results <- utils::read.table(paste0(output_dir, "/ldsc_final_output.cell_type_results.txt"), header = T, stringsAsFactors = FALSE)
  colnames(results)[4] <- "P_value"
  results$Name <- gobj$gsn$orig_names[match(results$Name, gobj$gsn$temp_names)]


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
#' @param window_size Size of window surrounding gene, for SNP-to-gene mapping purposes. Default = 100000 (100kb).
#' @param pop Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param ld_wind_cm Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM). Default = 1.
#' @param gene_nomenclature Type of gene IDs fed in from the expression data (default = hgnc; alternatives = entrez).
#' @param deps A named list, containing all of the dependencies for LDSC (e.g. Python script paths, SNP linkage disequilibrium data, etc.). Provided by default by another function, called gwascelltyper::dep_path().
#'
#' @import doParallel
#' @import data.table
#' @importFrom R.utils gzip
#' @importFrom R.utils compressFile
#' @importFrom utils write.table
#'
#' @return Dataframe
#'
#' @export
compute_LDSC_enrichment_discrete <- function(gene_sets, sumstats_file, output_dir = NULL, number_of_threads = 1, window_size = 100000, pop = "eur", ld_wind_cm = 1,
                                                       deps = dep_path(pop, "ldsc")) {
  # Checks
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("GWAS sumstats file does not exist at provided path.")}
  # dep_path() already checks the existence of the other required files!
  
  
  # A few other boilerplate checks:
  window_size <- format(window_size, scientific = FALSE)
  if (names(gene_sets)[length(gene_sets)] != "All_genes_control") stop("No <All_genes_control> contained within the last geneset from the gene_sets object.")
  gobj <- gene_object_renamer(gene_sets)
  gene_sets <- gobj$gobj


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
  doParallel::registerDoParallel(cores=number_of_threads)


  print("Writing lists of genes to files and making thin-annots for each gene set...")
  foreach::foreach(i=1:(length(gene_sets)-1), .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    geneset_file <- paste0(output_dir, "/genes_", gsub(" ", "_", names(gene_sets)[i]), "_list.txt") # The gsub command is the cell name
    con <- file(geneset_file)
    writeLines(gene_sets[[i]], con)
    close(con)
    for (j in seq(from = 1, to = 22, by = 1)) { # Compute annot file / each chromosome
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", deps$ldsc_scripts[2], " --gene-set-file ",
                        geneset_file, " --gene-coord-file ", deps$geneloc_ldsc,
                        " --windowsize ", window_size,
                        " --bimfile ", deps$bim_prefix, as.character(j), ".bim ",
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
  foreach::foreach(j=1:22, .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", deps$ldsc_scripts[2], " --gene-set-file ",
                      geneset_file, " --gene-coord-file ", deps$geneloc_ldsc,
                      " --windowsize ", window_size,
                      " --bimfile ", deps$bim_prefix, as.character(j), ".bim ",
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
  foreach::foreach(i=1:(length(gene_sets)-1), .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    for (j in seq(from = 1, to = 22, by = 1)) {
      command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", deps$ldsc_scripts[1], " --l2 --bfile ",
                        deps$bim_prefix, as.character(j), " --ld-wind-cm ", deps$ld_wind_cm,
                        " --annot ", output_dir, "/celltype.", as.character(i), ".", as.character(j), ".annot.gz",
                        " --thin-annot --out ", output_dir, "/celltype.", as.character(i), ".", as.character(j),
                        " --print-snps ", deps$print_snps, as.character(j), ".snp")
      print(command)
      system(command = command)
    }
  }


  print("Calling garbage collection now, to free memory after the LD/geneset parallel operations done by DoParallel.")
  gc()


  # Compute LD for the control too
  print("Compute LD for the control too")
  foreach::foreach(j=1:22, .options.snow=list(preschedule=TRUE)) %dopar% {
    gc()
    command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python ", deps$ldsc_scripts[1], " --l2 --bfile ",
                      deps$bim_prefix, as.character(j), " --ld-wind-cm ", deps$ld_wind_cm,
                      " --annot ", output_dir, "/control.", as.character(j), ".annot.gz",
                      " --thin-annot --out ", output_dir, "/control.", as.character(j),
                      " --print-snps ", deps$print_snps, as.character(j), ".snp")
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
  command <- paste0("source activate ldsc && OPENBLAS_NUM_THREADS=", number_of_threads, " MKL_NUM_THREADS=", number_of_threads, " python ", deps$ldsc_scripts[1], " --h2-cts ", sumstats_file,
                    " --ref-ld-chr ", deps$ref_ld_chr, " ",
                    "--out ", output_dir, "/ldsc_final_output ",
                    "--ref-ld-chr-cts ", output_dir, "/main.ldcts ",
                    "--w-ld-chr ", deps$w_ld_chr)
  print(command)
  system(command)
  print("Done. Returning results dataframe.")


  # Here we build the results dataframe for plotting.
  results <- utils::read.table(paste0(output_dir, "/ldsc_final_output.cell_type_results.txt"), header = T, stringsAsFactors = FALSE)
  colnames(results)[4] <- "P_value"
  results$Name <- gobj$gsn$orig_names[match(results$Name, gobj$gsn$temp_names)]


  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)


  return(results)
}


#' Compute MAGMA enrichment
#'
#' Compute GWAS enrichment in discrete lists of genes with MAGMA.
#'
#' @param gene_sets List of genesets. Each geneset must be a character vector, containing HGNC gene names.
#' @param genesets_conditioned_for Character vector of name(s) of the geneset(s) (from the above list) that should be conditioned for (if NULL, no geneset is conditioned for). Default is NULL.
#' @param sumstats_file Address of properly formatted MAGMA GWAS file.
#' @param output_dir Folder where you want the output to be stored. Default = NULL, which makes this function create a directory where results will be temporarily stored. The directory will be deleted afterwards, if this parameter is null. The only output will be the return. Otherwise, both the files and the return are kept.
#' @param upstream_kb Number of kb upstream of gene, SNP-to-gene mapping window. Default = 10.
#' @param downstream_kb Number of kb downstream of gene, SNP-to-gene mapping window. Default = 1.5.
#' @param gene_nomenclature Type of gene IDs fed in from the expression data (default = hgnc; alternatives = entrez).
#' @param pop Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param gwas_sample_size Sample size of GWAS.
#' @param deps A named list, containing all of the dependencies for MAGMA (e.g. executables, SNP linkage disequilibrium data, etc.). Provided by default by another function, called gwascelltyper::dep_path().
#'
#' @import doParallel
#' @import data.table
#' @importFrom R.utils gzip
#' @importFrom R.utils compressFile
#' @importFrom utils read.table
#' @importFrom utils write.table
#'
#' @return Dataframe
#'
#' @export
compute_MAGMA_enrichment_discrete <- function(gene_sets, genesets_conditioned_for = NULL, sumstats_file, output_dir = NULL, upstream_kb = 10, downstream_kb = 1.5,
                                              gene_nomenclature = "hgnc", pop = "eur", gwas_sample_size = NULL,
                                              deps = dep_path(pop, "magma")) {
  if (sum(gene_nomenclature %in% c("hgnc", "entrez")) != 1) {stop("Gene nomenclature of your input expression data must be either of: hgnc, entrez.")}
  if (sum(pop %in% c("eur", "eas")) != 1) {stop("Unsupported population param supplied. Possible entries: eur, eas.")} # You gotta enlarge this.
  if (is.null(gene_sets)) {stop("What gene sets?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  

  # Checking if mapping of SNPs to genes exists.
  prefix <- get_magma_paths(sumstats_file, upstream_kb, downstream_kb, pop, gene_nomenclature)
  if (!(file.exists(paste0(prefix, ".genes.annot")) & file.exists(paste0(prefix, ".genes.out")) & file.exists(paste0(prefix, ".genes.raw")))) {
    print(".genes.annot file containing the mapping of SNPs (from the supplied sumstats file) to genes not found. Generating it.")
    map_SNPs_to_genes_for_MAGMA(sumstats_file = sumstats_file,
                                upstream_kb = upstream_kb,
                                downstream_kb = downstream_kb,
                                N = gwas_sample_size,
                                pop = pop,
                                gene_nomenclature = gene_nomenclature,
                                deps = deps)
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
  gobj <- gene_object_renamer(gene_sets)
  gene_sets <- gobj$gobj
  if (!is.null(genesets_conditioned_for)) {
    cntrl_gs <- paste(gobj$gsn$temp_names[match(genesets_conditioned_for, gobj$gsn$orig_names)], collapse=",")
  }


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
  utils::write.table(geneCovarFile_contents, file = geneCovarFile, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


  cat("Making sure MAGMA is chmod +x.\n")
  system(paste0("chmod +x ", deps[["magma_binary"]]))


  # Now you run the analysis.
  if (!is.null(genesets_conditioned_for)) {
    magma_cmd = sprintf("%s --gene-results '%s.genes.raw' --set-annot '%s' --model direction=positive condition='%s' --out '%s/%s'",
                        deps[["magma_binary"]],         prefix,                    geneCovarFile,                           cntrl_gs,  output_dir, analysis_name)
  } else {
    magma_cmd = sprintf("%s --gene-results '%s.genes.raw' --set-annot '%s' --out '%s/%s'",
                        deps[["magma_binary"]],         prefix,                    geneCovarFile, output_dir, analysis_name)
  }
  print(magma_cmd)
  system(magma_cmd)

  # Read and properly format the results (e.g. missing name column)
  res = utils::read.table(paste0(output_dir,"/",analysis_name, ".gsa.out"), stringsAsFactors = FALSE, header = TRUE)
  if ("FULL_NAME" %in% colnames(res)) {
    colnames(res)[which(colnames(res)=="FULL_NAME")] <- "Name"
  } else {
    colnames(res)[which(colnames(res)=="VARIABLE")] <- "Name"
  }
  colnames(res)[which(colnames(res)=="P")] <- "P_value"

  # This returns the true gene-set names from the map.
  res$Name <- gobj$gsn$orig_names[match(res$Name, gobj$gsn$temp_names)]

  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)

  return(res)
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
#' @param pop Word describing ancestry of the GWAS data. Relevant so that the proper LD panel will be used (default = eur; alternatives = afr, amr, eas, sas, subpop)
#' @param gwas_sample_size Sample size of GWAS.
#' @param deps A named list, containing all of the dependencies for MAGMA (e.g. executables, SNP linkage disequilibrium data, etc.). Provided by default by another function, called gwascelltyper::dep_path().
#'
#' @importFrom utils write.table
#'
#' @return Dataframe
#'
#' @export
compute_MAGMA_enrichment_linear <- function(gene_score_matrix, sumstats_file, output_dir = NULL, upstream_kb = 10, downstream_kb = 1.5,
                                            gene_nomenclature = "hgnc", pop = "eur", gwas_sample_size,
                                            deps = dep_path(pop, "magma")) {
  if (sum(gene_nomenclature %in% c("hgnc", "entrez")) != 1) {stop("Gene nomenclature of your input expression data must be either of: hgnc, entrez.")}
  if (sum(pop %in% c("eur", "eas")) != 1) {stop("Unsupported population param supplied. Possible entries: eur, eas.")} # You gotta enlarge this.
  if (is.null(gene_score_matrix)) {stop("What gene score matrix?")}
  if (any(is.na(gene_score_matrix))) {stop("Gene score matrix contains NAs. This input is not supported.")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  
  
  # Checking if mapping of SNPs to genes exists.
  prefix <- get_magma_paths(sumstats_file, upstream_kb, downstream_kb, pop, gene_nomenclature)
  if (!(file.exists(paste0(prefix, ".genes.annot")) & file.exists(paste0(prefix, ".genes.out")) & file.exists(paste0(prefix, ".genes.raw")))) {
    print(".genes.annot file containing the mapping of SNPs (from the supplied sumstats file) to genes not found. Generating it.")
    map_SNPs_to_genes_for_MAGMA(sumstats_file = sumstats_file,
                                upstream_kb = upstream_kb,
                                downstream_kb = downstream_kb,
                                N = gwas_sample_size,
                                pop = pop,
                                gene_nomenclature = gene_nomenclature,
                                deps = deps)
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
  gobj <- gene_object_renamer(gene_score_matrix)
  gene_score_matrix <- gobj$gobj
  if (gene_nomenclature == "hgnc") {
    gene_score_matrix <- cbind(hgnc = rownames(gene_score_matrix), gene_score_matrix)
  } else {
    gene_score_matrix <- cbind(entrez = rownames(gene_score_matrix), gene_score_matrix)
  }
  
  
  cat("Generating gene_covars.\n")
  geneCovarFile = tempfile(pattern = "geneCovarFile", tmpdir = output_dir, fileext = ".txt")
  cat("Writing gene_covar file for MAGMA.\n")
  utils::write.table(gene_score_matrix, file = geneCovarFile, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  cat("Making sure MAGMA is chmod +x.\n")
  system(paste0("chmod +x ", deps[["magma_binary"]]))
  
  
  magma_cmd = sprintf("%s --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos --out '%s/linear_%s'",
                      deps[["magma_binary"]], prefix, geneCovarFile, output_dir, analysis_name)
  print(magma_cmd)
  system(magma_cmd)
  
  
  # Read and properly format the results (e.g. missing name column)
  res = utils::read.table(paste0(output_dir,"/linear_",analysis_name, ".gsa.out"), stringsAsFactors = FALSE, header = TRUE)
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
  results$Name <- gobj$gsn$orig_names[match(results$Name, gobj$gsn$temp_names)]
  
  
  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)
  
  
  return(res)
}


#' Map SNPs to genes for MAGMA
#'
#' Compute GWAS enrichment in discrete lists of genes with MAGMA.
#'
#' @param sumstats_file Path to summary statistics file.
#' @param upstream_kb Default = 10
#' @param downstream_kb Default = 1.5
#' @param N Default = NULL
#' @param pop Default = eur
#' @param gene_nomenclature Default = hgnc
#' @param deps A named list, containing all of the dependencies for MAGMA (e.g. executables, SNP linkage disequilibrium data, etc.). Provided by default by another function, called gwascelltyper::dep_path().
#'
#' @return Character
#'
#' @export
map_SNPs_to_genes_for_MAGMA <- function (sumstats_file, upstream_kb = 10, downstream_kb = 1.5, N = NULL,
                                         pop = "eur", gene_nomenclature = "hgnc",
                                         deps = dep_path(pop, "magma")) {

  sumstats_file = path.expand(sumstats_file)
  prefix = get_magma_paths(sumstats_file, upstream_kb, downstream_kb, pop, gene_nomenclature)

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
  system(paste0("chmod +x ", deps[["magma_binary"]]))

  cat("GWAS Sumstats appear to come from genome build:", genome_build, "\n")
  magma_cmd = sprintf("%s --annotate window=%s,%s --snp-loc '%s' --gene-loc '%s' --out '%s'",
                      deps[["magma_binary"]], upstream_kb, downstream_kb, sumstats_file, genomeLocFile, prefix)
  system(magma_cmd)
  magma_cmd = sprintf("%s --bfile '%s' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'",
                      deps[["magma_binary"]], path.expand(deps$genome_ref_path), sumstats_file, n_arg, prefix, prefix)
  system(magma_cmd)

  return(sprintf("%s.genes.out", prefix))
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
#' @param snpsea_path Path to the SNPsea executable
#'
#' @importFrom utils read.csv
#' @importFrom utils write.table
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
  if (is.null(gene_sets)) {stop("What gene score matrix?")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(gene_intervals)) {stop("What SNPsea gene_intervals param?")}
  if (is.null(snp_intervals)) {stop("What SNPsea snp_intervals param?")}
  if (is.null(null_snps)) {stop("What SNPsea null_snps param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  
  
  gobj <- gene_object_renamer(gene_sets)
  gene_sets <- gobj$gobj
  
  
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
  utils::write.table(gene_set_matrix, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, file = paste0(output_dir, "/gene_set_matrix.gct"))
  # The line underneath adds a necessary header to our file
  system(paste0("( echo -e '#1.2\n",paste(nrow(gene_set_matrix), ncol(gene_set_matrix)-2, sep = "\t"),"'; cat ",output_dir,"/gene_set_matrix.gct) > ",output_dir,"/tmp && mv ",output_dir,"/tmp ",output_dir,"/gene_set_matrix.gct"))
  R.utils::gzip(paste0(output_dir,"/gene_set_matrix.gct"), destname=paste0(output_dir,"/gene_set_matrix.gct.gz"), overwrite=TRUE, remove=FALSE)
  
  
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
  results$Name <- gobj$gsn$orig_names[match(results$Name, gobj$gsn$temp_names)]
  
  
  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)
  
  
  return(results)
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
#' @param snpsea_path Path to the SNPsea executable
#'
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @importFrom R.utils gzip
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
  if (is.null(gene_score_matrix)) {stop("What gene score matrix?")}
  if (any(is.na(gene_score_matrix))) {stop("Gene score matrix contains NAs. This input is not supported.")}
  if (is.null(sumstats_file)) {stop("What sumstats file?")}
  if (is.null(gene_intervals)) {stop("What SNPsea gene_intervals param?")}
  if (is.null(snp_intervals)) {stop("What SNPsea snp_intervals param?")}
  if (is.null(null_snps)) {stop("What SNPsea null_snps param?")}
  if (!file.exists(sumstats_file)) {stop("GWAS summary statistics file doesn't exist.")}
  
  gobj <- gene_object_renamer(gene_score_matrix)
  gene_score_matrix <- gobj$gobj
  
  
  
  
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
  utils::write.table(gene_score_matrix, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, file = paste0(output_dir,"/gene_score_matrix.gct"))
  system(paste0("( echo -e '#1.2\n",paste(nrow(gene_score_matrix), ncol(gene_score_matrix)-2, sep = "\t"),"'; cat ",output_dir,"/gene_score_matrix.gct) > ",output_dir,"/tmp && mv ",output_dir,"/tmp ",output_dir,"/gene_score_matrix.gct"))
  R.utils::gzip(paste0(output_dir,"/gene_score_matrix.gct"), destname=paste0(output_dir,"/gene_score_matrix.gct.gz"), overwrite=TRUE, remove=FALSE)
  
  
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
  results$Name <- gobj$gsn$orig_names[match(results$Name, gobj$gsn$temp_names)]
  
  
  
  # Now that the results are read into memory, we're deleting all of the intermediary files the UNIX R way.
  if (!store_results) unlink(output_dir, recursive = TRUE)
  
  
  
  return(results)
}


#' Get genome build from sumstats
#'
#' Tries to estimate what genome build a sumstats file comes from.
#'
#' @param path Path of summary statistics file.
#'
#' @importFrom utils data
#'
#' @return Character
#'
#' @export
get_genomebuild_for_sumstats <- function (path) {
  topLines = utils::read.table(path, nrows = 30000, header = TRUE, stringsAsFactors = FALSE)
  topSNPS = topLines$SNP
  topLOCs = sprintf("%s-%s-%s", topLines$SNP, topLines$CHR, topLines$BP)
  snp_loc <- gwascelltyper::snp_loc
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
#' @param pop Population (default = eur)
#' @param gene_nomenclature Gene nomenclature (default = hgnc)
#'
#' @return List
#'
#' @export
get_magma_paths <- function (sumstats_file, upstream_kb, downstream_kb, pop, gene_nomenclature) {
  prefix <- paste0(dirname(sumstats_file), "/MAGMA_Files/", basename(sumstats_file), ".", upstream_kb, "UP.", downstream_kb, "DOWN_", pop, "_", gene_nomenclature)
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
#' @importFrom stats quantile
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
                                                    breaks=unique(stats::quantile(input_matrix[input_matrix>0], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)),
                                                    include.lowest=TRUE))
    return(quantileValues)
  }
  specificity_quantiles = apply(input_matrix,2,FUN=bin.columns.into.quantiles,numberOfBins=100)
  rownames(specificity_quantiles) = rownames(input_matrix)
  return(specificity_quantiles)
}


#' Gene object renamer
#'
#' Takes a named-list (of gene sets) or a matrix (of gene ranking sets) and renames the sets (e.g. from "CD4+ T-cell" to "CD4_T_cell")
#'
#' @param gobj Either: 1. named list of character vectors (each vector is a gene-set, containing gene IDs), or 2. numeric matrix (each row is a gene, each column is a condition/group/set)
#'
#' @return List
#'
#' @export
gene_object_renamer <- function(gobj) {
  if (is.matrix(gobj) | is.list(gobj)) {
    if (is.list(gobj)) {
      # Then it's a list
      df <- data.frame(orig_names = names(gobj),
                       temp_names = paste0("GS", 1:length(gobj)))
      names(gobj) <- df$temp_names
      return(list(gobj = gobj, gsn = df)) # gsn stands for gene-set names
    } else {
      # Then it's a matrix
      df <- data.frame(orig_names = colnames(gobj),
                       temp_names = paste0("GS", 1:ncol(gobj)))
      colnames(gobj) <- df$temp_names
      return(list(gobj = gobj, gsn = df))
    }
  } else {
    stop("Gene expression object is not a list or matrix.")
  }
}


