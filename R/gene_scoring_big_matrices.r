#' Make gene-set from gene-by-celltype specificity matrix
#'
#' Returns list of the top N\% (largest values) genes for each cell-type. Output is a list of named character vectors, where each character vector represents a cell-type.
#'
#' @param gene_score_matrix Gene-by-cell-type specificity matrix.
#' @param gene_selection_cutoff Selection cutoff as a percentile. Default = 10.
#'
#' @return List
#'
#' @export
make_geneset_from_score_matrix <- function(gene_score_matrix, gene_selection_cutoff = 10) {
  # Gene_sets needs to be a named list of character vectors.
  cat("Extracting gene sets with top", gene_selection_cutoff, "% genes.\n")
  gene_sets = list()
  for (type_index in 1:ncol(gene_score_matrix)) {
      gene_sets[[type_index]] = as.vector(names(gene_score_matrix[order(gene_score_matrix[,type_index], decreasing = TRUE),type_index][1:((gene_selection_cutoff*nrow(gene_score_matrix))/100)]))
      }
  print("Adding the controls (All_genes_control) to the gene_sets object.")
  gene_sets[[ncol(gene_score_matrix)+1]] = rownames(gene_score_matrix)
  print("Naming the vectors.")
  names(gene_sets) = c(colnames(gene_score_matrix), "All_genes_control")
  return(gene_sets)
}


#' Scale matrix between 0 and 1.
#'
#' Scales matrix by min and max between 0 and 1. Code from https://stackoverflow.com/a/5665644. Adjuvant function.
#'
#' @param m Numeric matrix
#'
#' @return Matrix
#'
#' @export
range01 <- function(m) {
  return((m - min(m, na.rm = TRUE))/(max(m, na.rm = TRUE)-min(m, na.rm = TRUE)))
}


#' Format celltype names
#'
#' Properly formats common symbols for MAGMA, LDSC & SNPsea.
#'
#' @param character_to_process Vector of characters.
#'
#' @return Character vector
#'
#' @export
format_celltype_names <- function(character_to_process) {
  character_to_process <- gsub("+", "plus", character_to_process, fixed = TRUE)
  character_to_process <- gsub("-", "minus", character_to_process, fixed = TRUE)
  character_to_process <- gsub(" ", "_", character_to_process, fixed = TRUE)
  # Using a bit of simple regex here:
  character_to_process <- gsub("[^A-Za-z0-9_]", "_", character_to_process) # Replace all special characters (that are not {alphanumeric and _}) with underscore
  character_to_process <- gsub("_{2,}", "_", character_to_process) # Replace repeated underscores with a single underscore
  return(character_to_process)
}

#' Read exp and annot
#'
#' Read expression and annotation data into the environment, then check the sample/cell tags for correspondence between exp and annot.
#'
#' @param obj Environment to read exp and annotLevels to.
#' @param exp_rds Path to RDS file containing the expression matrix.
#' @param annot_rds Path to RDS file containing the annotation dataframe.
#'
#' @return None
#'
#' @export
read_exp_and_annot_environ <- function(obj, exp_rds, annot_rds) {
  if (!file.exists(exp_rds)) stop("The given expression matrix RDS file does not actually exist. Halting execution.")
  cat("Reading exp object...\n")
  obj$exp <- readRDS(exp_rds)
  cat("The exp has", ncol(obj$exp), "cells and", nrow(obj$exp), "genes.\n")
  obj$annotLevels <- NULL
  if (file.exists(annot_rds)) {
    cat("Clustering and cell-type labelling supplied in arguments. Reading annotLevels file...\n")
    obj$annotLevels <- readRDS(annot_rds)
    cat("The annotLevels contains", nrow(obj$annot), "cells, with", length(unique(obj$annot[,2])), "specific cell types and", length(unique(obj$annot[,1])), "main (broader) cell types.\nProceeding to QC.\n")
    # Do sanity checking of the input
    sanity_check_for_exp_and_annot_environ(obj)
    gc()
    # Some code for subsampling here, if necessary
  } else {
    stop("Annot file does not exist. Clustering inside Multicelltyper not supported yet. Halting execution.")
  }
}

#' Sanity checking exp-annot correspondence
#'
#' Checking that exp and annot have 1-to-1 correspondence. Eliminating excesive/needless annots.
#'
#' @param obj Environment to read exp and annotLevels to.
#' @param format_cell_names Prepare cell/sample names for MAGMA, LDSC, SNPsea (default = TRUE).
#' @param verbose Whether to print intermediate work (default = TRUE).
#'
#' @return None
#'
#' @export
sanity_check_for_exp_and_annot_environ <- function(obj, format_cell_names = TRUE, verbose = TRUE) {
  if (any(obj$exp<0, na.rm = TRUE)) {stop("Fatal error: Negative values inside expression matrix.")}
  if (any(is.na(obj$annot))) {stop("Error, NAs present in obj$annot. That is not manageable input.")}
  # This makes the names of the samples and celltypes correctly formatted for MAGMA, LDSC & SNPsea.
  if (format_cell_names) {
    colnames(obj$exp) <- format_celltype_names(colnames(obj$exp))
    rownames(obj$annot) <- format_celltype_names(rownames(obj$annot))
    for (i in 1:ncol(obj$annot)) {
      obj$annot[,i] <- format_celltype_names(obj$annot[,i])
    }
  }

  if (verbose) print("Verifying that there are equal numbers of cells in the expression matrix and in the annotation list.") # Which happens to not be the case for GTEX

  # There needs to be a 1-to-1 correspondence between obj$annot and obj$exp. If there isn't, then we're missing something out.
  if (!all(c(rownames(obj$annot) %in% colnames(obj$exp), colnames(obj$exp) %in% rownames(obj$annot)))) {

    orphan_exp <- sum(!colnames(obj$exp) %in% rownames(obj$annot)) # How many cells in exp don't have an annotation
    orphan_annot <- sum(!rownames(obj$annot) %in% colnames(obj$exp)) # How many annotated cells don't have a corresponding column in exp

    if (orphan_exp>0) {
      if (verbose) {cat("Cell IDs missing an annotation:", colnames(obj$exp)[!colnames(obj$exp) %in% rownames(obj$annot)], "\n\n")}
      stop(paste0("There are ", orphan_exp," cells with no annotation present in your object! Something is very wrong. Cell IDs shown above, if verbose=TRUE."))
    }
    else {
      # then logically, it must be that orphan_annot is > 0.
      if (verbose) {
        cat("There are", orphan_annot, "annotations with no corresponding cell in exp. Removing them from downstream analyses.\n")
        cat("Cell IDs present in annotation but missing from exp:", rownames(obj$annot)[!rownames(obj$annot) %in% colnames(obj$exp)], "\n\n")
      }
      obj$annot <- obj$annot[rownames(obj$annot) %in% colnames(obj$exp),]

      if (!all(c(rownames(obj$annot) %in% colnames(obj$exp), colnames(obj$exp) %in% rownames(obj$annot)))) {
        stop("Still no 1-to-1 mapping between colnames(obj$exp) and rownames(obj$annot). Something is very wrong.")
      } else {
        if (verbose) print("Problem fixed. Equal number of cells in annot and exp.")
      }
    }
  }
}

#' Expression matrix (in environ object) missingness evaluation
#'
#' Warning: this function modifies any object it is given. This function does two things at once:
#' 1. automatically removes completely-NA rows and columns.
#' 2. finds out if some rows or columns have more missingness than a threshold (default = 0.1 (otherwise written, 10\%)), and eliminates them (if requested).
#'
#' It then returns a list of two vectors, containing the bad row & column indexes respectively.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param col_max_miss Column maximum missingness (default = 0.1)
#' @param row_max_miss Row maximum missingness (default = 0.1)
#' @param elim Whether to eliminate rows & cols with too many NAs (default = TRUE), alternative is to just throw a warning.
#' @param verbose Print what function has done.
#'
#' @return List (of two vectors: bad row & column indexes)
#'
#' @export
exp_environ_missingness <- function(obj, col_max_miss = 0.1, row_max_miss = 0.1, elim = TRUE, verbose = TRUE) {
  # Sanity checks:
  if (is.null(obj$annot) || is.null(obj$exp)) stop("obj$exp or obj$annot does not exist. Check with is.null(obj$exp) and is.null(obj$annot).")
  if (any(c(col_max_miss, row_max_miss) < 0, c(col_max_miss, row_max_miss) >= 1)) stop("col_max_miss & row_max_miss need to be numerics between [0; 1). Fatal error.")
  # Check if there's any NA.
  if(!any(is.na(obj$exp))) {
    if (verbose) print("No NA present in obj$exp.")
    return(NULL)
  }
  # Remove the completely-NA ones first.
  keep_col_indexes <- !(colSums(is.na(obj$exp)) == nrow(obj$exp))
  keep_row_indexes <- !(rowSums(is.na(obj$exp)) == ncol(obj$exp))
  obj$exp <- obj$exp[keep_row_indexes, keep_col_indexes]
  sanity_check_for_exp_and_annot_environ(obj, format_cell_names = FALSE, verbose = FALSE)

  if (verbose) cat(sum(!keep_row_indexes), "rows and", sum(!keep_col_indexes), "columns had been completely NA-filled. They have been removed from your object.\n")
  rm(keep_col_indexes, keep_row_indexes); gc()

  # Continue evaluating missingness
  if (col_max_miss != 1 || row_max_miss != 1) {
    bad_rows <- which(  rowMeans(is.na(obj$exp)) > row_max_miss  )
    bad_cols <- which(  colMeans(is.na(obj$exp)) > col_max_miss  )
    if (elim == TRUE) {
      if (length(bad_cols) > 0) obj$exp <- obj$exp[,-bad_cols]
      if (length(bad_rows) > 0) obj$exp <- obj$exp[-bad_rows,]
      if (verbose) cat(length(bad_rows), "rows and", length(bad_cols), "passed the missingness thresholds,", col_max_miss, "and", row_max_miss, "and have thus been removed from your object.\n")
      sanity_check_for_exp_and_annot_environ(obj, format_cell_names = FALSE, verbose = FALSE)
      return(list(bad_rows, bad_cols))
    } else {
      if (verbose) cat(length(bad_rows), "rows and", length(bad_cols), "passed the missingness thresholds,", col_max_miss, "and", row_max_miss, "and they have not been removed, on your request.\n")
      return(list(bad_rows, bad_cols))
    }
  } else {
    cat("Nothing to return, given you did not supply missingness parameters smaller than 100%.\n")
    return(NULL)
  }
}

#' Matrix missingness evaluation
#'
#' This function does two things at once:
#' 1. automatically removes completely-NA rows and columns.
#' 2. finds out if some rows or columns have more missingness than a threshold (default = 0.1 (otherwise written, 10%)), and eliminates them (if requested).
#'
#' It then returns your fixed matrix.
#'
#' @param mat Matrix
#' @param col_max_miss Column maximum missingness (default = 0.1)
#' @param row_max_miss Row maximum missingness (default = 0.1)
#' @param elim Whether to eliminate rows & cols with too many NAs (default = TRUE), alternative is to just throw a warning
#' @param verbose Print what function has done
#'
#' @return Matrix
#'
#' @export
matrix_missingness <- function(mat, col_max_miss = 0.1, row_max_miss = 0.1, elim = TRUE, verbose = TRUE) {
  # Sanity checks:
  if (is.null(mat)) stop("Matrix does not exist.")
  if (any(c(col_max_miss, row_max_miss) < 0, c(col_max_miss, row_max_miss) >= 1)) stop("col_max_miss & row_max_miss need to be numerics between [0; 1). Fatal error.")

  # Remove the completely-NA ones first.
  keep_col_indexes <- !(colSums(is.na(mat)) == nrow(mat))
  keep_row_indexes <- !(rowSums(is.na(mat)) == ncol(mat))
  mat <- mat[keep_row_indexes, keep_col_indexes]

  if (verbose) cat(sum(!keep_row_indexes), "rows and", sum(!keep_col_indexes), "columns had been completely NA-filled. They have been removed from your object.\n")

  # Continue evaluating missingness
  if (col_max_miss != 1 && row_max_miss != 1) {
    bad_rows <- which(  rowMeans(is.na(mat)) > row_max_miss  )
    bad_cols <- which(  colMeans(is.na(mat)) > col_max_miss  )
    if (elim == TRUE && length(bad_cols) > 0 && length(bad_rows) > 0) {
      mat <- mat[-bad_rows, -bad_cols]
      if (verbose) cat(length(bad_rows), "rows and", length(bad_cols), "passed the missingness thresholds,", col_max_miss, "and", row_max_miss, "and have thus been removed from your object.\n")
      return(mat)
    } else {
      if (verbose) cat(length(bad_rows), "rows and", length(bad_cols), "passed the missingness thresholds,", col_max_miss, "and", row_max_miss, ". Their IDs are (cols then rows):\n", bad_cols, "\n\n", bad_rows, "\n")
      return(mat)
    }
  } else {
    return(mat)
  }
}

#' Compute t-stat matrix - from environment-stored expression data
#'
#' Compute gene-by-celltype specificity matrix with the t-statistic. The input must be an environment containing an exp & annot.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param col_max_miss Maximum missingness (is.na) for the column (samples), between 0 and 1. Default is 10\% (col_max_miss=0.1). If a column contains more than 10\% NAs, it is eliminated from further analyses.
#' @param row_max_miss Maximum missingness (is.na) for the row (genes), between 0 and 1. Default is 10\% (row_max_miss=0.1). If a row contains more than 10\% NAs, it is eliminated from further analyses.
#'
#' @return Dataframe
#'
#' @export
compute_t_stat_matrix_environ <- function(obj, number_of_threads = 1, col_max_miss = 0.1, row_max_miss = 0.1) {
  # Sanity checks:
  if (number_of_threads <= 0 | !(is.numeric(number_of_threads))) {stop("Number_of_threads param from compute_t_stat_matrix_environ function needs to be a positive integer. Fatal error.")} else {number_of_threads <- as.integer(number_of_threads)}
  if (any(c(col_max_miss, row_max_miss) < 0, c(col_max_miss, row_max_miss) >= 1)) stop("col_max_miss & row_max_miss need to be numerics between [0; 1). Fatal error.")
  if (is.null(obj$annot) || is.null(obj$exp)) stop("obj$exp or obj$annot does not exist.")
  if (!all(c(rownames(obj$annot) %in% colnames(obj$exp), colnames(obj$exp) %in% rownames(obj$annot)))) stop("colnames(obj$exp) and rownames(obj$annot) are not in 1 to 1 mapping. Something is wrong.")

  # Evaluate/correct missingness in expression matrix
  if (any(is.na(obj$exp))) {exp_environ_missingness(obj)}
  sanity_check_for_exp_and_annot_environ(obj)
  cat(dim(obj$exp), "dims\n")

  require(doParallel)
  registerDoParallel(cores=number_of_threads)
  print("Computing t statistic for each gene per each cell cluster...")
  nrow_exp <- nrow(obj$exp) # Speeds up the printing processes if matrix has dozens of thousands of rows
  t_stats <- matrix(NA, nrow = nrow_exp, ncol = length(unique(obj$annot[,2])))
  colnames(t_stats) <- unique(obj$annot[,2])
  rownames(t_stats) <- rownames(obj$exp)

  type2_classes <- unique(obj$annot[,2])
  gene_index <- 1:nrow_exp
  t_stats = foreach(gene_index <- 1:nrow_exp, .combine=rbind, .options.snow=list(preschedule=TRUE)) %dopar% {

    # Get the t-statistic of this gene for each cell type
    if (gene_index %% 1000 == 0) system(paste0("echo Processed ", gene_index, " genes out of ", nrow_exp))

    row <- rep(NA,length(type2_classes))

    for (k in 1:length(type2_classes)){
      type2_class <- type2_classes[k]
      main_type <- obj$annot[match(type2_class, obj$annot[,2]),1]
      #The target set of expressions for the gene:
      target_expressions <- obj$exp[gene_index, which(obj$annot[,2] == type2_class)] # All cells of the same L2 type
      #The control set of expressions for the gene:
      control_expressions <- obj$exp[gene_index, which(obj$annot[,1] != main_type)] # All cells that are not of the same main type
      if (sum(!is.na(target_expressions)) == 1 & sum(!is.na(control_expressions)) == 1) next # If there is just one non-NA replicate in both at the same time
      if (sum(!is.na(target_expressions)) == 0 | sum(!is.na(control_expressions)) == 0) next # If there is no non-NA replicate in any of them
      if (class(try(stats::t.test(target_expressions,control_expressions,var.equal=T)[["statistic"]], silent=TRUE))=="try-error") next # If there is any other error in performing the t_stat, continue to the next pair
      row[k] <- stats::t.test(target_expressions,control_expressions,var.equal=T)[["statistic"]] # You need to have at least two in each group
    }
    row
  }
  colnames(t_stats) <- unique(obj$annot[,2])
  rownames(t_stats) <- rownames(obj$exp)
  print("Finished computing matrix of t statistics.")
  registerDoSEQ() # This unregisters the paralel backend of the 'parallel' library
  return(t_stats)
}

#' Compute EWCE matrix - from environment-stored expression data
#'
#' Compute gene-by-celltype specificity matrix with the EWCE score. Code modified from the EWCE package by Skene et al. (2016, DOI 10.3389/fnins.2016.00016).
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel mclapply
#' @importFrom stats median
#'
#' @return Dataframe
#'
#' @export
generate_celltype_data_environ <- function(obj, number_of_threads = 1) {
  require("parallel")
  cat("Starting multithreading with", number_of_threads, "cores.\n")
  cl <- parallel::makeCluster(number_of_threads)

  if (nrow(obj$annot) != ncol(obj$exp)) {
    stop("Error: number of cells in annot must equal the number of columns in exp matrix")
  }
  storage.mode(obj$exp) <- "double"
  ctd = list()
  for (i in 1:length(obj$annot)) {
    ann <- as.character(obj$annot[,i])
    names(ann) <- rownames(obj$annot)
    ctd[[length(ctd) + 1]] = list(annot = ann)
  }
  aggregate.over.celltypes <- function(rowOfMeans, celltypes, func = "mean") {
    if (func == "mean") {
      exp_out = as.matrix(data.frame(stats::aggregate(rowOfMeans,
                                               by = list(celltypes), FUN = mean)))
    }
    else if (func == "median") {
      exp_out = as.matrix(data.frame(stats::aggregate(rowOfMeans,
                                               by = list(celltypes), FUN = stats::median)))
    }
    rownames(exp_out) = exp_out[, "Group.1"]
    exp_out = exp_out[, 2]
    exp_out2 = as.numeric(exp_out)
    names(exp_out2) = names(exp_out)
    return(exp_out2)
  }
  calculate.meanexp.for.level <- function(ctd_oneLevel, expMatrix) {
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
      print(dim(expMatrix)[2])
      print(length(ctd_oneLevel$annot))
      if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) !=
          0) {
        stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
      }
      ctd_oneLevel$mean_exp = expMatrix
    }
    else {
      mean_exp = apply(expMatrix, 1, aggregate.over.celltypes,
                       ctd_oneLevel$annot)
      ctd_oneLevel$mean_exp = t(mean_exp)
    }
    return(ctd_oneLevel)
  }
  calculate.medianexp.for.level <- function(ctd_oneLevel, expMatrix) {
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
      print(dim(expMatrix)[2])
      print(length(ctd_oneLevel$annot))
      if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) !=
          0) {
        stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
      }
      ctd_oneLevel$median_exp = expMatrix
    }
    else {
      median_exp = apply(expMatrix, 1, aggregate.over.celltypes,
                         ctd_oneLevel$annot, func = "median")
      ctd_oneLevel$median_exp = t(median_exp)
    }
    return(ctd_oneLevel)
  }
  calculate.specificity.for.level <- function(ctd_oneLevel) {
    normalised_meanExp = t(t(ctd_oneLevel$mean_exp) * (1/colSums(ctd_oneLevel$mean_exp)))
    normalised_medianExp = t(t(ctd_oneLevel$median_exp) *
                               (1/colSums(ctd_oneLevel$mean_exp)))
    ctd_oneLevel$specificity = normalised_meanExp/(apply(normalised_meanExp,
                                                         1, sum) + 1e-12)
    ctd_oneLevel$median_specificity = normalised_medianExp/(apply(normalised_meanExp,
                                                                  1, sum) + 1e-12)
    return(ctd_oneLevel)
  }
  ctd = parallel::mclapply(ctd, calculate.meanexp.for.level, obj$exp, mc.cores = number_of_threads)
  ctd = parallel::mclapply(ctd, calculate.medianexp.for.level, obj$exp, mc.cores = number_of_threads)
  ctd = parallel::mclapply(ctd, calculate.specificity.for.level, mc.cores = number_of_threads)
  parallel::stopCluster(cl)
  #ctd = lapply(ctd, bin.specificity.into.quantiles, numberOfBins = 100)
  # This is the line that was giving errors with the Blueprint_ENCODE dataset ^
  # And it really ISN'T a mandatory feature for the output of this function. Using this function for specificity scores, not for specificity quantiles.

  return(ctd)
}

#' Compute t-statistic scores - from environment-stored expression data
#'
#' Compute gene-by-celltype specificity matrix with the t-statistic, then replaces the NA values for each row and finally scales the matrix between 0 and 1.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param col_max_miss Column maximum missingness (default = 0.1)
#' @param row_max_miss Row maximum missingness (default = 0.1)
#'
#' @return Matrix
#'
#' @export
compute_t_stat_scores_environ <- function(obj, number_of_threads = 1, col_max_miss = 0.1, row_max_miss = 0.1) {
  if (number_of_threads <= 0 | !(is.numeric(number_of_threads))) {stop("Number_of_threads param from compute_t_stat_scores_environ function needs to be a positive integer. Fatal error.")} else {number_of_threads <- as.integer(number_of_threads)}
  if (!all(c(rownames(obj$annot) %in% colnames(obj$exp), colnames(obj$exp) %in% rownames(obj$annot)))) stop("colnames(obj$exp) and rownames(obj$annot) are not in 1 to 1 mapping. Something is wrong.")

  t_stats <- compute_t_stat_matrix_environ(obj, number_of_threads, col_max_miss, row_max_miss)

  # Diagnostic:
  cat("Ncol(t_stats): ", ncol(t_stats), "Nrow(t_stats): ", nrow(t_stats), "Length(is.na(t_stats)): ", sum(is.na(t_stats)), "out of Length(t_stats): ", length(t_stats), "\n")

  if (any(is.na(t_stats))) {
    cat("The t_stats matrix has NA values on some rows, which isn't proper input for linear enrichment tests. There are", sum(rowSums(is.na(t_stats)) > 0), "rows that contain NAs. Out of these,", sum(rowSums(is.na(t_stats)) == ncol(t_stats)), "contain(s) *only* NAs. The NA-only rows will be removed, and the entries remaining that contain NAs will be replaced with their row's average. Example below:\n\n0.1, 0.2, 0.3\n0.4, NA, 0.6\nNA, NA, NA\n\nWill become\n\n0.1, 0.2, 0.3\n0.4, 0.5, 0.6\n\n")

    # This removes rows with missingness above threshold (otherwise said, if more than e.g. 10% of data on the row is missing, the row is eliminated from further computations). If this is set to 0, then any row that contains NAs will be eliminated. If this is set to 1, then all rows containing any amount of NAs will be kept.
    t_stats <- matrix_missingness(t_stats, col_max_miss = 0.1, row_max_miss = 0.1, elim = TRUE, verbose = TRUE)

    # This averages out NAs over rows for the remaining rows which contain some still
    if (any(is.na(t_stats))) {
      for (i in 1:nrow(t_stats)) {
        t_stats[i, is.na(t_stats[i,])] <- mean(t_stats[i,], na.rm = TRUE)
      }
    }

    if (any(is.na(t_stats))) {stop("Even after processing, there still are NAs in t_stats. Something weird is going on? Fatal error. If you can't figure this out on your own, email the package author with a minimally-reproducible example leading to this error.")}
  } else {
    print("No NAs in t_stats. Awesome!")
  }

  print("Scaling t statistics matrix between 0 (least specific) and 1 (most specific).")
  # The minimum of c(1,2,3,NA) is NA, so define NA as the smallest NON-NA. Also, the max of c(1,2,3,NA) is also NA!
  norm_t_stats <- range01(t_stats)

  # Replace all NAs with zero specificity.
  #print("Replacing specificities marked as NA with zeros.")
  #norm_t_stats[is.na(norm_t_stats)] <- 0

  return(norm_t_stats) # Returns L2 specificity scores (numeric matrix, where cols=celltypes, rows=genes)
}

#' Compute EWCE scores - from environment-stored expression data
#'
#' Compute gene-by-celltype specificity matrix with the EWCE score. Code modified from the EWCE package by Skene et al. (2016, DOI 10.3389/fnins.2016.00016). If the expression matrix contains NAs, then this function removes NA-only rows and replaces the remaining NAs with the row average.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#' @param col_max_miss Column maximum missingness (default = 0.1)
#' @param row_max_miss Row maximum missingness (default = 0.1)
#'
#' @return Matrix
#'
#' @export
compute_ewce_scores_environ <- function(obj, number_of_threads = 1, col_max_miss = 0.1, row_max_miss = 0.1) {
  # EWCE is throwing errors when encountering NA values. This deals with it before launching EWCE.
  if (any(is.na(obj$exp))) {
    cat("The exp has NA values on some rows, which isn't proper input for EWCE. There are", sum(rowSums(is.na(obj$exp)) > 0), "rows that contain NAs. Out of these,", sum(rowSums(is.na(obj$exp)) == ncol(obj$exp)), "contain(s) *only* NAs. The NA-only rows will be removed, and the entries remaining that contain NAs will be replaced with their row's average. Example below:\n\n1, 2, 3\n4, NA, 6\nNA, NA, NA\n\nWill become\n\n1, 2, 3\n4, 5, 6\n\n")

    exp_environ_missingness(obj)

    # This averages out NAs over rows for the remaining rows which contain some still
    for(i in 1:nrow(obj$exp)) {
      obj$exp[i, is.na(obj$exp[i,])] <- mean(obj$exp[i,], na.rm = TRUE)
    }
    if (any(is.na(obj$exp))) {stop("Even after processing, there still are NAs in exp. Something weird is going on?")}
  } else {
    print("No NAs in exp. Awesome!")
  }

  print("Running generate.celltype.data().")
  ctd <- generate_celltype_data_environ(obj, number_of_threads)
  return(ctd[[2]][["specificity"]]) # Returns L2 specificity scores (numeric matrix, where cols=celltypes, rows=genes)
}

#' Compute top N\% t-statistic gene-sets - from environment-stored expression data
#'
#' Compute top N\% t-statistic ranked genes for each cell-type. Output is a list of named character vectors, where each character vector represents a cell-type.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param gene_selection_cutoff Selection cutoff as a percentile. Default = 10.
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#'
#' @return List
#'
#' @export
compute_t_stat_selection_environ <- function(obj, gene_selection_cutoff = 10, number_of_threads = 1) {
  gene_score_matrix <- compute_t_stat_matrix_environ(obj, args)
  return(make_geneset_from_score_matrix(gene_score_matrix, gene_selection_cutoff))
}

#' Compute top N\% EWCE gene-sets - from environment-stored expression data
#'
#' Compute top N\% EWCE ranked genes for each cell-type. Output is a list of named character vectors, where each character vector represents a cell-type.
#'
#' @param obj Environment object containing an exp (numeric matrix of genes-by-cell-ID) & an annot (dataframe of cell annotations).
#' @param gene_selection_cutoff Selection cutoff as a percentile. Default = 10.
#' @param number_of_threads Number of threads to parallelize calculation over. Default = 1.
#'
#' @return List
#'
#' @export
compute_ewce_selection_environ <- function(obj, gene_selection_cutoff = 10, number_of_threads = 1) {
  gene_score_matrix <- compute_ewce_scores_environ(obj, number_of_threads)
  return(make_geneset_from_score_matrix(gene_score_matrix, gene_selection_cutoff))
}

