#' Plot multiple results
#'
#' Plot heatmap of P-values per celltype from multiple methods (e.g. both MAGMA & LDSC results).
#'
#' @param results List of result dataframes. Each result dataframe needs at least two columns, one called Name and another called P_value)
#' @param methods Vector containing titles of all the result dataframes in the list (e.g. c("10pc_t_stat_ldsc", "ewce_magma_linear")).
#' @param main_title Title of the graph (default = "Untitled")
#' @param adj Method of correcting for multiple comparisons (default = "bonferroni", other options are NA or "BH" - for benj-hoch fdr).
#' @param top_n How many cells to plot (top N, default = FALSE, options are numbers below length(unique(celltypes)))
#'
#' @return List of plots
#'
#'
#' @export
plot_multiple <- function(results, methods, main_title = "Untitled", adj = "bonferroni", top_n = FALSE) {
  require(ggplot2)
  require(reshape2)

  merged_results <- merge(results[[1]][,c("Name", "P_value")],
                          results[[2]][,c("Name", "P_value")],
                          all = TRUE,
                          by = "Name",
                          suffixes = c("",2))
  for (i in 3:length(results)) {
    merged_results <- merge(merged_results,
                            results[[i]][,c("Name", "P_value")],
                            all = TRUE,
                            by = "Name",
                            suffixes = c("",i))
  }
  colnames(merged_results)[2:ncol(merged_results)] <- methods


  if (!isFALSE(top_n) && top_n >= nrow(merged_results)) {
    stop("The top_n parameter can't be equal or larger to the number of cell-types in your results.")
  }


  # Create var with cell names to eliminate (that aren't in top N)
  if (!isFALSE(top_n) && top_n < nrow(merged_results)) {
    avg_pval_per_row <- rowSums(merged_results[,2:ncol(merged_results)])/ncol(merged_results)
    cells_to_elim <- as.character(merged_results[order(avg_pval_per_row)[(top_n+1):nrow(merged_results)],1])
  }


  # This is the sort of dataframe that can be easily plotted with ggplot2
  processed_tidy <- melt(merged_results, id.vars = "Name")
  colnames(processed_tidy)[2:3] = c("Method", "P_value")


  if (!is.na(adj)) {
    for (i in 1:nrow(processed_tidy)) {
      if (!is.na(processed_tidy$P_value[i])) {
        if (p.adjust(processed_tidy$P_value, method=adj)[i]<0.05) {processed_tidy$signif[i] <- "••"}
        else {
          this_method = which(processed_tidy$Method == processed_tidy$Method[i])
          if (p.adjust(processed_tidy$P_value[this_method], method=adj)[which(this_method==i)] < 0.05) {
            processed_tidy$signif[i] <- "•"
          } else {
            processed_tidy$signif[i] <- ""
          }
        }
      } else {
        processed_tidy$signif[i] <- ""
      }
    }
  } else {
    processed_tidy$signif <- ""
  }
  processed_tidy[,3] <- -log10(processed_tidy[,3])
  colnames(processed_tidy)[3] <- "log10_P_value"


  # This is only required for the bonferroni plots.
  glo_pval = -log10(0.05/nrow(processed_tidy))
  loc_pval = -log10(0.05/nrow(merged_results))


  # Eliminate cells that are not in top N (if user requested that)
  if (!isFALSE(top_n) && top_n < nrow(merged_results)) {
    processed_tidy <- processed_tidy[!(processed_tidy$Name %in% cells_to_elim),]
    merged_results <- merged_results[!(merged_results$Name %in% cells_to_elim),]
  }



  p <- ggplot(data = processed_tidy, aes(x=Method, y=Name, fill=log10_P_value)) +
    scale_fill_gradientn(colours = c("skyblue4", "red3"), name= "Shades are:\n-log_10(P-value)") +
    geom_tile() + labs(title = main_title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.na(adj)) {
    if (adj == "bonferroni") {
      # This plots whether cell-types pass global Bonferroni
      p1 <- ggplot(data = processed_tidy, aes(x=Method, y=Name, fill=interaction(log10_P_value >= glo_pval))) +
        geom_tile() + labs(title = main_title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      # Same but shaded instead of binary
      p2 <- ggplot(data = processed_tidy, aes(x=Method, y=Name, fill=log10_P_value)) +
        scale_fill_gradientn(colours = c("skyblue4", "skyblue4", "skyblue4","red3", "red3", "red3", "red3", "red3"),
                             breaks=c(0,glo_pval-0.001,glo_pval-0.0001,glo_pval-0.000001,glo_pval,Inf),
                             name="-log_10(P-value)\nWhite = Bonferroni") +
        geom_tile() + labs(title = main_title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      # Shaded plot with stars for significance
      p3 <- ggplot(data = processed_tidy, aes(x=Method, y=Name, fill=log10_P_value)) +
        scale_fill_gradientn(colours = c("skyblue4", "red3"),
                             name=paste0("Shades are:\n-log_10(P-value)\nBonf_global=",round(glo_pval,digits=2)," (••)\nBonf_loc=",round(loc_pval,digits=2), " (•)")) +
        geom_tile() + geom_text(aes(label=signif), color="black", size=5) + labs(title = main_title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      return(list(merged_results, processed_tidy, p, p1, p2, p3))
    }
    if (adj != "bonferroni") {
      p1 <- ggplot(data = processed_tidy, aes(x=Method, y=Name, fill=log10_P_value)) +
        scale_fill_gradientn(colours = c("skyblue4", "red3"),
                             name=paste0("Shades are:\n-log_10(P-value)\n• loc_",adj,"\n•• glo_",adj)) +
        geom_tile() + geom_text(aes(label=signif), color="black", size=5) + labs(title = main_title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      return(list(merged_results, processed_tidy, p, p1))
    }
  }


  return(list(merged_results, processed_tidy, p))
}
