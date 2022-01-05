# gwascelltyper

**gwascelltyper** is an R package for testing the enrichment of GWAS risk signal in gene expression data. For any questions, email `alexandru.voda@seh.ox.ac.uk`

Current version: `0.9 beta` (24th of September 2019)

Changes log file to be added soon.

## The idea behind

![picture](https://i.imgur.com/FYG5u82.png)

**gwascelltyper** can compute cell-type-specificity matrices of genes or gene lists in several ways. Demos are included in the example folder.

Then, starting either from the lists computed above or those computed by yourself outside of this package, it can use 3 different methods (LDSC, MAGMA & SNPsea) to test the enrichment of GWAS risk loci.

## Aims & current state

![picture](https://i.imgur.com/2LOyFmF.png)

[Hilary K. Finucane et al. (2018)](https://www.nature.com/articles/s41588-018-0081-4) have used a t-statistic to score the specificity of genes to each cell type, then tested the enrichment of risk variants (from GWAS) in the top 10% genes for each cell type with LDSC.

[N. Skene et al. (2018)](https://www.nature.com/articles/s41588-018-0129-5) have developed [EWCE](https://github.com/NathanSkene/EWCE) to score the specificity of genes to each cell type, then tested the enrichment of risk variants in each cell-type with [MAGMA.Celltyping](https://ctg.cncr.nl/software/magma).


Right now, **gwascelltyper** users can easily compute the specificity scores with either `t-stat`/`EWCE`/`differential expression`, and can use the scores to test enrichment with either `LDSC top N%`, `MAGMA top N%`, `SNPsea top N%`, `MAGMA linear`, `SNPsea linear`.

Soon, we intend to add the following features:

* [`INRICH`](https://atgu.mgh.harvard.edu/inrich/) - the enrichment test
* More `differential-expression` features

## Pre-requisites & Installation

CRAN pre-requisites: `install.packages(c("data.table","doParallel", "R.utils"))`. `Git` and `conda` are needed to install the `LDSC` environment, but only if you want to use `LDSC`. `tar`, `gunzip` and `wget` are required to automatically download the `MAGMA` & `LDSC`-formatted 1000G & HapMap data, which are required for these analyses.

To install this package, run: `devtools::install_bitbucket("alexandru_ioan_voda/cell_typing_package")`

After installing the lightweight code through the previous command, you must run `download_dependencies(population = "eur")` (or `afr`, `eas`, etc. run `?download_dependencies()` for details) to download the required baseline data (LD reference panels & others, approx. 2 to 6GB on disk, depending on what populations and packages you need). After answering the questions this function asks, you will be able to run your analyses.

## Tutorial

### First interactive run

Before running any analyses, you need to open R, load the package (`library("gwascelltyper")`), and then download the reference LD panel data (run `download_dependencies()` and follow the instructions that it outputs).

After these steps, you can run the analysis either in batch mode or interactively with just the following lines:

```
library("gwascelltyper")

# load example pre-processed gene specificity data
data(gene_score_matrix_bpe_ewce)

# download pre-processed GWAS summary statistics for SNPsea
download.file("https://docs.google.com/uc?export=download&id=122wracczXXbbG34jE75m15lWNFPwilJf", destfile = "~/snpsea_CD.gwas")

# run enrichment test
results <- compute_SNPSEA_enrichment_linear(gene_score_matrix_bpe_ewce,
                                            sumstats_file = "~/snpsea_CD.gwas",
                                            number_of_threads = 2)

# view results
head(results)
tail(results)
```

This analysis should show that there is a strong enrichment of Crohn's disease GWAS risk variants in genes expressed specifically by immune cells (top hit: Monocytes).

### Breakdown

`library("gwascelltyper")` loads the package

`data(gene_score_matrix_bpe_ewce)` Loads a gene-to-celltype specificity score matrix (computed from the gene expression data from the Blueprint & ENCODE projects, taken from the SingleR package by Dvir Aran).

`download.file("https://docs.google.com/uc?export=download&id=122wracczXXbbG34jE75m15lWNFPwilJf", destfile = "~/snpsea_CD.gwas")` Downloads the preformatted file of Crohn's disease european GWAS summary statistics for SNPsea, which only uses genome-wide significant SNPs (also making the file a lot quicker for downloading purposes). MAGMA & LDSC use all the data, not just genome-wide significant SNPs. There are several functions in this package aiming to help you to process the summary statistics for each package (all of them start with `munge_sumstats_for` [...]).

```
results <- compute_SNPSEA_enrichment_linear(gene_score_matrix_bpe_ewce,
                                            sumstats_file = "~/snpsea_CD.gwas",
                                            number_of_threads = 2)
```

This function runs SNPsea with 2 threads, with the default dependencies.


### Preparing your own summary statistics

```
library("gwascelltyper")

munge_sumstats_for_MAGMA(sumstats_raw = "/path/to/your/raw_sumstats.txt", sumstats_munged = "/path/to/your/processed_for_MAGMA_sumstats.txt")

munge_sumstats_for_LDSC_from_MAGMA(sumstats_magma = "/path/to/your/processed_for_MAGMA_sumstats.txt", sumstats_munged = "/path/to/your/processed_for_LDSC_sumstats.txt",
                                  number_of_threads = 1)

munge_sumstats_for_SNPSEA_from_MAGMA(sumstats_magma = "/path/to/your/processed_for_MAGMA_sumstats.txt", sumstats_munged = "/path/to/your/processed_for_SNPsea_sumstats.txt")
```

### Computing your own gene specificity scores

To compute gene-to-tissue or gene-to-celltype specificity scores, you need two objects:

1. A numeric matrix containing normalized gene counts (e.g. TPM), with column names as sample IDs or cell IDs, and row names as gene IDs (HGNC).

2. A list of two objects, specifically called `annotLevels$level1class` and `annotLevels$level2class`. Each of these is a named character vector. The names corespond to sample/cell IDs, and the value coresponds to the cluster (cell-type) or sample origin.

Example data:

```

library("gwascelltyper")

data(exp_hpca)
data(annot_hpca)

gene_score_matrix_hpca <- compute_ewce_scores(exp = exp_hpca, annotLevels = annot_hpca, number_of_threads = 1)

# This will compute a gene-to-celltype specificity matrix just like the gene_score_matrix_bpe_ewce one, but for the HPCA expression dataset instead of the Blueprint/ENCODE (BPE) dataset.

?compute_t_stat_matrix()
# This function computes t-statistics for specificity of expression, as described in H Finucane et al. 2018 Nat. Genet.

?compute_t_stat_scores()
# This function is basically compute_t_stat_matrix() but with NA filters and normalization of t-statistics between 0 and 1 (min max scaling).

```


More vignettes to be added soon.

## Citation
**If you are publishing results from this package**, please cite the following works:

* If you are using the LDSC enrichment test: [Hilary K. Finucane et al. Nature Genetics, 2018.](https://www.nature.com/articles/s41588-018-0081-4) & [Bulik-Sullivan, et al. Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)
* If you are using the MAGMA enrichment test: [C.A. de Leeuw et al. PLOS Comp Bio, 2015.](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004219) & [N. Skene et al. Nature Genetics, 2018.](https://www.nature.com/articles/s41588-018-0129-5)
* If you are using the SNPsea enrichment test: [K. Slowikowski et al. Bioinformatics, 2014.](https://doi.org/10.1093/bioinformatics/btu326)
* If you are using the EWCE gene-scoring method: [N. Skene et al. Front. Neurosci., 2016.](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full)
* If you are using the differential expression module: [A. Butler et al. Nature Biotechnology, 2018](https://doi.org/10.1038/nbt.4096) (the single cell DE module, which uses Seurat) or [M. Love et al. BMC Genome Biology, 2014](https://doi.org/10.1186/s13059-014-0550-8) (the bulk tissue DE module, which uses DESeq2)
