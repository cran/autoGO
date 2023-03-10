% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering_DE.R
\name{filtering_DE}
\alias{filtering_DE}
\title{Filtering DESeq2 results}
\usage{
filtering_DE(
  padj_threshold = 0.05,
  log2FC_threshold = 1,
  outfolder = "./results",
  save_excel = FALSE
)
}
\arguments{
\item{padj_threshold}{(Default = 0.05) Threshold value for adjusted p-value filtering.}

\item{log2FC_threshold}{(Default = 0) Threshold value for log2(Fold Change) filtering.}

\item{outfolder}{(Default = "./results") Name of the folder in which the new output is written.}

\item{save_excel}{(Default = FALSE) Write output in MS Excel file format (.xlsx).}
}
\value{
No return value. Files will be produced as part of normal execution.
}
\description{
We could be in a position to carry out multiple filters on the results of the differential analysis, in order not to repeat all the deseq_analysis.R code which provides for the actual computation of the differential analysis, the filtering_DE.R function has been implemented to be able to filter the file(s) " * _allres.tsv " and generate all the folders and files associated with the specific filters applied.

The function automatically searches inside the folders where_results and outfolder the file(s) (See ?deseq_analysis()) "_allres.tsv" and generates folders and files in the same folders with the new filters for foldchange and pvalue respectively.
}
\examples{
\dontrun{
filtering_DE(
  padj_threshold = 0.05,
  log2FC_threshold = 1,
  outfolder = "./results",
  save_excel = F
)
}
}
