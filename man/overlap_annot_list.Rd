% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDSC_Creatingannot_Functions.R
\name{overlap_annot_list}
\alias{overlap_annot_list}
\title{Find the overlapping SNPs between a list of dataframes containing genomic
regions and a genomic ranges object.}
\usage{
overlap_annot_list(list, query_GR, seqname_col, start_col, end_col, cores = 2)
}
\arguments{
\item{list}{List with dataframes containing genomic regions to be overlapped
with the query_GR.}

\item{query_GR}{Genomic ranges object which user wants to query -- this will
typically be the SNPs in the baseline model.}

\item{seqname_col}{Column name for column in inputted dataframes referencing
chromosome.}

\item{start_col}{Column name for column in inputted dataframes referencing
start BP for search.}

\item{end_col}{Column name for column in inputted dataframes referencing end
BP.}

\item{cores}{integer. Number of cores to parallelise across. Default = 2.}
}
\value{
List of dataframes with SNPs overlapping between input dataframes and
  query genomic ranges object.
}
\description{
Function that can, in principle, be used to overlap any list of dataframes
containing CHR, BP with a single genomic ranges object. In practice, this
function is used to find the overlap between a list of dataframes containing
genomic regions within an annotation (e.g. could be genomic coorgindates for
SNPs or genes found within the annotation) and the baseline model. Any
regions that are found to overlap will be assigned a value of 1 in the newly
created 'Binary' column.
}
