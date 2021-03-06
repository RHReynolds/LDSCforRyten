% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDSC_Pipeline_Functions.R
\name{Calculate_H2}
\alias{Calculate_H2}
\title{Function to run LDSC heritability script for a number of annotations and
GWASs.}
\usage{
Calculate_H2(
  Annotation_Basedir = NULL,
  Annot_name = NULL,
  Annotation_Subcategories = NULL,
  Fixed_Arguments = NULL,
  GWAS_df = NULL,
  cores = 1
)
}
\arguments{
\item{Annotation_Basedir}{Common directory wherein specific annotation is
stored.}

\item{Annot_name}{Name of annotation. This should be written exactly as
written in the annotation directory name.}

\item{Annotation_Subcategories}{Annoation subcategories. If there is more
than 1, these should be written precisely as in the folder, and should be
entered with a comma separator without any spaces.}

\item{Fixed_Arguments}{List of fixed arguments. This can be created using the
get_LDSC_fixed_args() function.}

\item{GWAS_df}{Dataframe of GWAS to be run in the H2 estimation, as generated
using the Create_GWAS_df() function. Columns should include: Full.paths
(full paths to GWAS), Original.name (original name of sumstat.gz file),
Output.prefix (alternative output name).}

\item{cores}{integer. Number of cores to parallelise across. Default = 1.}
}
\value{
Will run heritability script for all annotations and GWASs.
}
\description{
#' Note: For this to be able to run, annot.gz files should be stored in:
$ANNOTATION_BASEDIR/$ANNOT_NAME/$ANNOTATION_SUBCATEGORY/ where
$ANNOTATION_BASEDIR refers to the common directory, $ANNOT_NAME refers to the
author of article from which the annotation is derived and
$ANNOTATION_SUBCATEGORY refers to the subcategories within the annotation.
Function should be run in the folder $ANNOT_NAME/
}
