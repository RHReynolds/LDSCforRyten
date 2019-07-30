#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(R.utils)
library(tidyverse)
library(optparse)
library(stringr)
library(LDSCforRyten)

#---Parse_args---------------------------------------------------------------------------------------------------------------------------####

arguments <- parse_args(OptionParser(usage = "%prog",
                                     prog = "LDSC",
                                     description="Script for running LDSC with multiple annotations and multiple GWASs. \n For script to run, annot.gz files should be stored in: $ANNOTATION_BASEDIR/$ANNOTATION_NAME/$ANNOTATION_SUBCATEGORY/. Function should be run in the directory $ANNOTATION_NAME. \n Required inputs:\n <ANNOTATION_BASEDIR>: Common directory wherein specific annotation is stored.\n <ANNOTATION_NAME>: Annotation name.\n <ANNOTATION_SUBCATEGORY>: Subcategories within the annotation. Each subcategory should be the same name as the directory in which annot.gz files are stored. Subcategories should be separated by a comma and no spaces.",
                                     option_list=list(
                                       make_option(c("-g","--GWAS"), default = "all", help = "If select GWASs to be run, use this flag. GWAS names supplied should be the same as their sumstat.gz name (excluding 'sumstat.gz'). Names should be separated by a comma and no spaces. E.g. To run with the schizophrenia GWAS from 2018, write SCZ2018. [default: all GWASs run"),
                                       make_option(c("-b","--baseline_model"), default = "53", help = "If a different baseline is desired for estimation of heritability, call this flag. Choose between: 53, 75, 76, 86 or 97. [default: Empty string]"),
                                       make_option(c("-l", "--calculate_LD_alone"), default = FALSE, help = "If user only wants to calculate LD scores for an annotation, argument 'TRUE' must be supplied [default: FALSE]"),
                                       make_option(c("-c", "--calculate_h2_alone"), default = FALSE, help = "If user only wants to calculate heritability estimates, argument 'TRUE' must be supplied [default: FALSE]")
                                     )),
                        positional_arguments = 3)

# # Comment in if want to test run script
# arguments <- list()
# arguments$args[1] <- "/home/rreynolds/data/Extra_Annotations/"
# arguments$args[2] <- "Dystonia2019"
# arguments$args[3] <- "blue,cyan"
# arguments$opt$GWAS <- "ADHD2017_EUR,Anxiety2016.cc"
# arguments$opt$baseline_model <- "53"
# arguments$opt$calculate_LD_alone <- FALSE
# arguments$opt$calculate_h2_alone <- FALSE


#---Defining arguments-------------------------------------------------------------------------------------------------------------------####

#---Positional parse_args arguments---
annot_basedir <- arguments$args[1] %>% as.character()
annot_name <- arguments$args[2] %>% as.character
annot_subcategories <- arguments$args[3] %>% as.character %>% str_split(",") %>% unlist()

#---Optional parse_args arguments---
opt <- arguments$opt

# If statement to set baseline_model if optional argument called.
if(!is.null(opt$baseline_model)){
  baseline_model <- opt$baseline_model
} else {
  baseline_model <- "53"
}

# Create dataframe with file paths and output names of GWAS.sumstats.gz
gwas_df <- Create_GWAS_df()

# Filter out necessary GWAS
if (is.null(opt$GWAS)) {
  gwas_df
} else {
  GWAS_selection <- unlist(as.character(opt$GWAS) %>% str_split(pattern = ","))

  gwas_df <- gwas_df %>%
    filter(Original.name %in% GWAS_selection)
}

#---Fixed arguments---
fixed_args <- get_LDSC_fixed_args(Baseline_model = baseline_model)

#---Main---------------------------------------------------------------------------------------------------------------------------------####
if(opt$calculate_LD_alone == FALSE & opt$calculate_h2_alone == FALSE | opt$calculate_LD_alone == TRUE & opt$calculate_h2_alone == TRUE) {

  Calculate_LDscore(Annotation_Basedir = annot_basedir, Annot_name = annot_name, Annotation_Subcategories = annot_subcategories, Fixed_Arguments = fixed_args)
  Calculate_H2(Annotation_Basedir = annot_basedir, Annot_name = annot_name, Annotation_Subcategories = annot_subcategories, Fixed_Arguments = fixed_args, GWAS_df = gwas_df)

}

if(opt$calculate_LD_alone == TRUE & opt$calculate_h2_alone == FALSE) {

  Calculate_LDscore(Annotation_Basedir = annot_basedir, Annot_name = annot_name, Annotation_Subcategories = annot_subcategories, Fixed_Arguments = fixed_args)

}

if(opt$calculate_LD_alone == FALSE & opt$calculate_h2_alone == TRUE) {

  Calculate_H2(Annotation_Basedir = annot_basedir, Annot_name = annot_name, Annotation_Subcategories = annot_subcategories, Fixed_Arguments = fixed_args, GWAS_df = gwas_df)

}


