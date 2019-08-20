#' Function to calculate LD score for a number of annotations.
#'
#' Note: For this to be able to run, annot.gz files should be stored in:
#' $ANNOTATION_BASEDIR/$ANNOT_NAME/$ANNOTATION_SUBCATEGORY/ where
#' $ANNOTATION_BASEDIR refers to the common directory, $ANNOT_NAME refers to the
#' author of article from which the annotation is derived and
#' $ANNOTATION_SUBCATEGORY refers to the subcategories within the annotation.
#' Function should be run in the folder $ANNOT_NAME/
#'
#' @param Annotation_Basedir Common directory wherein specific annotation is
#'   stored.
#' @param Annot_name  Name of annotation. This should be written exactly as
#'   written in the annotation directory name.
#' @param Annotation_Subcategories Annotation subcategories. If there is more
#'   than 1, these should be written precisely as in the folder, and should be
#'   entered with a comma separator without any spaces.
#' @param Fixed_Arguments List of fixed arguments. This can be created using the
#'   get_LDSC_fixed_args() function.
#'
#' @return Will run LDscore script for all annotations and chromosomes.
#' @export
#'

Calculate_LDscore <- function(Annotation_Basedir = NULL, Annot_name = NULL, Annotation_Subcategories = NULL, Fixed_Arguments = NULL){

  # Arguments
  fixed_args <- Fixed_Arguments
  current_dir <- paste0(Annotation_Basedir, Annot_name, "/")

  # Loop for annotation subcategories
  for(i in seq_along(Annotation_Subcategories)){

    annot_prefix <- paste0(Annotation_Basedir, Annot_name, "/", Annotation_Subcategories[i], "/", Annotation_Subcategories[i])
    out_prefix <- paste0(Annotation_Subcategories[i], ".")

    # Loop for chromsomes within annotation subcategories
    for(CHR in 1:22){

      # Print annotation and chromosome for tracking of progress
      print(str_c("Annotation: ", Annotation_Subcategories[i], ", CHR: ", CHR))

      # Creating necessary sub arguments with chromosome number attached
      refLD_chr <- paste0(fixed_args$refLD_basedir, "1000G.EUR.QC.", CHR)
      annot_chr <- paste0(annot_prefix, ".", CHR,".annot.gz")
      out_chr <- paste0(out_prefix, CHR)
      printsnps_chr <- paste0(fixed_args$printsnps_basedir, "SNPsinBaselinel2ld.", CHR, ".snp")

      # Creating entire argument
      LDscoreARG <- paste0(" ", fixed_args$ldsc,
                           " --l2 --bfile ", refLD_chr,
                           " --ld-wind-cm 1",
                           " --annot ", annot_chr,
                           " --out ", out_chr,
                           " --print-snps ", printsnps_chr)

      print(LDscoreARG)

      # Running command
      system2(command = "python", args = LDscoreARG)

    }

  }

  # Move results to appropriate directories
  for(i in seq_along(Annotation_Subcategories)){

    list.of.files <- list.files(current_dir, paste0(Annotation_Subcategories[i], "."), full.names = T)
    new.location <- paste0(Annotation_Basedir, Annot_name, "/", Annotation_Subcategories[i], "/")

    file.copy(list.of.files, new.location)

    file.remove(list.of.files)

  }


}

#' Function to run LDSC heritability script for a number of annotations and
#' GWASs.
#'
#' #' Note: For this to be able to run, annot.gz files should be stored in:
#' $ANNOTATION_BASEDIR/$ANNOT_NAME/$ANNOTATION_SUBCATEGORY/ where
#' $ANNOTATION_BASEDIR refers to the common directory, $ANNOT_NAME refers to the
#' author of article from which the annotation is derived and
#' $ANNOTATION_SUBCATEGORY refers to the subcategories within the annotation.
#' Function should be run in the folder $ANNOT_NAME/
#'
#' @param Annotation_Basedir Common directory wherein specific annotation is
#'   stored.
#' @param Annot_name  Name of annotation. This should be written exactly as
#'   written in the annotation directory name.
#' @param Annotation_Subcategories Annoation subcategories. If there is more
#'   than 1, these should be written precisely as in the folder, and should be
#'   entered with a comma separator without any spaces.
#' @param Fixed_Arguments List of fixed arguments. This can be created using the
#'   get_LDSC_fixed_args() function.
#' @param GWAS Dataframe of GWAS to be run in the H2 estimation, as generated
#'   using the Create_GWAS_df() function. Columns should include: Full.paths
#'   (full paths to GWAS), Original.name (original name of sumstat.gz file),
#'   Output.prefix (alternative output name).
#'
#' @return Will run heritability script for all annotations and GWASs.
#' @export
#'

Calculate_H2 <- function(Annotation_Basedir = NULL, Annot_name = NULL, Annotation_Subcategories = NULL, Fixed_Arguments = NULL, GWAS_df = NULL){

  # Arguments
  fixed_args <- Fixed_Arguments
  current_dir <- paste0(Annotation_Basedir, Annot_name, "/")

  # Make directory for output
  directory <- paste0(Annotation_Basedir, Annot_name, "/", "Output")

  if(dir.exists(directory)){
    cat("Directory exists")
  } else {
    cat("Directory does not exist -- creating")
    system2(command = "mkdir", args = directory)
  }

  # Loop for annotation subcategories
  for(i in seq_along(Annotation_Subcategories)){

    annot_prefix <- paste0(Annotation_Basedir, Annot_name, "/", Annotation_Subcategories[i], "/", Annotation_Subcategories[i], ".")
    out_suffix <- paste0(Annotation_Subcategories[i])

    # Loop for GWASs
    for(j in 1:nrow(GWAS_df)){

      # Creating necessary sub arguments
      h2 <- as.character(GWAS_df$Full.paths[j])
      out_prefix <- as.character(GWAS_df$Output.prefix[j])

      # Creating entire argument
      H2ARG <- paste0(" ", fixed_args$ldsc,
                      " --h2 ", h2,
                      " --w-ld-chr ", fixed_args$regression_weights,
                      " --ref-ld-chr ", annot_prefix, ",", fixed_args$baseline_annot_path,
                      " --overlap-annot ",
                      " --frqfile-chr ", fixed_args$freqfile,
                      " --out ", out_prefix, "_", Annot_name, "_", fixed_args$baseline_model_name, "baseline_", out_suffix, # baseline_model_name to indicate which one is in use
                      " --print-coefficients ")


      # Print annotation and chromosome for tracking of progress
      print(str_c("Annotation: ", Annotation_Subcategories[i], ", GWAS: ", out_prefix))

      print(H2ARG)

      # Running command
      system2(command = "python", args = H2ARG)

    }

  }

  # Move results to appropriate directories
  list.of.files.results <- list.files(current_dir, ".results", full.names = T)
  list.of.files.log <- list.files(current_dir, ".log", full.names = T)

  file.copy(list.of.files.results, directory)
  file.copy(list.of.files.log, directory)

  file.remove(list.of.files.results)
  file.remove(list.of.files.log)

}

#' Create GWAS dataframe.
#'
#' Function that creates a GWAS dataframe based on all of the .sumstats.gz files available in the LDscore GWAS directory.
#'
#' @return Dataframe of availabe GWASs in the LDscore GWAS directory.
#' @export
#'

Create_GWAS_df <- function(){

  GWAS.paths <- list.files(path = "/data/LDScore/GWAS/", pattern = ".sumstats.gz", recursive = TRUE)

  Full.paths <- paste0("/data/LDScore/GWAS/", GWAS.paths)
  Original.name <- GWAS.paths %>%
    str_replace(".*/", "") %>%
    str_replace(".sumstats.gz", "")
  Output.prefix <- GWAS.paths %>%
    str_replace(".*/", "") %>%
    str_replace(".sumstats.gz", "")

  # Change names of those .sumstats.gz files that include '_', so that output name works with Assimilate_H2_results() function
  # (i.e. provides full name of the GWAS, otherwise anything after an '_' would be removed)
  Output.prefix[Output.prefix == "ILAE_All_Epi_11.8.14"] <- "EPI2014"
  Output.prefix[Output.prefix == "Epilepsy2018"] <- "EPI2018"
  Output.prefix[Output.prefix == "epipgx_dre_hc"] <- "EPI.Dre.Hc"
  Output.prefix[Output.prefix == "epipgx_dre_res"] <- "EPI.Dre.Res"
  Output.prefix[Output.prefix == "PD"] <- "PD2017"
  Output.prefix[Output.prefix == "SCZ"] <- "SCZ2014"
  Output.prefix[Output.prefix == "ALS"] <- "ALS2014"
  Output.prefix[Output.prefix == "PD2018_meta5_ex23andMe"] <- "PD2018.ex23andMe"
  Output.prefix[Output.prefix == "PD2018_AOO"] <- "PD2018.AOO"



  GWAS_df <- data.frame(Full.paths) %>%
    bind_cols(data.frame(Original.name)) %>%
    bind_cols(data.frame(Output.prefix))

  return(GWAS_df)

}


#' Get fixed LDSC arguments.
#'
#' Function to return a list of the fixed arguments. To be used in combination
#' with the Calculate_LDscore() and Calculate_H2() functions.
#'
#' @param Baseline_model Choose baseline model to run heritability estimates
#'   with. Choose by specifying number of annotations: 53, 75, 76, 86, or 97.
#'
#' @return List of fixed LDSC arguments.
#' @export
#'

get_LDSC_fixed_args <- function(Baseline_model = NULL){

  fixed_args <- list()

  # Fixed paths to python programme and reference files
  fixed_args$ldsc <- "/tools/LDScore-master/ldsc.py"
  fixed_args$refLD_basedir <- "/data/LDScore/Reference_Files/1000G_EUR_Phase3_plink/"
  fixed_args$regression_weights <- "/data/LDScore/Reference_Files/weights_hm3_no_hla/weights."
  fixed_args$freqfile <- "/data/LDScore/Reference_Files/1000G_Phase3_frq/1000G.EUR.QC."

  if (!is.null(Baseline_model)) {

    # Directory path to baseline models
    # Baseline model argument
    if(!Baseline_model %in% c("53", "75", "76", "86", "97")){
      stop("Invalid baseline model specified. Please use either denote which baseline model to use: 53, 75, 76, 86, 97.")
    }

    # Choose baseline model
    fixed_args$baseline_model_name <- as.character(Baseline_model)

    # Read in correct baseline model depending on model specified in argument baseline.model
    if (fixed_args$baseline_model_name == "53") {
      # Need to overwrite refLD path and weights path, as v1.2 is aligned to GRCh38 and therefore cannot use original reference files.
      fixed_args$baseline_annot_path <- "/data/LDScore/Reference_Files/GRCh38/baseline_v1.2/baseline."
      fixed_args$refLD_basedir <- "/data/LDScore/Reference_Files/GRCh38/plink_files/"
      fixed_args$printsnps_basedir <- "/data/LDScore/Reference_Files/GRCh38/baseline_v1.2/SNPsinBaselinel2ldscore/"
      fixed_args$regression_weights <- "/data/LDScore/Reference_Files/GRCh38/weights/weights.hm3_noMHC."
    }

    if (fixed_args$baseline_model_name == "75") {
      fixed_args$baseline_annot_path <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v1.1_ldscores_75annot/baselineLD."
      fixed_args$printsnps_basedir <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v1.1_ldscores_75annot/SNPsinBaselinel2ldscore/"
    }

    if (fixed_args$baseline_model_name == "76") {
      fixed_args$baseline_annot_path <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.0/baselineLD."
      fixed_args$printsnps_basedir <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.0/SNPsinBaselinel2ldscore/"
    }

    if (fixed_args$baseline_model_name == "86") {
      fixed_args$baseline_annot_path <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.1_ldscores/baselineLD."
      fixed_args$printsnps_basedir <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.1_ldscores/SNPsinBaselinel2ldscore/"
    }

    if (fixed_args$baseline_model_name == "97") {
      fixed_args$baseline_annot_path <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."
      fixed_args$printsnps_basedir <- "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.2_ldscores/SNPsinBaselinel2ldscore/"
    }

  }

  return(fixed_args)

}
