#---First level code---------------------------------------------------------------------------------------------------------------------####
#' Function assembling list of SNPs in baseline LD score model
#'
#' @param baseline_model Specify the baseline model to be used i.e. 53, 75, 76,
#'   86 or 97 annotation model. Recommendations for use in cell-type specific
#'   analyses from Alkes Price group: 1. We recommend that for identifying
#'   critical tissues/cell-types via P-value of tau, it is best to use the
#'   baseline model, specifically baseline v1.2. 2. We recommend that for
#'   estimating heritability enrichment (i.e., %h2/%SNPs) of any annotation,
#'   including tissue-specific annotations, it is best to use baselineLD v2.2.
#'
#' @return Dataframe of baseline SNPs with BP and CHR.
#' @export
#'
creating_baseline_df <- function(baseline_model = c("53", "75", "76", "86", "97")){

  # Read in correct baseline model depending on model specified in argument baseline_model
  if (baseline_model == "53") {
    file.paths <- list.files(path = "/data/LDScore/Reference_Files/GRCh38/baseline_v1.2/",
                                 pattern = "annot.gz", full.names = T)
  }

  if (baseline_model == "75") {
    file.paths <- list.files(path = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v1.1_ldscores_75annot",
                             pattern = "annot.gz", full.names = T)
  }

  if (baseline_model == "76") {
    file.paths <- list.files(path = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.0",
                             pattern = "annot.gz", full.names = T)
  }

  if (baseline_model == "86") {
    file.paths <- list.files(path = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.1_ldscores",
                             pattern = "annot.gz", full.names = T)
  }

  if (baseline_model == "97") {
    file.paths <- list.files(path = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.2_ldscores",
                             pattern = "annot.gz", full.names = T)
  }


  for(i in seq_along(file.paths)){

    # Read in file with only CHR, BP, SNP and CM columns
    file <- readr::read_delim(file.paths[i],
                              delim = "\t",
                              col_types = readr::cols_only(CHR = readr::col_integer(),
                                                           BP = readr::col_integer(),
                                                           SNP = readr::col_character(),
                                                           CM = readr::col_double()))

    if(i == 1){
      Master <- file
    } else{
      Master <- dplyr::bind_rows(Master, file)
    }

  }

  return(Master)

}

#' Find the overlapping SNPs between a list of dataframes containing genomic
#' regions and a genomic ranges object.
#'
#' Function that can, in principle, be used to overlap any list of dataframes
#' containing CHR, BP with a single genomic ranges object. In practice, this
#' function is used to find the overlap between a list of dataframes containing
#' genomic regions within an annotation (e.g. could be genomic coorgindates for
#' SNPs or genes found within the annotation) and the baseline model. Any
#' regions that are found to overlap will be assigned a value of 1 in the newly
#' created 'Binary' column.
#'
#' @param list List with dataframes containing genomic regions to be overlapped
#'   with the query_GR.
#' @param query_GR Genomic ranges object which user wants to query -- this will
#'   typically be the SNPs in the baseline model.
#' @param seqname_col Column name for column in inputted dataframes referencing
#'   chromosome.
#' @param start_col Column name for column in inputted dataframes referencing
#'   start BP for search.
#' @param end_col Column name for column in inputted dataframes referencing end
#'   BP.
#' @param cores integer. Number of cores to parallelise across. Default = 2.
#'
#' @return List of dataframes with SNPs overlapping between input dataframes and
#'   query genomic ranges object.
#' @export
#'

overlap_annot_list <- function (list, query_GR, seqname_col, start_col, end_col, cores = 2) {

  # Run in parallel
  cl <- parallel::makeCluster(cores)

  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()

  # Run the cluster
  master_list <- foreach::foreach(i = 1:length(list),
                         .verbose = TRUE,
                         .packages = c("LDSCforRyten", "tidyverse", "stringr")) %dopar% {
                           annotation <- names(list)[i]
                           Subject_GR <-
                             list[[i]] %>%
                             LDSCforRyten::df2GRanges(seqname_col = seqname_col,
                                                      start_col = start_col,
                                                      end_col = end_col)
                           overlap <- as.data.frame(GenomicRanges::findOverlaps(query_GR, Subject_GR))
                           hits <- cbind(as.data.frame(query_GR)[overlap$queryHits, 6],
                                         as.data.frame(Subject_GR)[overlap$subjectHits, ])
                           colnames(hits)[1] <- c("SNP")

                           hits %>%
                             dplyr::mutate(Binary = 1)
                         }

  # Stop cluster
  parallel::stopCluster(cl)

  names(master_list) <- names(list)
  return(master_list)

}

#' Map annotation SNPs overlapping the baseline model back to the full baseline
#' model.
#'
#' Function is used to map annotation SNPs overlapping the baseline model back
#' to the full baseline model, with a column of 1s and 0s to distinguish between
#' SNPs overlapping the baseline model (1) and those that do not overlap (0).
#' Function is intended for use after using the overlap_annot_list() function.
#'
#' @param list_of_annotations List of annotation dataframes containing only
#'   those SNPs that were found to overlap the baseline model using the
#'   overlap_annot_list() function.
#' @param BM Baseline model with columns: CHR, BP, SNP, CM.
#' @param cores integer. Number of cores to parallelise across. Default = 2.
#'
#' @return List of annotations with all SNPs present in the baseline model, in
#'   addition to a column distinguishing between annotation SNPs present in the
#'   baseline model (1) and those not found within the baseline model (0).
#' @export
#'

overlap_annot_hits_w_baseline <- function(list_of_annotations, BM, cores = 2){

  # Run in parallel
  cl <- parallel::makeCluster(cores)

  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()

  master_list <- foreach::foreach(i = 1:length(list_of_annotations),
                         .verbose = TRUE,
                         .packages = c("LDSCforRyten", "tidyverse", "stringr")) %dopar% {
                           df <- list_of_annotations[[i]] %>%
                             dplyr::select(SNP, Binary) %>%
                             dplyr::distinct(SNP, .keep_all = TRUE)
                           df <- BM %>%
                             dplyr::left_join(df[,c("SNP", "Binary")], by = c("SNP"))
                           df[is.na(df)] <- 0

                           df
                         }

  # Stop cluster
  parallel::stopCluster(cl)

  names(master_list) <- names(list_of_annotations)
  return(master_list)

}

#' Function to create 22 .annot.gz files for each annotation.
#'
#' @param list_of_annotations A named list with all annotations that have been
#'   overlapped with the baseline model. It is important that annotation names
#'   do not include an underscore ('_'). If it does, when using the
#'   Assimilate_H2_results() function anything prior to the last word in the
#'   annotation name will be removed, due to the way the Assimilate_H2_results()
#'   function works. Consider using a different separate, such as a full stop or
#'   colon.
#' @param annot_basedir Base directory in which files should be created. Do not
#'   include "/" at the end of the directory name.
#'
#' @return Create directories representing each of the named annotations, with
#'   .annot.gz files for each of 22 chromosomes.
#' @export
#'

create_annot_file_and_export <- function(list_of_annotations, annot_basedir){

  annot_basedir <- as.character(annot_basedir)

  for(i in 1:length(list_of_annotations)){

    annotation.name <- names(list_of_annotations)[i]
    annotation <- list_of_annotations[[i]]

    print(stringr::str_c("Annotation:", annotation.name))

    # Make directory for output
    directory <- paste0(annot_basedir, "/", annotation.name)

    if(dir.exists(directory)){
      cat("Directory exists")
    } else {
      cat("Directory does not exist -- creating")
      system2(command = "mkdir", args = directory)
    }

    # Split annotation into chromosomes
    annotation.split <- annotation %>%
      split(annotation$CHR)

    # For each chromosome, export chromosome as .annot file and
    for(j in 1:length(annotation.split)){

      CHR.name <- names(annotation.split)[j]
      CHR.list <- annotation.split[[j]]

      print(stringr::str_c("CHR:", CHR.name))

      file.name <- paste0(annot_basedir, "/", annotation.name, "/", annotation.name, ".", CHR.name, ".annot")

      readr::write_delim(as.data.frame(CHR.list), path= file.name, delim = "\t")

      R.utils::gzip(file.name, overwrite = TRUE, remove = TRUE)

    }

  }

}

#---Second level code--------------------------------------------------------------------------------------------------------------------####


#' Remove duplicates
#'
#' Function for removing duplicates within defined column.
#'
#' @param Dataframe Dataframe with a column containing duplicated values.
#' @param Column_w_duplicates Column name in quotation marks.
#'
#' @return Original dataframe without duplicated rows.
#' @export
#'

remove_duplicates <- function(Dataframe, Column_w_duplicates){
  Dataframe[!(duplicated(Dataframe[[Column_w_duplicates]])),]
}

#' Keep duplicates
#'
#' Function for extracting duplicated values in a defined column.
#'
#' @param Dataframe Dataframe with a column containing duplicated values.
#' @param Column_w_duplicates Column name in quotation marks.
#'
#' @return Dataframe with only duplicated rows.
#' @export
#'

keep_duplicates <- function(Dataframe, Column_w_duplicates){
  Dataframe[duplicated(Dataframe[[Column_w_duplicates]], fromLast = TRUE) | duplicated(Dataframe[[Column_w_duplicates]], fromLast = FALSE),]
}

#' Convert from dataframe to GRanges
#'
#' @param df Dataframe with columns as sequence name, start, and end. Extra
#'   columns kept as metadata.
#' @param seqname_col Name of column containing chromsome name. Must be in
#'   quotation marks.
#' @param start_col Name of column containing start position. Must be in
#'   quotation marks.
#' @param end_col Name of column containing end position. Must be in quotation
#'   marks.
#'
#' @return GRanges object with original data from dataframe.
#' @export
#'

df2GRanges <- function(df, seqname_col, start_col, end_col ) {

  gr <- GenomicRanges:: makeGRangesFromDataFrame(df,
                                                 keep.extra.columns = TRUE,
                                                 ignore.strand = TRUE,
                                                 seqinfo = NULL,
                                                 seqnames.field = seqname_col,
                                                 start.field = start_col,
                                                 end.field = end_col)
  return(gr)
}

#' Add window to genes
#'
#' Function for adding a window to genes.
#'
#' @param dataset Dataframe with columns with gene start_position and
#'   end_position. Check columns names are start_position and end_position for
#'   human queries.
#' @param windowsize Size of window to be added to annotation. Enter as integer.
#' @param mouse If original genes are from mouse then select "true". Necessary
#'   as output column names from biomart queries differs between species.
#'
#' @return Dataframe with x bp added to gene start/end.
#' @export
#'

AddBPWindow <- function(dataset, windowsize, mouse = NULL){

  if (!is.null(mouse)) {

    N <- transform(dataset, Gene.Start.MinusKB = (dataset$Gene.start..bp. - windowsize),
                   Gene.end.PlusKB = dataset$Gene.end..bp. + windowsize)

  } else {

    N <- transform(dataset, Gene.Start.MinusKB = (dataset$start_position - windowsize),
                   Gene.end.PlusKB = dataset$end_position + windowsize)

  }

  return(N)

}

#' Add gene positional details
#'
#' Function for extracting human gene positional information and removing X, Y
#' and MT chromosomal genes.
#'
#' @param dataframe Dataframe with gene names.
#' @param columnToFilter Name of column in dataframe, which contains gene names.
#' @param mart Specify genome build. Baseline from LDSC based on hg19/37.
#' @param attributes Vector of attributes to extract from BioMart.
#' @param filter Vector of filter to be used for BioMart query.
#'
#' @return Original dataframe together with gene positional information.
#' @export
#'

AddGenePosDetails_RemoveXandYandMT <- function(dataframe, columnToFilter, mart = 38, attributes, filter){

  if(mart != 38 && mart != 37) stop("Mart must be 38 or 37...")

  if(mart == 38){

    ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }else if(mart == 37){

    ensembl_mart <-
      biomaRt::useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }

  # Query genes as a vector
  genes <- dataframe %>% .[[columnToFilter]] %>% unique()
  print(stringr::str_c("Number of unique genes to search: ", length(genes)))

  # BioMart search
  biomart_query <- biomaRt::getBM(attributes = attributes, filters = filter, values = genes , mart = ensembl_mart)
  print(stringr::str_c("Number of matches found:", nrow(biomart_query)))

  # Create new data frame with positional information + remainder of the original dataframe
  # First requires creating join vector for the by argument in inner_join
  join_vector <- filter
  names(join_vector) <- columnToFilter
  Merged <- dplyr::inner_join(dataframe, biomart_query, by = join_vector)

  # Remove any rows with "CHR" or "HG" in chomosome_name
  CHR_grep <- Merged[grep("CHR", Merged$chromosome_name),]
  Merged <- dplyr::anti_join(Merged, CHR_grep)
  HG_grep <- Merged[grep("HG", Merged$chromosome_name),]
  Merged <- dplyr::anti_join(Merged, HG_grep)

  # Remove X, Y and MT
  Merged <- Merged[!(Merged$chromosome_name=="X" | Merged$chromosome_name=="Y" |  Merged$chromosome_name=="MT" ),]

  return(Merged)

}
