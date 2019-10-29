#' Liftover GWAS from hg19 to hg38
#'
#' @param GWAS GWAS with a minimum of 2 columns CHR and BP. If
#'   additional columns included these will be preserved.
#' @param path_to_chain Path to UCSC chain file for transformation from hg19 to
#'   hg38 coordinates.
#'
#' @return GWAS with hg38 coordinates
#' @export

liftover_hg19_to_hg38 <- function(GWAS, path_to_chain){

  library(rtracklayer)
  library(tidyverse)

  # Import chain file
  hg19_to_hg38 <- rtracklayer::import.chain(path_to_chain)

  # If GWAS CHR column does not have "chr" in name, add to allow liftover
  if(!str_detect(GWAS$CHR[1], "chr")){

    GWAS <- GWAS %>%
      dplyr::mutate(CHR = str_c("chr", CHR))

  }

  # Convert GWAS to GRanges object
  GWAS_GR <- makeGRangesFromDataFrame(GWAS,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqinfo = NULL,
                                      seqnames.field = "CHR",
                                      start.field = "BP",
                                      end.field = "BP",
                                      starts.in.df.are.0based = FALSE)

  GWAS_hg38 <- liftOver(GWAS_GR, hg19_to_hg38) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::rename(CHR = seqnames,
                  BP = start) %>%
    dplyr::select(-end, -width, -strand) %>%
    dplyr::mutate(CHR = str_replace(CHR, "chr", ""))

  return(GWAS_hg38)

}

#' Find overlapping RS ids in GWAS using genomic co-ordinates.
#'
#' Function to find the RS ids in a GWAS where SNPs have been provided as
#' genomic coordinates i.e. some combination of SNP and BP.
#'
#' @param GWAS GWAS with a minimum of 2 columns labelled CHR and BP. If
#'   additional columns included these will be preserved.
#' @param dbSNPref BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build)
#'
#' @return GWAS with RS ids.
#' @export

add_RS_to_GWAS <- function(GWAS, dbSNPref){

  library(BSgenome)
  library(GenomicRanges)

  # If GWAS CHR column has "chr" in name, remove
  if(str_detect(GWAS$CHR[1], "chr")){

    GWAS <- GWAS %>%
      dplyr::mutate(CHR = str_replace(CHR, "chr", ""))

  }

  # Convert GWAS to GRanges object
  GWAS_GR <- makeGRangesFromDataFrame(GWAS,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqinfo = NULL,
                                      seqnames.field = "CHR",
                                      start.field = "BP",
                                      end.field = "BP",
                                      starts.in.df.are.0based = FALSE)

  # Genomic position object as dataframe with SNP locations converted to RS id.
  GWAS_GR <-
    snpsByOverlaps(dbSNPref, GWAS_GR, minoverlap = 1L) %>%
    # Note that the default value for minoverlap is 0 which means that, by default, in addition to the SNPs that are
    # located within the genomic regions specified thru the ranges argument, snpsByOverlaps also returns SNPs that are
    # adjacent to these regions. Use minoverlap=1L to omit these SNPs.
    as.data.frame()

  # Inner join of GWAS_RS with GWAS
  Combined <-
    GWAS_GR %>%
    dplyr::rename(SNP = RefSNP_id, CHR = seqnames, BP = pos) %>%
    # Remove all positions with more than one rs id attached
    dplyr::mutate(CHR_BP = str_c(CHR, ":", BP)) %>%
    dplyr::group_by(CHR_BP) %>%
    dplyr::filter(!any(row_number() > 1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-CHR_BP) %>%
    # Inner join with GWAS
    dplyr::inner_join(GWAS, by = c("CHR", "BP")) %>%
    dplyr::select(-strand, -alleles_as_ambig)

  return(Combined)

}
