#----Function----

#' Create a file containing all SNPs within the baseline model
#'
#' This function will create a .snp file (which can be used in the --print-snps
#' argument of ldsc) from the baseline files supplied. One .snp file will be
#' created per chromosome and each will contain one column detailing the SNPs in
#' the inputted baseline files.
#'
#' @param path_to_BM Path to directory containing baseline .l2.ldscore.gz files.
#' @param output_dir Path to directory where output will be stored.
#'
#' @return A .snp file for each chromsome, containing one column with all the
#'   SNPs that are used in the inputted baseline files.

SNPsinBM <- function(path_to_BM, output_dir){

  BM <- list.files(path = path_to_BM, pattern = ".l2.ldscore.gz", full.names = TRUE)

  for(i in 1:length(BM)){

    chr_number <- BM[i] %>%
      str_replace("/.*/", "") %>%
      str_replace(".l2.ldscore.gz", "") %>%
      str_replace(".*\\.", "")

    chr <- fread(BM[i]) %>%
      dplyr::select(SNP)

    write_delim(chr, path = str_c(output_dir, "/SNPsinBaselinel2ld.", chr_number, ".snp"), delim = "\t", col_names = FALSE)

  }

}

#----Main-----
SNPsinBM(path_to_BM = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v1.1_ldscores_75annot/",
         output_dir = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v1.1_ldscores_75annot/SNPsinBaselinel2ldscore/")

SNPsinBM(path_to_BM = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.0/",
         output_dir = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.0/SNPsinBaselinel2ldscore/")

SNPsinBM(path_to_BM = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.1_ldscores/",
         output_dir = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.1_ldscores/SNPsinBaselinel2ldscore/")

SNPsinBM(path_to_BM = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.2_ldscores/",
         output_dir = "/data/LDScore/Reference_Files/1000G_Phase3_baselineLD_v2.2_ldscores/SNPsinBaselinel2ldscore/")

SNPsinBM(path_to_BM = "/data/LDScore/Reference_Files/GRCh38/baseline_v1.2/",
         output_dir = "/data/LDScore/Reference_Files/GRCh38/baseline_v1.2/SNPsinBaselinel2ldscore/")
