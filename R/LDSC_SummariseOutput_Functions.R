#---First level code---------------------------------------------------------------------------------------------------------------------####
#' Function assembling herirtability results from LD score regression
#'
#' @param path_to_results List of file paths where .results output from LDSC
#'   stored. For the function to work (i.e. provide the full name of the
#'   annotation and the name of the GWAS), it is important that neither the GWAS
#'   name or the annotation name includes '_'. If it does, anything after the
#'   first word in the GWAS name will be removed, or anything before the last
#'   word in the annotation will be removed. Consider using a different
#'   separator in annotation/GWAS names e.g. a full stop, or colon.
#'
#' @return Dataframe of results together with annotation name and GWAS.
#' @export
#'

Assimilate_H2_results <- function(path_to_results){

  assimilated.results <- data.frame()

  for(i in seq_along(path_to_results)){

    path.name <- path_to_results[i]

    GWAS <-
      path.name %>%
      str_replace("/.*/", "") %>%
      str_replace("\\.results", "") %>%
      str_replace("_.*", "")

    annot_name <-
      path.name %>%
      str_replace("/.*/", "") %>%
      str_replace("\\.results", "") %>%
      str_replace(".*_", "")

    results <- read_delim(file = path.name, delim = "\t")

    # Select relevant row referring to category of interest
    results <- results[1,] %>%
      dplyr::select(-Category) %>%
      mutate(annot_name = annot_name, GWAS = GWAS)

    # If statement to ensure rows only bound after first iteration of loop
    if(i == 1){
      assimilated.results <- results
    } else{
      assimilated.results <- bind_rows(assimilated.results, results)
    }

  }

  return(assimilated.results)

}

#' Assimilate LDSC output from multiple annotations run together.
#'
#' Function assembling heritability results from LD score regression run with
#' multiple annotations (i.e. conditional upon one another). In other words, in
#' the --ref-ld-chr of 'ldsc' multiple annotations (in addition to baseline) are
#' used.
#'
#' @param path_to_results List of file paths where .results output from LDSC
#'   stored. For the function to work (i.e. provide the full name of the
#'   annotation and the name of the GWAS), it is important that neither the GWAS
#'   name or the annotation name includes '_'. If it does, anything after the
#'   first word in the GWAS name will be removed, or anything before the last
#'   word in the annotation will be removed. Consider using a different
#'   separator in annotation/GWAS names e.g. a full stop, or colon.
#'
#' @return Dataframe of results together with associated annotation names and
#'   GWAS.
#' @export
#'
Assimilate_H2_results_multipleannot <- function(path_to_results){

  assimilated.results <- data.frame()

  for(i in seq_along(path_to_results)){

    path.name <- path_to_results[i]

    GWAS <-
      path.name %>%
      str_replace("/.*/", "") %>%
      str_replace("\\.results", "") %>%
      str_replace("_.*", "")

    annot_name <-
      path.name %>%
      str_replace("/.*/", "") %>%
      str_replace("\\.results", "") %>%
      str_replace(".*_", "")

    annot_name_split <-
      annot_name %>%
      str_split(":") %>%
      .[[1]]

    Annot1 <- annot_name_split[1]

    Annot2 <- annot_name_split[2]

    results <- read_delim(file = path.name, delim = "\t")

    # Select relevant row referring to category of interest
    results.annot1 <- results[1,] %>%
      dplyr::select(-Category) %>%
      mutate(Model = annot_name, Annotation.Subset = Annot1, GWAS = GWAS)

    results.annot2 <- results[2,] %>%
      dplyr::select(-Category) %>%
      mutate(Model = annot_name, Annotation.Subset = Annot2, GWAS = GWAS)

    # If statement to ensure rows only bound after first iteration of loop
    if(i == 1){
      assimilated.results <- results.annot1
      assimilated.results <- bind_rows(assimilated.results, results.annot2)
    } else{
      assimilated.results <- bind_rows(assimilated.results, results.annot1)
      assimilated.results <- bind_rows(assimilated.results, results.annot2)
    }

  }

  return(assimilated.results)

}

#' Calculate enrichment SE, logP and z-score p-value.
#'
#' Function calculating upper and lower bounds of enrichment SE, log
#' P(enrichment) and z-score p-value.
#'
#' @param df Dataframe with .results output of LDSC, as created using the
#'   Assimilate_H2_results() function.
#'
#' @return Dataframe with upper and lower bounds of enrichment SE, log
#'   P(enrichment) and z-score p-value
#' @export
#'

Calculate_enrichment_SE_and_logP <- function(df, one_sided){

  # Calculations
  Enrichment.Lower.SE <- df$Enrichment - df$Enrichment_std_error
  Enrichment.Upper.SE <- df$Enrichment + df$Enrichment_std_error
  Log_P <- -log10(df$Enrichment_p)
  Coefficient.Lower.SE <- df$Coefficient - df$Coefficient_std_error
  Coefficient.Upper.SE <- df$Coefficient + df$Coefficient_std_error
  # One-sided as only interested in testing whether significantly positive
  Z_score_P <- convert_z_score(df$`Coefficient_z-score`, one_sided=one_sided)
  Z_score_logP <- -log10(Z_score_P)

  # Bind to df
  df <- cbind(df, Enrichment.Lower.SE, Enrichment.Upper.SE, Log_P,
              Coefficient.Lower.SE, Coefficient.Upper.SE, Z_score_P, Z_score_logP)

  return(df)

}

#' Plot S-LDSC H2 enrichment & coefficient per annotation
#'
#' Function plotting enrichment, coefficient and coefficient p-value of
#' annotation, with GWASs facetted.
#'
#' @param h2_results Dataframe with overall heritability.
#' @param x_axis character. Variable to be plotted on the x-axis.
#' @param xlab character. X-axis name in quotation marks.
#' @param fill_variable character. Variable to be used to determine the fill
#'   colour of bars/points. Default is NULL.
#' @param colour character vector. If fill_variable not NULL then supply colours
#'   to be used.
#' @param pvalue_cutoff numeric. P-value cutoff after multiple test correction
#'   for number of annotations and GWASs.
#'
#' @return Plot
#' @export
#'

Plot_H2_enrichment_coefficient <- function(h2_results, x_axis, xlab, fill_variable = NULL, colour, pvalue_cutoff){

  fill_variable <- vars(!!!fill_variable)

  h2_results <- h2_results %>%
    dplyr::mutate(GWAS = str_replace(GWAS, "PD2018.ex23andMe", "PD2018 ex23andMe"))

  Enrichment <- ggplot(data = h2_results, aes(x = h2_results %>% .[[x_axis]],
                                              y = Enrichment)) +
    geom_col(aes(fill = if(!is.null(fill_variable)) h2_results %>% .[[fill_variable]] else "arb_colour_cat"),
             colour = "black") +
    geom_hline(aes(yintercept = 1),
               linetype = "22",
               size = 0.60) + #No enrichment
    geom_errorbar(aes(ymin = Enrichment.Lower.SE,
                      ymax = Enrichment.Upper.SE),
                  width = 0.2) +
    facet_grid(facets = GWAS ~ .) +
    labs(x = xlab, y = expression("Enrichment"), title = "") +
    theme_bw() +
    theme(axis.text.x = if(!is.null(fill_variable)) element_blank() else element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 8),
          legend.position = "bottom") +
    scale_fill_manual(guide = guide_legend(
      title.theme = element_text(size = 8, face = "bold", colour = "black", angle = 0),
      label.theme = element_text(size = 7, angle = 0),
      nrow = 1,
      byrow = TRUE),
      name = element_blank(),
      values = if(!is.null(fill_variable)) colour else "#8C8C8C")

  Coefficient_plot <- ggplot(data = h2_results,
                             aes(x = h2_results %>% .[[x_axis]],
                                 y = Coefficient)) +
    geom_col(aes(fill = if(!is.null(fill_variable)) h2_results %>% .[[fill_variable]] else "arb_colour_cat"),
             colour = "black") +
    # geom_hline(aes(yintercept = 1), linetype = "22", size = 0.60) + #No enrichment
    geom_errorbar(aes(ymin = Coefficient.Lower.SE,
                      ymax = Coefficient.Upper.SE),
                  width = 0.2) +
    facet_grid(facets = GWAS ~ .) +
    labs(x = xlab, y = expression("Regression coefficient,"~tau*"c"), title = "") +
    theme_bw() +
    theme(axis.text.x = if(!is.null(fill_variable)) element_blank() else element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 8),
          legend.position = "bottom") +
    scale_fill_manual(guide = guide_legend(
      title.theme = element_text(size = 8, face = "bold", colour = "black", angle = 0),
      label.theme = element_text(size = 7, angle = 0),
      nrow = 1,
      byrow = TRUE),
      name = element_blank(),
      values = if(!is.null(fill_variable)) colour else "#8C8C8C")

  Coefficient_pvalue_plot <- ggplot(data = h2_results, aes(x = h2_results %>% .[[x_axis]] ,
                                                           y = Z_score_logP)) +
    geom_point(aes(colour = if(!is.null(fill_variable)) h2_results %>% .[[fill_variable]] else "arb_colour_cat"),
               shape = 19,
               size = 1.5) +
    geom_point(data = subset(h2_results, Z_score_logP >= pvalue_cutoff),
               aes(x = subset(h2_results, Z_score_logP >= pvalue_cutoff) %>% .[[x_axis]]),
               shape = 1,
               size = 1.5,
               stroke = 0.75) +
    geom_hline(aes(yintercept = pvalue_cutoff), linetype = "77", size = 0.6) + # Bonferroni correct for number of tests
    facet_grid(facets = GWAS ~ .) +
    labs(x = xlab, y = expression("Coefficient P-value (-log"[1][0]~"p)"), title = "") +
    theme_bw() +
    theme(axis.text.x = if(!is.null(fill_variable)) element_blank() else element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 8),
          legend.position = "bottom") +
    scale_colour_manual(guide = FALSE,
                        name = element_blank(),
                        values = if(!is.null(fill_variable)) colour else "#8C8C8C")

  ggarrange(Enrichment, Coefficient_plot, Coefficient_pvalue_plot,
            widths = c(7.5,7.5,10), ncol = 3, nrow = 1,
            common.legend = if(!is.null(fill_variable)) TRUE else FALSE,
            legend = if(!is.null(fill_variable)) "bottom" else "none"
  )

}

#' Plot S-LDSC coefficient p-value per annotation
#'
#' Function plotting coefficient p-value of annotation, with GWASs facetted.
#'
#' @param h2_results Dataframe with overall heritability.
#' @param x_axis character. Variable to be plotted on the x-axis.
#' @param xlab character. X-axis name in quotation marks.
#' @param fill_variable character. Variable to be used to determine the fill
#'   colour of bars/points. Default is NULL.
#' @param colour character vector. If fill_variable not NULL then supply colours
#'   to be used.
#' @param pvalue_cutoff numeric. P-value cutoff after multiple test correction
#'   for number of annotations and GWASs.
#' @param show_legend logical. Default = FALSE.
#'
#'
#' @return Ggplot of co-efficient p-values.
#' @export
#'

Plot_H2_coefficient <- function(h2_results, x_axis, xlab, fill_variable = NULL, colour, pvalue_cutoff, show_legend = FALSE){

    Coefficient_pvalue_plot <- ggplot(data = h2_results, aes(x = h2_results %>% .[[x_axis]] ,
                                                             y = Z_score_logP)) +
      geom_point(aes(colour = if(!is.null(fill_variable)) h2_results %>% .[[fill_variable]] else "arb_colour_cat"),
                 shape = 19,
                 size = 1.5) +
      geom_point(data = subset(h2_results, Z_score_logP >= pvalue_cutoff),
                 aes(x = subset(h2_results, Z_score_logP >= pvalue_cutoff) %>% .[[x_axis]]),
                 shape = 1,
                 size = 1.5,
                 stroke = 0.75) +
      geom_hline(aes(yintercept = pvalue_cutoff), linetype = "77", size = 0.6) + # Bonferroni correct for number of tests
      facet_grid(facets = GWAS ~ .) +
      labs(x = xlab, y = expression("Coefficient P-value (-log"[1][0]~"p)"), title = "") +
      theme_bw() +
      theme(axis.title = element_text(size = 8, face = "bold"),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            strip.text = element_text(size = 8),
            legend.position = "bottom",
            legend.key.size = unit(1,"line")) +
      scale_colour_manual(guide = if(show_legend) guide_legend(label.theme = element_text(size = 6, angle = 0), ncol = 2, byrow = TRUE) else FALSE,
                          name = element_blank(),
                          values = if(!is.null(fill_variable)) colour else "#8C8C8C")

      return(Coefficient_pvalue_plot)

}

#---Second level code--------------------------------------------------------------------------------------------------------------------####
#' Convert z-score to p-value
#'
#' Function converting z-score to p-value.
#'
#' @param z Z-score
#' @param one_sided Specify whether one-sided or two-sided p-value required. If
#'   NULL, two-sided test will be used. If one-sided preferred then specify
#'   direction of one-sided test. "-" for negative, and "+" for positive.
#'
#' @return P-value.
#' @export
#'

convert_z_score<-function(z, one_sided=NULL) {
  if(is.null(one_sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one_sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);
  }
  return(pval);
}
