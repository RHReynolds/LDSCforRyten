# S-LDSC for Ryten Lab 
1. [Installing the package](#install)
2. [Running S-LDSC](#running)
3. [Scripts](#scripts)
4. [Reference files](#reference_files)

S-LDSC scripts for use on MR server. **This will only run on the MR server**, as some of the arguments are hard-coded directory paths. If user wishes to use the package locally, clone the package into your own server and change the appropriate arguments listed in `get_LDSC_fixed_args` and `Create_GWAS_df` functions found in [LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R), as well as the `creating_baseline_df` function found in [LDSC_Creatingannot_Functions.R](R/LDSC_Creatingannot_Functions.R).

Scripts make use of the command line tool `ldsc`. For more information on S-LDSC, please refer to: 
- https://github.com/bulik/ldsc/wiki
- https://www.ncbi.nlm.nih.gov/pubmed/26414678

## Installing the package <a name="install"></a>
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/LDSCforRyten")
```

## Running S-LDSC <a name="running"></a>
Running S-LDSC can be divided into the following steps.
1. Creating an .annot file. This is a file consisting of CHR, BP, SNP, and CM columns, followed by a column for your annotation, with the value of the annotation for each SNP (0/1 for binary annotations or arbitrary numbers for continuous annotations). It is important that SNPs are provided in the same order as the .bim file used for the computation LD scores. To ensure this is the case, the easiest thing to do is to find those SNPs from an annotation (which may be genes, or genomic regions) that overlap with the baseline model that is used for computation of LD scores. 
2. Computing LD scores for the annotation.
3. Partitioning heritability by annotation, using the baseline model.
4. Collating and summarising the output of S-LDSC.

For an example of this process run from end to end, please refer to this [tutorial](vignettes/tutorial.Rmd).

## Scripts <a name="scripts"></a>

Script | Description | Author(s)
------ | ----------- | ---------
[GWAS_formatting_functions.R](R/GWAS_formatting_functions.R) | Functions that can be used to format and liftover GWAS from hg19/GRCh37 to GRCh38. | RHR
[LDSC_Creatingannot_Functions.R](R/LDSC_Creatingannot_Functions.R) | Functions that can be used for creating binary .annot.gz files prior to running [LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R). | RHR
[LDSC_Pipeline_Entire.R](pipelines/LDSC_Pipeline_Entire.R) | Pipeline for running stratified LDSC with a binary annotation and it's subcategories. This requires that the user has created the appropriate .annot.gz files. Note that this pipeline is divided into two steps: 1) calculating LD scores for an annotation and 2) partitioning heritability in the annotation. If necessary, these two steps can be run separately. Call the script using: `Rscript /path/to/script/LDSC_Pipeline_Functions.R -h`. The `-h` flag will list the required inputs and optional arguments. | RHR
[LDSC_SummariseOutput_Functions.R](R/LDSC_SummariseOutput_Functions.R) | Functions that can be used to summarise the output of S-LDSC once pipeline has been run. | RHR


## Reference files <a name="reference_files"></a>
- To run S-LDSC requires a number of reference files. On our server these are located in the following directory: `/data/LDScore/Reference_Files/`. They are also available via the Alkes lab [repository](https://data.broadinstitute.org/alkesgroup/LDSCORE/), and descriptions of necessary reference files can be found described in the `ldsc` [wiki page](https://github.com/bulik/ldsc/wiki).

### Baseline models
- An important decision to make when running S-LDSC is the choice of baseline model. As of August 2019, the Price group recommend the following: 
    1. We recommend that for identifying critical tissues/cell-types via P-value of tau, it is best to use the baseline model, specifically baseline v1.2.
    2. We recommend that for estimating heritability enrichment (i.e., %h2/%SNPs) of any annotation, including tissue-specific annotations, it is best to use baselineLD v2.2.
- For a short overview of the various baseline models, please refer to: https://data.broadinstitute.org/alkesgroup/LDSCORE/readme_baseline_versions
- For a detailed overview of the various baselines, please refer to [LDSC_Baseline_Models.xlsx](misc/LDSC_Baseline_Models.xlsx).
- **IMPORTANT:** All baseline models, except baseline v1.2, are aligned to GRCh37. Baseline v1.2, however, is aligned to GRCh38. This is important, as most sumstat.gz files are generated from summary statistic files based on GRCh37 co-ordinates. This is also important when making your .annot files, as the annotations are created by looking for overlaps between your gene list and the SNPs in the baseline version. **Thus, it is super important that gene positions are based on the same genome build as the baseline version being used.**
    - According to the following [thread](https://groups.google.com/forum/#!topic/ldsc_users/_wIQrqK57Nc), running GRCh37-based summary statistics with the GRCh38-based baseline v1.2 model should not make much of a difference to the outputted estimates. There is, however, a big difference between running with v1.1 (GRCh37-based) or v1.2. 
    - Until conversion of summary stats files from GRCh37 to GRCh38 has occurred, it is recommended to run with baseline v1.2 (despite it being GRCh38 based).  
    - If users to encounter a GRCh37-based GWAS, [functions](R/GWAS_formatting_functions.R) are available for GWAS liftover to GRCh38.

### Available GWASs
- Available GWASs on our server can be found in the following directory: `/data/LDScore/GWAS/`
- For details of these GWASs (including samples numbers and references), please refer to [LDSC_GWAS_details.xlsx](misc/LDSC_GWAS_details.xlsx).
- Additional GWASs can be downloaded, but must be prepared. Please refer to the `ldsc` [wiki page](https://github.com/bulik/ldsc/wiki) for details on how to prepare GWAS summary statistics for use with `ldsc`.
- **If new GWASs are downloaded, please update [LDSC_GWAS_details.xlsx](misc/LDSC_GWAS_details.xlsx) with details.**
  
