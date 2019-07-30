---WORK IN PROGRESS---

# S-LDSC for Ryten Lab
S-LDSC scripts for use on MR server. This will only run on the MR server, as some of the arguments are hard-coded directory paths. If user wishes to use the package locally, change the appropriate arguments listed in `get_LDSC_fixed_args` and `Create_GWAS_df` functions found in [LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R).

Scripts make use of the command line tool `ldsc`. For more information on S-LDSC, please refer to: 
- https://github.com/bulik/ldsc/wiki
- https://www.ncbi.nlm.nih.gov/pubmed/26414678

## Installing the package
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/LDSCforRyten", auth_token = "")
```

As this is currently a private repository, you will have to generate a [personal access token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line), and insert this into the ```auth_token``` argument. **Remember to save this token, as you may need it to access other private repositories.**

## Running S-LDSC
Running S-LDSC can be divided into the following steps.
1. Creating an .annot file. This is a file consisting of CHR, BP, SNP, and CM columns, followed by a column for your annotation, with the value of the annotation for each SNP (0/1 for binary annotations or arbitrary numbers for continuous annotations). It is important that SNPs are provided in the same order as the .bim file used for the computation LD scores. To ensure this is the case, the easiest thing to do is to find those SNPs from an annotation (which may be genes, or genomic regions) that overlap with the baseline model that is used for computation of LD scores. 
2. Computing LD scores for the annotation.
3. Partitioning heritability by annotation, using the baseline model.
4. Collating and summarising the output of S-LDSC.

For an example of this process run from end to end, please see the following markdown: (insert link to html).

## Scripts

Script | Description | Author(s)
------ | ----------- | ---------
[LDSC_Creatingannot_Functions.R](R/LDSC_Creatingannot_Functions.R) | Functions that can be used for creating binary .annot.gz files prior to running [LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R). | RHR
[LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R) | Pipeline for running stratified LDSC with a binary annotation and it's subcategories. This requires that the user has created the appropriate .annot.gz files. Note that this pipeline is divided into two steps: 1) calculating LD scores for an annotation and 2) partitioning heritability in the annotation. If necessary, these two steps can be run separately. Call the script using: `Rscript /path/to/script/LDSC_Pipeline_Functions.R -h`. The `-h` flag will list the required inputs and optional arguments. | RHR
[LDSC_SummariseOutput_Functions.R](R/LDSC_SummariseOutput_Functions.R) | Functions that can be used to summarise the output of S-LDSC once pipeline has been run. | RHR
