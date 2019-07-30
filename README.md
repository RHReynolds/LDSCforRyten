# LDSCforRyten
LDSC scripts for use on MR server. This will only run on the MR server, as some of the arguments are hard-coded directory paths. If user wishes to use the package locally, change the appropriate arguments listed in `get_LDSC_fixed_args` and `Create_GWAS_df` functions found in [LDSC_Pipeline_Functions.R](R/LDSC_Pipeline_Functions.R).

Scripts make use of the command line tool `ldsc`. For more information on S-LDSC, please refer to: https://github.com/bulik/ldsc/wiki.

## Installing the package
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/LDSCforRyten", auth_token = "")
```

As this is a private repository, you will have to generate a [personal access token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line), and insert this into the ```auth_token``` argument. **Remember to save this token, as you may need it to access other private repositories.**

## Scripts

Script | Description | Author(s)
------ | ----------- | ---------
[]() | Script for running LDSC with an annotation at it's subcategories. | RHR