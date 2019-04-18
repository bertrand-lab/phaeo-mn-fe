## Phaeocystis Mn/Fe Proteomics Analysis

This is the github page for a project submitted recently for review. These scripts reproduce the analysis presented in the paper:

`irs_normalization_plotting_de.R` does the differential expression analyses, produces plots, and produces differential expression spreadsheets as an output.
`phys_anova.R` does the plots, ANOVAs, and Tukeys Post Hoc testing for the physiological data.
`antarctica_map.R` makes the map shown in Figure 2.
`peptide_transition_plots.R` is for producing the figures of ion chromatograms.
'irs_normalization_plotting_de-imputed-values.R' does the same as the differential expression analysis, except it imputes data where proteins are missing from only one channel of 20, and the imputed value is equal to the 
0.5*minimum value of that channel.

Here is the R session info that details required packages:

```R
R version 3.4.3 (2017-11-30)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                    LC_TIME=English_Canada.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.0.0       reshape2_1.4.3     cowplot_0.9.2      multcompView_0.1-7 broom_0.4.4        forcats_0.3.0      stringr_1.3.1      dplyr_0.7.4        purrr_0.2.4       
[10] readr_1.1.1        tidyr_0.8.0        tibble_1.4.2       tidyverse_1.2.1    testthat_2.0.0     ggfortify_0.4.4    gridExtra_2.3      readxl_1.1.0       xlsx_0.5.7        
[19] xlsxjars_0.6.1     rJava_0.9-9        edgeR_3.20.9       limma_3.34.9       ggplot2_2.2.1.9000

loaded via a namespace (and not attached):
 [1] locfit_1.5-9.1   haven_1.1.1      lattice_0.20-35  colorspace_1.3-2 rlang_0.2.0      pillar_1.2.2     foreign_0.8-69   glue_1.2.0       withr_2.1.2      modelr_0.1.2    
[11] bindrcpp_0.2.2   bindr_0.1.1      plyr_1.8.4       munsell_0.5.0    gtable_0.2.0     cellranger_1.1.0 rvest_0.3.2      psych_1.8.4      parallel_3.4.3   Rcpp_0.12.16    
[21] jsonlite_1.5     mnormt_1.5-5     hms_0.4.2        stringi_1.2.2    grid_3.4.3       cli_1.0.0        tools_3.4.3      magrittr_1.5     lazyeval_0.2.1   crayon_1.3.4    
[31] pkgconfig_2.0.1  xml2_1.2.0       lubridate_1.7.4  rstudioapi_0.7   assertthat_0.2.0 httr_1.3.1       R6_2.2.2         nlme_3.1-131     compiler_3.4.3  
```

In order to run these scripts, you need to have the following directories (~ denotes main directory):

`~/data/` which contains the two supplementary files (`SupplementalFileS1.xlsx` and `SupplementalFileS2.xlsx`).

A directory for figures is required `~/figures/`. A subdirectory in data is required for the intermediate data products: `~/data/intermediate-data/`, and a subdirectory is required for the MS data `~/data/ms_data`.
