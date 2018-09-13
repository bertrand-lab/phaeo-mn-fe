## Phaeocystis Mn/Fe Proteomics Analysis

This is the github page for a project submitted recently for review. These scripts reproduce the analysis presented in the paper:

`irs_normalization_plotting_dge.R` does the differential expression analyses, produces plots, and produces differential expression spreadsheets as an output.
`phys_anova.R` does the plots, ANOVAs, and Tukeys Post Hoc testing for the physiological data.
`antarctica_map.R` makes the map shown in Figure 2.
`peptide_transition_plots.R` is for producing the figures of ion chromatograms.

Required packages (shown as loading them in R):

```R
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(tidyverse)
library(edgeR)
library(xlsx)
library(readxl)
library(gridExtra)
library(ggfortify)
library(testthat)
library(tidyverse)
library(edgeR)
library(xlsx)
library(readxl)
library(ggfortify)
library(testthat)
library(broom)
library(multcompView)
library(gridExtra)
library(cowplot)
library(maps)
library(mapdata)
library(bitops)
library(RCurl)
library(png)
library(RJSONIO)
library(RgoogleMaps)
library(TeachingDemos)
library(dplyr)
library(mapproj)
library(rgdal)
library(maptools)
library(ggplot2)
```

In order to run these scripts, you need to have the following directories (~ denotes main directory):

`~/data/` which contains the files `Dataset 1_1F2F_.xlsx` and `Dataset 2_1F2F_.xlsx`. A directory for figures is required `~/figures/`. A subdirectory in data is required for the intermediate data products: `~/data/intermediate-data/`, and a subdirectory is required for the MS data `~/data/ms_data`.
