READ ME

This code is associated with the paper:
Wauchope, H. S., Jones, J. P. G., Geldmann, J., Simmons, B. I., Amano, T., Blanco, D., Fuller, R.A., Johnston, A., Langendoen, T., Mundkur T., Nagy, S., Sutherland, W.J. (2022) "Protected areas have a mixed impact on waterbirds but management helps" Nature. (Volume and Issue not assigned at time of writing)
It was written by Hannah Wauchope (current affiliation: The University of Exeter). Please do not hesitate to direct any queries to hannah.wauchope@gmail.com

***

This code is all in the language R, and provides the full workflow. It is divided into four sequential sections:

1: CollateCleanRawData.R - This script takes data as received from Wetlands International and Audubon and collates and cleans the data, including harmonising taxonomy (according to Birdlife international), managing cases of multiple counts, removing any suspicious sites, and imputing zeroes
2: ExtractCovariates.R - This script takes data output from Script 1, and extracts site and species specific covariates (see Extended Data Tables 1 & 2 for those covariates)
3: Matching.R - This script takes data output from Script 2, subsets data to the requirements for analysis, identifies sites that fall within protected areas, and then matches these to unprotected sites, according to the 21 analysis scenarios (focal analysis plus full-parameter analyses)
4: AnalysisResults.R - This script takes data output from Script 3, assesses the quality of matches obtained in Script 3, calculates protected area impact on populations under BA, CI and BACI frameworks, and, for BACI, runs cumulative link mixed models to correlate protected area impact with predictors. It also produces all data figures found in the paper.

To fully replicate analysis, data should first be acquired from Wetlands International and Audubon:

	Data from Wetlands International can be acquired via http://iwc.wetlands.org/index.php/requestingdata. Once an account has been created, it's easiest to get in touch with the Wetlands International data managers to request the full dataset, rather than obtaining from individual countries. This involves explaining the purpose for the data usage, and this is then sent out to National Coordinators for approval. Russian National Coordinators declined for Russian data to be used in this paper.
	The dataset used in this publication was received from Wetlands international on the 5th March, 2020.

	Data from Audubon can be acquired via http://netapp.audubon.org/cbcobservation/ . Here, individual species or site data can be downloaded, but in the case of this study it's easier to email Audubon (using the email on the website) to request the full dataset. This again involves explaining the purpose for the study and acquiring permissions. 
	The dataset used in this publication was received from Audubon on 9th August, 2019.

	Due to covariate restrictions, all data in the paper is restricted to 2018 and earlier, so any replicates of analysis that have data up until 2018 will have the data we had access to. 

Covariate data should then be acquired, according to the various sources given in Extended Data Tables 2 and 3. Script 2 (and a small part of Script 4) perform any required data transformations on these.

***

We thank the coordinators, thousands of volunteer counters, and funders of the International Waterbird Census. 
CBC Data is provided by National Audubon Society and through the generous efforts of Bird Studies Canada and countless volunteers across the western hemisphere to whom we are most grateful

***

The full analysis and code were run in the following R environment:

	R version 4.0.3 (2020-10-10)
	Platform: x86_64-apple-darwin17.0 (64-bit)
	Running under: macOS Big Sur 10.16

	Matrix products: default
	LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

	locale:
	en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

	attached base packages:
	grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

	other attached packages:
	gridExtra_2.3      RColorBrewer_1.1-2 cowplot_1.1.0      ggeffects_1.0.1   ggalluvial_0.12.3  glmmTMB_1.1.2.2    ggbeeswarm_0.6.0  ordinal_2019.12-10  lme4_1.1-27.1      nlme_3.1-149       lmtest_0.9-38      zoo_1.8-8  DHARMa_0.4.3       car_3.0-11         carData_3.0-4      fields_11.6    spam_2.6-0         dotCall64_1.0-0    DOS_1.0.0          maptools_1.0-2   scales_1.1.1       foreign_0.8-80     resample_0.4       StatMatch_1.4.0   lpSolve_5.6.15     survey_4.0         survival_3.2-7     Matrix_1.2-18  proxy_0.4-24       usdm_1.1-18        MASS_7.3-53        wdpar_1.0.5    sf_0.9-6           ClusterR_1.2.2     gtools_3.8.2       abind_1.4-5   forcats_0.5.0      purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
	tibble_3.0.4       tidyverse_1.3.0    rlist_0.4.6.1      chron_2.3-56   dplyr_1.0.2        rredlist_0.7.0     stringr_1.4.0      taxize_0.9.99   ggalt_0.6.2        maps_3.3.0         pbapply_1.4-3      ncdf4_1.17    plyr_1.8.6         rgeos_0.5-5        reshape2_1.4.4     pbmcapply_1.5.0  data.table_1.13.4  raster_3.4-5       ggplot2_3.3.5      rgdal_1.5-18    sp_1.4-4          

	loaded via a namespace (and not attached):
	readxl_1.3.1        uuid_0.1-4          backports_1.2.1   TMB_1.7.18          splines_4.0.3       gmp_0.6-1     foreach_1.5.1       fansi_0.4.1         magrittr_2.0.1   openxlsx_4.2.4      modelr_0.1.8        extrafont_0.17   extrafontdb_1.0     colorspace_2.0-0    rvest_0.3.6    mitools_2.4         haven_2.3.1         xfun_0.19   crayon_1.3.4        jsonlite_1.7.2      iterators_1.0.13   ape_5.4-1           glue_1.4.2          gtable_0.3.0    emmeans_1.5.3       proj4_1.0-10        Rttf2pt1_1.3.8    mvtnorm_1.1-1       DBI_1.1.0           Rcpp_1.0.5   xtable_1.8-4        units_0.6-7         bold_1.1.0    httr_1.4.2          ellipsis_0.3.1      pkgconfig_2.0.3   reshape_0.8.8       dbplyr_2.0.0        conditionz_0.1.0   crul_1.0.0          tidyselect_1.1.0    rlang_0.4.10   munsell_0.5.0       cellranger_1.1.0    tools_4.0.3   cli_2.2.0           generics_0.1.0      sjlabelled_1.1.7   broom_0.7.3         fs_1.5.0            zip_2.1.1    ash_1.0-15          xml2_1.3.2          compiler_4.0.3   rstudioapi_0.13     beeswarm_0.2.3      curl_4.3           
	e1071_1.7-4         reprex_0.3.0        stringi_1.5.3    lattice_0.20-41     classInt_0.4-3      nloptr_1.2.2.2   vctrs_0.3.5         pillar_1.4.7        lifecycle_0.2.0   ucminf_1.1-4        estimability_1.3    insight_0.11.1   R6_2.5.0            KernSmooth_2.23-17  rio_0.5.27      vipor_0.4.5         codetools_0.2-16    boot_1.3-25      assertthat_0.2.1    withr_2.3.0         httpcode_0.3.0   hms_0.5.3           coda_0.19-4         class_7.3-17      minqa_1.2.4         numDeriv_2016.8-1.1   lubridate_1.7.9.2     tinytex_0.28  e