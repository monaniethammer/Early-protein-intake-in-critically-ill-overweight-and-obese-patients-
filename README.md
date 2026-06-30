# Project - Early Protein Intake

## Outcomes of Overweight, Critically Ill Patients: The Impact of Early Protein Intake

## Overview

This repository contains the code used for the analyses presented in the manuscript:

**"Prognostic relevance of protein intake during the early acute phase in critically ill overweight or obese patients - analysis of a large international database"**

This project investigates the association between early protein intake and time-to-event outcomes in critically ill patients with overweight or obesity. The analysis is based on a large international intensive care database. Competing risks models were used to assess the association of early protein intake with in-ICU death and ICU discharge alive.

The aim of the project was to examine whether protein intake during the early acute phase of critical illness is associated with clinically relevant ICU outcomes in patients with overweight or obesity.

## Contributors

This project was conducted as a supervised statistical practice project.

The manuscript was led by Dr. med. Michael Neuberger and Prof. Dr. med. Wolfgang H. Hartl.
Statistical supervision was carried out by Mona Niethammer.

The practical implementation of the statistical analysis, including code development, validation, and visualization, was carried out jointly by:

- [Alexander Göstl](https://github.com/agoestl)
- [Daniel Müller](https://github.com/daniel-mueller17)
- [Julian Ehrenberg](https://github.com/ju-eberg)
- [Regina Pichler](https://github.com/reginapichler)
- [Rudolf Schwaab](https://github.com/rschwaab)

These five student contributors contributed equally to the practical implementation of the analysis.

## Project Structure

The repository is divided into several main folders:

* [`Code/`](./Code/) - Source code, scripts, and functions developed for the project

  * [`Analysis/`](./Code/Analyse/)
  * [`Descriptive/`](./Code/Deskriptiv/)
  * [`Preprocessing/`](./Code/Preprocessing/)
  * [`Subgroup_Analysis/`](./Code/Subgruppenanalyse/)
* [`Figures/`](./Grafiken/) - Figures, diagrams, and visualizations
  * [`Analysis/`](./Grafiken/Analyse/)
  * [`Descriptive/`](./Grafiken/Deskriptiv/)
  * [`Background/`](./Grafiken/Hintergrund/)
  * [`Subgroup_Analysis/`](./Grafiken/Subgruppenanalyse/)

## Computational Environment

The analysis was run using R version 4.4.2 on Windows 11. The full session information is provided below.

```r
sessionInfo()

R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape2_1.4.4    survival_3.7-0    broom_1.0.7       gridExtra_2.3     data.table_1.16.4 lubridate_1.9.4  
 [7] haven_2.5.4       checkmate_2.3.2   tm_0.7-16         NLP_0.3-2         viridis_0.6.5     viridisLite_0.4.2
[13] dplyr_1.1.4       tidyr_1.3.1       purrr_1.0.4       here_1.0.1        ggplot2_3.5.1     mgcv_1.9-1       
[19] nlme_3.1-166      pammtools_0.7.3  

loaded via a namespace (and not attached):
 [1] gtable_0.3.6        xfun_0.51           lattice_0.22-6      numDeriv_2016.8-1.1 vctrs_0.6.5         tools_4.4.2        
 [7] generics_0.1.3      parallel_4.4.2      tibble_3.2.1        pkgconfig_2.0.3     Matrix_1.7-1        lifecycle_1.0.5    
[13] scam_1.2-18         stringr_1.5.1       compiler_4.4.2      munsell_0.5.1       codetools_0.2-20    htmltools_0.5.8.1  
[19] yaml_2.3.10         lazyeval_0.2.2      prodlim_2024.06.25  Formula_1.2-5       pillar_1.10.2       iterators_1.0.14   
[25] foreach_1.5.2       parallelly_1.45.0   lava_1.8.1          timereg_2.0.6       tidyselect_1.2.1    digest_0.6.37      
[31] stringi_1.8.7       mvtnorm_1.3-3       slam_0.1-55         future_1.67.0       listenv_0.9.1       forcats_1.0.0      
[37] splines_4.4.2       rprojroot_2.1.1     fastmap_1.2.0       colorspace_2.1-1    cli_3.6.3           magrittr_2.0.3     
[43] future.apply_1.20.0 withr_3.0.2         scales_1.3.0        backports_1.5.0     timechange_0.3.0    rmarkdown_2.29     
[49] globals_0.18.0      hms_1.1.3           evaluate_1.0.4      knitr_1.50          pec_2023.04.12      rlang_1.1.4        
[55] Rcpp_1.0.14         glue_1.8.0          xml2_1.3.7          rstudioapi_0.17.1   plyr_1.8.9          R6_2.6.1
```

