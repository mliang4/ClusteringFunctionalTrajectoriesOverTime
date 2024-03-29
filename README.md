# ClusteringFunctionalTrajectoriesOverTime

This folder contains the code of the Rcpp package hDFPmcmc, DFPmcmc, and other R files associated with experiments in the manuscript:

A Bayesian Nonparametric Approach for Clustering Functional Trajectories over Time

by Mingrui Liang, Matthew D. Koslovsky, Emily T. Hebert, Darla E. Kendzor and Marina Vannucci

The main functions operate in C++ via the R package Rcpp. 
This package relies on various R packages that need to be installed in the R environment before running. 
To install, use the install.packages("") command for the following packages:  

'Rcpp'   
'RcppArmadillo'  
'mcclust'   
'mvtnorm'  
'ggplot2'
and also 
'devtools'

To install the package, download the ‘hDFPmcmc’ folder and set the working directory in Rstudio to its path, then run  
library( devtools )  
devtools::build(vignettes = F)  
devtools::install()  
And do the same for folder 'DFPmcmc'.

In the 'sample code' folder, we provide a code demo (Rmd file and html/PDF output), as well as R files to run simulations in the manuscript. 
