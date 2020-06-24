# Long_read_meth_assessment
Scripts from my internship at QIMR assessing the ONT PromethION platform for methylation calling

This repo contains scripts from my internship with QIMR Berghofer's medical genomics group investigating the Oxford Nanopore PromethION platform for methylation calling. 

A description of the scripts in this repo are as follows:

 * minfi.R is the R script used to call methylation and generate beta values in the EPIC array data using the [minfi](https://www.bioconductor.org/packages/release/bioc/html/minfi.html) R library
 * shared_loci.R was used to identify common loci between the PromethION platform and EPIC array data
 * Plots_for_report.R was used to generate the plots and do the statistical analysis for my report.
 * The f5c folder contains the job submission scripts used by PBS to generate the methylation data for the PromethION platform. Methylation calling was done using [f5c](https://github.com/hasindu2008/f5c)
 
 
 
 The companion project investigating genomic dark regions in short- and long- read data undertaken concurrently is located [here](https://github.com/brookshenry3/simpler_DRF)
