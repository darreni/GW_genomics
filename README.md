[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

# GW_genomics

##Greenish Warbler genomics

This repository contains the R code used in this manuscript:

Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in a ring species. In revision, *Molecular Ecology*.

An earlier version of the manuscript was posted on bioRxiv:

Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. 2016. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in greenish warblers. *bioRxiv*, doi: http://dx.doi.org/10.1101/041467

There are two key R files included in this repository:

**GW_GBS_R_analysis_script_for_Dryad.R**   This is the primary script that will conduct the analyses and produce the figures. See notes within the file for guidance in using the script.

**genomics_R_functions.R**   This contains the custom-designed functions that are called by the primary script above. These functions are written to be general, such that they can be applied to datasets other than just the greenish warblers.

The above two files have also been deposited in a Dryad package that contains data and other useful files needed to reproduce the entire analysis. (The R files only contain scripts for the analysis from the "012NA" genotype file stage; for earlier processing and data files, see the Dryad package)

Both R files were written by Darren Irwin. You are welcome to use these scripts; if you do, please site the paper above, the Dryad package, and/or this GitHub repository.

The purpose of posting to this GitHub repositore as well as the Dryad site is to provide a way to provide updates to these processing scripts. The Dryad site archives the scripts at the time of publication; improvements to the code for future analyses can be posted here.

I welcome questions and comments: irwin@zoology.ubc.ca