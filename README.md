# GB4GT

This repository contains R scripts for the manuscript "Gradient boosting for group testing data" by Erica M. Porter, Christopher S. McMahan, Joshua M. Tebbs, and Christopher R. Bilder.  

1. <b>GB_functions.R</b> contains the necessary functions to simulate a group-testing data set and implement the proposed gradient boosting method in our manuscript.

2. <b>CV_functions.R</b> contains functions that can be used to perform cross-validation for 3 types of weak learners to be used for gradient boosting: (i) regression trees, (ii) kernel smoothing, and (iii) splines.  These functions can easily be adapted to select other tuning parameters of interest or to accommodate other types of weak learners.

3. <b>GB_example.R</b> creates a simulated group testing data set according to Dorfman testing, with one individual-level predictor variable.  Given values for tuning parameters of the weak learners, gradient boosting is used to estimate the probability of infection from the one predictor variable using three types of weak learners studied in the manuscript.
