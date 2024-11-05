# GB4GT

This repository contains R scripts for the manuscript "Gradient boosting for group testing data" by Erica M. Porter, Christopher S. McMahan, Joshua M. Tebbs, and Christopher R. Bilder.  

1. <b>GT_functions.R</b> contains all of the necessary functions to simulate a group-testing data example and implement the proposed gradient boosting method in our manuscript.

2. <b>GB_example.R</b> creates a simulated group testing data set according to Dorfman testing, with one individual-level predictor variable.  Performs cross-validation and gradient boosting to estimate the probability of infection from the one predictor variable using three types of weak learners:
	(i) Regression trees
	(ii) Kernel smoothing
	(iii) Splines
