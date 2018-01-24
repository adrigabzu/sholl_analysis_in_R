# Sholl analysis in R

Initially described in the 50s [[1](#1)], Sholl analysis has been a standard method to quantitatively assess the morphological complexity of neurons. It is performed by counting dendrite intersections in concentric circles emanating from the soma of a cell at a fixed radius interval. Diverse programs can be used for processing the images [[2](#2)] but also different statistical approaches have been implemented to solve hypotheses based on these analyses.

This [implementation in R](sholl_analysis_git.R) of a Sholl profile analysis based on mixed-effect models was adapted from a SAS as described in:

> Wilson, M.D., Sethi, S., Lein, P.J., Keil, K.P., 2017. Valid statistical approaches for analyzing sholl data: Mixed effects versus simple linear models. Journal of Neuroscience Methods 279, 33–43. https://doi.org/10.1016/j.jneumeth.2017.01.003

In this case, mixed-effect models allow accounting for the variability per mice and neurons as a random effect in order to evaluate differences between the conditions to test, described as fixed effects. For more information about linear mixed effects models and their implementation in R:

* https://pagepiccinini.com/r-course/lesson-6-part-1-linear-mixed-effects-models/

* http://rpsychologist.com/r-guide-longitudinal-lme-lmer

### Features implemented:
* Sholl profile plot with error bars (ggplot)
* Overall difference between conditions considering autoregressive covariance
* Statistical test of differences between conditions per radius
* Area under the curve (AUC) calculation per neuron and statistical comparison based also on linear mixed-effect models.
* Boxplots of AUCs per condition (ggplot2)

Cite this code: [![DOI](https://zenodo.org/badge/118108207.svg)](https://zenodo.org/badge/latestdoi/118108207)
____

#1 Sholl, D.A., 1953. Dendritic organization in the neurons of the visual and motor cortices of the cat. J. Anat. 87, 387–406.

#2 https://imagej.net/Sholl_Analysis