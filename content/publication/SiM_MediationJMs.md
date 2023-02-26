+++
abstract = "Many clinical trials repeatedly measure several longitudinal outcomes on patients. Patient follow-up can discontinue due to an outcome-dependent event, such as clinical diagnosis, death, or dropout. Joint modeling is a popular choice for the analysis of this type of data. Using example data from a prodromal Alzheimer's disease trial, we propose a new type of multivariate joint model in which longitudinal brain imaging outcomes and memory impairment ratings are allowed to be associated both with time to open-label medication and dropout, and where the brain imaging outcomes may also directly affect the memory impairment ratings. Existing joint models for multivariate longitudinal outcomes account for the correlation between the longitudinal outcomes through the random effects, often by assuming a multivariate normal distribution. However, for these models, it is difficult to interpret how the longitudinal outcomes affect each other. We model the dependence between the longitudinal outcomes differently so that a first longitudinal outcome affects a second one. Specifically, for each longitudinal outcome, we use a linear mixed-effects model to estimate its trajectory, where, for the second longitudinal outcome, we include the linear predictor of the first outcome as a time-varying covariate. This facilitates an easy and direct interpretation of the association between the longitudinal outcomes and provides a framework for latent mediation analysis to understand the underlying biological processes. For the trial considered here, we found that part of the intervention effect is mediated through hippocampal brain atrophy. The proposed joint models are fitted using a Bayesian framework via MCMC simulation."
abstract_short = "Statistics in Medicine 39, 4120-4132"
authors = ["F van Oudenhoven", "SHN Swinkels", "T Hartmann", "D Rizopoulos"]
date = "2022-07-30"
math = true
publication_types = ["2"]
publication = "Statistics in Medicine 41, 3421-3433"
publication_short = "Statistics in Medicine 41, 3421-3433"
selected = true
title = "Modeling the underlying biological processes in Alzheimer's disease using a multivariate competing risk joint model"
url_pdf = "https://doi.org/10.1002/sim.9425"
+++