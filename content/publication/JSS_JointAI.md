+++
abstract = "Missing data occur in many types of studies and typically complicate the analysis. Multiple imputation, either using joint modeling or the more flexible fully conditional specification approach, are popular and work well in standard settings. In settings involving nonlinear associations or interactions, however, incompatibility of the imputation model with the analysis model is an issue often resulting in bias. Similarly, complex outcomes such as longitudinal or survival outcomes cannot be adequately handled by standard implementations. In this paper, we introduce the R package JointAI, which utilizes the Bayesian framework to perform simultaneous analysis and imputation in regression models with incomplete covariates. Using a fully Bayesian joint modeling approach it overcomes the issue of uncongeniality while retaining the attractive flexibility of fully conditional specification multiple imputation by specifying the joint distribution of analysis and imputation models as a sequence of univariate models that can be adapted to the type of variable. JointAI provides functions for Bayesian inference with generalized linear and generalized linear mixed models and extensions thereof as well as survival models and joint models for longitudinal and survival data, that take arguments analogous to the corresponding well known functions for the analysis of complete data from base R and other packages. Usage and features of JointAI are described and illustrated using various examples and the theoretical background is outlined."
abstract_short = "Journal of Statistical Software (2021)"
authors = ["NS Erler", "D Rizopoulos", "E Lesaffre"]
date = "2021-11-30"
math = true
publication_types = ["2"]
publication = "Journal of Statistical Software 100(20), 1-56"
publication_short = "Journal of Statistical Software (2021)"
selected = true
title = "JointAI: Joint analysis and imputation of incomplete data in R "
url_pdf = "https://doi.org/10.18637/jss.v100.i20"
+++