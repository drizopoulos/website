+++
abstract = "Intensive longitudinal data from ecological momentary assessment (EMA) is widely used in clinical research but often suffers from dropout, leading to reduced statistical power, invalid results, and poor treatment outcomes. Predicting dropout could help with its prevention. While existing methods utilise baseline covariates, few studies account for the temporal dynamics of EMA data or identify the exact timing of dropout. Joint models (JM) enable simultaneous modelling of longitudinal processes and time-to-event data, offering dynamic predictions. However, conventional JMs assume limited measurement occasions and do not account for the autocorrelation inherent in EMA data. We extended the standard JM by incorporating an autoregressive submodel, capturing temporal dependencies in EMA measurements. We validated our approach through simulation studies, demonstrating good parameter recovery across different missingness mechanisms (MCAR, MAR, MNAR) and high dropout prediction accuracy. We applied the JM to an existing empirical EMA dataset, using baseline (e.g., depression) and time-varying (affect, intermittent missingness) predictors of dropout. The extended JM outperformed a baseline-only survival model in predicting dropout. The sensitivity analysis of the missingness mechanism revealed that fixed effect estimates remained stable across different missing data mechanisms, whereas random effect estimates for autocorrelation were sensitive to these assumptions. By integrating autoregressive components, the extended JM accommodates temporal dependencies and dynamically updates predictions of dropout risk. This approach improves dropout prediction in EMA studies and highlights the importance of utilising JMs for predicting clinically relevant outcomes while integrating EMA data."
abstract_short = "Psychological Assessment, to appear"
authors = ["F Petersen", "L Bringmann", "D Rizopoulos"]
date = "2025-05-25"
math = true
publication_types = ["2"]
publication = "Psychological Assessment, to appear"
publication_short = "Psychological Assessment, to appear"
selected = true
title = "Predicting dropout in intensive longitudinal data: Extending the joint model for auto-correlated data"
url_pdf = "https://doi.org/10.1037/pas0001397"
+++