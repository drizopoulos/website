+++
abstract = "Joint models for longitudinal and survival data have become a popular framework for studying the association between repeatedly measured biomarkers and clinical events. Nevertheless, addressing complex survival data structures, especially handling both recurrent and competing event times within a single model, remains a challenge. This causes important information to be disregarded. Moreover, existing frameworks rely on a Gaussian distribution for continuous markers, which may be unsuitable for bounded biomarkers, resulting in biased estimates of associations. To address these limitations, we propose a Bayesian shared-parameter joint model that simultaneously accommodates multiple (possibly bounded) longitudinal markers, a recurrent event process, and competing risks. We use the beta distribution to model responses bounded within any interval (a, b) without sacrificing the interpretability of the association. The model offers various forms of association, discontinuous risk intervals, and both gap and calendar timescales. A simulation study shows that it outperforms simpler joint models. We utilize the US Cystic Fibrosis Foundation Patient Registry to study the associations between changes in lung function and body mass index, and the risk of recurrent pulmonary exacerbations, while accounting for the competing risks of death and lung transplantation. Our efficient implementation allows fast fitting of the model despite its complexity and the large sample size from this patient registry. Our comprehensive approach provides new insights into cystic fibrosis disease progression by quantifying the relationship between the most important clinical markers and events more precisely than has been possible before. The model implementation is available in the R package JMbayes2."
abstract_short = "Statistics in Medicine, to appear"
authors = ["P Afonso", "D Rizopoulos", "A Palipana", "E Gecili", "C Brokamp", "J Clancy", "R Szczesniak", "ER Andrinopoulou"]
date = "2025-04-25"
math = true
publication_types = ["2"]
publication = "Statistics in Medicine, to appear"
publication_short = "Statistics in Medicine, to appear"
selected = true
title = "A joint model for (un)bounded longitudinal markers, competing risks, and recurrent events using patient registry data"
url_pdf = "https://doi.org/10.1002/sim.70057"
+++