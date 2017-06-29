+++
# Date this page was created.
date = "2017-06-20"

# Project title.
title = "Multivariate Joint Models"

# Project summary to display on homepage.
summary = "Computational methods for multivariate joint models."

# Optional image to display on homepage (relative to `static/img/` folder).
image_preview = "multivariate.jpeg"

# Tags: can be used for filtering projects.
# Example: `tags = ["machine-learning", "deep-learning"]`
tags = ["mvJM"]

# Optional external URL for project (replaces project detail page).
external_link = ""

# Does the project detail page use math formatting?
math = false

+++

Joint models for longitudinal and survival data have gained a lot of attention the recent years. There have been extended to handle among others multivariate longitudinal data, competing risks and recurrent events, and nowadays there also exist several freely available software packages for their implementation. 

From the aforementioned extensions, the one that is most practically relevant is the multivariate longitudinal data one. Even though this extension is mathematically straightforward, from a computational viewpoint joint models with multiple longitudinal outcomes remain difficult to fit in practice due to the high number of random effects they require. This difficulty has also hampered to a degree their practical application. 

In this project we work on novel approaches that enable fitting such joint models in realistic computing times. The idea behind our methods is to split the estimation in two steps, first to estimate a multivariate mixed model for the longitudinal outcomes, and then use the output of this model to fit the survival submodel. Such two-stage approaches have been previously proposed in the literature and have been shown to be biased. What is different in our approach is a correction we apply in the resulting estimates that transform them to the estimates we would expect to obtain if we would fit the multivariate joint model. The type of corrections we use are based on importance sampling ideas.

