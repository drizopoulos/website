+++
abstract = "Joint models for longitudinal and survival data have garnered a lot of attention in recent years, with the development of myriad extensions to the basic model, including those which allow for multivariate longitudinal data, competing risks and recurrent events. Several software packages are now also available for their implementation. Although mathematically straightforward, the inclusion of multiple longitudinal outcomes in the joint model remains computationally difficult due to the large number of random effects required, which hampers the practical application of this extension. We present a novel approach that enables the fitting of such models with more realistic computational times. The idea behind the approach is to split the estimation of the joint model in two steps; estimating a multivariate mixed model for the longitudinal outcomes, and then using the output from this model to fit the survival submodel. So called two-stage approaches have previously been proposed, and shown to be biased. Our approach differs from the standard version, in that we additionally propose the application of a correction factor, adjusting the estimates obtained such that they more closely resemble those we would expect to find with the multivariate joint model. This correction is based on importance sampling ideas. Simulation studies show that this corrected-two-stage approach works satisfactorily, eliminating the bias while maintaining substantial improvement in computational time, even in more difficult settings."
abstract_short = "Statistics and Computing, under review"
authors = ["K Mauff", "EW Steyerberg", "I Kardys", "H Boersma", "D Rizopoulos"]
date = "2020-03-04"
math = true
publication_types = ["2"]
publication = "Statistics and Computing 30, 999–1014"
publication_short = "Statistics and Computing (2020)"
selected = true
title = "Joint models with multiple longitudinal outcomes and a time-to-event outcome: a corrected two-stage approach"
url_pdf = "https://link.springer.com/article/10.1007/s11222-020-09927-9"
url_slides = "https://emcbiostatistics.shinyapps.io/IBC2018_Dimitris/"
url_project = "project/mvJM/"
+++