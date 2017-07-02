###########################################################################################
# Description: This R script contains the code included in the slides (+ extra code) for  #
#   the short course 'An Introduction to Joint Models for Longitudinal & Survival Data,   #
#    with Applications in R' and illustrates the basic capabilities of package JM         #
# Author: Dimitris Rizopoulos                                                             #
# Last update: 2015-09-30                                                                 #
###########################################################################################


# load first package JM to make everything available
library("JM")

########################################################
# Section 2.3: Fitting linear mixed models using lme() #
########################################################

# random intercepts + random slopes model for the AIDS dataset
lmeFit <- lme(CD4 ~ obstime + obstime:drug, data = aids,
    random = ~ obstime | patient)

summary(lmeFit)

# random intercepts alone
lmeFit.int <- lme(CD4 ~ obstime + obstime:drug, data = aids,
    random = ~ 1 | patient)

summary(lmeFit.int)

# random intercepts + random slopes with diagonal covariance matrix
lmeFit.diag <- lme(CD4 ~ obstime + obstime:drug, data = aids, 
    random = list(patient = pdDiag(form = ~ obstime)))

summary(lmeFit.diag)


###########################################################
# Section 3.3: Fitting relative risk models using coxph() #
###########################################################

# Cox PH model for the PBC dataset (define first the composite event 
# 'status2')
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")
coxFit <- coxph(Surv(years, status2) ~ drug + sex + age, data = pbc2.id)

summary(coxFit)


########################################################
# Section 3.3: Fitting joint models using jointModel() #
########################################################

# joint model for the AIDS dataset: First we fit the linear mixed model as above
lmeFit.aids <- lme(CD4 ~ obstime + obstime:drug, 
    random = ~ obstime | patient, data = aids)

# then we fit a Cox model for the time-to-death, including the argument 'x = TRUE'
coxFit.aids <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# jointModel() takes the above fitted models as arguments, and fits the
# joint model; below we fit a joint model with a relative risk submodel
# for the event time outcome, in which the baseline risk function is assumed
# piecewise-constant
jointFit.aids <- jointModel(lmeFit.aids, coxFit.aids, timeVar = "obstime", 
    method = "piecewise-PH-aGH")

summary(jointFit.aids)

# time-dependent Cox model for the AIDS dataset
td.Cox <- coxph(Surv(start, stop, event) ~ drug + CD4, data = aids)

summary(td.Cox)


#######################################
# Section 5.1: Extending joint models #
#######################################

# lagged effect for the PBC dataset
lmeFit.pbc <- lme(log(serBilir) ~ year, random = ~ year | id, data = pbc2)
coxFit.pbc <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)
jointFit.pbc <- jointModel(lmeFit.pbc, coxFit.pbc, timeVar = "year", 
    method = "piecewise-PH-aGH", lag = 1)

summary(jointFit.pbc)

# Time-dependent slopes parameterization: We refit the standard joint model 
# for the PBC dataset include also drug as a covariate

# first we fit the linear mixed effects and Cox models separately
lmeFit.pbc2 <- lme(log(serBilir) ~ drug * year, random = ~ year | id, data = pbc2)
coxFit.pbc2 <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# the default parameterization assumes that the risk for an event at
# time t depends on the true value of log serum Bilirubin at the same
# time point (by true value we mean the value without measurement error).
# For the survival submodel we postulate a relative risk model, with a 
# piecewise-constant baseline risk function
jointFit.pbc2 <- jointModel(lmeFit.pbc2, coxFit.pbc2, timeVar = "year", 
    method = "piecewise-PH-aGH")

summary(jointFit.pbc2)

# we can extend the above standard parameterization by assuming that
# the risk for an event at time t depends on both the true value of
# log serum Bilirubin *and* the true value of the slope of the
# patient-specific trajectory at the same time point.
# To fit this joint model which first need to specify the 'derivForm'
# argument of jointModel(), which is a list with first component the
# derivative of the 'fixed' argument of 'lmeFit' wrt 'year', second
# component the index of values of the fixed effects corresponding to
# the derivative, third and fourth components the same but for the 
# random effects structure. For the above model we have:
dform <- list(fixed = ~ drug, indFixed = 3:4, random = ~ 1, indRandom = 2)
jointFit.pbc2.slp <- jointModel(lmeFit.pbc2, coxFit.pbc2, timeVar = "year", 
    method = "piecewise-PH-aGH", parameterization = "both", derivForm = dform)

# the output of the summary() method is augmented with the 
# coefficient 'Assoct.s' that corresponds to the association parameter
# for the slope
summary(jointFit.pbc2.slp)

# an LRT can be performed with anova()
anova(jointFit.pbc2, jointFit.pbc2.slp)


#################################
# Section 5.5: AFT joint models #
#################################

# a joint model with an accelerated failure time model for the survival
# outcome
lmeFit.pbc3 <- lme(log(serBilir) ~ year, random = ~ year | id, data = pbc2)
coxFit.pbc3 <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)
jointFit.pbc3 <- jointModel(lmeFit.pbc3, coxFit.pbc3, timeVar = "year", 
    method = "weibull-AFT-aGH")

summary(jointFit.pbc3)


#############################################################
# Section 6.2: Dynamic predictions for the survival outcome #
#############################################################

# Dynamic predictions of survival probabilities: we see how the survival probabilities
# for Patient 2 are updated as more longitudinal information is recorded

# first we fit the joint model with linear + quadratic slopes
lmeFit.quad.pbc <- lme(log(serBilir) ~ drug * (year + I(year^2)), data = pbc2,
    random = list(id = pdDiag(form = ~ year + I(year^2))))
coxFit.pbc2 <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
jointFit.quad.pbc <- jointModel(lmeFit.quad.pbc, coxFit.pbc2, timeVar = "year", 
    method = "piecewise-PH-aGH")

ND <- pbc2[pbc2$id %in% "2", ] # data of Patient 2

# calculate dynamic predictions using survfitJM()
slis <- vector("list", nrow(ND))
for (i in seq_along(slis)) {
    set.seed(123)
    ND.i <- ND[1:i, ]
    slis[[i]] <- survfitJM(jointFit.quad.pbc, newdata = ND.i, idVar = "id", M = 500)
}

# plot the resulting estimates
labs <- c("baseline available", "2nd measurement available", "3rd measurement available",
    paste(4:9, "th measurement available", sep = ""))
op <- par(mfrow = c(3, 3), oma = c(0, 2, 2, 0))
for (i in seq_along(labs)) {
    plot(slis[[i]], estimator = "median", conf.int = TRUE, lwd = 2, xlab = "Time", ylab = "", 
        main = labs[i], col = 1)
    grid()
}
mtext(expression(paste("Pr(", T[i]^symbol("*")  >= u, " | ", T[i]^symbol("*") > t, ", ", Y[i](t), ", ", 
    D[n], ")", sep = " ")), 2, cex = 1.5, line = -1, outer = TRUE)
par(op)


#################################################################
# Section 6.3: Dynamic predictions for the longitudinal outcome #
#################################################################

# we will again use Patient 2, and the same joint model as in Section 6.2

# we calculate dynamic predictions of future longitudinal responses using
# the first 5 measurements of Patient 2
lfit <- predict(jointFit.quad.pbc, newdata = ND[1:5, ], 
    type = "Subject", interval = "conf", returnData = TRUE)
lfit

library(lattice)
xyplot(pred + low + upp ~ year, data = lfit, type = "l", 
    lty = c(1,2,2), col = 1, lwd = 2)


###################################
# Section 6.6: Dynamic ROC curves #
###################################

# we will again use the same joint model as in Section 6.2

# we create a dataset with the time points at which we expect
# longitudinal measurements + covariate information
NewData <- expand.grid(
    year = c(0, 0.5, 2, 3, 5),
    drug = c("placebo", "D-penicil")
)
NewData$id <- rep(1:2, each = 5)

roc <- rocJM(jointFit.quad.pbc, data = NewData, dt = c(1, 2, 4))
roc

plot(roc)
