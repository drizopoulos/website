################################################################################################
# Sample joint modelling analysis for the Prothrobin Data                                      #
# Description and source: Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993).     #
#   Statistical Models Based on Counting Processes. New York: Springer.                        #
################################################################################################


# install the latest version of the 'JM' package
# (if you do not already have it)
install.packages("JM", dependencies = TRUE)

# load packages 'JM' & 'lattice'
library("JM")
library("lattice")


########################
# Descriptive Analysis #
########################


# for the longitudinal process

# number of measurements per patient
ni <- with(prothro, tapply(pro, id, length))
summary(ni)

# Visiting Process Plots
# first we split the dead patients in 4 groups
prothro$DeadGroups <- rep(NA, nrow(prothro))
prothro$DeadGroups[prothro$death == 1] <- rep(1:4, 
    sum(prothro$death == 1)%/%4, length.out = sum(prothro$death == 1))

# we plot per group the visit times
xyplot(id ~ time | DeadGroups, data = prothro, groups = id, subset = death == 1,
    type = "b", xlab = "Visit Time (months)", ylab = "Patients", 
    main = "Visiting Process")

# we plot per group the death time - visit times
# (to check if more measurements were taken closer to
# the event time)
xyplot(id ~ I(Time - time) | DeadGroups, data = prothro, groups = id, subset = death == 1,
    type = "b", xlab = "Death Time - Visit Times (months)", ylab = "Patients", 
    main = "Visiting Process")


# subject-specific longitudinal profiles for a random sample of
# 10 patients (each time you run the code a new sample of 10 patients
# is used)
xyplot(pro ~ time, group = id, data = prothro[prothro$id %in% sample(prothros$id, 10), ], 
    xlab = "Months", type = "l", ylim = range(prothro$pro))

# smooth average longitudinal evolutions per treatment
xyplot(pro ~ time | treat, data = prothro, xlab = "Months", 
    type = c("p", "smooth"), lwd = 2)

# smooth average longitudinal evolutions per treatment
# and death status
prothro$deathF <- factor(prothro$death, labels = c("censored", "dead"))
xyplot(pro ~ time | deathF , data = prothro, group = treat, xlab = "Months", 
    type = c("p", "smooth"), auto.key = TRUE, lwd = 2)


# smooth average longitudinal evolutions per treatment
# and dropout time
prothro$binTime <- cut(prothro$Time, quantile(prothros$Time), 
    include.lowest = TRUE)
xyplot(pro ~ time | binTime*deathF, data = prothro, group = treat, xlab = "Months", 
    type = c("p", "smooth"), auto.key = TRUE, lwd = 2, cex = 0.3)


# for the survival process
sf <- survfit(Surv(Time, death) ~ treat, data = prothros)
sf

# Kaplan-Meier estimate for each treatment group
plot(sf, mark.time = FALSE, col = 1:2, lty = 1:2, lwd = 2, 
    ylab = "Survival", xlab = "Months") 
legend("topright", c("placebo", "prednisone"), lty = 1:2, lwd = 2, 
    col = 1:2, bty = "n")


##########################################
# Cox analysis treating prothrobin as an #
# exogenous time-dependent covariate     #
##########################################

timdepnt.Cox <- coxph(Surv(start, stop, event) ~ treat + pro, data = prothro)
summary(timdepnt.Cox)


############################
# Joint Modelling Analysis #
############################

# first we fit the linear mixed model to describe the prothrobin
# evolutions; we include an indicator for the baseline measurement
prothro$t0 <- as.numeric(prothro$time == 0)
lmeFit <- lme(pro ~ treat * (time + t0), random = ~ time | id, data = prothro)

# then we fit a Cox model for the time-to-death
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)

# jointModel() takes the above fitted models as arguments, and fits the
# joint model; below we fit two joint models with a relative risk submodel
# for the event time outcome. In the first the baseline risk function is assumed
# piecewise-constant, whereas in the second it is approximated with B-splines
# (the following requires some time to run)
fitJoint.pw <- jointModel(lmeFit, survFit, timeVar = "time", method = "piecewise-PH-aGH")
fitJoint.sp <- jointModel(lmeFit, survFit, timeVar = "time", method = "spline-PH-aGH")

# we compare the two joint models; these are not nested and therefore,
# we do not perform a likelihood ratio test -- the joint model with the
# piecewise-constant baseline risk function seems a bit better according
# to both AIC and BIC
anova(fitJoint.pw, fitJoint.sp, test = FALSE)

# detailed results are returned by the summary() method
summary(fitJoint.pw)

# asymptotic 95% CIs for the parameters of the event process;
# we also compare with the time-dependent Cox model
confint(fitJoint.pw, "Event")
confint(timdepnt.Cox)


# check the fit for the event time outcome in the marginal survival function
par(mfrow = c(1, 2))
plot(fitJoint.pw, 3, add.KM = TRUE, col = 2, lwd = 2, main = "Piecewise-constant Baseline Risk Function")
plot(fitJoint.sp, 3, add.KM = TRUE, col = 2, lwd = 2, main = "Spline-approximated Baseline Risk Function")

# check the fit for the event time outcome using the Cox-Snell residuals
par(mfrow = c(1, 2))

res.surv <- residuals(fitJoint.pw, process = "Event", type = "Cox")
sfit <- survfit(Surv(res.surv, fitJoint.pw$y$d) ~ 1)
plot(sfit, mark.time = FALSE, conf.int = TRUE, lty = 1:2, 
    xlab = "Cox-Snell Residuals", ylab = "Survival Probability", 
    main = "Piecewise-constant Baseline Risk Function")
curve(pexp(x, lower.tail = FALSE), from = min(res.surv), to = max(res.surv), 
    add = TRUE, col = "red", lwd = 2)
legend("topright", c("Survival function of unit Exponential", 
    "Survival function of Cox-Snell Residuals"), lty = 1, lwd = 2, col = 2:1, bty = "n")

res.surv <- residuals(fitJoint.sp, process = "Event", type = "Cox")
sfit <- survfit(Surv(res.surv, fitJoint.sp$y$d) ~ 1)
plot(sfit, mark.time = FALSE, conf.int = TRUE, lty = 1:2, 
    xlab = "Cox-Snell Residuals", ylab = "Survival Probability", 
    main = "Spline-approximated Baseline Risk Function")
curve(pexp(x, lower.tail = FALSE), from = min(res.surv), to = max(res.surv), 
    add = TRUE, col = "red", lwd = 2)
legend("topright", c("Survival function of unit Exponential", 
    "Survival function of Cox-Snell Residuals"), lty = 1, lwd = 2, col = 2:1, bty = "n")


# check the residuals for longitudinal outcome: subject-specific residuals vs 
# fitted values
par(mfrow = c(1, 2))
plot(fitJoint.pw, 1, main = "Piecewise-constant Baseline Risk Function")
plot(fitJoint.sp, 1, main = "Spline-approximated Baseline Risk Function")

# check the residuals for longitudinal outcome: QQ plot of subject-specific residuals
res.Ypw <- residuals(fitJoint.pw, type = "stand-Subject")
res.Ysp <- residuals(fitJoint.sp, type = "stand-Subject")

par(mfrow = c(1, 2))
qqnorm(res.Ypw, main = "Piecewise-constant Baseline Risk Function"); qqline(res.Ypw)
qqnorm(res.Ysp, main = "Spline-approximated Baseline Risk Function"); qqline(res.Ysp)



# Conditional survival probabilities

# a patient who died
ND.dead <- prothro[prothro$id == 155, ]

# probability of surviving at least one month after the last
# longitudinal measurement was recorded
survProbs <- lapply(seq_len(nrow(ND.dead)), function (i) {
    Di <- ND.dead[1:i, ]
    survfitJM(fitJoint.pw, newdata = Di, survTimes = max(Di$time) + 1, simulate = FALSE)
})

vals <- sapply(survProbs, function (x) x$summaries[[1]][, 2])

layout(cbind(c(1,2)), heights = c(1.6, 1))
plot(ND.dead$time, vals, col = 2, cex = 1.2, pch = 16, ylim = c(0, 1), type = "b", 
    xlab = "", ylab = "Survival", main = "Extra 1 Month Survival")
abline(h = c(0.5, 0.7, 0.9), lty = 2, col = "grey")
plot(pro ~ time, data = ND.dead, pch = 8, xlab = "Time (months)", ylim = range(prothro$pro))


# the same but with 95% pointwise CIs for the extra one month
# survival probabilities
set.seed(1234)
survProbs.ci <- lapply(seq_len(nrow(ND.dead)), function (i) {
    Di <- ND.dead[1:i, ]
    survfitJM(fitJoint.pw, newdata = Di, survTimes = max(Di$time) + 1, M = 500)
})

vals.ci <- t(sapply(survProbs.ci, function (x) x$summaries[[1]][, 3:5]))

library(plotrix)
layout(cbind(c(1,2)), heights = c(1.6, 1))
plotCI(ND.dead$time, vals.ci[, 1], li = vals.ci[, 2], ui = vals.ci[, 3], scol = 4, lwd = 2, 
    col = 2, cex = 1.2, pch = 16, ylim = c(0, 1), xlab = "", ylab = "Survival",
    main = "Extra 1 Month Survival")
lines(ND.dead$time, vals.ci[, 1], col = 2)
abline(h = c(0.5, 0.7, 0.9), lty = 2, col = "grey")
plot(pro ~ time, data = ND.dead, pch = 8, xlab = "Time (months)", ylim = range(prothro$pro))

