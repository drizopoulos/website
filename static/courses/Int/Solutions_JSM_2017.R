# load package JM
library("JM")

###############
# Practical 1 #
###############

# the linear mixed model
lmeFit.p1 <- lme(log(serBilir) ~ year + drug:year, data = pbc2,
    random = ~ year | id)

# create the indicator for the composite event
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# the Cox model
survFit.p1 <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# the joint model
jointFit.p1 <- jointModel(lmeFit.p1, survFit.p1, timeVar = "year",
    method = "piecewise-PH-aGH")
    
summary(jointFit.p1)

confint(jointFit.p1, parm = "Longitudinal")
exp(confint(jointFit.p1, parm = "Event"))


# interaction term between serBilir and drug
intFact <- list(value = ~ drug, data = pbc2.id)
jointFit2.p1 <- update(jointFit.p1, interFact = intFact)

summary(jointFit2.p1)

# the joint model under the null H_0: beta_2 = 0
lmeFit2.p1 <- lme(log(serBilir) ~ year, data = pbc2,
    random = ~ year | id)
jointFit3.p1 <- update(jointFit2.p1, lmeObject = lmeFit2.p1)

anova(jointFit3.p1, jointFit2.p1)

# the joint model under the null H_0: gamma = alpha_2 = 0
survFit2.p1 <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)
jointFit4.p1 <- update(jointFit2.p1, survObject = survFit2.p1, interFact = NULL)

anova(jointFit4.p1, jointFit2.p1)

# the joint model under the null H_0: beta_2 = gamma = alpha_2 = 0
jointFit5.p1 <- update(jointFit2.p1, lmeObject = lmeFit2.p1, 
    survObject = survFit2.p1, interFact = NULL)

anova(jointFit5.p1, jointFit2.p1)


###############
# Practical 2 #
###############

# the linear mixed model
lmeFit.p4 <- lme(pro ~ time + time:treat, data = prothro,
    random = ~ time | id)

# the Cox model
survFit.p4 <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)

# the joint model
jointFit.p4 <- jointModel(lmeFit.p4, survFit.p4, timeVar = "time",
    method = "piecewise-PH-aGH")

summary(jointFit.p4)


# the data of Patient 155
dataP155 <- prothro[prothro$id == 155, ]

# survival probabilities using only the 1st measurement
sfit <- survfitJM(jointFit.p4, newdata = dataP155[1, ])

sfit
plot(sfit)
plot(sfit, include.y = TRUE)

# we will use a for-loop to update predictions
pdf()
for (i in 1:nrow(dataP155)) {
    data.i <- dataP155[1:i, ]
    sfit.i <- survfitJM(jointFit.p4, newdata = data.i)
    plot(sfit.i, estimator = "mean", include.y = TRUE,
        conf.int = TRUE, fill.area = TRUE, col.area = "lightgrey", 
        ylab2 = "Prothrobin")
}
dev.off()


# predictions of future longitudinal responses, using only the 1st measurement
lfit <- predict(jointFit.p4, newdata = dataP155[1, ], type = "Subject",
    interval = "conf", returnData = TRUE)
    
library("lattice")
xyplot(pred + low + upp ~ time, data = lfit, type = "l",
    lty = c(1,2,2), col = 1, lwd = 2)

# we will use a for-loop to update predictions
pdf()
for (i in 1:nrow(dataP155)) {
    data.i <- dataP155[1:i, ]
    lfit.i <- predict(jointFit.p4, newdata = data.i, type = "Subject",
        interval = "conf", returnData = TRUE)
    print(xyplot(pred + low + upp ~ time, data = lfit.i, type = "l",
        lty = c(1,2,2), col = 1, lwd = 2), 
        xlab = "Time", ylab = "Predicted Prothrobin Level")
}
dev.off()

# AUC at 2 years with 0.5 year window
aucJM(jointFit.p4, newdata = prothro, Tstart = 2, Dt = 0.5)

# prediction error at 2 years with 0.5 year window
prederrJM(jointFit.p4, newdata = prothro, Tstart = 2, Thoriz = 2.5)




