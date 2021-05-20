# load package JM
library("JMbayes2")

###############
# Practical 1 #
###############

# the linear mixed model
lmeFit.p1 <- lme(log(serBilir) ~ year + drug:year,
                 data = pbc2, random = ~ year | id)

# create the indicator for the composite event
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# the Cox model
survFit.p1 <- coxph(Surv(years, status2) ~ drug, data = pbc2.id)

# the joint model
jointFit.p1 <- jm(survFit.p1, lmeFit.p1, time_var = "year")

summary(jointFit.p1)

# Hazard Ratios
stab <- summary(jointFit.p1)$Survival
exp(stab[c(1,3,4)])


# interaction term between serBilir and drug
jointFit2.p1 <- update(jointFit.p1,
                       functional_forms = ~ value(log(serBilir)):drug)

summary(jointFit2.p1)


###############
# Practical 3 #
###############

# the linear mixed model
lmeFit.p3 <- lme(log(serBilir) ~ year + I(year^2), data = pbc2,
    random = list(id = pdDiag(form = ~ year + I(year^2))))

# create the indicator for the composite event
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# the Cox model
survFit.p3 <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)

# the joint model
jointFit.p3 <- jm(survFit.p3, lmeFit.p3, time_var = "year")

summary(jointFit.p3)


# the joint model that include both terms
jointFit2.p3 <- update(jointFit.p3,
    functional_forms = ~ value(log(serBilir)) + slope(log(serBilir)))

summary(jointFit2.p3)


# the linear mixed model with natural cubic splines
lmeFit2.p3 <- lme(log(serBilir) ~ ns(year, 3), data = pbc2,
    random = list(id = pdDiag(form = ~ ns(year, 3))))

# the joint model with splines that include both terms
jointFit3.p3 <- jm(survFit.p3, lmeFit2.p3, time_var = "year",
    functional_forms = ~ value(log(serBilir)) + slope(log(serBilir)),
    n_iter = 6500L, n_burnin = 2500L)

summary(jointFit3.p3)


###############
# Practical 4 #
###############

# the linear mixed model
lmeFit.p4 <- lme(pro ~ time + time:treat, data = prothro,
    random = ~ time | id)

# the Cox model
survFit.p4 <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)

# the joint model
jointFit.p4 <- jm(survFit.p4, lmeFit.p4, time_var = "time")

summary(jointFit.p4)


# the data of Patient 155
dataP155 <- prothro[prothro$id == 155, ]

# survival probabilities using only the 1st measurement
sfit <- predict(jointFit.p4, newdata = dataP155[1:2, ],
                process = "event",
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

sfit
plot(sfit)

roc <- tvROC(jointFit.p4, newdata = prothro, Tstart = 3, Dt = 1,
             n_samples = 600L, n_mcmc = 85L)
plot(roc)
# AUC at 2 years with 0.5 year window
tvAUC(roc)





