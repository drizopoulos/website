################################################################################
# Title: Solutions to the practicals of the Joint Models in R Course           #
# Author: Dimitris Rizopoulos                                                  #
# Requirements: the latest version of R and R package JMbayes2;                #
# R can be installed from CRAN https://cran.r-project.org/ and package         #
# JMbayes2 by running the command: install.packages("JMbayes2")                #
################################################################################


# load package JMbayes2
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
jointFit1.p1 <- jm(survFit.p1, lmeFit.p1, time_var = "year")

summary(jointFit1.p1)

# Hazard Ratios
stab <- summary(jointFit1.p1)$Survival
exp(stab[c(1,3,4)])


# interaction term between serBilir and drug
jointFit2.p1 <- update(jointFit1.p1,
                       functional_forms = ~ value(log(serBilir)) * drug)

summary(jointFit2.p1)

compare_jm(jointFit1.p1, jointFit2.p1)


###############
# Practical 2 #
###############

# the linear mixed model
lmeFit.p2 <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 14.4)), data = pbc2,
    random = ~ ns(year, 3, B = c(0, 14.4)) | id,
    control = lmeControl(opt = "optim"))

# create the indicator for the composite event
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# the Cox model
survFit.p2 <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)

# the joint model
jointFit1.p2 <- jm(survFit.p2, lmeFit.p2, time_var = "year")

summary(jointFit1.p2)


# the joint model that include both terms
jointFit2.p2 <- update(jointFit1.p2, n_iter = 8500L, n_burnin = 3500L,
    functional_forms = ~ value(log(serBilir)) + slope(log(serBilir)))

summary(jointFit2.p2)

jointFit3.p2 <- update(jointFit2.p2,
functional_forms = ~ value(log(serBilir)) +
    slope(log(serBilir), direction = "back", eps = 1))

summary(jointFit3.p2)


jointFit4.p2 <- update(jointFit2.p2,
                       functional_forms = ~ area(log(serBilir)))

summary(jointFit4.p2)


###############
# Practical 3 #
###############

# the linear mixed model
lmeFit.p3 <- lme(pro ~ time + time:treat, data = prothro,
    random = ~ time | id)

# the Cox model
survFit.p3 <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)

# the joint model
jointFit.p3 <- jm(survFit.p3, lmeFit.p3, time_var = "time")

summary(jointFit.p3)


# the data of Patient 155
dataP155 <- prothro[prothro$id == 155, ]
dataP155$Time <- dataP155$death <- NULL


# survival probabilities using the first measurement
Spred <- predict(jointFit.p3, newdata = dataP155[1, ],
                process = "event",
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Spred
plot(Spred)

# longitudinal predictions using the first measurement
Lpred <- predict(jointFit.p3, newdata = dataP155[1, ],
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Lpred
plot(Lpred)

# combine in one plot
plot(Lpred, Spred)


# dynamic predictions, each time an extra measurement
n <- nrow(dataP155)
for (i in seq_len(n)) {
    Spred <- predict(jointFit.p3, newdata = dataP155[1:i, ],
                    process = "event",
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)
    Lpred <- predict(jointFit.p3, newdata = dataP155[1:i, ],
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)

    plot(Lpred, Spred)
}


# ROC at 3 years with 1 year window
roc <- tvROC(jointFit.p3, newdata = prothro, Tstart = 3, Dt = 1)
roc
plot(roc)
# AUC at 3 years with 1 year window
tvAUC(roc)

# Calibration plot at 3 years with 1 year window
calibration_plot(jointFit.p3, newdata = prothro, Tstart = 3, Dt = 1)

# Brier score at 3 years with 1 year window
tvBrier(jointFit.p3, newdata = prothro, Tstart = 3, Dt = 1)

