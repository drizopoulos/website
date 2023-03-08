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
lmeFit <- lme(pro ~ time + time:treat, data = prothro,
                 random = ~ time | id)

# the Cox model
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros)

# the joint model
jointFit <- jm(survFit, lmeFit, time_var = "time")

summary(jointFit)


# the data of Patient 155
dataP155 <- prothro[prothro$id == 155, ]
dataP155$Time <- dataP155$death <- NULL


# survival probabilities using the first measurement
Spred <- predict(jointFit, newdata = dataP155[1, ],
                process = "event",
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Spred
plot(Spred)

# longitudinal predictions using the first measurement
Lpred <- predict(jointFit, newdata = dataP155[1, ],
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Lpred
plot(Lpred)

# combine in one plot
plot(Lpred, Spred)


# dynamic predictions, each time an extra measurement
n <- nrow(dataP155)
for (i in seq_len(n)) {
    Spred <- predict(jointFit, newdata = dataP155[1:i, ],
                    process = "event",
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)
    Lpred <- predict(jointFit, newdata = dataP155[1:i, ],
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)

    plot(Lpred, Spred)
}


# ROC at 3 years with 1 year window
roc <- tvROC(jointFit, newdata = prothro, Tstart = 3, Dt = 1)
roc
plot(roc)
# AUC at 3 years with 1 year window
tvAUC(roc)

# Calibration plot at 3 years with 1 year window
calibration_plot(jointFit, newdata = prothro, Tstart = 3, Dt = 1)

# Brier score at 3 years with 1 year window
tvBrier(jointFit, newdata = prothro, Tstart = 3, Dt = 1)

