################################################################################
# Title: Solutions to the practicals of the Dynamic Risk Predictions course    #
# Authors: Dimitris Rizopoulos and Christos Thomadakis                         #
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

# the first longitudinal measurement
first <- dataP155[1, ]
first$Time <- max(first$time)
first$death <- 0

# survival probabilities using the first measurement
Spred <- predict(jointFit, newdata = first,
                process = "event",
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Spred
plot(Spred)

# longitudinal predictions using the first measurement
Lpred <- predict(jointFit, newdata = first,
                times = seq(0, 10, length.out = 51),
                return_newdata = TRUE)

Lpred
plot(Lpred)

# combine in one plot
plot(Lpred, Spred)


# dynamic predictions, each time an extra measurement
n <- nrow(dataP155)
for (i in seq_len(n)) {
    data_i <- dataP155[1:i, ]
    data_i$Time <- max(data_i$time)
    data_i$death <- 0
    Spred <- predict(jointFit, newdata = data_i,
                    process = "event",
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)
    Lpred <- predict(jointFit, newdata = data_i,
                    times = seq(0, 10, length.out = 51),
                    return_newdata = TRUE)

    plot(Lpred, Spred)
}


# ROC at 3 years with 1 year window using model-based censoring weights
roc <- tvROC(jointFit, newdata = prothro, Tstart = 3, Dt = 1)
roc
plot(roc)
# AUC at 3 years with 1 year window
tvAUC(roc)

# ROC at 3 years with 1 year window using IPCW
roc2 <- tvROC(jointFit, newdata = prothro, Tstart = 3, Dt = 1,
              type_weights = "IPCW")
roc2
plot(roc2)
# AUC at 3 years with 1 year window
tvAUC(roc2)


# Brier score at 3 years with 1 year window using model-based censoring weights
tvBrier(jointFit, newdata = prothro, Tstart = 3, Dt = 1)

# Brier score at 3 years with 1 year window using IPCW
tvBrier(jointFit, newdata = prothro, Tstart = 3, Dt = 1,
        type_weights = "IPCW")


###############
# Practical 2 #
###############

##################################
### Competing risk joint model ###
##################################

# Have a look at individuals 1,2 , and 5
pbc2.id[pbc2.id$id %in% c(1, 2, 5), c("id", "years", "status")]

# Transform data to a competing risk form
pbc2.idCR <- crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive",
                         nameStrata = "CR")

pbc2.idCR[pbc2.idCR$id %in% c(1, 2, 5),
          c("id", "years", "status", "status2", "CR")]

# Fit the competing risk cox model
CoxFit_CR <- coxph(Surv(years, status2) ~ (age + drug) * strata(CR),
                   data = pbc2.idCR)
summary(CoxFit_CR)

# Fit the linear mixed models for log(serBilir) and prothrombin
fm1 <- lme(log(serBilir) ~ poly(year, 2) * drug, data = pbc2,
           random = ~ poly(year, 2) | id)
summary(fm1)

fm2 <- lme(prothrombin ~ year * drug, data = pbc2, random = ~ year | id)
summary(fm2)

# Functional forms, interaction of current marker value with event cause
CR_forms <- list(
    "log(serBilir)" = ~ value(log(serBilir)):CR,
    "prothrombin" = ~ value(prothrombin):CR
)

# Fit the competing risk joint model
jFit_CR <- jm(CoxFit_CR, list(fm1, fm2), time_var = "year",
              functional_forms = CR_forms,
              n_iter = 5000L, n_burnin = 1000L, n_thin = 5L)

summary(jFit_CR)

#############################
### The data of patient 2 ###
#############################

# Longitudinal dataset
ND_long <- pbc2[pbc2$id == 2, ]
ND_long <- ND_long[1,]
ND_long$years <- 0.2
ND_long$status2 <- 0

# Competing risk dataset
ND_event <- pbc2.idCR[pbc2.idCR$id == 2, ]
ND_event$status2 <- 0
ND_event$years <- 0.2

# Combine in a list()
ND <- list(newdataL = ND_long, newdataE = ND_event)

# Predictions about longitudinal outcomes
predLong <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                    times = seq(2, 10, length = 25))

# Predictions about cumulative risks
predEvent <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                     process = "event",times = seq(2, 10, length = 25))

# Combine results in a single plot
plot(predLong, predEvent, outcomes = 1:2)
legend(x = 2, y = 0.89, legend = levels(pbc2.idCR$CR),
       lty = 1, lwd = 2, col = c("#03BF3D", "#FF0000"), bty = "n", cex = 1)

# Plot for future longitudinal outcomes
par(mfrow = c(1, 2))
plot(predLong, outcomes = 1)
plot(predLong, outcomes = 2)

# Censor data at the following times
censTimes <- c(0.2, 0.5, 1, 5, 8)

for (i in seq_along(censTimes)) {
    # Longitudinal dataset
    ND_long <- pbc2[pbc2$id == 2 & pbc2$year <= censTimes[i], ]
    ND_long$years <- censTimes[i]

    # Competing risk dataset
    ND_event <- pbc2.idCR[pbc2.idCR$id == 2, ]
    ND_event$status2 <- 0
    ND_event$years <- censTimes[i]
    # Combine the two datasets in a list.
    # prerequisite from the predict() method for competing risk joint models
    ND <- list(newdataL = ND_long, newdataE = ND_event)

    predLong <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                        times = seq(censTimes[i], 10, length = 25))

    predEvent <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                         process = "event",
                         times = seq(censTimes[i], 10, length = 25))

    # One figure
    par(mfrow = c(1,1))

    # Combined plot: Evolution of markers + Cumulative risks for the competing events
    plot(predLong, predEvent, outcomes = 1:2, ylim_long_outcome_range = T,
         col_line_event = c("#03BF3D", "#FF0000"),
         fill_CI_event = c("#03BF3D4D", "#FF00004D"))
    legend(x = 2, y = 0.89, legend = levels(pbc2.idCR$CR),
           lty = 1, lwd = 2, col = c("#03BF3D", "#FF0000"), bty = "n", cex = 1)

    # Predictions about future values of longitudinal markers
    par(mfrow = c(1, 2))
    plot(predLong, outcomes = 1, ylim_long_outcome_range = TRUE)
    plot(predLong, outcomes = 2, ylim_long_outcome_range = TRUE)
}

