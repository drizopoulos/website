##########################################################################################
# Title: Biostatistics Methods II: Classical Regression Models -                         #
#        Solutions for the Software Practicals on Likelihood & Survival Analysis         #
#                                                                                        #
# Required R packages: survival, splines, lattice, JM                                    #
#                                                                                        #
# Author: Dimitris Rizopoulos                                                            #
##########################################################################################

##########################
# Load Required Packages #
##########################

library("survival")
library("splines")
library("lattice")
library("JM")

lung$sex <- factor(lung$sex, levels = 1:2, labels = c("male", "female"))
lung <- with(lung, lung[complete.cases(time, status, sex, age, ph.karno), ])

###########################################################
# Likelihood - Practical 1: Linear Regression Models in R #
###########################################################

# Q1
fit1 <- lm(log(serBilir) ~  drug + age + sex + serChol + prothrombin, 
           data = pbc2.id)

summary(fit1)

# Q2
par(mfrow = c(2, 2))
plot(fit1)

# Q3
fit2 <- lm(log(serBilir) ~  drug + age + sex + serChol + prothrombin + 
               sex:age + sex:serChol + sex:prothrombin, 
           data = pbc2.id)

fit2 <- lm(log(serBilir) ~  drug + sex * (age + serChol + prothrombin), 
           data = pbc2.id)

summary(fit2)

# Q4
anova(fit1, fit2)

# Q5
summ_fit2 <- summary(fit2)

coef(summ_fit2)[7:9, 4]

p.adjust(coef(summ_fit2)[7:9, 4])

# Q6
fit1 <- lm(log(serBilir) ~ poly(serChol, 3) + age + sex, 
           data = pbc2.id[complete.cases(pbc2.id$serChol), ])

summary(fit1)

# Q7
fit2 <- lm(log(serBilir) ~ serChol + age + sex, 
           data = pbc2.id[complete.cases(pbc2.id$serChol), ])

anova(fit2, fit1)

##########################################################################################
##########################################################################################

###########################################################################
# Survival - Practical 1: Features Survival Data, Basic Functions & Tests #
###########################################################################

# Q1
fitKM <- survfit(Surv(Time, death) ~ drug, data = aids.id)

fitKM

plot(fitKM, lty = 1:2, col = 1:2)

# Q2
fitB <- survfit(Surv(Time, death) ~ drug, data = aids.id, type = "fleming-harrington")

fitB

plot(fitB, lty = 1:2, col = 1:2)

quantile(fitB, 1 - c(0.5, 0.6, 0.7))

# Q3
summary(fitB, times = c(8, 10))

# Q4
survdiff(Surv(Time, death) ~ drug, data = aids.id)

survdiff(Surv(Time, death) ~ drug, data = aids.id, rho = 1)

# Q5
fitKM_gender <- survfit(Surv(Time, death) ~ gender, data = aids.id)

fitKM_gender

plot(fitKM_gender)

survdiff(Surv(Time, death) ~ gender, data = aids.id)

survdiff(Surv(Time, death) ~ gender, data = aids.id, rho = 1)

##########################################################################################
##########################################################################################

###########################################################
# Survival - Practical 2: Accelerated Failure Time Models #
###########################################################

# Q1
fit1 <- survreg(Surv(time, status == 2) ~ sex * (age + ph.karno), data = lung,
                dist = "weibull")

summary(fit1)

# Q2
fit2 <- survreg(Surv(time, status == 2) ~ age + ph.karno, data = lung)

anova(fit2, fit1)

fit3 <- survreg(Surv(time, status == 2) ~ sex + age + ph.karno, data = lung)

anova(fit3, fit1)

# Q3
with(lung, median(age))

with(lung, range(ph.karno))

ND <- expand.grid(sex = c("male", "female"), age = 63, 
                  ph.karno = seq(50, 100, length.out = 50))

prs <- predict(fit3, ND, se.fit = TRUE, type = "lp")
ND$pred <- prs[[1]]
ND$se <- prs[[2]]
ND$lo <- exp(ND$pred - 1.96 * ND$se)
ND$up <- exp(ND$pred + 1.96 * ND$se)
ND$pred <- exp(ND$pred)

xyplot(pred + lo + up ~ ph.karno | sex, data = ND, type = "l", 
       lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 3, 
       xlab = "Karnofsky Performance Score", ylab = "Survival Time")

# Q4
fits <- fit3$linear.predictors

resids <- (fit3$y[, 1] - fits) / fit3$scale

resKM <- survfit(Surv(resids, status) ~ 1, data = lung)

plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", 
     ylab = "Survival Probability", main = "AIDS Data Set")
xx <- seq(min(resids), max(resids), length.out = 35)
yy <- exp(- exp(xx))
lines(xx, yy, col = "red", lwd = 2)

##########################################################################################
##########################################################################################

###########################################################
# Survival - Practical 3: Cox Proportional Hazards Models #
###########################################################

# Q1
fit1 <- coxph(Surv(Time, death) ~ ns(CD4, 3) * (drug + AZT), data = aids.id)

summary(fit1)

# Q2
fit2 <- coxph(Surv(Time, death) ~ ns(CD4, 3) + drug + AZT, data = aids.id)

anova(fit2, fit1)

# Q3
summary(fit2)

# Q4
ND <- expand.grid(CD4 = seq(0, 10, length.out = 20), drug = c("ddC", "ddI"), 
                  AZT = c("intolerance", "failure"))

prs <- predict(fit1, newdata = ND, type = "lp", se.fit = TRUE)
ND$pred <- prs[[1]]
ND$se <- prs[[2]]
ND$lo <- ND$pred - 1.96 * ND$se
ND$up <- ND$pred + 1.96 * ND$se

xyplot(pred + lo + up ~ CD4 | AZT * drug, data = ND, 
       col = c("red", "black", "black"), lty = c(1, 2, 2), lwd = 2, type = "l",
       abline = list(h = 0, lty = 2), xlab = "CD4", ylab = "log Hazard Ratio")

# Q5
fit <- survfit(Surv(Time, death) ~ AZT, data = aids.id)

par(mfrow = c(1, 2))

plot(fit, xlab = "Months", ylab = "Survival", col = 1:2)

plot(fit, fun = function (s) -log(-log(s)), xlab = "Months", 
     ylab = "-log(- log(Survival))", col = 1:2)

