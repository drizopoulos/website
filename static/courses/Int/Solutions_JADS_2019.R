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


