library(CNVtools)
source('R/CNVtools_class.R')
source('R/model_fitting.R')
dyn.load('src/CNVtools.so')

#library(survival); source("CNVtools.r"); dyn.load("../src/CNVtools.so"); load("../../CNVtools/data/A112.RData")

set.seed(10)
n <- 1000
CNV <- rbinom(size = 2, n = 1000, prob = 0.3)

beta <- 0  ##log odds parameter
case.control <- rbinom(prob = -0.3 + exp(beta*CNV)/(1 +  exp(beta*CNV) ), n = n, size = 1)
                       
signal <- rnorm(n = n, mean = CNV, sd = 0.05)

test <- new('CNVtools', signal = signal, 
            trait = case.control)


test <- fit(test, ncomp = 3, hyp = 'H0', EM.starting.point = 'kmeans')
test <- fit(test, ncomp = 3, hyp = 'H1')
stat.CNVtools <- 2*(test@mle.H1$lnL - test@mle.H0$lnL)


pdf('test.pdf')
plot.cnv(test, hyp = 'H0', col = 'red')
plot.cnv(test, hyp = 'H1', col = 'red')
dev.off()

#################

source('R/CNVtools_class.R')
source('R/model_fitting.R')
dyn.load('src/CNVtools.so')

n <- 1000
beta <- log(2)
association.test.strata <- rbinom( n = n, size = 1, prob = 0.6)
CNV <- rbinom(size = 2, n = 1000, prob = 0.3)
case.control <- rbinom(prob = -0.15 - association.test.strata*0.15 + exp(beta*CNV)/(1 +  exp(beta*CNV) ), n = n, size = 1)
signal <- rnorm(n = n, mean = CNV, sd = 0.05)
my.covar <- data.frame(cc = case.control, region = factor(association.test.strata), cn = CNV)

test.with.strata <- new('CNVtools', signal = signal, 
                        trait = case.control,
                        model.association.H0 = formula(' ~ factor(region)'),
                        model.association.H1 = formula(' ~ factor(region) + cn'),
                        covariates = my.covar)

test.with.strata <- fit(test.with.strata, ncomp = 3, hyp = 'H0', EM.starting.point = 'kmeans')
test.with.strata <- fit(test.with.strata, ncomp = 3, hyp = 'H1', EM.starting.point = 'H1.from.H0')
stat.CNVtools.strata <- 2*(test.with.strata@mle.H1$lnL - test.with.strata@mle.H0$lnL)
message(stat.CNVtools.strata)



test.without.strata <- new('CNVtools',
                           signal = signal, 
                           trait = case.control,
                           model.association.H0 = formula(' ~ 1 '),
                           model.association.H1 = formula(' ~ cn '))
                           
test.without.strata <- fit(test.without.strata, ncomp = 3, hyp = 'H0', EM.starting.point = 'kmeans')
test.without.strata <- fit(test.without.strata, ncomp = 3, hyp = 'H1')
stat.CNVtools.no.strata <- 2*(test.without.strata@mle.H1$lnL - test.without.strata@mle.H0$lnL)

mod0 <- glm(data = my.covar, family = binomial, formula = 'cc ~ factor(region) ')
mod1 <- glm(data = my.covar, family = binomial, formula = 'cc ~ factor(region) + CNV')
stat.glm <- diff(anova(mod1, mod0)[[2]])

message(stat.CNVtools.no.strata, ' ', stat.CNVtools.strata, ' ', stat.glm)
