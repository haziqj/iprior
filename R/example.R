## classic Fisher cats data from package MASS
data(cats, package="MASS")

###
### Comparison of I-prior and OLS
### 
mod.iprior <- iprior(Hwt~Bwt, data=cats)
mod.lm <- lm(Hwt~Bwt, data=cats)

### Compare value of sigma
sigma.iprior <-  1/sqrt(mod.iprior$psi)
sigma.lm <- summary(mod.lm)$sigma
sigma.iprior; sigma.lm

### classical simple regression model
fit.iprior <- fitted(mod.iprior)
fit.lm <- fitted(mod.lm)

### comparing fitted values
plot(fit.lm, fit.iprior, type="n", xlab="Classical regression model estimates", ylab="I-prior estimates", main="Comparison between I-prior and classical regression predicted values")
text(fit.lm, fit.iprior, pch=as.character(1:length(fit.lm)), col=1:length(fit.lm), cex=0.7)
abline(a=0, b=1)
