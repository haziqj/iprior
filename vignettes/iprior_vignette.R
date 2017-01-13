## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 6.5, fig.height = 4.33, 
                      fig.align = "center")
library(iprior)

## ----syntax, eval=FALSE--------------------------------------------------
#  mod <- iprior(y ~ x, data = dat)  # formula based input
#  mod <- iprior(y = dat$y, x = dat$x)  # non-formula based input

## ----data1---------------------------------------------------------------
str(stackloss)

## ----mod1----------------------------------------------------------------
mod.iprior <- iprior(stack.loss ~ ., data = stackloss)

