# Example usage for I-prior modelling with interactions

>**WARNING: `iprior()` function may appear slow or frozen with large datasets. The following examples have been tested and found to be working, for `n=600` roughly.**

There are two ways we can model interactions using I-priors. The first, is the parsimonious method. Here, the scale parameter for the interaction term between `x1` and `x2` is the product between the two scale parameters, i.e. `lambda1*lambda2`. For variables `x1, ..., xp` with all two-way interactions, the number of scale parameters is still `p`.

The second method is to assign another scale parameter to for the interaction terms. Thus, the scale for `x1:x2` is given by `lambda12`, say. The maximum number of parameters is `p + p(p-1)/2`. There are more parameters to be estimated in this method.

The `iprior()` function can now handle both of these methods. The formula call is exactly the same as for `lm()`:
```r
iprior(Hwt ~ Bwt + Sex + Bwt:Sex, data=cats)  #This model includes the interaction between Bwt and Sex
iprior(Hwt ~ (Bwt + Sex)^2, data=cats)        
iprior(Hwt ~ .^2, data=cats)  
```
These are all equivalent calls to the same model. The formula call `^2` models all two-way interactions between the terms. It is useful to use `.^2` as a shorthand, or when it is not known how many variables the dataset contains. The parsimonious interactions method is called by default, but if one wishes to do the second non-parsimonious method, one simply needs to add the option `parsm=F`.

## Example 1
