# TMMP (Transformed mixed model prediction of general finite population parameters)
This package implements transformed mixed model prediction of general finite population parameters, as proposed by the following papers.

Sugasawa, S. and Kubokawa, T.  (2019). Adaptively transformed mixed model prediction of general finite population parameters. *Scandinavian Journal of Statistics*,  to appear
https://doi.org/10.1111/sjos.12380

Functions are implemented in `TMMP-function.R` available in the repository.
```{r}
source("TMMP-function.R")   # require "sae" package
```

Load demo dataset
```{r}
load("data.RData")
m=length(ni)   # number of areas
z=0.6*median(sY)   # poverty level
```

Fit the transformed mixed model with dual power transformation

Input of `TMMP`

- `y`: response vector
- `X`: matrix of sampled covariates (the first column should be 1's)
- `rX`: matrix of non-sampled covariates (the first column should be 1's)
- `m`: number of areas
- `Ni`: vector of the numbers of population in each area
- `ni`: vector of the numbers of samples in each area
- `z`: poverty line
- `la.set`: set of lambda (the default is 0, 0.05,..., 0.95, 1)
- `mc`: number of Monte Carlo samples (the default is 1000)

Output of `TMMP`

- `Pred`: estimated area-wise poverty measure
- `Estimates`: parameter estimates 
- `RE`: Estimates of random effects 
- `MC`: All Monte calro samples  
```{r}
fit=TMMP(y=sY,X=sX,rX=rX,m=m,Ni=Ni,ni=ni,z=z)
fit$Pred   # estimated poverty rate
fit$Estimates   # parameter estimates
fit$RE   # estimated random effects
```

Compute empirical Bayes confidence intervals with dual power transformation

Input of `TMMP.CI`  (Additionally to the input in `TMMP`)

- `alpha`: significance level (the default is 0.05)
- `B`: number of boostrap replications (the default is 100)

Output of `TMMP.CI`

- area-wise confidence intervals
```{r}
TMMP.CI(y=sY,X=sX,rX=rX,m=m,Ni=Ni,ni=ni,z=z,mc=100,B=100,print=T)
```



