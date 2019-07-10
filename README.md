# CWGEE
Cluster-weighted generalized estimating equations for clustered longitudinal data with informative cluster size

#### R installation Instructions
```
install.packages("devtools")
library(devtools)
devtools::install_github("AyaMitani/CWGEE")
library(CWGEE)
```
#### Example for ordinal clustered longitudinal outcome with informative cluster size
```
data(perio)
fitmod <- ordCWGEE(formula = cal ~ mets + edu + age + smoking, data = perio,
id = subject, cluster.var = tooth, time.var = visit, time.str = "ind")
summary(fitmod)
```

#### Example for clustered multivariate binary outcomes with informative cluster size
```
data(perio_base)
fitmod <- mvoCWGEE(formula = y ~ smoking + age + edu, data = perio_base,
cluster = subject, resp.ind = outcome, unit = tooth,
common.slope = c("smoking", "edu"), corr.str = "exch")
summary(fitmod)
```