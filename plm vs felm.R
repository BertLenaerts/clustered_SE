
library(plm)
library(lfe) 
library(lmtest)
library(dplyr)
library(readr)
source("vcov_plm.R")

######################################################

### COMPPARISON PLM VS FELM

d = read_csv("d.csv") %>%
  mutate(idyear = as.numeric(paste0(id, year)))

# n = 41
# 9 id groups
# 5 year groups
# 3 cluster groups (gid)
# 2 x-vars (x & u)

## GROUP-FIXED EFFECTS

plm2a=plm(y~x +u,
          data=d,
          model="within", 
          effect="individual", 
          index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm2a) 
summary(felm(y ~ u + x| id | 0 | 0, d))

# HC CORRECTION
coeftest(plm2a, vcov=(41/(41-9-2))*vcovHC(plm2a, type="HC0", method = "white1"))
coeftest(plm2a, vcovHR(plm2a))

# GROUP CLUSTER
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=d$id))
summary(felm(y ~ u + x | id | 0 | id, d))

vcovHC(plm2a, type="HC0", cluster="group", method="arellano") 
# clustering can also be done using vcovHC but without Stata-like df adjustment
# default: cluster = group & method = arellano
vcovG(plm2a, type = "HC0", l = 0, inner = "cluster", cluster="group")
# this is the workhorse behind vcovHC

# TIME CLUSTER
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=d$year))
summary(felm(y ~ u + x| id | 0 | year, d))

# HIGHER-LEVEL CLUSTERING
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=d$gid))
summary(felm(y ~ u + x| id | 0 | gid, d))

# TWOWAY CLUSTERING
coeftest(plm2a, vcov=vcovDC(plm2a, type="HC0"))
coeftest(plm2a, vcov=vcovTC(x=plm2a, d$year, d$id))
summary(felm(y ~ u + x| id | 0 | (id+year), d))
# felm uses the intersection of both the group and time cluster to generate
# a third cluster variable; if the id-year combos are unique then this is equivalent
# to White's heteroskedasticity-robust correction
# since felm allows non-unique year-id combos (and also allows for multi-way clustering), 
# the third variable method is implemented, whereas plm uses the White method
# the df correction for the substraction term takes the number of clusters for the
# third variable for felm whereas plm uses the standard White HC df correction
vcovCL(x=plm2a, cluster=d$idyear)
vcovHC(plm2a, type="HC0", method = "white1")


## TIME-FIXED EFFECTS

plm2b=plm(y~x +u,
          data=d,
          model="within", 
          effect="time", 
          index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm2b)
summary(felm(y ~ u + x| year | 0 | 0, d))

# HC CORRECTION
coeftest(plm2b, vcov=(41/(41-5-2))*vcovHC(plm2b, type="HC0", method = "white1"))
coeftest(plm2b, vcovHR(plm2b))

# GROUP CLUSTER
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$id))
summary(felm(y ~ u + x| year | 0 | id, d))

# TIME CLUSTER
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$year))
summary(felm(y ~ u + x| year | 0 | year, d))

# HIGHER-LEVEL CLUSTERING
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$gid))
summary(felm(y ~ u + x| year | 0 | gid, d))

# TWOWAY CLUSTERING
coeftest(plm2b, vcov=vcovDC(plm2b, type="HC0"))
coeftest(plm2b, vcov=vcovTC(x=plm2b, d$year, d$id))
summary(felm(y ~ u + x| year | 0 | (id+year), d))


## TWOWAY FIXED EFFECTS

plm3=plm(y~x +u,
         data=d,
         model="within", 
         effect="twoway", 
         index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm3) 
summary(felm(y ~ u + x| id+year | 0 | 0, d))

# HC CORRECTION
coeftest(plm3, vcov=(41/(41-14+1-2))*vcovHC(plm3, type="HC0", method = "white1"))
coeftest(plm3, vcovHR(plm3))

# GROUP CLUSTER
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$id))
summary(felm(y ~ u + x| id+year | 0 | id, d))

# TIME CLUSTER
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$year))
summary(felm(y ~ u + x| id+year | 0 | year, d))

# HIGHER-LEVEL CLUSTERING
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$gid))
summary(felm(y ~ u + x| id+year | 0 | gid, d))

# TWOWAY CLUSTERING
coeftest(plm3, vcov=vcovDC(plm3, type="HC0"))
coeftest(plm3, vcov=vcovTC(x=plm3, d$year, d$id))
summary(felm(y ~ u + x| id+year | 0 | (id+year), d))
