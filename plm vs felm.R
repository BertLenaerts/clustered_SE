
library(plm)
library(lfe) 
library(lmtest)
library(dplyr)
library(readr)
source("vcovCL.R")

######################################################

### COMPPARISON PLM VS FELM

d = read_csv("d.csv") %>%
  mutate(idyear = as.numeric(paste0(id, year)))

# n = 41
# 9 id groups
# 5 year groups
# 3 cluster groups (gid)
# 2 x-vars (x & u)
# 1 end. var (x)
# 1 IV-vars (e)

## GROUP-FIXED EFFECTS

plm2a=plm(y~x +u | u+e,
          data=d,
          model="within", 
          effect="individual", 
          index=c("id", "year"))

# GROUP CLUSTER
coeftest(plm2a, vcov=(40/(41-9-2))*(9/8)*vcovHC(plm2a, type="HC0", cluster="group"))
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=d$id, stata = T))
summary(felm(y ~ u | id | (x ~ e) | id, d))

vcovHC(plm2a, type="HC0") # default cluster = group (method = arellano)
vcovHC(plm2a, type="HC0", cluster="group")
vcovHC(plm2a, type="HC0", cluster="group", method="arellano")
vcovG(plm2a, type = "HC0", l = 0, inner = "cluster", cluster="group")

# TIME CLUSTER
coeftest(plm2a, vcov=(40/(41-9-2))*(5/4)*vcovHC(plm2a, type="HC0", cluster="time"))
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=d$year, stata = T))
summary(felm(y ~ u | id | (x ~ e) | year, d))

# NO CLUSTERING (no HC & no df correction)
summary(plm2a) # plm2a$vcov
summary(felm(y ~ u | id | (x ~ e) | 0, d))

# HC COREECTION
coeftest(plm2a, vcov=(40/(41-9-2))*vcovG(plm2a, type = "HC0", l = 0, inner = "white"))
coeftest(plm2a, vcov=(40/(41-9-2))*vcovHC(plm2a, type="HC0", method = "white1"))

# HIGHER-LEVEL CLUSTERING
coeftest(plm2a, vcov=(40/(41-9-2))*(3/2)*vcovCL(x=plm2a, cluster=d$gid))
summary(felm(y ~ u | id | (x ~ e) | gid, d))

# TWOWAY CLUSTERING
coeftest(plm2a, vcov=vcovDC(plm2a, type="HC0"))
coeftest(plm2a, vcov=vcovTC(x=plm2a, d$year, d$id))
coeftest(plm2a, vcov=vcovTC(x=plm2a, d$year, d$id, stata=T))
coeftest(plm2a, vcov=( (40/(41-9-2))*(9/8)*vcovCL(x=plm2a, cluster=d$id)+
                         (40/(41-9-2))*(5/4)* vcovCL(x=plm2a, cluster=d$year)-
                         (40/(41-9-2))*(41/40)*vcovHC(plm2a, type="HC0", method = "white1") ) )
summary(felm(y ~ u | id | (x ~ e) | (id+year), d))
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

plm2b=plm(y~x +u | u+e,
          data=d,
          model="within", 
          effect="time", 
          index=c("id", "year"))

# GROUP CLUSTER
coeftest(plm2b, vcov=(40/(41-5-2))*(9/8)*vcovHC(plm2b, type="HC0", cluster="group"))
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$id, stata=T))
summary(felm(y ~ u | year | (x ~ e) | id, d))

# TIME CLUSTER
coeftest(plm2b, vcov=(40/(41-5-2))*(5/4)*vcovHC(plm2b, type="HC0", cluster="time"))
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$year, stata=T))
summary(felm(y ~ u | year | (x ~ e) | year, d))

# NO CLUSTER (no HC too)
summary(plm2b)
summary(felm(y ~ u | year | (x ~ e) | 0, d))

# HC COREECTION
coeftest(plm2b, vcov=(40/(41-5-2))*vcovG(plm2b, type = "HC0", l = 0, inner = "white"))
coeftest(plm2b, vcov=(40/(41-5-2))*vcovHC(plm2b, type="HC0", method = "white1"))

# HIGHER-LEVEL CLUSTER
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=d$gid, stata=T))
summary(felm(y ~ u | year | (x ~ e) | gid, d))

# TWOWAY CLUSTERING
coeftest(plm2b, vcov=vcovTC(x=plm2b, d$year, d$id, stata=T))
summary(felm(y ~ u | year | (x ~ e) | (id+year), d))


## TWOWAY FIXED EFFECTS

plm3=plm(y~x +u | u+e,
         data=d,
         model="within", 
         effect="twoway", 
         index=c("id", "year"))

# GROUP CLUSTER
coeftest(plm3, vcov=(40/(41-14+1-2))*(9/8)*vcovHC(plm3, type="HC0", cluster="group"))
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$id, stata=T))
# cannot subtract intercept twice
summary(felm(y ~ u | id+year | (x ~ e) | id, d))

# TIME CLUSTER
coeftest(plm3, vcov=(40/(41-14+1-2))*(5/4)*vcovHC(plm3, type="HC0", cluster="time"))
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$year, stata=T))
# cannot subtract intercept twice
summary(felm(y ~ u | id+year | (x ~ e) | year, d))

# NO CLUSTER (no HC too)
summary(plm3)
summary(felm(y ~ u | id+year | (x ~ e) | 0, d))

# HC COREECTION
coeftest(plm3, vcov=(40/(41-14+1-2))*vcovG(plm3, type = "HC0", l = 0, inner = "white"))
coeftest(plm3, vcov=(40/(41-14+1-2))*vcovHC(plm3, type="HC0", method = "white1"))

# HIGHER-LEVEL CLUSTER
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=d$gid, stata=T))
summary(felm(y ~ u | id+year | (x ~ e) | gid, d))

# TWOWAY CLUSTERING
coeftest(plm3, vcov=vcovDC(plm3, type="HC0"))
coeftest(plm3, vcov=vcovTC(x=plm3, d$year, d$id, stata=F))
coeftest(plm3, vcov=( (40/(41-14+1-2))*(9/8)*vcovCL(x=plm3, cluster=d$id)+
                         (40/(41+1-14-2))*(5/4)* vcovCL(x=plm3, cluster=d$year)-
                         (40/(41-14+1-2))*(41/40)*vcovHC(plm3, type="HC0", method = "white1") ) )
coeftest(plm3, vcov=vcovTC(x=plm3, d$year, d$id, stata=T))
summary(felm(y ~ u | id+year | (x ~ e) | (id+year), d))

vcovCL(x=plm3, cluster=d$idyear)
vcovHC(plm3, type="HC0", method = "white1")
