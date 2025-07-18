###########################################################################
## RDHTE R Package
## Code for Empirical Illustration
## Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Filippo Palomba and Rocio Titiunik 
###########################################################################

### Install R library
### NOTE: depending on your system, you may need to do it as root
#install.packages('rdhte')
#install.packages('rdrobust')


### Clear R environment
rm(list=ls(all=TRUE))


##########################################################################
# Set up packages, etc
##########################################################################
library("rdhte")
library("rdrobust")


##########################################################################
# Load RDHTE illustration data (from Granzier, Pons and Tricaud (2023, AEJ)
##########################################################################
attach(rd.data <- read.csv("rdhte_dataset.csv", stringsAsFactors = TRUE))



##########################################################################
# Single binary variable
##########################################################################
summary(rd_left <- rdhte(y = y, x = x, covs.hte = factor(w_left), cluster = cluster_var))


### Post-estimation illustration -- rdhte_lincom
# How large is the advantage for left-of-center candidates?
rdhte_lincom(model=rd_left, linfct = c("`factor(w_left)1` - `factor(w_left)0` = 0"))
#Note: must specify the hypothesis as A-B=0 and not as A=B. The command below will fail
#rdhte_lincom(model=rd_left, linfct = c("`factor(w_left)1` = `factor(w_left)0`"))

# Extra illustration:
# Note that forcing a common bandwidth makes the effect for `center or right' 
# (incorrectly) statistically significant
summary(rdhte(y = y, x = x, covs.hte = w_left, cluster = cluster_var, bw.joint = TRUE))

# Extra illustration:
#rdhte will automatically treat a variable coded 0/1 as a factor.
summary(rd_left <- rdhte(y = y, x = x, covs.hte = w_left, cluster = cluster_var))
rdhte_lincom(model=rd_left, linfct = c("`w_left1` - `w_left0` = 0"))

##########################################################################
# Single categorical variable -- unordered
##########################################################################
summary(rd_ideology <- rdhte(y = y, x = x, covs.hte = factor(w_ideology), cluster = cluster_var))


### Post-estimation illustration
# All the non-left groups are statistically 
# indistinguishable from each other and from zero
rdhte_lincom(model=rd_ideology, linfct = c("`factor(w_ideology)4` - `factor(w_ideology)3` = 0", "`factor(w_ideology)4` - `factor(w_ideology)2` = 0", "`factor(w_ideology)4` = 0"))
rdhte_lincom(model=rd_ideology, linfct = c("`factor(w_ideology)2` = 0", "`factor(w_ideology)3` = 0", "`factor(w_ideology)4` = 0"))



##########################################################################
# Single categorical variable -- ordered
##########################################################################
summary(rdhte(y = y, x = x, covs.hte = factor(w_strength_qrt), cluster = cluster_var))

# -> result: advantage increases with strength 
#candidates with higher average strength at the national level have higher effect, monotonic!


# Extra illustration:
#Will be treated as continuous by default
summary(rdhte(y = y, x = x, covs.hte = w_strength_qrt, cluster = cluster_var))



##########################################################################
# Two binary variables - interaction
##########################################################################

summary(rdhte(y = y, x = x, covs.hte = factor(w_left):factor(w_strong), cluster = cluster_var))

# Extra illustration:
# Above syntax is the same as creating a new factor variable flagging each combination
interactions <- 1*(w_left==0)*(w_strong==1) + 2*(w_left==0)*(w_strong==2) + 3*(w_left==1)*(w_strong==1) + 4*(w_left==1)*(w_strong==2)
summary(rdhte(y = y, x = x, covs.hte = factor(interactions), cluster = cluster_var))


##########################################################################
# Single continuous variable
##########################################################################

summary(rd_continuous <- rdhte(y = y, x = x, covs.hte = w_strength, kernel="uni", cluster = cluster_var))

# to aid with interpretation, notice that the coefficient on T#W
# is _precisely_ a slope coefficient in a linear model, and can be interpreted the same.
# To see this, use the uniform kernel and a set bandwidth to fit local unweighted least squares
# then match this with the base command -lm()- 

trt <- (x>0)
new.data <- data.frame(y,x,w_strength,trt)
using.lm <- coef(lm(y ~ trt*x*w_strength, data=new.data[abs(new.data$x)<rd_continuous$h[1],]))
rd_continuous$Estimate[1] - using.lm["trtTRUE"]
rd_continuous$Estimate[2] - using.lm["trtTRUE:w_strength"]

# IMPORTANT: inference requires robust bias correction and cannot be obtained from this regression


##########################################################################
# Interaction effect: binary#continuous
##########################################################################

#Full interaction is most natural (different intercept and slope for each category)
summary(rd_interaction <- rdhte(y = y, x = x, covs.hte = "w_left*w_strength", cluster = cluster_var))

# Extra illustration:
# Using the factor() syntax is available, but will lead to different baseline categories for the factor 
# variable, leading to different intercepts.
summary(rdhte(y = y, x = x, covs.hte = "factor(w_left)*w_strength", cluster = cluster_var))

# Extra illustration:
# Each effect is insignificant, but the joint test shows there is information
rdhte_lincom(model=rd_interaction, linfct = c("`T` = 0", "`T:w_left` = 0", "`T:w_strength` = 0", "`T:w_left.w_strength` = 0"))

# Extra illustration:
# to aid interpretation, the fully interacted model will match results from 
# category-specific estimation. Fix the bandwidth to make results match exactly.
summary(rd_interaction <- rdhte(y = y, x = x, covs.hte = "w_left*w_strength", h=0.1, cluster = cluster_var))
summary(rdhte(y = y[w_left==0], x = x[w_left==0], covs.hte = w_strength[w_left==0], h=0.1, cluster = cluster_var[w_left==0]))
summary(rdhte(y = y[w_left==1], x = x[w_left==1], covs.hte = w_strength[w_left==1], h=0.1, cluster = cluster_var[w_left==1]))
rdhte_lincom(model=rd_interaction, linfct = c("`T` + `T:w_left` = 0", "`T:w_strength` + `T:w_left.w_strength` = 0"))



##########################################################################
# Replicating rdrobust
##########################################################################

### Average effects

#Using default settings, packages will not match, because of different settings
summary(rdhte(y = y, x = x))
summary(rdrobust(y = y, x = x))

# to replicate exactly with RDROBUST, set rho=1 and specific vce() option
# bandwidth selection is also different, so enforce the same bandwidth
summary(rdhte(y = y, x = x, h=0.1, vce="hc3"))
summary(rdrobust(y = y, x = x, h=0.1, rho=1, vce="hc3"))


### subgroup analysis

# keeping the settings the same as above, rdhte will match rdrobust if the latter
# is run for each subgroup separately
summary(rdhte(y = y, x = x, covs.hte=w_left, h=c(0.078,0.116), vce="hc3"))
summary(rdrobust(y = y[w_left==1], x = x[w_left==1], h=0.116, rho=1, vce="hc3"))

