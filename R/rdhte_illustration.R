################################################################################
## RDHTE Package
## Empirical Illustration
################################################################################

##########################################################################
## Data and design
##########################################################################
## This script illustrates `rdhte` on a real-data extract from Granzier,
## Pons, and Tricaud (2023, AEJ: Applied), "Coordination and Bandwagon
## Effects: How Past Rankings Shape the Behavior of Voters and
## Candidates." The authors study French two-round elections, where
## candidates must clear a qualifying-vote threshold in the first round
## to advance to the runoff. The institutional rule creates a sharp RD
## design on every candidate's first-round margin against that threshold:
## candidates just above the cutoff advance, those just below do not.
##
## The bundled extract `rdhte_dataset.csv` has 39,534 candidate-race
## observations with the following variables:
##
##   y                outcome: 1 if the candidate advances to the runoff
##   x                running variable: first-round margin against the
##                    qualifying threshold (cutoff at zero)
##   cluster_var      district identifier (cluster-robust inference)
##   w_left           1 if the candidate's party is left of center
##   w_ideology       unordered party-ideology bucket (4 levels)
##   w_strength       continuous proxy for ex-ante candidate strength
##   w_strong         1 if above-median strength
##   w_strength_qrt   ordered quartile bucket of w_strength

library(rdhte)
library(rdrobust)

rd.data <- read.csv("rdhte_dataset.csv")


##########################################################################
## Single binary variable
##########################################################################

rd_left <- rdhte(y = y, x = x, covs.hte = w_left,
                 cluster = cluster_var, data = rd.data)
summary(rd_left)

rdhte_lincom(rd_left, "covs.hte1 - covs.hte0 = 0")

rdhte(y = y, x = x, covs.hte = w_left,
      cluster = cluster_var, bw.joint = TRUE, data = rd.data)

rdhte(y = y, x = x, covs.hte = factor(w_left),
      cluster = cluster_var, data = rd.data)


##########################################################################
## Single categorical variable -- unordered
##########################################################################

rd_ideology <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                     cluster = cluster_var, data = rd.data)
summary(rd_ideology)

rdhte_lincom(rd_ideology,
             c("covs.hte4 - covs.hte3 = 0",
               "covs.hte4 - covs.hte2 = 0",
               "covs.hte4 = 0"))


##########################################################################
## Single categorical variable -- ordered
##########################################################################

rdhte(y = y, x = x, covs.hte = factor(w_strength_qrt),
      cluster = cluster_var, data = rd.data)

rdhte(y = y, x = x, covs.hte = w_strength_qrt,
      cluster = cluster_var, data = rd.data)


##########################################################################
## Two binary variables -- interaction
##########################################################################

rdhte(y = y, x = x, covs.hte = factor(w_left):factor(w_strong),
      cluster = cluster_var, data = rd.data)


##########################################################################
## Single continuous variable
##########################################################################

rd_strength <- rdhte(y = y, x = x, covs.hte = w_strength,
                     kernel = "uni", cluster = cluster_var,
                     data = rd.data)
summary(rd_strength)

bw <- rd_strength$h[1, 1]
lm(y ~ I(x > 0) * x * w_strength,
   data = subset(rd.data, abs(x) <= bw))


##########################################################################
## Interaction effect: binary x continuous
##########################################################################

rd_interaction <- rdhte(y = y, x = x,
                        covs.hte = "w_left*w_strength",
                        cluster = cluster_var, data = rd.data)
summary(rd_interaction)

rdhte_lincom(rd_interaction,
             c("T = 0",
               "T:w_left = 0",
               "T:w_strength = 0",
               "T:w_left:w_strength = 0"))

rd_interaction <- rdhte(y = y, x = x,
                        covs.hte = "w_left*w_strength",
                        h = 0.1, cluster = cluster_var, data = rd.data)
summary(rd_interaction)

rdhte(y = y, x = x, covs.hte = w_strength, h = 0.1,
      cluster = cluster_var, subset = w_left == 1, data = rd.data)

rdhte(y = y, x = x, covs.hte = w_strength, h = 0.1,
      cluster = cluster_var, subset = w_left == 0, data = rd.data)


##########################################################################
## Standalone bandwidth selection (rdbwhte)
##########################################################################

rd_bw <- rdbwhte(y = y, x = x, covs.hte = factor(w_ideology),
                 cluster = cluster_var, data = rd.data)
summary(rd_bw)
rd_bw$h

rdbwhte(y = y, x = x, covs.hte = factor(w_ideology),
        cluster = cluster_var, bw.joint = TRUE, data = rd.data)


##########################################################################
## Efficiency-improving covariates (covs.eff)
##########################################################################

rd_eff_off <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                    cluster = cluster_var, data = rd.data)
rd_eff_on <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                   covs.eff = w_strength, cluster = cluster_var,
                   data = rd.data)

cbind(no_covs_eff = rd_eff_off$se.rb, covs_eff = rd_eff_on$se.rb)


##########################################################################
## Plotting
##########################################################################

plot(rd_ideology)
plot(rd_ideology, sort = TRUE)
plot(rd_ideology, sort = TRUE,
     title = "Heterogeneity by ideology bucket",
     ylab = "Sharp RD ITT")


##########################################################################
## Building publication-ready tables
##########################################################################

broom::tidy(rd_ideology)
broom::glance(rd_ideology)

rd_hc1 <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                vce = "hc1", data = rd.data)
rd_hc2 <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                vce = "hc2", data = rd.data)
rd_hc3 <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                vce = "hc3", data = rd.data)
rd_cr1 <- rdhte(y = y, x = x, covs.hte = factor(w_ideology),
                vce = "cr1", cluster = cluster_var, data = rd.data)

cbind(HC1 = rd_hc1$Estimate, HC2 = rd_hc2$Estimate,
      HC3 = rd_hc3$Estimate, CR1 = rd_cr1$Estimate)
cbind(HC1 = rd_hc1$se.rb, HC2 = rd_hc2$se.rb,
      HC3 = rd_hc3$se.rb, CR1 = rd_cr1$se.rb)


##########################################################################
## Replicating rdrobust
##########################################################################

rdhte(y = y, x = x, data = rd.data)
rdrobust(y = rd.data$y, x = rd.data$x)

rdhte(y = y, x = x, h = 0.1, vce = "hc3", data = rd.data)
rdrobust(y = rd.data$y, x = rd.data$x, h = 0.1, rho = 1, vce = "hc3")

rdhte(y = y, x = x, covs.hte = w_left, h = c(0.078, 0.116),
      vce = "hc3", data = rd.data)
rdrobust(y = rd.data$y[rd.data$w_left == 1],
         x = rd.data$x[rd.data$w_left == 1],
         h = 0.116, rho = 1, vce = "hc3")
