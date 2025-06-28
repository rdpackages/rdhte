******************************************************************************** 
** RDHTE Stata Package
** Do-file for Empirical Illustration 
** Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Filippo Palomba and Rocio Titiunik 
********************************************************************************
** hlp2winpdf, cdn(rdhte) replace
** hlp2winpdf, cdn(rdbwhte) replace
** hlp2winpdf, cdn(rdhte_lincom) replace
********************************************************************************
** net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
** net install rdhte, from(https://raw.githubusercontent.com/rdpackages/rdhte/main/stata) replace
********************************************************************************

clear all

use "rdhte_dataset.dta", clear


**************************************************************************
* Single binary variable
**************************************************************************

rdhte y x, covs_hte(w_left) vce(hc2 cluster_var)

*** Post-estimation illustration -- rdhte_lincom
* How large is the advantage for left-of-center candidates?

rdhte_lincom 1.w_left - 0.w_left


* Extra illustration:
* Note that forcing a common bandwidth makes the effect for `center or right' 
* (incorrectly) statistically significant

rdhte y x, covs_hte(w_left) vce(hc2 cluster_var) bwjoint

* Extra illustration:

*rdhte supports many, but not all, factor variable specifications in Stata
rdhte y x, covs_hte(w_left) vce(hc2 cluster_var)
rdhte y x, covs_hte(i.w_left) vce(hc2 cluster_var)
rdhte y x, covs_hte(ibn.w_left) vce(hc2 cluster_var)


**************************************************************************
* Single categorical variable -- unordered
**************************************************************************

rdhte y x, covs_hte(i.w_ideology) vce(hc2 cluster_var) labels

*** Post-estimation illustration -- Stata's test command
* All the non-left groups are statistically 
* indistinguishable from each other and from zero

test 4.w_ideology = 3.w_ideology = 2.w_ideology = 0


**************************************************************************
* Single categorical variable -- ordered
**************************************************************************

rdhte y x, covs_hte(i.w_strength_qrt) vce(hc2 cluster_var)

* -> result: advantage increases with strength 
             *candidates with higher average strength at the national level have higher effect, monotonic!

			 
* Extra illustration:
*Will be treated as continuous by default

rdhte y x, covs_hte(w_strength_qrt) vce(hc2 cluster_var)



**************************************************************************
* Two binary variables - interaction
**************************************************************************

rdhte y x, covs_hte(i.w_left#i.w_strong) vce(hc2 cluster_var)


**************************************************************************
* Single continuous variable
**************************************************************************

rdhte y x, covs_hte(w_strength) kernel(uni) vce(hc2 cluster_var)

* to aid with interpretation, notice that the coefficient on T#c.mean_strength_nat1
* is _precisely_ a slope coefficient in a linear model, and can be interpreted the same.
* To see this, use the uniform kernel and a set bandwidth to fit local unweighted least squares
* then match this with the base command -regress- 

local bw = e(h)[1,1]
gen T = (x>0)
reg y T##c.x##c.w_strength if abs(x)<=`bw'

* IMPORTANT: inference requires robust bias correction and cannot be obtained from this regression


**************************************************************************
* Interaction effect: binary*continuous
**************************************************************************

*Full interaction is most natural (different intercept and slope for each category)

rdhte y x, covs_hte(i.w_left##c.w_strength) vce(hc2 cluster_var)

* Extra illustration:
* Each effect is insignificant, but the joint test shows there is information
test T  T#1.w_left  T#c.w_strength  T#1.w_left#c.w_strength
test (T = T#1.w_left = T#c.w_strength = T#1.w_left#c.w_strength = 0)

* Extra illustration:
* to aid interpretation, the fully interacted model will match results from 
* category-specific estimation. Fix the bandwidth to make results match exactly
*Match using a fixed bandwidth

rdhte y x, covs_hte(i.w_left##c.w_strength) h(0.1) vce(hc2 cluster_var)
rdhte_lincom T#c.w_strength + T#1.w_left#c.w_strength
rdhte y x if w_left==1, covs_hte(c.w_strength) h(0.1) vce(hc2 cluster_var)


**************************************************************************
* Replicating rdrobust
**************************************************************************

*** Average effects

*Using default settings, packages will not match, because of different settings

rdhte y x
rdrobust y x


* to replicate exactly with RDROBUST, set rho=1 and specific vce() option
* bandwidth selection is also different, so enforce the same bandwidth

rdhte y x, h(0.1) vce(hc3)
rdrobust y x, h(0.1) rho(1) vce(hc3)


*** subgroup analysis

* keeping the settings the same as above, rdhte will match rdrobust if the latter
* is run for each subgroup separately

rdhte y x, covs_hte(w_left) h(0.078 0.116) 
rdrobust y x if w_left==1, h(0.116) rho(1) vce(hc3)

