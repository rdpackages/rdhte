********************************************************************************
** RDHTE Package
** Empirical Illustration
********************************************************************************

**************************************************************************
* Data and design
**************************************************************************
* This do-file illustrates `rdhte` on a real-data extract from Granzier,
* Pons, and Tricaud (2023, AEJ: Applied), "Coordination and Bandwagon
* Effects: How Past Rankings Shape the Behavior of Voters and
* Candidates." The authors study French two-round elections, where
* candidates must clear a qualifying-vote threshold in the first round
* to advance to the runoff. The institutional rule creates a sharp RD
* design on every candidate's first-round margin against that threshold:
* candidates just above the cutoff advance, those just below do not.
* The authors use this design to ask whether being just-qualified
* causally affects subsequent voter and candidate behavior, and
* document substantive heterogeneity across ideology and
* candidate-strength dimensions.
*
* The bundled extract `rdhte_dataset.dta` has 39,534 candidate-race
* observations with the following variables:
*
*   y                outcome: 1 if the candidate advances to the runoff
*   x                running variable: first-round margin against the
*                    qualifying threshold (cutoff at zero)
*   cluster_var      district identifier (cluster-robust inference)
*   w_left           1 if the candidate's party is left of center
*   w_ideology       unordered party-ideology bucket (4 levels)
*   w_strength       continuous proxy for ex-ante candidate strength
*   w_strong         1 if above-median strength
*   w_strength_qrt   ordered quartile bucket of w_strength
*
* These covariates support every covariate-incorporation pattern that
* rdhte exposes: binary cells, multi-level (unordered and ordered)
* factor cells, factor-by-factor interactions, continuous slopes, and
* binary x continuous interactions.

clear all

use "rdhte_dataset.dta", clear


**************************************************************************
* Single binary variable
**************************************************************************

rdhte y x, covs_hte(w_left) vce(cluster cluster_var)

*** Post-estimation illustration -- rdhte_lincom
* How large is the advantage for left-of-center candidates?

rdhte_lincom 1.w_left - 0.w_left


* Extra illustration:
* Note that forcing a common bandwidth makes the effect for `center or right'
* (incorrectly) statistically significant

rdhte y x, covs_hte(w_left) vce(cluster cluster_var) bwjoint

* Extra illustration:

*rdhte supports many, but not all, factor variable specifications in Stata
rdhte y x, covs_hte(w_left) vce(cluster cluster_var)
rdhte y x, covs_hte(i.w_left) vce(cluster cluster_var)
rdhte y x, covs_hte(ibn.w_left) vce(cluster cluster_var)


**************************************************************************
* Single categorical variable -- unordered
**************************************************************************

rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var) labels

*** Post-estimation illustration -- Stata's test command
* All the non-left groups are statistically
* indistinguishable from each other and from zero

test 4.w_ideology = 3.w_ideology = 2.w_ideology = 0


**************************************************************************
* Single categorical variable -- ordered
**************************************************************************

rdhte y x, covs_hte(i.w_strength_qrt) vce(cluster cluster_var)

* -> result: advantage increases with strength
             *candidates with higher average strength at the national level have higher effect, monotonic!


* Extra illustration:
*Will be treated as continuous by default
rdhte y x, covs_hte(w_strength_qrt) vce(cluster cluster_var)



**************************************************************************
* Two binary variables - interaction
**************************************************************************

rdhte y x, covs_hte(i.w_left#i.w_strong) vce(cluster cluster_var)


**************************************************************************
* Single continuous variable
**************************************************************************

rdhte y x, covs_hte(w_strength) kernel(uni) vce(cluster cluster_var)

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

rdhte y x, covs_hte(i.w_left##c.w_strength) vce(cluster cluster_var)

* Extra illustration:
* Each effect is insignificant, but the joint test shows there is information
test T  T#1.w_left  T#c.w_strength  T#1.w_left#c.w_strength
test (T = T#1.w_left = T#c.w_strength = T#1.w_left#c.w_strength = 0)

* Extra illustration:
* to aid interpretation, the fully interacted model will match results from
* category-specific estimation. Fix the bandwidth to make results match exactly

rdhte y x, covs_hte(i.w_left##c.w_strength) h(0.1) vce(cluster cluster_var)
rdhte_lincom T#c.w_strength + T#1.w_left#c.w_strength
rdhte y x if w_left==1, covs_hte(c.w_strength) h(0.1) vce(cluster cluster_var)
rdhte y x if w_left==0, covs_hte(c.w_strength) h(0.1) vce(cluster cluster_var)


**************************************************************************
* Standalone bandwidth selection (rdbwhte)
**************************************************************************

* `rdbwhte` runs the same data-driven bandwidth selectors as `rdhte`
* but returns the bandwidths without estimating the CATEs. Useful when
* you want to inspect or fix `h` and then explore alternative variance
* estimators or compare specifications at a common bandwidth.

rdbwhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var)
matrix list e(h)

* `bwjoint` forces a single shared bandwidth across cells.
rdbwhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var) bwjoint
matrix list e(h)


**************************************************************************
* Efficiency-improving covariates (covs_eff)
**************************************************************************

* `covs_eff` adds covariates to the local-polynomial regression to
* shrink standard errors WITHOUT changing the identification of the
* CATE. They enter additively (and as covs_eff x W interactions in the
* heterogeneity-aware paths) but never with the treatment indicator or
* the running-variable polynomial. Useful when you have strong
* pretreatment predictors of the outcome.

quietly rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var)
estimates store m_eff_off
quietly rdhte y x, covs_hte(i.w_ideology) covs_eff(w_strength) ///
        vce(cluster cluster_var)
estimates store m_eff_on
estimates table m_eff_off m_eff_on, b(%9.4f) se(%9.4f) ///
    title("Effect of adding covs_eff = w_strength on SEs")


**************************************************************************
* Plotting (rdhte_plot)
**************************************************************************

* `rdhte_plot` is a post-estimation command for categorical covs_hte:
* one point per cell at the conventional point estimate, with the
* robust bias-corrected CI shown as an error bar.

rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var)
rdhte_plot

* Sort cells by point estimate.
rdhte_plot, sort

* Customize titles and pass extra twoway options.
rdhte_plot, sort title("Heterogeneity by ideology bucket") ///
            ytitle("Sharp RD ITT") ///
            graph_options(scheme(s2color))


**************************************************************************
* Building publication-ready tables
**************************************************************************
* rdhte stores all per-cell results in e() matrices; the summary panel
* printed at estimation time already covers the typical use case, but
* for paper-ready tables it is useful to extract the same numbers into
* a stand-alone display or LaTeX file. Three patterns are illustrated:
*
*   (A) per-cell estimate / SE / CI / N / h table from a SINGLE rdhte call
*   (B) cell x specification comparison (varying vce across columns)
*   (C) LaTeX export of (A) ready for inclusion in a paper
*
* All three pull directly from e() after estimation, so they work with
* any covs_hte spec that produces a categorical group structure.

* ---- (A) per-cell table from a single call ----
quietly rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var) labels
local I        = colsof(e(tau_hat))
local lab_list "`e(lab_list)'"

di _n "Table A: per-cell estimates  (covs_hte = i.w_ideology, vce(cluster))"
di in gr "{hline 92}"
di in gr  %-18s "Cell"          ///
       _col(20) %9s "Estimate"  ///
       _col(30) %9s "SE"        ///
       _col(40) %9s "z"         ///
       _col(50) %9s "Pr>|z|"    ///
       _col(60) %18s "[95% Conf. Int.]" ///
       _col(80) %6s "Nh-"       ///
       _col(87) %6s "Nh+"
di in gr "{hline 92}"
forval i = 1/`I' {
    local lab : word `i' of `lab_list'
    di in ye %-18s "`lab'" ///
       _col(20) %9.3f e(tau_hat)[1, `i']                   ///
       _col(30) %9.3f e(tau_se)[`i', 1]                    ///
       _col(40) %9.3f e(tau_t)[`i', 1]                     ///
       _col(50) %9.3f e(tau_pv)[`i', 1]                    ///
       _col(60) %9.3f e(tau_ci_lb)[`i', 1]                 ///
       _col(70) %9.3f e(tau_ci_ub)[`i', 1]                 ///
       _col(80) %6.0f e(tau_N)[`i', 1]                     ///
       _col(87) %6.0f e(tau_N)[`i', 2]
}
di in gr "{hline 92}"

* The same numbers can be turned into a temporary dataset via
* `svmat` so they round-trip through `outsheet`, `frames`, `putexcel`,
* etc. Skipped here -- the manual `display` is enough for inline use.


* ---- (B) cell x specification comparison ----
* Fix the covs_hte spec and vary vce across columns. Each column stores
* the per-cell conventional point estimate and its robust BC standard
* error in parentheses.

local specs hc1 hc2 hc3 cr1
local n_specs : word count `specs'

* One pass to get the cell count and labels.
quietly rdhte y x, covs_hte(i.w_ideology) vce(hc1) labels
local I        = colsof(e(tau_hat))
local lab_list "`e(lab_list)'"

* Accumulate I x n_specs matrices for points and SEs.
matrix POINT = J(`I', `n_specs', .)
matrix SE_   = J(`I', `n_specs', .)
local j = 1
foreach v of local specs {
    * Cluster-robust specs need the cluster variable as the second token.
    if inlist("`v'", "cr1", "cr2", "cr3") {
        quietly rdhte y x, covs_hte(i.w_ideology) vce(`v' cluster_var) labels
    }
    else {
        quietly rdhte y x, covs_hte(i.w_ideology) vce(`v') labels
    }
    forval i = 1/`I' {
        matrix POINT[`i', `j'] = e(tau_hat)[1, `i']
        matrix SE_[`i', `j']   = e(tau_se)[`i', 1]
    }
    local ++j
}

di _n "Table B: per-cell point estimates by variance option"
di       "         (SE under each estimate in parentheses)"
di in gr "{hline 72}"
di in gr %-18s "Cell" _col(20) %12s "HC1" %12s "HC2" %12s "HC3" %12s "CR1"
di in gr "{hline 72}"
forval i = 1/`I' {
    local lab : word `i' of `lab_list'
    di in ye %-18s "`lab'" _col(20) ///
       %12.3f POINT[`i', 1] %12.3f POINT[`i', 2] ///
       %12.3f POINT[`i', 3] %12.3f POINT[`i', 4]
    di in ye %-18s ""      _col(20) ///
       "(" %5.3f SE_[`i', 1] ")  (" %5.3f SE_[`i', 2] ")  (" ///
           %5.3f SE_[`i', 3] ")  (" %5.3f SE_[`i', 4] ")"
}
di in gr "{hline 72}"


* ---- (C) LaTeX export of Table A ----
* Manual file-write (esttab works off e(b) which rdhte posts as the BC
* point estimate; the convention here is to report the conventional
* point with the RBC CI).

quietly rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var) labels
local I        = colsof(e(tau_hat))
local lab_list "`e(lab_list)'"

tempname fh
file open `fh' using "rdhte_table_A.tex", write replace text
file write `fh' "\begin{tabular}{lccccc}" _n
file write `fh' "\toprule" _n
* `$h$' would be eaten by Stata global-macro expansion; write the
* dollar signs as `\$` so the .tex file gets a literal `$h$`.
file write `fh' " & Estimate & SE & 95\% RBC CI & Eff.~N & \$h\$ \\" _n
file write `fh' "\midrule" _n
forval i = 1/`I' {
    local lab : word `i' of `lab_list'
    local pt  : di %5.3f e(tau_hat)[1, `i']
    local se  : di %5.3f e(tau_se)[`i', 1]
    local lo  : di %5.3f e(tau_ci_lb)[`i', 1]
    local hi  : di %5.3f e(tau_ci_ub)[`i', 1]
    local nl  : di %6.0f e(tau_N)[`i', 1]
    local nr  : di %6.0f e(tau_N)[`i', 2]
    local hl  : di %5.3f e(h)[`i', 1]
    local hr  : di %5.3f e(h)[`i', 2]
    file write `fh' "`lab' & " "`=trim("`pt'")'" " & " "`=trim("`se'")'" ///
                 " & [`=trim("`lo'")', `=trim("`hi'")'] " ///
                 " & `=trim("`nl'")'/`=trim("`nr'")' " ///
                 " & `=trim("`hl'")'/`=trim("`hr'")' \\" _n
}
file write `fh' "\bottomrule" _n
file write `fh' "\end{tabular}" _n
file close `fh'
di _n "Wrote rdhte_table_A.tex"


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
