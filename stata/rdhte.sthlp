{smcl}
{* *!version 0.1.1 2025-08-08}{...}
{viewerjumpto "Syntax" "rdhte##syntax"}{...}
{viewerjumpto "Description" "rdhte##description"}{...}
{viewerjumpto "Options" "rdhte##options"}{...}
{viewerjumpto "Examples" "rdhte##examples"}{...}
{viewerjumpto "Stored results" "rdhte##stored_results"}{...}
{viewerjumpto "References" "rdhte##references"}{...}
{viewerjumpto "Authors" "rdhte##authors"}{...}

{title:Title}

{p 4 8}{cmd:rdhte} {hline 2} RD Heterogeneous Treatment Effects Estimation and Inference.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdhte} {it:depvar} {it:runvar} {ifin} 
[{cmd:,} 
{cmd:covs_hte(}{it:covars}{cmd:)}
{cmd:c(}{it:#}{cmd:)} 
{cmd:p(}{it:#}{cmd:)} 
{cmd:q(}{it:#}{cmd:)} 
{cmd:h(}{it:#}{cmd:)} 
{cmd:h_l(}{it:#}{cmd:)} 
{cmd:h_r(}{it:#}{cmd:)} 
{cmd:kernel(}{it:kernelfn}{cmd:)}
{cmd:vce(}{it:vcetype}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:covs_eff(}{it:covars}{cmd:)}
{cmd:bwjoint}
{cmd:labels}
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdhte} provides estimation and inference for heterogeneous treatment effects in RD designs using local polynomial regressions, allowing for interactions with pretreatment covariates ({browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_HTERD.pdf":Calonico, Cattaneo, Farrell, Palomba and Titiunik, 2025a}). Inference is implemented using robust bias-correction methods ({browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf":Calonico, Cattaneo, and Titiunik, 2014}){p_end}

{p 8 8} Companion commands are: {help rdbwhte:rdbwhte} for data-driven bandwidth selection and {help rdhte_lincom:rdhte_lincom} for testing linear restrictions of paramaters. More general post-estimation linear hypotheses can be tested with the Stata function {help test}.{p_end}

{p 8 8}A detailed introduction to {cmd:rdhte} in Stata is given in
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_Stata.pdf":Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025b)}.{p_end}

{p 4 8}Related software packages for analysis and interpretation of RD designs and related methods are available in:{p_end}

{p 8 8}{browse "https://rdpackages.github.io/":https://rdpackages.github.io/}{p_end}

{p 4 8}For background methodology, see {browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf":Calonico, Cattaneo, Farrell, and Titiunik (2019}), {browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf":Calonico, Cattaneo and Farrell (2020)}, 
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf":Cattaneo and Titiunik (2022)}.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8}{cmd:c(}{it:#}{cmd:)} specifies the RD cutoff for {it:indepvar}.
Default is {cmd:c(0)}.{p_end}

{p 4 8}{cmd:covs_hte(}{it:covars}{cmd:)} specifies covariate(s) for heterogeneous treatment effects. 
{help fvvarlist:Factor variables} notation can be used to distinguish between continuous and categorical variables, select reference categories, specify interactions between variables, and include polynomials of continuous variables.
If not specified, the RD Average Treatment Effect is computed.{p_end}

{p 4 8}{cmd:labels} displays the final RD estimates using variable labels from {cmd:covs_hte(}{it:covars}{cmd:)}.{p_end}


{dlgtab:Local Polynomial Regression}

{p 4 8}{cmd:p(}{it:#}{cmd:)} specifies the order of the local polynomial used to construct the point estimator.
Default is {cmd:p(1)} (local linear regression).{p_end}

{p 4 8}{cmd:q(}{it:#}{cmd:)} specifies the order of the local polynomial used to construct the bias correction. 
Default is q(2) (local quadratic regression).{p_end}

{p 4 8}{cmd:h(}{it:#}{cmd:)}, {cmd:h_l(}{it:#}{cmd:)} and {cmd:h_r(}{it:#}{cmd:)} set the bandwidths to construct the RD estimator. The same choice could be used on each side of the cutoff (via {cmd:h(}{it:#}{cmd:)}), or different to the left and right (using {cmd:h_l(}{it:#}{cmd:)} and {cmd:h_r(}{it:#}{cmd:)}). More than one bandwidth can be specified for categorical covariates.
If not specified, bandwidths are computed by the companion command {help rdbwhte:rdbwhte}.
{p_end}

{p 4 8}{cmd:kernel(}{it:kernelfn}{cmd:)} specifies the kernel function used to construct the local-polynomial estimator(s). Options are: {opt tri:angular}, {opt epa:nechnikov}, and {opt uni:form}.
Default is {cmd:kernel(triangular)}.{p_end}

{p 4 8}{cmd:covs_eff(}{it:covars}{cmd:)} specifies additional covariates to be used for efficiency improvements.{p_end}


{dlgtab:Data-Driven Bandwidth Selection}

{p 4 8}{cmd:bwselect(}{it:bwmethod}{cmd:)} specifies the bandwidth selection procedure to be used.
Options are:{p_end}
{p 8 12}{opt mserd} one common MSE-optimal bandwidth selector for the RD treatment effect estimator.{p_end}
{p 8 12}{opt msetwo} two different MSE-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.{p_end}
{p 8 12}{opt msesum} one common MSE-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).{p_end}
{p 8 12}{opt msecomb1} for min({opt mserd},{opt msesum}).{p_end}
{p 8 12}{opt msecomb2} for median({opt msetwo},{opt mserd},{opt msesum}), for each side of the cutoff separately.{p_end}
{p 8 12}{opt cerrd} one common CER-optimal bandwidth selector for the RD treatment effect estimator.{p_end}
{p 8 12}{opt certwo} two different CER-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.{p_end}
{p 8 12}{opt cersum} one common CER-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).{p_end}
{p 8 12}{opt cercomb1} for min({opt cerrd},{opt cersum}).{p_end}
{p 8 12}{opt cercomb2} for median({opt certwo},{opt cerrd},{opt cersum}), for each side of the cutoff separately.{p_end}
{p 8 12}Note: MSE = Mean Square Error; CER = Coverage Error Rate.{p_end}
{p 8 12}Default is {cmd:bwselect(mserd)}. {p_end}

{p 4 8}{cmd:bwjoint} forces all bandwidths to be the same across groups.
Default is {cmd:level(95)}.{p_end}


{dlgtab:Variance-Covariance Estimation}

{p 4 8}{cmd:vce(}{it:{help regress##vcetype:vcetype}}{cmd:)} vcetype may be robust, cluster clustvarlist, hc2 [clustvar], or hc3.{p_end}
{p 8 12}Default is {cmd:vce(hc3)}.{p_end}

{p 4 8}{cmd:level(}{it:#}{cmd:)} specifies confidence level for confidence intervals.{p_end}


   {hline}

   {marker examples}{...}
{title:Example:}

{p 4 8}Setup using {browse "https://www.aeaweb.org/articles?id=10.1257/app.20210840":Granzier, Pons, and Tricaud (2023)} Data{p_end}
{p 8 8}{cmd:. use rdhte_dataset.dta}{p_end}

{p 4 8}RD-HTE Estimation by left/right groups{p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(i.left)}{p_end}

{p 4 8}RD-HTE Estimation by left/right groups with common bandwidth{p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(i.left) bwjoint}{p_end}

{p 4 8}RD-HTE Estimation by left/right groups and strong {p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(i.left#i.strong)}{p_end}

{p 4 8}RD-HTE Estimation using a continuous variable{p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(c.mean_strength_nat1) }{p_end}

{p 4 8}RD-HTE Estimation using a continuous variable with clustered standard errors {p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(c.mean_strength_nat1)  vce(cluster id_district)}{p_end}


   {hline}
   
{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rdhte} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}original number of observations{p_end}
{synopt:{cmd:e(c)}}cutoff value{p_end}
{synopt:{cmd:e(p)}}order of the polynomial used for estimation of the regression function{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(runningvar)}}name of running variable{p_end}
{synopt:{cmd:e(outcomevar)}}name of outcome variable{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(covs)}}name of covariates{p_end}
{synopt:{cmd:e(vce_select)}}vcetype specified in vce(){p_end}
{synopt:{cmd:e(kernel)}}kernel choice{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(h)}}bandwidth{p_end}
{synopt:{cmd:e(tau_hat)}}p-order local-polynomial estimates for the outcome variable{p_end}
{synopt:{cmd:e(tau_bc)}}bias-corrected local-polynomial estimates for the outcome variable{p_end}
{synopt:{cmd:e(tau_se)}}robust standard errors{p_end}
{synopt:{cmd:e(tau_V)}}robust variance-covariance matrix{p_end}
{synopt:{cmd:e(tau_t)}}robust t-statistics{p_end}
{synopt:{cmd:e(tau_pv)}}robust p-values{p_end}
{synopt:{cmd:e(tau_N)}}sample size{p_end}
{synopt:{cmd:e(tau_ci_lb)}}robust lower bound confidence interval{p_end}
{synopt:{cmd:e(tau_ci_ub)}}robust upper bound confidence interval{p_end}
	
	   {hline}
	
{marker references}{...}
{title:References}

{p 4 8}Calonico, Cattaneo, Farrell, Palomba and Titiunik. 2025a.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_HTERD.pdf":Treatment Effect Heterogeneity in Regression Discontinuity Designs}.
{it:Working Paper}.{p_end}

{p 4 8}Calonico, Cattaneo, Farrell, Palomba and Titiunik. 2025b.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_Stata.pdf":rdhte: Conditional Average Treatment Effects in RD Designs}.
{it:Working Paper}.{p_end}

{p 4 8}Granzier, Pons, and Tricaud. 2023.
{browse "https://www.aeaweb.org/articles?id=10.1257/app.20210840":Coordination and Bandwagon Effects: How Past Rankings Shape the Behavior of Voters and Candidates}.
{it:American Economic Journal: Applied Economics}, 15(4): 177–217.{p_end}

{p 4 8}Cattaneo and Titiunik. 2022.
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf":Regression Discontinuity Designs}.
{it:Annual Review of Economics}, 14: 821-851.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf":Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs}.
{it:Econometrics Journal}, 23(2): 192-210.{p_end}

{p 4 8}Calonico, Cattaneo, Farrell, and Titiunik. 2019.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf":Regression Discontinuity Designs using Covariates}.
{it:Review of Economics and Statistics}, 101(3): 442-451.{p_end}

{p 4 8}Calonico, Cattaneo, and Titiunik. 2014.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf":Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs}.
{it:Econometrica}, 82(6): 2295-2326.{p_end}


   {hline}

{marker authors}{...}
{title:Authors}

{p 4 8}Sebastian Calonico, University of California, Davis, CA.
{browse "mailto:scalonico@ucdavis.edu":scalonico@ucdavis.edu}.{p_end}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Max H. Farrell, University of California, Santa Barbara, CA.
{browse "mailto:maxhfarrell@ucsb.edu":maxhfarrell@ucsb.edu}.{p_end}

{p 4 8}Filippo Palomba, Princeton University, Princeton, NJ.
{browse "mailto:fpalomba@princeton.edu":fpalomba@princeton.edu}.{p_end}

{p 4 8}Rocio Titiunik, Princeton University, Princeton, NJ.
{browse "mailto:titiunik@princeton.edu":titiunik@princeton.edu}.{p_end}


