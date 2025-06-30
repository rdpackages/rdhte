{smcl}
{* *!version 0.1.0 2025-06-30}{...}
{viewerjumpto "Syntax" "rdhte_lincom##syntax"}{...}
{viewerjumpto "Description" "rdhte_lincom##description"}{...}
{viewerjumpto "Options" "rdhte_lincom##options"}{...}
{viewerjumpto "Examples" "rdhte_lincom##examples"}{...}
{viewerjumpto "Stored results" "rdhte_lincom##stored_results"}{...}
{viewerjumpto "References" "rdhte_lincom##references"}{...}
{viewerjumpto "Authors" "rdhte_lincom##authors"}{...}

{title:Title}

{p 4 8}{cmd:rdhte_lincom} {hline 2} RD Heterogeneous Treatment Effects. Linear combinations of parameters.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdhte_lincom} {it:exp} 
[{cmd:, options} 
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdhte_lincom} computes point estimates, p-values, and robust bias-corrected confidence intervals for linear combinations of parameters after any estimation using {cmd:rdhte}. 
It is based on the Stata function {help lincom}. More general post-estimation linear hypotheses can be tested with the Stata function {help test}.{p_end}

{p 8 8} Companion commands are: {help rdhte:rdhte} for estimation and inference of RD-HTE, and 
{help rdbwhte:rdbwhte} for data-driven bandwidth selection.{p_end}

{p 8 8}A detailed introduction to {cmd:rdhte} in Stata is given in
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2025_Stata.pdf":Calonico, Cattaneo, Farrell, Palomba and Titiunik (2025b)}.{p_end}

{p 4 8}Related software packages for analysis and interpretation of RD designs and related methods are available in:{p_end}

{p 8 8}{browse "https://rdpackages.github.io/":https://rdpackages.github.io/}{p_end}

{p 4 8}For background methodology, see 
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf":Calonico, Cattaneo, Farrell, and Titiunik (2019}), 
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf":Calonico, Cattaneo and Farrell (2020)}, 
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf":Cattaneo and Titiunik (2022)}.{p_end}

{marker options}{...}
{title:Options}
				
{p 4 8}{cmd:level(}{it:#}{cmd:)} specifies the confidence level, as a percentage, for confidence intervals.
The default is {cmd:level(95)} or as set by set level.{p_end}

{p 4 8}{cmd:display_options}{cmd:} cformat(%fmt), pformat(%fmt), and sformat(%fmt).{p_end}

{hline}


   {marker examples}{...}
{title:Example:}

{p 4 8}Setup using {browse "https://www.aeaweb.org/articles?id=10.1257/app.20210840":Granzier, Pons, and Tricaud (2023)} Data{p_end}
{p 8 8}{cmd:. use rdhte_dataset.dta}{p_end}

{p 4 8}RD-HTE Estimation by left/right groups{p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(i.w_ideology) vce(cluster cluster_var)}{p_end}

{p 4 8}Robust RD Estimation of HTE {p_end}
{p 8 8}{cmd:. rdhte_lincom 4.w_ideology - 3.w_ideology}{p_end}

{p 4 8}Testing for equality of the effects{p_end}
{p 8 8}{cmd:. test 4.w_ideology = 3.w_ideology = 2.w_ideology}{p_end}

{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rdhte_lincom} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(rdhte_lincom_est)}}point estimate{p_end}
{synopt:{cmd:e(rdhte_lincom_pv)}}p-value{p_end}
{synopt:{cmd:e(rdhte_lincom_lb)}}lower bound of confidence interval{p_end}
{synopt:{cmd:e(rdhte_lincom_ub)}}upper bound of confidence interval{p_end}
{synopt:{cmd:e(rdhte_lincom_level)}}confidence level{p_end}
  
	  
	  
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

{p 4 8}Calonico, Cattaneo, and Farrell. 2020.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf":Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs}.
{it:Econometrics Journal}, 23(2): 192-210.{p_end}

{p 4 8}Calonico, Cattaneo, Farrell, and Titiunik. 2019.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf":Regression Discontinuity Designs using Covariates}.
{it:Review of Economics and Statistics}, 101(3): 442-451.{p_end}

{p 4 8}Calonico, Cattaneo, and Titiunik. 2014.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf":Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs}.
{it:Econometrica}, 82(6): 2295-2326.{p_end}




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



