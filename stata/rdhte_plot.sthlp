{smcl}
{* *! version 0.2.0  2026-05-15}{...}
{viewerjumpto "Syntax" "rdhte_plot##syntax"}{...}
{viewerjumpto "Description" "rdhte_plot##description"}{...}
{viewerjumpto "Options" "rdhte_plot##options"}{...}
{viewerjumpto "Examples" "rdhte_plot##examples"}{...}
{title:Title}

{p 4 8}{cmd:rdhte_plot} {hline 2} Plot RD heterogeneous treatment effects across groups.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:rdhte_plot}
[{cmd:,}
{opt sort}
{opt nozero}
{opt ti:tle(string)}
{opt xti:tle(string)}
{opt yti:tle(string)}
{opt gr:aph_options(string)}
]

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdhte_plot} is a post-estimation command for {help rdhte}. It reads
the {help ereturn:e()} results left by the most recent {cmd:rdhte} run and draws
one point per heterogeneity group at the conventional point estimate
({cmd:e(tau_hat)}) with the robust bias-corrected confidence interval
({cmd:e(tau_ci_lb)}, {cmd:e(tau_ci_ub)}). A dashed horizontal line at
{it:y = 0} provides a visual reference for the null effect.{p_end}

{p 4 8}The command requires that the preceding {cmd:rdhte} call used a
categorical {cmd:covs_hte()} (i.e. a factor variable). Continuous
{cmd:covs_hte()} or no {cmd:covs_hte()} are refused with an explanatory
error.{p_end}

{p 4 8}The intervals shown are the same robust bias-corrected CIs reported in
the {cmd:rdhte} table: they are centered on {cmd:e(tau_bc)} (not
{cmd:e(tau_hat)}) and use the robust standard error. Because the point and
the CI center can differ slightly, the point may sit just inside or just
outside the bar; this is the rdrobust convention and is not a plotting bug.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{opt sort} reorders groups along the x-axis by point estimate
(ascending). Default keeps the original group order.{p_end}

{p 4 8}{opt nozero} suppresses the dashed reference line at {it:y = 0}.{p_end}

{p 4 8}{opt title(string)} overrides the plot title. Default is derived from
the {cmd:rdhte} model and a "(Cluster-Robust)" suffix when {cmd:e(clustvar)}
is set.{p_end}

{p 4 8}{opt xtitle(string)} overrides the x-axis title. Default is the name of
the heterogeneity variable ({cmd:e(covs_hte)}).{p_end}

{p 4 8}{opt ytitle(string)} overrides the y-axis title. Default is
"Treatment effect".{p_end}

{p 4 8}{opt graph_options(string)} is a free-form string passed verbatim to
{help twoway} after the built-in options. Use to set scheme, axis ranges,
custom legends, etc.{p_end}

{marker examples}{...}
{title:Examples}

{p 4 8}Categorical heterogeneity with cluster-robust inference{p_end}
{p 8 8}{cmd:. rdhte y x, covs_hte(i.region) vce(cluster id)}{p_end}
{p 8 8}{cmd:. rdhte_plot}{p_end}

{p 4 8}Sort groups by effect size and rename the y-axis{p_end}
{p 8 8}{cmd:. rdhte_plot, sort ytitle("ATE on outcome")}{p_end}

{p 4 8}Pass extra twoway options{p_end}
{p 8 8}{cmd:. rdhte_plot, graph_options(scheme(s2color) ylabel(0(0.5)2.5))}{p_end}

{title:Author}

{p 4 8}Sebastian Calonico, UC Davis. {browse "mailto:scalonico@ucdavis.edu":scalonico@ucdavis.edu}.{p_end}

{title:See also}
{p 4 8}{help rdhte}, {help rdbwhte}, {help rdhte_lincom}{p_end}
