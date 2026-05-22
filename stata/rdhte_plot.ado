*!version 0.2.0  2026-05-15
* rdhte_plot -- post-estimation plot of group-specific RD treatment effects.
*
* Reads e() left by rdhte and draws one point per group (conventional
* Estimate from e(tau_hat)) with the robust bias-corrected confidence
* interval (e(tau_ci_lb), e(tau_ci_ub)). A dashed reference line at
* y = 0 visualizes the null effect.
*
* Categorical covs_hte only. Continuous covs_hte (or no covs_hte) is
* refused with a clear error.
*
* 2026-05-14: bumped version directive 14 -> 16 to match the rest of the
* package (frames migration). Wrapped the temp-frame work in a
* nobreak/capture-noisily block so the user is restored to their
* original frame on error -- mirrors the M1 retrofit applied to
* rdhte.ado and rdbwhte.ado on 2026-05-12.

capture program drop rdhte_plot
program define rdhte_plot
	version 16.0
	syntax [, sort NOZero TItle(string) XTItle(string) ///
	          YTItle(string) GRaph_options(string asis)]

	if "`e(cmd)'" != "rdhte" {
		di as error "rdhte_plot: last estimation command is not -rdhte-"
		exit 301
	}
	if "`e(is_factor)'" != "true" {
		di as error "rdhte_plot: only categorical (factor) {bf:covs_hte()} is supported."
		di as error "  The fitted rdhte model has continuous (or no) {bf:covs_hte()}."
		exit 198
	}

	tempname tau cilb ciub
	matrix `tau'  = e(tau_hat)
	matrix `cilb' = e(tau_ci_lb)
	matrix `ciub' = e(tau_ci_ub)

	* tau_hat is a 1xI row vector with colnames = group labels.
	local I = colsof(`tau')
	if `I' < 1 {
		di as error "rdhte_plot: e(tau_hat) is empty."
		exit 459
	}

	* Group labels: prefer e(lab_list) (prettified by the -labels- option),
	* fall back to the colnames of e(tau_hat).
	local labs : colfullnames `tau'
	if "`e(lab_list)'" != "" {
		local n_labs : word count `e(lab_list)'
		if `n_labs' == `I' local labs `e(lab_list)'
	}

	* Build a temporary frame with one row per group, then plot from it.
	* Using a frame keeps the user's data untouched and avoids
	* preserve/restore overhead. The nobreak/capture-noisily wrap
	* guarantees the user is returned to their original frame even if
	* the twoway call (or any setup step) errors out -- mirrors the
	* M1 retrofit pattern in rdhte.ado and rdbwhte.ado.
	local _orig = c(frame)
	tempname _plotfr
	capture frame drop `_plotfr'
	frame create `_plotfr' double(gpos est cilo cihi) str244 glab

	nobreak {
	cwf `_plotfr'
	capture noisily {

		quietly {
			set obs `I'
			forvalues i = 1/`I' {
				local lab : word `i' of `labs'
				replace gpos = `i'                in `i'
				replace est  = `tau'[1, `i']      in `i'
				replace cilo = `cilb'[`i', 1]     in `i'
				replace cihi = `ciub'[`i', 1]     in `i'
				replace glab = `"`lab'"'          in `i'
			}

			if "`sort'" != "" {
				sort est
				replace gpos = _n
			}
		}

		* xlabel string: 1 "A" 2 "B" 3 "C" ...
		* Auto-rotate x labels when there are many groups or long labels.
		local xlab ""
		local maxlen = 0
		forvalues i = 1/`I' {
			local lab = glab[`i']
			local xlab `xlab' `i' `"`lab'"'
			local len = strlen(`"`lab'"')
			if `len' > `maxlen' local maxlen = `len'
		}
		local xrotate ""
		if `I' > 6 | `maxlen' > 6 local xrotate "angle(45)"

		* Titles: defaults derived from e().
		if `"`title'"' == "" {
			local ttl "RD heterogeneous treatment effects"
			if "`e(clustvar)'" != "" local ttl "`ttl' (Cluster-Robust)"
			local title `"`ttl'"'
		}
		if `"`xtitle'"' == "" local xtitle `"`e(covs_hte)'"'
		if `"`xtitle'"' == "" local xtitle "Group"
		if `"`ytitle'"' == "" local ytitle "Treatment effect"

		local zline ""
		if "`nozero'" == "" local zline "yline(0, lpattern(dash) lcolor(gs8))"

		local xhi = `I' + 0.5

		twoway ///
			(rcap cilo cihi gpos, lcolor(black) lwidth(medthin)) ///
			(scatter est gpos, mcolor(black) msize(medium) msymbol(O)) ///
			, ///
			`zline' ///
			xlabel(`xlab', `xrotate' noticks) ///
			xscale(range(0.5 `xhi')) ///
			xtitle(`"`xtitle'"') ///
			ytitle(`"`ytitle'"') ///
			title(`"`title'"') ///
			legend(off) ///
			graphregion(color(white)) ///
			plotregion(color(white)) ///
			`graph_options'

	}   /* close capture noisily */
	local _rc = _rc
	cwf `_orig'
	capture frame drop `_plotfr'
	if (`_rc' != 0) exit `_rc'
	}   /* close nobreak */
end
