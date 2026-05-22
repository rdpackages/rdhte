************************************************************************************************
* RDHTE STATA PACKAGE -- rdbwhte
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
* 2026-05-15: Extended the missing-value drop scan to cover covs_eff /
*             covs_hte (via fvrevar ... list) / clustvar / weights, so
*             e(N) matches the actual analytic sample. Mirrors the
*             companion fix in rdhte.ado.
*             M6 cosmetic: summary-header magic offsets hoisted into
*             col_* locals (mirrors rdhte.ado). Pure formatting.
* 2026-05-15 (returns audit): e(vce_select) now holds the canonical raw
*             name (matching rdrobust.ado:1271); added e(vce_type) for the
*             display label. Added e(N_h_l) / e(N_h_r) (sum of tau_N rows),
*             mirroring rdhte.ado and rdrobust's per-side counts.
*!version 0.2.0  2026-05-22

capture program drop rdbwhte
program define rdbwhte, eclass
    version 16.0
    syntax varlist(min=2 max=2)  [if] [in] [ , c(real 0) p(integer 1) q(real 0) covs_hte(string) bwselect(string) covs_eff(varlist) kernel(string) weights(string) vce(string) cluster(string) precision(string) bwjoint labels]

	marksample touse, novarlist

	* Extract main variables from varlist
	local y : word 1 of `varlist'
	local x : word 2 of `varlist'

	* M1: do all work in an isolated temporary frame instead of preserve
	* on the user's data. See rdhte.ado for the full rationale. The
	* nobreak/capture wrap guarantees the user is returned to their
	* original frame even if the work block errors (mirrors rdrobust
	* M1 pattern 2026-05-12).
	local _orig_frame `c(frame)'
	tempname _rdbwhte_frame
	capture frame drop `_rdbwhte_frame'
	frame put * if `touse', into(`_rdbwhte_frame')

	nobreak {
	cwf `_rdbwhte_frame'
	capture noisily {

	qui drop if mi(`y') | mi(`x')

	if ("`bwselect'"=="") local bwselect = "mserd"

	******************** Set Kernel ***************************
	local kernel   = lower("`kernel'")
	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
		local kernel_type = "Epanechnikov"
	}
	else if ("`kernel'"=="uniform" | "`kernel'"=="uni") {
		local kernel_type = "Uniform"
	}
	else {
		local kernel = "tri"
		local kernel_type = "Triangular"
	}

	******************** Set VCE ***************************
	tokenize `vce'
	local count : word count `vce'
	if `count' == 0 {
		local vce_select ""
	}
	else if `count' == 1 {
		local vce_select `"`1'"'
	}
	else if `count' == 2 {
		local vce_select `"`1'"'
		if !inlist("`vce_select'", "cluster", "hc0", "hc1", "hc2", "hc3", "cr0", "cr1", "cr2", "cr3") {
			di as error "{err}{cmd:vce()} incorrectly specified: only " ///
			            "{cmd:vce(cluster} {it:clustvar}{cmd:)}, " ///
			            "{cmd:vce(cr1/cr2/cr3} {it:clustvar}{cmd:)}, or " ///
			            "{cmd:vce(hc0/hc1/hc2/hc3} {it:clustvar}{cmd:)} accept a 2nd argument"
			exit 198
		}
		local clustvar `"`2'"'
	}
	else if `count' > 2 {
		di as error "{err}{cmd:vce()} incorrectly specified: too many tokens"
		exit 198
	}

	* With a cluster variable, map hc0/hc1/hc2/hc3 to cr1/cr1/cr2/cr3.
	if ("`clustvar'"!="") {
		if ("`vce_select'"=="hc0" | "`vce_select'"=="hc1" | "`vce_select'"=="cr0") {
			di as text "Warning: vce(`vce_select' `clustvar') is not a cluster option. Switching to vce(cr1 `clustvar')."
			local vce_select = "cr1"
		}
		if ("`vce_select'"=="hc2") {
			di as text "Warning: vce(hc2 `clustvar') is not a cluster option. Switching to vce(cr2 `clustvar')."
			local vce_select = "cr2"
		}
		if ("`vce_select'"=="hc3") {
			di as text "Warning: vce(hc3 `clustvar') is not a cluster option. Switching to vce(cr3 `clustvar')."
			local vce_select = "cr3"
		}
		if ("`vce_select'"=="cluster") local vce_select = "cr1"
	}

	if inlist("`vce_select'","cr1","cr2","cr3") & "`clustvar'"=="" {
		di as error "{err}{cmd:vce(`vce_select' clustervar)} requires a cluster variable"
		exit 125
	}
	if "`vce_select'"=="cr0" & "`clustvar'"=="" {
		di as text "Warning: vce(cr0) requires a cluster variable. Falling back to vce(hc0)."
		local vce_select = "hc0"
	}

	if ("`vce_select'"=="")              local vce_select = "hc3"

	local vce_type = "HC3"
	if ("`vce_select'"=="hc0")      	 local vce_type = "HC0"
	if ("`vce_select'"=="hc1")      	 local vce_type = "HC1"
	if ("`vce_select'"=="robust")      	 local vce_type = "HC1"
	if ("`vce_select'"=="hc2")      	 local vce_type = "HC2"
	if ("`vce_select'"=="hc3")      	 local vce_type = "HC3"
	if ("`vce_select'"=="cr1")  	 	 local vce_type = "CR1"
	if ("`vce_select'"=="cr2")  	 	 local vce_type = "CR2"
	if ("`vce_select'"=="cr3")  	 	 local vce_type = "CR3"

	if ("`clustvar'"!="") local cluster = "cluster"

	* Build an rdbwselect-compatible vce string. rdbwselect (rdrobust 4.0.0+)
	* accepts cr1/cr2/cr3 directly; older versions only know hc0/hc1/hc2/hc3
	* and cluster <var>, so we keep the legacy path too.
	if ("`vce_select'"=="robust")      	 local vce_rdbw = "hc1"
	if ("`vce_select'"=="hc0")      	 local vce_rdbw = "hc0"
	if ("`vce_select'"=="hc1")      	 local vce_rdbw = "hc1"
	if ("`vce_select'"=="hc2")      	 local vce_rdbw = "hc2"
	if ("`vce_select'"=="hc3")      	 local vce_rdbw = "hc3"
	if ("`vce_select'"=="cr1")           local vce_rdbw = "cr1 `clustvar'"
	if ("`vce_select'"=="cr2")           local vce_rdbw = "cr2 `clustvar'"
	if ("`vce_select'"=="cr3")           local vce_rdbw = "cr3 `clustvar'"


	******************** COVS EFF ****************************************************
	local covseff_count : word count `covs_eff'

	******************** Numeric precision for generated work variables **************
	local precision = lower("`precision'")
	if ("`precision'"=="") local precision "double"
	if !inlist("`precision'", "double", "single") {
		di as error "{err}{cmd:precision()} incorrectly specified: use {cmd:precision(double)} or {cmd:precision(single)}"
		exit 198
	}
	local _rdbwhte_dtype "double"
	if ("`precision'"=="single") local _rdbwhte_dtype "float"

    * Create temporary variables for centered running variable and treatment indicator
    tempvar Xc T
	qui gen `_rdbwhte_dtype' `Xc' = `x' - `c'
    qui gen `_rdbwhte_dtype' `T' = `Xc' >= 0
	lab var `T' "T"
	lab var `Xc' "Xc"

	* Drop rows missing on any variable that participates in the downstream
	* rdbwselect call, so e(N) / N_l / N_r reflect the same analytic sample
	* that rdhte uses. Mirrors rdrobust.ado's drop_cond pattern.
	local _drop_extra ""
	if ("`covs_eff'"!="") {
		foreach _v of varlist `covs_eff' {
			local _drop_extra "`_drop_extra' | mi(`_v')"
		}
	}
	if ("`clustvar'"!="") local _drop_extra "`_drop_extra' | mi(`clustvar')"
	if ("`weights'"!="")  local _drop_extra "`_drop_extra' | mi(`weights') | `weights'<=0"
	if ("`covs_hte'"!="") {
		fvrevar `covs_hte', list
		local _hte_bases "`r(varlist)'"
		foreach _v of local _hte_bases {
			local _drop_extra "`_drop_extra' | mi(`_v')"
		}
	}
	if (`"`_drop_extra'"' != "") {
		local _drop_extra = substr(`"`_drop_extra'"', 4, .)
		qui drop if `_drop_extra'
	}

	* Per-side post-NA-drop sample sizes.
	qui count
	local N = r(N)
	qui count if `T'==1
	local N_r = r(N)
	qui count if `T'==0
	local N_l = r(N)

	if ("`q'"=="0") local q = `p'+1

	* Input validation (mirror rdhte.ado)
	if (`p' < 0) {
		di as error "{err}{cmd:p()} must be a non-negative integer"
		exit 125
	}
	if (`q' <= `p' | `q' != round(`q')) {
		di as error "{err}{cmd:q()} must be an integer with q > p"
		exit 125
	}

	***********************************************************
	***** Check for factor variables in covs_hte:
	**********************************************************
	if "`covs_hte'" != "" {

	local w = "`covs_hte'"
	local w_count : word count `w'

	if (`w_count'==1) {
		capture assert `w' == 0 | `w' == 1
		if _rc == 0 {
			* Warn if degenerate -- mirrors the rdhte.ado constant-W check.
			qui sum `w', meanonly
			if (r(min) == r(max)) {
				di as txt "{txt}note: covs_hte() variable `w' is constant ({txt}all == " r(min) "{txt}); heterogeneity analysis collapses to overall RD ATE."
			}
			local w i.`w'
		}
	}

	local d_hash = (length("`w'") - length(subinstr("`w'", "##", "", .)))/2
	local s_hash = (length("`w'") - length(subinstr("`w'", "#", "", .))) - (length("`w'") - length(subinstr("`w'", "##", "", .)))

	rdhte_countvars `w'

	local quad  = 0
	if (`w_count'==1 & `s_hash'==1 & `e(n_ic)'==1) | (`d_hash'==1 & `e(n_c)'==1 & `e(n_ic)'==1) | (`s_hash'==1 & `e(n_c)'==1 & `e(n_ic)'==1)  {
		local quad = 1
	}
	if (`d_hash'==1 & `e(n_f)'==2) {
		di as error  "{err}{cmd:covs_hte()} does not allow factorial specifications (##) between factor variables. Use one (#) instead"
		exit 125
	}
	if (`w_count'==2 & `e(n_c)'==1 & `e(n_f)'==1) {
		local quad = 1
	}


	rdhte_fvexpand `w'

	local is_factor = "`r(fvops)'"

	if (`quad'==1) {
		local is_factor = ""
		local f_list = r(varlist)
		fvexpand i.`T'#c.(`f_list')
		foreach wrd in `r(varlist)' {
			if strpos("`wrd'", "b.") == 0 local result `result' `wrd'
		}
		local f_list =  "1.`T' `result'"
		local f_list2: subinstr local f_list "1.`T'" "T", all
	}
	else if "`is_factor'" == "true" {
			local f_list = r(varlist)
	}
	else {
		fvexpand i.`T'#c.(`w')
		foreach wrd in `r(varlist)' {
			if strpos("`wrd'", "b.") == 0 local result `result' `wrd'
		}
		local f_list =  "1.`T' `result'"
		local f_list2: subinstr local f_list "1.`T'" "T", all
	}

	}   /* close section with covs_hte*/
	else {
		local f_list =  "1.`T' `result'"
		local f_list2: subinstr local f_list "1.`T'" "T", all
	}

	local I : word count `f_list'

	* Tempnames so working matrices don't leak into the caller's
	* namespace. Promoted to user-visible e() at the end.
	tempname tau_h tau_N
	matrix `tau_h' = J(`I', 2, .)
	matrix `tau_N' = J(`I', 2, `N')


	*** Extract labels ***************************
	local lab_list = "`f_list'"

	if ("`w'" != "") & ("`labels'" != "") {

		if (`e(n_f)'==1 & `s_hash'==0 & `d_hash'==0) {

			qui {
			local rawvar `covs_hte'
			if regexm("`covs_hte'", "^(i[bn]?[0-9]*\.)?(.*)$") {
				local rawvar `=regexs(2)'
			}

			local lblname : value label `rawvar'
			if "`lblname'" != "" {
				levelsof `rawvar', local(levels)

				local lab_list
				foreach l of local levels {
					local lab : label `lblname' `l'
					local lab_list `lab_list'  `lab'
				}
			}
			else {
				local lab_list = "`f_list'"
			}
			}

		}
		else if  (`s_hash'>0 & `d_hash'==0) {

			local fvvar = "`covs_hte'"
			// 3. Normalize "##" to "#" and split into part1...part`n'
			local cleaned = subinstr("`fvvar'", "##", "#", .)
			local n = 0
			while strpos("`cleaned'", "#") {
				local ++n
				local part`n' = substr("`cleaned'", 1, strpos("`cleaned'", "#") - 1)
				local cleaned   = substr("`cleaned'", strpos("`cleaned'", "#") + 1, .)
			}
			local ++n
			local part`n' = "`cleaned'"

			// 4. For each part, strip off i./ib./ibn., get its value-label name,
			//    then pull its text-labels (unquoted) into lablist`i'
			forvalues i = 1/`n' {
				if regexm("`part`i''", "^(i[bn]?[0-9]*\.)?(.*)$") {
					local rawvar`i' = "`=regexs(2)'"
				}
				local lblname`i' : value label `rawvar`i''
				quietly levelsof `rawvar`i'', local(levels`i')
				local lablist`i' ""
				foreach l of local levels`i' {
					local lab : label `lblname`i'' `l'
					local lablist`i' `lablist`i'' `lab'
				}
			}

			// Build the cross-product of value-labels into `combos'.
			local combos `lablist1'
			forvalues j = 2/`n' {
				local newcombos ""
				foreach prefix of local combos {
					foreach lbl of local lablist`j' {
						local newcombos "`newcombos' `prefix'#`lbl'"
					}
				}
				local combos "`newcombos'"
			}

			// 6. Save into `interaction_labels' (unquoted tokens)
			local interaction_labels `combos'


			local lab_list = "`interaction_labels'"

		}

	}

		if ("`bwjoint'" ~= "" & "`is_factor'" == "true") {
		qui rdbwselect `y' `x', c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce_rdbw') kernel(`kernel') covs(`covs_eff') weights(`weights')
			tempname _rdbwhte_mat_h
			matrix `_rdbwhte_mat_h' = e(mat_h)
			local _h_l = el(`_rdbwhte_mat_h', 1, 1)
			local _h_r = el(`_rdbwhte_mat_h', 1, 2)

			local i = 1
			foreach fvar in `f_list' {

				matrix `tau_h'[`i', 1] = `_h_l'
				matrix `tau_h'[`i', 2] = `_h_r'

				qui count if `fvar'==1 & abs(`Xc') <= `_h_l' & `T'==0
				matrix `tau_N' [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= `_h_r' & `T'==1
				matrix `tau_N' [`i', 2] = r(N)

				local ++i
            }
		}


	if ("`is_factor'" ~= "true") {
		qui rdbwselect `y' `x', c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce_rdbw') kernel(`kernel') covs(`covs_eff') weights(`weights')
				tempname _rdbwhte_mat_h
				matrix `_rdbwhte_mat_h' = e(mat_h)
				local _h_l = el(`_rdbwhte_mat_h', 1, 1)
				local _h_r = el(`_rdbwhte_mat_h', 1, 2)

				matrix `tau_h'[1, 1] = `_h_l'
				matrix `tau_h'[1, 2] = `_h_r'

				qui count if (abs(`Xc') <= `_h_l' & `T'==0)
				matrix `tau_N' [1, 1] = r(N)
				qui count if (abs(`Xc') <= `_h_r' & `T'==1)
				matrix `tau_N' [1, 2] = r(N)

		}

        if ("`bwjoint'" == "" & "`is_factor'" == "true") {
			local i = 1
			foreach fvar in `f_list' {
				qui rdbwselect `y' `x' if `fvar'==1, c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce_rdbw') kernel(`kernel') covs(`covs_eff') weights(`weights')
				tempname _rdbwhte_mat_h_i
				matrix `_rdbwhte_mat_h_i' = e(mat_h)
				local _h_l_i = el(`_rdbwhte_mat_h_i', 1, 1)
				local _h_r_i = el(`_rdbwhte_mat_h_i', 1, 2)

				matrix `tau_h'[`i', 1] = `_h_l_i'
				matrix `tau_h'[`i', 2] = `_h_r_i'

				qui count if `fvar'==1 & abs(`Xc') <= `_h_l_i' & `T'==0
				matrix `tau_N' [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= `_h_r_i' & `T'==1
				matrix `tau_N' [`i', 2] = r(N)

				local ++i
            }
        }



	local ml = 0
	foreach word of local f_list {
    local len = strlen("`word'")
    if `len' > `ml' {
        local ml = `len'
		}
	}

	local ml = `len' + 3




	************************************************
	********* OUTPUT TABLE *************************
	************************************************

	disp ""
	if "`is_factor'"=="true" {
		disp "Bandwidth estimators for Sharp RD Heterogeneous Treatment Effects: Subgroups."
	}
	else {
		disp "Bandwidth estimators for Sharp RD Heterogeneous Treatment Effects: Continuous."
	}
	disp ""


	* Summary-header column layout (mirrors rdhte.ado). See there for the
	* col_lab / col_bar / col_L / col_Rh / col_R / col_info / info_w
	* semantics.
	local col_lab  = 18
	local col_bar  = 19
	local col_L    = 21
	local col_Rh   = 33
	local col_R    = 34
	local col_info = 55
	local info_w   = 10
	disp in smcl in gr "{ralign `col_lab': Cutoff c = `c'}"      _col(`col_bar') " {c |} " _col(`col_L') in gr "Left of " in yellow "c"  _col(`col_Rh') in gr "Right of " in yellow "c"  _col(`col_info') in gr "Number of obs = "  in yellow %10.0f `N'
	disp in smcl in gr "{hline `col_bar'}{c +}{hline 22}"                                                                                                          _col(`col_info') in gr "BW type       = "  in yellow "{ralign `info_w':`bwselect'}"
	disp in smcl in gr "{ralign `col_lab':Number of obs}"        _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `N_l'  _col(`col_R') %9.0f `N_r'  _col(`col_info') in gr "Kernel        = "  in yellow "{ralign `info_w':`kernel_type'}"
	disp in smcl in gr "{ralign `col_lab':Order est. (p)}"       _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `p'    _col(`col_R') %9.0f `p'
	disp in smcl in gr "{ralign `col_lab':Order bias (q)}"       _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `q'    _col(`col_R') %9.0f `q'


	local ml1 = `ml'  + 1
	local ml2 = `ml1' + 8      /* Est */
	local ml3 = `ml2' + 8     /* Z   */
	local ml4 = `ml3' + 8     /* pv  */
	local ml5 = `ml4' + 10     /* CI  */


	local ml1b = `ml'   + 1
	local ml2b = `ml1b' + 3     /* Est */
	local ml3b = `ml2b' + 6	    /* Z   */
	local ml4b = `ml3b' + 14	/* pv  */
	local ml5b = `ml4b' + 10	/* CIl */


	disp ""
	disp           "Outcome: `y'. Running variable: `x'."

	if ("`is_factor'"=="true") {
		disp in smcl in gr "{hline `ml1'}{c TT}{hline 50}"
		disp in smcl in gr "{ralign `ml':Group}"      _col(`ml1') " {c |} " _col(`ml2')  "Nh-" _col(`ml3')   "Nh+"  _col(`ml4') "h-"  _col(`ml5') "h+"
		disp in smcl in gr "{hline `ml1'}{c +}{hline 50}"
		forval i = 1/`I'  {
			local lab : word `i' of `lab_list'
			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') %8.0f `tau_N'[`i',1] _col(`ml3b') %8.0f `tau_N'[`i',2] _col(`ml4b') %5.3f `tau_h'[`i',1] _col(`ml5b') %5.3f `tau_h'[`i',2]
		}
		disp in smcl in gr "{hline `ml1'}{c BT}{hline 50}"
	}
	else {
			disp in smcl in gr "{hline `ml1'}{c TT}{hline 50}"
			disp in smcl in gr "{ralign `ml':}" _col(`ml1') " {c |} " _col(`ml2')  "Nh-"           _col(`ml3')   "Nh+"           _col(`ml4')  "h-"           _col(`ml5')  "h+"
			disp in smcl in gr "{hline `ml1'}{c +}{hline 50}"
			disp in smcl in gr "{ralign `ml':`covs_hte'}" _col(`ml1b') " {c |} " _col(`ml2b') %8.0f `tau_N'[1,1] _col(`ml3b') %8.0f `tau_N'[1,2] _col(`ml4b') %5.3f `tau_h'[1,1] _col(`ml5b') %5.3f `tau_h'[1,2]
			disp in smcl in gr "{hline `ml1'}{c BT}{hline 50}"
	}

	if ("`covs_eff'"!="")   {
		disp "Covariate-adjusted estimates. Additional covariates included: `covseff_count' "
	}

	if ("`cluster'"=="cluster")   {
		* Bandwidth-selection runs multiple internal rdbwselect calls on
		* different bw-filtered subsets, so there is no single meaningful
		* per-regression cluster count to display (matches rdbwselect.ado
		* in rdrobust, which prints the same shape without a count).
		disp "	 (Std. err. adjusted for clusters in `clustvar')"
	}
	disp ""

	* M1: switch back to caller's frame and drop the work frame.
	* Always-run cleanup wrapped by nobreak/capture above.
	}
	local _rc = _rc
	cwf `_orig_frame'
	frame drop `_rdbwhte_frame'
	if `_rc' error `_rc'
	}

	* Return results
	ereturn clear

	ereturn local runningvar "`x'"
	ereturn local depvar     "`y'"
	ereturn local cmd        "rdbwhte"

	* Per-side effective-sample totals (sum across groups). Mirrors
	* rdhte.ado's e(N_h_l) / e(N_h_r); convenient by analogy with
	* rdrobust's e(N_h_l) / e(N_h_r).
	tempname Nh_l_total Nh_r_total
	matrix `Nh_l_total' = J(1, 1, 0)
	matrix `Nh_r_total' = J(1, 1, 0)
	local _I = rowsof(`tau_N')
	forval _i = 1/`_I' {
		matrix `Nh_l_total' = `Nh_l_total' + `tau_N'[`_i', 1]
		matrix `Nh_r_total' = `Nh_r_total' + `tau_N'[`_i', 2]
	}

	ereturn scalar N         = `N'
	ereturn scalar N_h_l     = `Nh_l_total'[1, 1]
	ereturn scalar N_h_r     = `Nh_r_total'[1, 1]
	ereturn scalar c         = `c'
	ereturn scalar p         = `p'
	ereturn scalar q         = `q'
	ereturn local  kernel    "`kernel_type'"
	ereturn local  bwselect  "`bwselect'"
	ereturn local  precision "`precision'"
	* e(vce_select) is the canonical normalized name (cr1/cr2/cr3/hc0/...),
	* e(vce_type) is the display label. Mirrors rdrobust convention.
	ereturn local  vce_select "`vce_select'"
	ereturn local  vce_type   "`vce_type'"
	if ("`covs_eff'"!="")    ereturn local covs     "`covs_eff'"
	if ("`clustvar'"!="")    ereturn local clustvar "`clustvar'"

	* User-visible LHS names; tempname tau_h / tau_N on the RHS.
	ereturn matrix h     = `tau_h'
	ereturn matrix tau_N = `tau_N'

end
