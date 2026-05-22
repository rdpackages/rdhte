************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
*
* 2026-05-14: M4 clarity (hash counting comments), vce-mapping consolidation
*             (single canonical-name -> (display, regress, rdbw) lookup),
*             and M7 dedup (six bwselect Nh-counting blocks routed through
*             new internal helper rdhte_fill_Nh.ado). Bit-for-bit identical
*             numerical outputs.
* 2026-05-15: Extended the missing-value drop scan so e(N) and per-side
*             N_l/N_r counts reflect the actual analytic sample. Now drops
*             rows missing on covs_eff / covs_hte (via fvrevar ... list) /
*             clustvar / weights, mirroring rdrobust's drop_cond pattern.
*             No-op on data with no missings in those vars (bit-for-bit
*             identical numerical outputs preserved).
*             M6 cosmetic: summary-header magic offsets (18/19/21/33/34/55)
*             hoisted into named col_* locals + ml-ladder header comment.
*             Pure formatting -- no numerical change.
* 2026-05-15 (returns audit): e(vce_select) now holds the canonical raw
*             name (was display label); added e(vce_type) for the display
*             label. Mirrors rdrobust.ado:1271-1272.
* 2026-05-15 (vce(hc0) fix): Stata's -regress- doesn't accept vce(hc0);
*             the previous mapping passed "hc0" through and hit a runtime
*             error at fit time. Now: vce(hc0) without cluster auto-remaps
*             to vce(hc1) with a warning (HC0 and HC1 differ only by the
*             (N-1)/(N-K) df correction). The cr0-without-cluster
*             fallback also remapped from hc0 to hc1.
************************************************************************************************
*!version 0.2.0  2026-05-22

capture program drop rdhte
program define rdhte, eclass
    version 16.0
    syntax varlist(min=2 max=2)  [if] [in] [ , c(real 0) p(integer 1) q(real 0) h(numlist) h_l(numlist) h_r(numlist) covs_hte(string) covs_eff(varlist) kernel(string) weights(string) bwselect(string) vce(string) level(real 95) precision(string) bwjoint labels]


    * Verify rdrobust is installed (anywhere on the adopath -- not just
    * sysdir_plus) and at version 11 or later. Uses findfile + a regex
    * scan of the *! version stamp rather than assuming a fixed line layout.
    *
    capture findfile rdrobust.ado
    if _rc {
        di as error "rdhte requires rdrobust to be installed."
        di as error "Install via: {stata net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace}"
        exit 9
    }
    local rdrobust_path "`r(fn)'"

    tempname fh
    file open `fh' using "`rdrobust_path'", read
    local rdversion ""
    local rd_major = .
    file read `fh' line
    while r(eof) == 0 & "`rdversion'" == "" {
        if regexm(`"`macval(line)'"', "^\*!.*v([0-9]+)\.([0-9]+)\.([0-9]+)") {
            local rdversion = regexs(0)
            local rd_major = real(regexs(1))
        }
        file read `fh' line
    }
    file close `fh'

    if "`rdversion'" == "" | `rd_major' < 11 {
        di as error "rdhte requires rdrobust v11.0.0 or newer (found: `rdversion')."
        di as error "Update via: {stata net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace}"
        exit 9
    }



	marksample touse, novarlist

	* Extract main variables from varlist
	local y : word 1 of `varlist'
	local x : word 2 of `varlist'

	* M1: do all work in an isolated temporary frame instead of
	* preserve+keep+drop on the user's data. Pros vs preserve:
	*   (a) no O(N) disk round-trip;
	*   (b) only the touse=1 subset is copied (typically a fraction);
	*   (c) the user's data + frame state is never mutated, even on error.
	* The nobreak/capture wrap guarantees the user is returned to their
	* original frame even if the work block errors out (mirrors the
	* pattern used in rdrobust.ado as of 2026-05-12).
	local _orig_frame `c(frame)'
	tempname _rdhte_frame
	capture frame drop `_rdhte_frame'
	frame put * if `touse', into(`_rdhte_frame')

	nobreak {
	cwf `_rdhte_frame'
	capture noisily {

	* Initial y/x-completeness drop. Extended below (after vce parse and
	* covs_hte base-var extraction) to also drop on covs_eff / covs_hte /
	* clustvar / weights missings, so e(N) and the per-side counts reflect
	* the actual analytic sample passed to -regress-. Without this, regress
	* silently drops missing rows at fit time but e(N) overstates the count.
	qui drop if mi(`y') | mi(`x')


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
	local user_vce = ("`vce'"!="")

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
		* bare vce(cluster clustvar) is equivalent to cr1.
		if ("`vce_select'"=="cluster") local vce_select = "cr1"
	}

	* cr1/cr2/cr3 require a cluster variable; otherwise error.
	if inlist("`vce_select'","cr1","cr2","cr3") & "`clustvar'"=="" {
		di as error "{err}{cmd:vce(`vce_select' clustervar)} requires a cluster variable"
		exit 125
	}
	* cr0 without cluster: warn + fall back to hc1 (NOT hc0 -- Stata's
	* -regress- does not expose vce(hc0); see hc0 handling below).
	if "`vce_select'"=="cr0" & "`clustvar'"=="" {
		di as text "Warning: vce(cr0) requires a cluster variable. Falling back to vce(hc1)."
		local vce_select = "hc1"
	}

	* Stata's -regress- does not accept vce(hc0): it exposes only
	* vce(robust)=HC1, vce(hc2), vce(hc3) (plus cluster/bootstrap/jackknife).
	* HC0 and HC1 differ only by the (N-1)/(N-K) finite-sample correction,
	* so HC1 is the closest available alternative. Warn and remap so the
	* downstream `regress vce(hc0)' call doesn't error out at fit time.
	if "`vce_select'"=="hc0" & "`clustvar'"=="" {
		di as text "Warning: Stata's regress does not support vce(hc0). Switching to vce(hc1) (HC0 with the (N-1)/(N-K) df correction)."
		local vce_select = "hc1"
	}

	if ("`vce_select'"=="")              local vce_select = "hc3"

	if ("`clustvar'"!="") local cluster = "cluster"

	* Single source of truth for the canonical-name -> (display, regress,
	* rdbwselect) triple. Each row pairs vce_select with three downstream
	* views:
	*   vce_type : user-facing display label (HC1/HC2/HC3/CR1/CR2/CR3)
	*   vce_reg  : argument passed to -regress, vce(...)-; note that
	*              regress does NOT accept "hc0" or "hc1" -- the only
	*              HC variants it exposes are vce(robust)=HC1, vce(hc2)
	*              and vce(hc3). hc0 is therefore remapped to hc1
	*              upstream (see Warning above). cr1/cr2/cr3 are
	*              rdhte-canonical names; -regress- has no direct CR2/CR3
	*              equivalent, so clustered regression fits route through
	*              vce(cluster clustvar).
	*   vce_rdbw : argument forwarded to rdbwselect (rdrobust 4.0.0+
	*              accepts cr1/cr2/cr3 directly).
	* Adding a new option = adding one row here.
	if      ("`vce_select'"=="hc1" | "`vce_select'"=="robust") {
		local vce_type = "HC1"
		local vce_reg  = "robust"
		local vce_rdbw = "hc1"
	}
	else if ("`vce_select'"=="hc2") {
		local vce_type = "HC2"
		local vce_reg  = cond("`clustvar'"!="", "hc2 `clustvar'", "hc2")
		local vce_rdbw = "hc2"
	}
	else if ("`vce_select'"=="hc3") {
		local vce_type = "HC3"
		local vce_reg  = cond("`clustvar'"!="", "hc3 `clustvar'", "hc3")
		local vce_rdbw = "hc3"
	}
	else if ("`vce_select'"=="cr1") {
		local vce_type = "CR1"
		local vce_reg  = "cluster `clustvar'"
		local vce_rdbw = "cr1 `clustvar'"
	}
	else if ("`vce_select'"=="cr2") {
		local vce_type = "CR2"
		local vce_reg  = "cluster `clustvar'"
		local vce_rdbw = "cr2 `clustvar'"
	}
	else if ("`vce_select'"=="cr3") {
		local vce_type = "CR3"
		local vce_reg  = "cluster `clustvar'"
		local vce_rdbw = "cr3 `clustvar'"
	}
	else {
		* Should not be reachable -- upstream validation already rejected
		* unknown tokens. Defensive fallback to the prior default behavior.
		local vce_type = "HC3"
		local vce_reg  = "`vce_select'"
		local vce_rdbw = "hc3"
	}


	******************** COVS EFF ****************************************************
	local covseff_count : word count `covs_eff'

	******************** Numeric precision for generated work variables **************
	local precision = lower("`precision'")
	if ("`precision'"=="") local precision "double"
	if !inlist("`precision'", "double", "single") {
		di as error "{err}{cmd:precision()} incorrectly specified: use {cmd:precision(double)} or {cmd:precision(single)}"
		exit 198
	}
	local _rdhte_dtype "double"
	if ("`precision'"=="single") local _rdhte_dtype "float"

    * Create temporary variables for centered running variable and treatment indicator
    tempvar Xc T
	qui gen `_rdhte_dtype' `Xc' = `x' - `c'
    qui gen `_rdhte_dtype' `T' = `Xc' >= 0
	lab var `T' "T"
	lab var `Xc' "Xc"

    * Construct polynomial terms up to order p/p+1
    local Xp "`Xc'"
    forval i = 2/`p' {
        tempvar Xc`i'
        gen `_rdhte_dtype' `Xc`i'' = `Xc'^`i'
        local Xp "`Xp' `Xc`i''"
    }

	if ("`q'"=="0") local q = `p'+1

	* Input validation
	if (`level' >= 100 | `level' <= 0) {
		di as error "{err}{cmd:level()} should be a number in (0, 100)"
		exit 125
	}
	if (`p' < 0) {
		di as error "{err}{cmd:p()} must be a non-negative integer"
		exit 125
	}
	if (`q' <= `p' | `q' != round(`q')) {
		di as error "{err}{cmd:q()} must be an integer with q > p"
		exit 125
	}

	local Xq "`Xc'"
    forval i = 2/`q' {
        tempvar Xc`i'
        qui gen `_rdhte_dtype' `Xc`i'' = `Xc'^`i'
        local Xq "`Xq' `Xc`i''"
    }

	* Drop rows missing on any variable that will participate in the
	* downstream -regress- call, so e(N) and the per-side N_l / N_r counts
	* below reflect the actual analytic sample. Mirrors rdrobust.ado's
	* unified drop_cond pattern (rdrobust.ado:194-202). Without this,
	* regress silently drops missings on covs_eff / covs_hte / clustvar /
	* weights at fit time but e(N) (computed off the y/x markout) overstates
	* the count.
	local _drop_extra ""
	if ("`covs_eff'"!="") {
		foreach _v of varlist `covs_eff' {
			local _drop_extra "`_drop_extra' | mi(`_v')"
		}
	}
	if ("`clustvar'"!="") local _drop_extra "`_drop_extra' | mi(`clustvar')"
	if ("`weights'"!="")  local _drop_extra "`_drop_extra' | mi(`weights') | `weights'<=0"
	if ("`covs_hte'"!="") {
		* fvrevar ... list returns the underlying base varnames, stripping
		* any i./ib./c./# factor-var decoration. Works for all syntactic
		* forms accepted by covs_hte (bare var, i.W, c.W, W#Z, W##Z, etc.).
		fvrevar `covs_hte', list
		local _hte_bases "`r(varlist)'"
		foreach _v of local _hte_bases {
			local _drop_extra "`_drop_extra' | mi(`_v')"
		}
	}
	if (`"`_drop_extra'"' != "") {
		* Strip the leading " | " and apply.
		local _drop_extra = substr(`"`_drop_extra'"', 4, .)
		qui drop if `_drop_extra'
	}

	* Per-side post-NA-drop sample sizes. _N is the total in-memory count
	* after all markouts. Combined into 2 passes (left + right);
	* the third was redundant since N = N_l + N_r.
	qui count if `T' == 1
	local N_r = r(N)
	qui count if `T' == 0
	local N_l = r(N)
	local N = `N_l' + `N_r'

	***********************************************************
	***** Check for factor variables in covs_hte:
	**********************************************************
	if "`covs_hte'" != "" {

	local w = "`covs_hte'"
	* Strip a leading c. prefix from a single-term continuous covs_hte
	* spec. Otherwise the continuous-W branch builds
	*   fvexpand i.T#c.(c.var)
	* which malforms (double c. prefix) and levelsof gets fv input it
	* does not handle. Factor (i.) prefixes, interactions (# / ##), and
	* bare names are unaffected.
	if (`: word count `w'' == 1) {
		if (regexm("`w'", "^c\.(.*)$")) local w = regexs(1)
	}
	local w_count : word count `w'

	if (`w_count'==1) {
		capture assert `w' == 0 | `w' == 1
		if _rc == 0 {
			* Warn if the var is degenerate (all 0 or all 1) -- the i.
			* prefix expands to a single level and the heterogeneity
			* analysis silently collapses to an overall RD ATE.
			qui sum `w', meanonly
			if (r(min) == r(max)) {
				di as txt "{txt}note: covs_hte() variable `w' is constant ({txt}all == " r(min) "{txt}); heterogeneity analysis collapses to overall RD ATE."
			}
			local w i.`w'
		}
	}

	* Count `##` (factorial interactions) and bare `#` (one-way
	* interactions) occurrences in the covs_hte spec.
	*   total_hash = total `#` characters; double_hash counts pairs.
	*   d_hash = number of `##` operators;
	*   s_hash = number of bare `#` operators = total_hash - 2*d_hash.
	local total_hash  = length("`w'") - length(subinstr("`w'", "#",  "", .))
	local double_hash = length("`w'") - length(subinstr("`w'", "##", "", .))
	local d_hash      = `double_hash' / 2
	local s_hash      = `total_hash' - `double_hash'

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


	if (`e(n_f)'>1 & `s_hash'==0 & `d_hash'==0) {
			di as error  "{err}{cmd:covs_hte()} does not allow more than one group specification without interactions. Try removing the i.() or include proper interactions"
			exit 125
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
		* L6: continuous-W branch. Warn if W has very few unique values --
		* the user probably wanted i.`w' (subgroup CATE) but is getting
		* the slope-coefficient model interpreted on a discrete variable.
		if (`w_count' == 1) {
			qui levelsof `w', missing
			local _nlev : word count `r(levels)'
			if (`_nlev' < 10) {
				di as txt "{txt}note: covs_hte() variable {bf:`w'} is " ///
				          "being treated as continuous but has only " ///
				          "`_nlev' unique values. " ///
				          "Use {bf:i.`w'} for subgroup-CATE estimation " ///
				          "instead of a slope coefficient."
			}
		}
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
	* Use tempnames so the working matrices don't leak into the caller's
	* namespace. They are promoted to permanent ereturn matrices at the end.
	tempname tau_h tau_hat tau_bc tau_se tau_V tau_t tau_pv tau_N
	tempname tau_h_l tau_h_r tau_ci_l tau_ci_r CV b_eret V_eret
	matrix `tau_h'   = J(`I', 2, .)
	matrix `tau_hat' = J(`I', 1, .)
	matrix `tau_bc'  = J(`I', 1, .)
	matrix `tau_se'  = J(`I', 1, .)
	matrix `tau_V'   = J(`I', `I', .)
	matrix `tau_t'   = J(`I', 1, .)
	matrix `tau_pv'  = J(`I', 1, .)
	matrix `tau_N'   = J(`I', 2, `N')

	*********************************** Extract labels ***************************
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

			local interaction_labels `combos'
			local lab_list = "`interaction_labels'"

		}

	}



************************************************************************************************************
****** Bandwidth Selection  ****
************************************************************************************************************

	****** Bwselect
	if ("`h'" == "" & "`h_l'" == "" & "`h_r'" == "") {
		tempvar h
		qui gen `_rdhte_dtype' `h' = .

		if ("`bwselect'"=="") local bwselect = "mserd"

		if ("`bwjoint'" ~= "" | "`is_factor'" ~= "true") {
			qui rdbwselect `y' `x', c(`c') p(`p') q(`q') vce(`vce_rdbw') kernel(`kernel') bwselect(`bwselect') covs(`covs_eff') weights(`weights')
			tempname _rdhte_mat_h
			matrix `_rdhte_mat_h' = e(mat_h)
			local _h_l = el(`_rdhte_mat_h', 1, 1)
			local _h_r = el(`_rdhte_mat_h', 1, 2)
			qui replace `h' = `_h_l' if `T'==0
			qui replace `h' = `_h_r' if `T'==1

			* Fill tau_h (rows = groups, cols = L/R). Joint case: every
			* row gets the same h. Counting of tau_N is delegated to
			* the helper.
			if ("`is_factor'" == "true") {
				local i = 1
				foreach fvar in `f_list' {
					matrix `tau_h'[`i', 1] = `_h_l'
					matrix `tau_h'[`i', 2] = `_h_r'
					local ++i
				}
			}
			else {
				* Continuous case: all parameters share one bandwidth.
				* Replicate across every row so e(tau_h) is well-formed.
				forvalues _i = 1/`I' {
					matrix `tau_h'[`_i', 1] = `_h_l'
					matrix `tau_h'[`_i', 2] = `_h_r'
				}
			}
			* When is_factor is "true" we pass f_list to the helper for
			* per-group counting. The continuous case omits f_list because
			* it contains factor-variable terms (e.g. 1.T#c.var) that
			* Stata's option-value parser chokes on when threaded through
			* `f_list(...)` in syntax(string). The helper does not use
			* f_list in the continuous (is_factor != "true") branch --
			* it counts once with row 1's h and replicates across all
			* rows of tau_N (all parameters share one bandwidth).
			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
		}


        if ("`bwjoint'" == "" & "`is_factor'" == "true") {
			* Per-cell bandwidth selection: one rdbwselect call per
			* factor level. Fill tau_h and h.vec in the loop; defer
			* tau_N counting to the helper at the end.
			local i = 1
			foreach fvar in `f_list' {
				qui rdbwselect `y' `x' if `fvar'==1, c(`c') p(`p') q(`q') vce(`vce_rdbw') kernel(`kernel') bwselect(`bwselect') covs(`covs_eff') weights(`weights')
				tempname _rdhte_mat_h_i
				matrix `_rdhte_mat_h_i' = e(mat_h)
				local _h_l_i = el(`_rdhte_mat_h_i', 1, 1)
				local _h_r_i = el(`_rdhte_mat_h_i', 1, 2)
				qui replace `h' = `_h_l_i' if `fvar'==1 & `T'==0
				qui replace `h' = `_h_r_i' if `fvar'==1 & `T'==1

				matrix `tau_h'[`i', 1] = `_h_l_i'
				matrix `tau_h'[`i', 2] = `_h_r_i'

				local ++i
			}
			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
        }
	}
	else {
		local bwselect = "Manual"		/* Manual bwselect*/

		if ("`h_l'"=="" & "`h_r'"=="") { /* same h on each side */

		local count: word count `h'

		* one common h
		if (`count'==1) {
			matrix `tau_h' = J(`I', 2, `h')
			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
		}
		else {   		/* different h by groups*/
			if (`count'!=`I') {
				di as error "{err}{cmd:h()} incorrectly specified"
				exit 125
			}

			local i = 1
			foreach v of local h {
				matrix `tau_h'[`i', 1] = `v'
				matrix `tau_h'[`i', 2] = `v'
				local ++i
			}

			tempvar h
			qui gen `_rdhte_dtype' `h' = .
			if ("`is_factor'" == "true") {
				local i = 1
				foreach fvar in `f_list' {
					local _h_i = el(`tau_h', `i', 1)
					qui replace `h' = `_h_i' if `fvar'==1
					local ++i
				}
			}
			else {
				local _h_1 = el(`tau_h', 1, 1)
				qui replace `h' = `_h_1'
			}
			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
		}
	}
	else {  /* different h on each side */

		local count_l: word count `h_l'
		local count_r: word count `h_r'

		if (`count_l'+`count_r' ==2) {
			matrix `tau_h_l' = J(`I', 1, `h_l')
			matrix `tau_h_r' = J(`I', 1, `h_r')
			matrix `tau_h' = `tau_h_l', `tau_h_r'

			tempvar h
			qui gen `_rdhte_dtype' `h' = .
			qui replace `h' = `h_l' if `T'==0
			qui replace `h' = `h_r' if `T'==1

			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
		}

		else {
			if (`count_l'!=`I' | `count_r'!=`I') {
				di as error "{err}{cmd:h_l()} and/or {cmd:h_r()} incorrectly specified"
				exit 125
			}

			local i = 1
			foreach v of local h_l {
				matrix `tau_h'[`i', 1] = `v'
				local ++i
			}
			local i = 1
			foreach v of local h_r {
				matrix `tau_h'[`i', 2] = `v'
				local ++i
			}

			tempvar h
			qui gen `_rdhte_dtype' `h' = .
			local i = 1
			foreach fvar in `f_list' {
				local _h_l_i = el(`tau_h', `i', 1)
				local _h_r_i = el(`tau_h', `i', 2)
				qui replace `h' = `_h_l_i' if `fvar'==1 & `T'==0
				qui replace `h' = `_h_r_i' if `fvar'==1 & `T'==1
				local ++i
			}
			if ("`is_factor'" == "true") {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(true) ///
				               xc(`Xc') t(`T') f_list(`f_list')
			}
			else {
				rdhte_fill_Nh, tau_h(`tau_h') tau_N(`tau_N') is_factor(false) ///
				               xc(`Xc') t(`T')
			}
		}


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

	// Weighted Regression
	tempvar kw
	qui gen `_rdhte_dtype' `kw' = .

	// Assign kernel weights based on the selected kernel type
	if ("`kernel_type'" == "Uniform") {
		qui replace `kw' = 1 if abs(`Xc') <= `h'
	}
	if ("`kernel_type'" == "Triangular") {
	    qui replace `kw' = (1 - abs(`Xc') / `h') if abs(`Xc') <= `h'
	}
	if ("`kernel_type'" == "Epanechnikov") {
		qui replace `kw' = 1 - (`Xc' / `h')^2 if abs(`Xc') <= `h'
	}


	if ("`weights'"~="") {
			qui replace `kw' = `kw'*`weights'
		}

	local covs_eff_add ""
	local covs_eff_factor ""
	local covs_eff_cont ""
	if ("`covs_eff'"!="") {
		local covs_eff_add "c.(`covs_eff')"
		local covs_eff_factor "c.(`covs_eff')##(`w')"
		local covs_eff_cont "c.(`covs_eff')##c.(`w')"
	}

	* Create the regression formula depending on the presence of W
    if "`w'" != "" {
		if "`is_factor'"=="true" {

			*di "Running regression with heterogeneity by factor `w'..."

			qui reg `y' c.(`Xp')##`T'##(`w') `covs_eff_factor' [aw = `kw'], vce(`vce_reg') level(`level')
			local i = 1
			foreach fvar in `f_list' {
				local tmp =  _b[1.`T'] + _b[1.`T'#`fvar']
				matrix `tau_hat' [`i', 1] = `tmp'
				local ++i
			}

			qui reg `y' c.(`Xq')##`T'##(`w') `covs_eff_factor' [aw = `kw'], vce(`vce_reg') level(`level')
			mat `CV' = e(V)

			local i = 1
			foreach fvar in `f_list' {
				local ind1 = "1.`T'"
				local ind2 = "1.`T'#`fvar'"

				local tmp1 =  _b[`ind1'] + _b[`ind2']
				local tmp2 = `CV'["1.`T'","1.`T'"] + `CV'["1.`T'#`fvar'","1.`T'#`fvar'"] + 2*`CV'["1.`T'", "1.`T'#`fvar'"]
				local tmp3 = `CV'["1.`T'","1.`T'"] + `CV'["1.`T'", "1.`T'#`fvar'"]

				matrix `tau_bc'[`i', 1] = `tmp1'
				matrix `tau_se'[`i', 1] = sqrt(`tmp2')
				matrix `tau_V'[`i', `i'] = `tmp2'
				matrix `tau_t'[`i', 1]  = `tmp1'/sqrt(`tmp2')

				matrix `tau_pv'[`i', 1] = 2*normal(-abs(`tau_t'[`i', 1]))

				if (`i'>1 & `i' <= `I') {
					matrix `tau_V' [1, `i'] = `tmp3'
					matrix `tau_V' [`i', 1] = `tmp3'
				}
				local ++i
			}


			*** Covariances
			if (`I' > 2) {
			tokenize "`f_list'"
			forvalues i = 2/`I' {
				local k = `i' + 1
				forvalues j = `k'/`I' {
					local lev_i : word `i' of `f_list'
					local lev_j : word `j' of `f_list'
					local tmp3 = `CV'["1.`T'","1.`T'"] + `CV'["1.`T'", "1.`T'#`lev_i'"] +  `CV'["1.`T'", "1.`T'#`lev_j'"] + `CV'["1.`T'#`lev_i'", "1.`T'#`lev_j'"]
					matrix `tau_V' [`i', `j'] = `tmp3'
					matrix `tau_V' [`j', `i'] = `tmp3'
				}
			}
			}

		}
		else{
			*di "Running regression with heterogeneity by `w'..."

			qui reg `y' c.(`Xp')##`T'##c.(`w') `covs_eff_cont' [aw = `kw'], vce(`vce_reg') level(`level')

			local i = 1
			foreach fvar in `f_list' {
				local tmp =  _b[`fvar']
				matrix `tau_hat' [`i', 1] = `tmp'
				local ++i
			}

			qui reg `y' c.(`Xq')##`T'##c.(`w') `covs_eff_cont' [aw = `kw'], vce(`vce_reg') level(`level')
			mat `CV' = e(V)

			local i = 1
			foreach fvar in `f_list' {
				local tmp1 = _b[`fvar']
				local tmp2 = _se[`fvar']
				matrix `tau_bc'[`i', 1] = `tmp1'
				matrix `tau_se'[`i', 1] = `tmp2'
				matrix `tau_t' [`i', 1] = `tmp1' / `tmp2'
				matrix `tau_pv'[`i', 1] = 2 * normal(-abs(`tau_t'[`i', 1]))
				local ++i
			}

			*** Covariances
			tokenize "`f_list'"
			forvalues i = 1/`I' {
				forvalues j = `i'/`I' {
					local lev_i : word `i' of `f_list'
					local lev_j : word `j' of `f_list'
					local tmp3 = `CV'["`lev_i'", "`lev_j'"]
					matrix `tau_V'[`i', `j'] = `tmp3'
					matrix `tau_V'[`j', `i'] = `tmp3'
				}
			}
		}
    }
    else {
        * No-heterogeneity branch (rdhte y x [, ...] with no covs_hte).
        qui reg `y' c.(`Xp')##`T' `covs_eff_add' [aw = `kw'], vce(`vce_reg') level(`level')
        local tmp = _b[1.`T']
        matrix `tau_hat'[1, 1] = `tmp'

        qui reg `y' c.(`Xq')##`T' `covs_eff_add' [aw = `kw'], vce(`vce_reg') level(`level')
        local tmp1 = _b[1.`T']
        local tmp2 = _se[1.`T']
        matrix `tau_bc' [1, 1] = `tmp1'
        matrix `tau_se' [1, 1] = `tmp2'
        matrix `tau_t'  [1, 1] = `tmp1' / `tmp2'
        matrix `tau_pv' [1, 1] = 2 * normal(-abs(`tau_t'[1, 1]))
        * Populate tau_V (1x1) -- previously left as . which made
        * `ereturn post' fail when its `cap' was removed by H4.
        matrix `tau_V'  [1, 1] = `tmp2'^2
    }

	************************************************
	********* OUTPUT TABLE *************************
	************************************************
	* Use the user-requested confidence level rather than hardcoding 1.96.
	local zq = invnormal(1 - (1 - `level'/100)/2)
	mat `tau_ci_l' = `tau_bc' - `zq'*`tau_se'
	mat `tau_ci_r' = `tau_bc' + `zq'*`tau_se'


	disp ""
	if "`covs_hte'" != "" {
		if "`is_factor'"=="true" {
			disp "Sharp RD Heterogeneous Treatment Effects: Subgroups."
		}
		else {
			disp "Sharp RD Heterogeneous Treatment Effects: Continuous."
		}
	}
	else {
			disp "Sharp RD Average Treatment Effect."
	}
	disp ""
	* Summary-header column layout. Pulled out of the `disp` chain so the
	* magic offsets only have to be tweaked in one place.
	*   col_lab  : right-align width of the row-label column (left of "|")
	*   col_bar  : column where the vertical "|" separator sits
	*   col_L    : column where the LEFT-of-cutoff value begins
	*   col_Rh   : column where the RIGHT-of-cutoff header label begins
	*   col_R    : column where the RIGHT-of-cutoff value begins
	*   col_info : column where the right-side info pairs begin
	*   info_w   : right-side `ralign` width for kernel / bw / vce labels
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

	if "`is_factor'"!="true" {
		disp in smcl in gr "{ralign `col_lab':Eff. Number of obs}" _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `tau_N'[1,1]  _col(`col_R') %9.0f `tau_N'[1,2]  _col(`col_info') in gr "VCE method    = "  in yellow "{ralign `info_w':`vce_type'}"
		disp in smcl in gr "{ralign `col_lab':Order est. (p)}"     _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `p'           _col(`col_R') %9.0f `p'
		disp in smcl in gr "{ralign `col_lab':Order bias (q)}"     _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `q'           _col(`col_R') %9.0f `q'
		disp in smcl in gr "{ralign `col_lab':BW est. (h)}"        _col(`col_bar') " {c |} " _col(`col_L') as result %9.3f `tau_h'[1,1]  _col(`col_R') %9.3f `tau_h'[1,2]
	}
	else {
		disp in smcl in gr "{ralign `col_lab':Order est. (p)}"     _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `p'           _col(`col_R') %9.0f `p'  _col(`col_info') in gr "VCE method    = "  in yellow "{ralign `info_w':`vce_type'}"
		disp in smcl in gr "{ralign `col_lab':Order bias (q)}"     _col(`col_bar') " {c |} " _col(`col_L') as result %9.0f `q'           _col(`col_R') %9.0f `q'
	}
	disp ""

	* Outcome-table column ladders. Each local is a cumulative `_col()`
	* offset from `ml` (the longest label-token width plus padding).
	*  - `ml1..ml9`   : header row (combined "[level% Conf. Interval]"
	*                   spans one column at ml5 -> ml6 jumps by 25)
	*  - `ml1b..ml10b`: data rows (CI lower / CI upper split into two
	*                   columns at ml5b and ml6b; extra Nh / h columns
	*                   shifted accordingly)
	* Adjusting a column width = bumping the SINGLE delta on that row.
	local ml1 = `ml'  + 2
	local ml2 = `ml1' + 3      /* Est */
	local ml3 = `ml2' + 10     /* Z   */
	local ml4 = `ml3' + 14     /* pv  */
	local ml5 = `ml4' + 10     /* CI  */
	local ml6 = `ml5' + 25     /* Nh- */
	local ml7 = `ml6' + 8      /* Nh+ */
	local ml8 = `ml7' + 10     /* h-  */
	local ml9 = `ml8' + 9      /* h+  */

	local ml1b = `ml'  + 2
	local ml2b = `ml1b' + 3     /* Est */
	local ml3b = `ml2b' + 10    /* Z   */
	local ml4b = `ml3b' + 14    /* pv  */
	local ml5b = `ml4b' + 12    /* CIl */
	local ml6b = `ml5b' + 10    /* CIr */
	local ml7b = `ml6b' + 8     /* Nh- */
	local ml8b = `ml7b' + 8     /* Nh+ */
	local ml9b = `ml8b' + 14    /* h-  */
	local ml10b = `ml9b' + 9    /* h+  */

	disp ""
	disp           "Outcome: `y'. Running variable: `x'."

	if "`is_factor'"=="true" {
		disp in smcl in gr "{hline `ml1'}{c TT}{hline 100}"

		disp in smcl in gr "{ralign `ml':}"      _col(`ml1') " {c |} " _col(`ml2') "Point"    _col(`ml3') " {c |} " "Robust Inference"
		disp in smcl in gr "{ralign `ml': `w'}"  _col(`ml1') " {c |} " _col(`ml2') "Estimate" _col(`ml3') " {c |} " "z-stat"  _col(`ml4') "P>|z|"  _col(`ml5')  `"[`level'% Conf. Interval]"' _col(`ml6')  "Nh-"   _col(`ml7')  "Nh+" _col(`ml8') "h-" _col(`ml9') "h+"

		disp in smcl in gr "{hline `ml1'}{c +}{hline 100}"

		forval i = 1/`I'  {
			local lab : word `i' of `lab_list'

			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') in ye %5.3f scalar(`tau_hat'[`i',1]) _col(`ml3b')  " {c |} " %5.3f scalar(`tau_t'[`i',1]) _col(`ml4b')  %5.3f scalar(`tau_pv'[`i',1]) _col(`ml5b') %5.3f scalar(`tau_ci_l'[`i',1]) _col(`ml6b') %5.3f scalar(`tau_ci_r'[`i',1]) _col(`ml7b') %8.0f `tau_N'[`i',1] _col(`ml8b') %8.0f `tau_N'[`i',2] _col(`ml9b') %5.3f `tau_h'[`i',1] _col(`ml10b') %5.3f `tau_h'[`i',2]

		}
		disp in smcl in gr "{hline `ml1'}{c BT}{hline 100}"

	}
	else {
		disp in smcl in gr "{hline `ml1'}{c TT}{hline 60}"
		disp in smcl in gr "{ralign `ml':}"      _col(`ml1') " {c |} " _col(`ml2') "Point"    _col(`ml3') " {c |} " "Robust Inference"
		disp in smcl in gr "{ralign `ml': `w'}"  _col(`ml1') " {c |} " _col(`ml2') "Estimate" _col(`ml3') " {c |} " "z-stat"  _col(`ml4') "P>|z|"  _col(`ml5')  `"[`level'% Conf. Interval]"' _col(`ml6')
		disp in smcl in gr "{hline `ml1'}{c +}{hline 60}"

		forval i = 1/`I'  {
			local lab : word `i' of `f_list2'
			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') in ye %5.3f scalar(`tau_hat'[`i',1]) _col(`ml3b')  " {c |} " %5.3f scalar(`tau_t'[`i',1]) _col(`ml4b')  %5.3f scalar(`tau_pv'[`i',1]) _col(`ml5b') %5.3f scalar(`tau_ci_l'[`i',1]) _col(`ml6b') %5.3f scalar(`tau_ci_r'[`i',1])
		}
		disp in smcl in gr "{hline `ml1'}{c BT}{hline 60}"
	}

	if ("`covs_eff'"!="")   {
		disp "Covariate-adjusted estimates. Additional covariates included: `covseff_count' "
	}

	if ("`cluster'"=="cluster")   {
		disp "	 (Std. err. adjusted for `e(N_clust)' clusters in `clustvar')"
	}
	disp ""

	* M1: switch back to the caller's frame and drop the work frame
	* before populating e() against the user's sample. Always-run
	* cleanup wrapped by nobreak/capture above.
	}
	local _rc = _rc
	cwf `_orig_frame'
	frame drop `_rdhte_frame'
	if `_rc' error `_rc'
	}

    * Return results
    ereturn clear

	* Build the e(b) / e(V) tuple from BC point estimate + robust V.
	matrix `b_eret' = `tau_bc''
	matrix `V_eret' = `tau_V'

	* Transpose tau_hat and tau_bc to row vectors for ereturn (consistent
	* with how callers index e(tau_hat)[1, k]).
	matrix `tau_hat' = `tau_hat''
	matrix `tau_bc'  = `tau_bc''

	if "`is_factor'"=="true" {
		matrix colnames `b_eret' = `f_list'
		matrix rownames `V_eret' = `f_list'
		matrix colnames `V_eret' = `f_list'

		matrix rownames `tau_se' = `f_list'
		matrix rownames `tau_V'  = `f_list'
		matrix rownames `tau_t'  = `f_list'
		matrix rownames `tau_pv' = `f_list'
		matrix rownames `tau_N'  = `f_list'
	}
	else {
		matrix colnames `b_eret' = `f_list2'
		matrix rownames `V_eret' = `f_list2'
		matrix colnames `V_eret' = `f_list2'
	}

	* `cap` was previously masking real conformability errors; with the
	* tempname migration the post should always succeed. Let it surface.
	ereturn post `b_eret' `V_eret', esample(`touse')

	ereturn local runningvar "`x'"
	ereturn local depvar     "`y'"
	ereturn local cmd        "rdhte"

	* Per-side effective sample sizes (sum over groups, within bandwidth).
	* Convenient for users by analogy with rdrobust's e(N_h_l)/e(N_h_r).
	tempname Nh_l_total Nh_r_total
	matrix `Nh_l_total' = J(1, 1, 0)
	matrix `Nh_r_total' = J(1, 1, 0)
	forval _i = 1/`I' {
		matrix `Nh_l_total' = `Nh_l_total' + `tau_N'[`_i', 1]
		matrix `Nh_r_total' = `Nh_r_total' + `tau_N'[`_i', 2]
	}

	ereturn scalar N         = `N'
	ereturn scalar N_h_l     = `Nh_l_total'[1, 1]
	ereturn scalar N_h_r     = `Nh_r_total'[1, 1]
	ereturn scalar c         = `c'
	ereturn scalar p         = `p'
	ereturn scalar q         = `q'
	ereturn scalar level     = `level'
	ereturn local  kernel    "`kernel_type'"
	ereturn local  bwselect  "`bwselect'"
	ereturn local  precision "`precision'"
	* e(vce_select) is the canonical normalized name (cr1/cr2/cr3/hc0/...),
	* e(vce_type) is the display label (CR1/CR2/CR3/HC0/...). Mirrors
	* rdrobust.ado:1271-1272.
	ereturn local  vce_select "`vce_select'"
	ereturn local  vce_type   "`vce_type'"
	if ("`covs_eff'"!="")    ereturn local covs     "`covs_eff'"
	if ("`clustvar'"!="")    ereturn local clustvar "`clustvar'"
	* Plotter-friendly extras: discriminator for categorical vs continuous
	* covs_hte and the prettified group labels (when -labels- was specified).
	ereturn local  is_factor "`is_factor'"
	ereturn local  lab_list  "`lab_list'"
	if ("`covs_hte'"!="")    ereturn local covs_hte "`covs_hte'"

	* User-visible e() matrix names (LHS literal) bound to the local
	* tempname matrices (RHS).
	ereturn matrix h         = `tau_h'
	ereturn matrix tau_hat   = `tau_hat'
	ereturn matrix tau_bc    = `tau_bc'
	ereturn matrix tau_se    = `tau_se'
	ereturn matrix tau_V     = `tau_V'
	ereturn matrix tau_t     = `tau_t'
	ereturn matrix tau_pv    = `tau_pv'
	ereturn matrix tau_N     = `tau_N'
	ereturn matrix tau_ci_lb = `tau_ci_l'
	ereturn matrix tau_ci_ub = `tau_ci_r'

end
