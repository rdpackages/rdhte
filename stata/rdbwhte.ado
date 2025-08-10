************************************************************************************************
* RDHTE STATA PACKAGE -- rdbwhte
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.1.1  2025-08-08

capture program drop rdbwhte
program define rdbwhte, eclass
    syntax varlist(min=2 max=2)  [if] [in] [ , c(real 0) p(integer 1) q(real 0) covs_hte(string) bwselect(string) covs_eff(varlist) kernel(string) vce(string) cluster(string) bwjoint labels]
	
	marksample touse, novarlist
	preserve
    qui keep if `touse'
	
    * Extract main variables from varlist
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'		
	
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
		if ("`vce'"=="") {
			local vce = "hc3"
			local vce_type = "HC3"
		}
		
	tokenize `vce'	
	local count : word count `vce'
	if `count' == 1 {
		local vce_select `"`1'"'
	}
	if `count' == 2 {
		local vce_select `"`1'"'
		if ("`vce_select'"=="cluster") local clustvar `"`2'"'	
		if ("`vce_select'"!="cluster") di as error  "{err}{cmd:vce()} incorrectly specified" 
	}
	if `count' > 2 {
		di as error "{err}{cmd:vce()} incorrectly specified"  
		*exit 125
	}
	
	if ("`vce_select'"=="hc0")     		 local vce_type = "HC0"
	if ("`vce_select'"=="hc1")      	 local vce_type = "HC1"
	if ("`vce_select'"=="hc2")      	 local vce_type = "HC2"
	if ("`vce_select'"=="hc3")      	 local vce_type = "HC3"
	if ("`vce_select'"=="cluster")  	 local vce_type = "Cluster"
	
	if ("`vce_select'"=="cluster")       local cluster = "cluster"
	if ("`vce_select'"=="cluster")       local vce_select = "hc0"
	if ("`vce_select'"=="")              local vce_select = "hc3"

	
	******************** COVS EFF ****************************************************
	local covseff_count : word count `covs_eff'
	
    * Create temporary variables for centered running variable and treatment indicator
    tempvar Xc T
	qui gen `Xc' = `x' - `c'  
    qui gen `T' = `Xc' >= 0   
	lab var `T' "T"
	lab var `Xc' "Xc"
	
	* Handle missing values
	qui count
	local N = r(N)
	qui count if `T'==1
	local N_r = r(N)
	qui count if `T'==0
	local N_l = r(N)
	
	if ("`q'"=="0") local q = `p'+1
	
	***********************************************************
	***** Check for factor variables in covs_hte:
	**********************************************************
	if "`covs_hte'" != "" {
	
	local w = "`covs_hte'"
	local w_count : word count `w'

	if (`w_count'==1) {
		capture assert `w' == 0 | `w' == 1
		if _rc == 0 {
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

	matrix tau_h   = J(`I', 2, .)
	matrix tau_N   = J(`I', 2, `N')
	
	
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
			// 3. Normalize "##"→"#" and split into part1…part`n'
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

				// ★ check that it exists
				*if "`lblname`i''" == "" {
				*	di as error "Variable `rawvar`i'' has NO value label -- aborting."
				*	exit 198
				*}

			// 5. Build the cross-product into `combos' (still unquoted)
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
		qui rdbwselect `y' `x', c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce') kernel(`kernel') covs(`covs_eff')
				
			local i = 1
			foreach fvar in `f_list' {
								
				matrix tau_h[`i', 1] = e(mat_h)[1,1]
				matrix tau_h[`i', 2] = e(mat_h)[1,2]
							
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,1] & `T'==0
				matrix tau_N [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,2] & `T'==1
				matrix tau_N [`i', 2] = r(N)				
				
				local ++i
            }				
		}
		

	if ("`is_factor'" ~= "true") {
		qui rdbwselect `y' `x', c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce') kernel(`kernel') covs(`covs_eff')
				
				matrix tau_h[1, 1] = e(mat_h)[1,1]
				matrix tau_h[1, 2] = e(mat_h)[1,2]
				
				qui count if (abs(`Xc') <= e(mat_h)[1,1] & `T'==0)
				matrix tau_N [1, 1] = r(N)
				qui count if (abs(`Xc') <= e(mat_h)[1,2] & `T'==1)
				matrix tau_N [1, 2] = r(N)	
							
		}
				
        if ("`bwjoint'" == "" & "`is_factor'" == "true") {
			local i = 1
			foreach fvar in `f_list' {
				qui rdbwselect `y' `x' if `fvar'==1, c(`c') p(`p') q(`q') bwselect(`bwselect') vce(`vce') kernel(`kernel') covs(`covs_eff')
				
				matrix tau_h[`i', 1] = e(mat_h)[1,1]
				matrix tau_h[`i', 2] = e(mat_h)[1,2]
				
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,1] & `T'==0
				matrix tau_N [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,2] & `T'==1
				matrix tau_N [`i', 2] = r(N)
				
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

	
	disp in smcl in gr "{ralign 18: Cutoff c = `c'}"        _col(19) " {c |} " _col(21) in gr "Left of " in yellow "c"  _col(33) in gr "Right of " in yellow "c"         _col(55) in gr "Number of obs = "  in yellow %10.0f `N'
	disp in smcl in gr "{hline 19}{c +}{hline 22}"                                                                                                                       _col(55) in gr "BW type       = "  in yellow "{ralign 10:`bwselect'}" 
	disp in smcl in gr "{ralign 18:Number of obs}"          _col(19) " {c |} " _col(21) as result %9.0f `N_l'           _col(34) %9.0f  `N_r'                            _col(55) in gr "Kernel        = "  in yellow "{ralign 10:`kernel_type'}" 
	disp in smcl in gr "{ralign 18:Order est. (p)}"         _col(19) " {c |} " _col(21) as result %9.0f `p'            _col(34) %9.0f  `p'         
	disp in smcl in gr "{ralign 18:Order bias (q)}"         _col(19) " {c |} " _col(21) as result %9.0f `q'            _col(34) %9.0f  `q'         
	
   
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
			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') %8.0f tau_N[`i',1] _col(`ml3b') %8.0f tau_N[`i',2] _col(`ml4b') %5.3f tau_h[`i',1] _col(`ml5b') %5.3f tau_h[`i',2]
		}
		disp in smcl in gr "{hline `ml1'}{c BT}{hline 50}"		
	}
	else {
			disp in smcl in gr "{hline `ml1'}{c TT}{hline 50}"
			disp in smcl in gr "{ralign `ml':}" _col(`ml1') " {c |} " _col(`ml2')  "Nh-"           _col(`ml3')   "Nh+"           _col(`ml4')  "h-"           _col(`ml5')  "h+" 
			disp in smcl in gr "{hline `ml1'}{c +}{hline 50}"
			disp in smcl in gr "{ralign `ml':`covs_hte'}" _col(`ml1b') " {c |} " _col(`ml2b') %8.0f tau_N[1,1] _col(`ml3b') %8.0f tau_N[1,2] _col(`ml4b') %5.3f tau_h[1,1] _col(`ml5b') %5.3f tau_h[1,2]
			disp in smcl in gr "{hline `ml1'}{c BT}{hline 50}"		
	}	
	
	if ("`covs_eff'"!="")   {
		disp "Covariate-adjusted estimates. Additional covariates included: `covseff_count' "
	}
	
	if ("`cluster'"=="cluster")   {
		disp "	 (Std. err. adjusted for `e(N_clust)' clusters in `clustvar')" 
	}
	disp ""
	
	*/
		
	restore
		
    * Return results
    ereturn clear
	
	ereturn local outcomevar "`y'"
	ereturn local runningvar "`x'"
	ereturn local depvar "`y'"
	ereturn local cmd "rdhte"
	
	*ereturn scalar N = `N'
	*ereturn scalar c = `c'
	*ereturn scalar p = `p'
	*ereturn scalar level   = `level'
	*ereturn local kernel     = "`Kernel'"
	*ereturn local bwselect   = "`bwselect'"
	*ereturn local vce_select = "`vce'"
	*if ("`covs_eff'"!="")    ereturn local covs "`covs_list'"
	*if ("`cluster'"!="")     ereturn local clustvar "`clustvar'"
		
	ereturn matrix h = tau_h
	ereturn matrix tau_N   = tau_N			

		
end




