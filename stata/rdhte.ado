************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.1.1  2025-08-08

capture program drop rdhte
program define rdhte, eclass
    syntax varlist(min=2 max=2)  [if] [in] [ , c(real 0) p(integer 1) q(real 0) h(numlist) h_l(numlist) h_r(numlist) covs_hte(string) covs_eff(varlist) kernel(string) weights(string) bwselect(string) vce(string) level(real 95) bwjoint labels]
	
	
    * ---- Build expected rdrobust ado path ----
    local ado_path "`c(sysdir_plus)'r/rdrobust.ado"
    local ado_path = subinstr("`ado_path'", "\", "/", .)
 
    * ---- Confirm file exists ----
    capture confirm file "`ado_path'"
    if _rc {
        di as error "rdhte error: rdrobust not installed at `ado_path'."
        exit 9
    }

	capture file close f
	file open f using "`ado_path'", read
	file read f line    // 1
	file read f line    // 2
	file read f line    // 3
	file read f line    // 4
	file read f line    // 5
	file close f

	* ---- Extract version assuming exact format "*!version <ver>  <date>" ----
    tokenize "`line'"
    local rdversion = "`2'"

	// extract major parts
	tokenize "`rdversion'", parse(".")
	local maj1 = real("`1'")
	
	* ---- Check extracted version ----
    if "`rdversion'" == "" {
        di as error "rdhte error: Could not find rdrobust."
        exit 9
    }
	
	if (`maj1' < 10) {
		di as error "error: rdhte requires rdrobust version 10.0.0 or newer, your version is `rdversion'"
		di as error "update from: net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace"
        exit 9		
	}
	
    * ---- Continue rdhte code  ----

	
	
	marksample touse, novarlist
	preserve
    qui keep if `touse'
	
    * Extract main variables from varlist
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'		
	
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
		if ("`vce_select'"=="cluster" | "`vce_select'"=="hc2") local clustvar `"`2'"'			
		*if ("`vce_select'"!="cluster") di as error  "{err}{cmd:vce()} incorrectly specified" 
	}
	*if `count' > 2 {
		*di as error "{err}{cmd:vce()} incorrectly specified"  
		*exit 125
	*}
	
	*if ("`vce_select'"=="hc0")     		 local vce_type = "HC0"
	if ("`vce_select'"=="robust")      	 local vce_type = "HC1"
	if ("`vce_select'"=="hc2")      	 local vce_type = "HC2"
	if ("`vce_select'"=="hc3")      	 local vce_type = "HC3"
	if ("`vce_select'"=="cluster")  	 local vce_type = "Cluster"
	
	if ("`vce_select'"=="cluster")       local cluster = "cluster"
	if ("`clustvar'"!="")                local cluster = "cluster"
	*if ("`vce_select'"=="cluster")       local vce_select = "hc0"
	if ("`vce_select'"=="")              local vce_select = "hc3"
	
	if ("`vce_select'"=="robust")      	 local vce_rdbw = "hc1"
	if ("`vce_select'"=="hc3")      	 local vce_rdbw = "hc3"

	if ("`vce_select'"=="hc2" & "`clustvar'"=="")      	 local vce_rdbw = "hc2"
	if ("`vce_select'"=="cluster")       local vce_rdbw = "cluster"
	
	if ("`vce_select'"=="hc2" & "`clustvar'"!="")      	 local vce_rdbw = "cluster"
	
	if ("`vce_rdbw'"=="cluster")  	 local vce_rdbw = "`vce_rdbw' `clustvar'"
	*if ("`clustvar'"!="")                local vce_rdbw = "`vce_rdbw' `clustvar'"

	
	******************** COVS EFF ****************************************************
	local covseff_count : word count `covs_eff'
	
    * Create temporary variables for centered running variable and treatment indicator
    tempvar Xc T
	qui gen `Xc' = `x' - `c'  
    qui gen `T' = `Xc' >= 0   
	lab var `T' "T"
	lab var `Xc' "Xc"
	
    * Construct polynomial terms up to order p/p+1
    local Xp "c.`Xc'"
    forval i = 2/`p' {
        tempvar Xc`i'
        gen `Xc`i'' = `Xc'^`i'
        local Xp "`Xp' c.`Xc`i''"
    }
	
	if ("`q'"=="0") local q = `p'+1
	
	local Xq "c.`Xc'"
    forval i = 2/`q' {
        tempvar Xc`i'
        qui gen `Xc`i'' = `Xc'^`i'
        local Xq "`Xq' c.`Xc`i''"
    }
	
    * Handle missing values
	qui count
	local N = r(N)
	qui count if `T'==1
	local N_r = r(N)
	qui count if `T'==0
	local N_l = r(N)
	
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
	matrix tau_hat = J(`I', 1, .)
	matrix tau_bc  = J(`I', 1, .)
	matrix tau_se  = J(`I', 1, .)
	matrix tau_V   = J(`I', `I', .)
	matrix tau_t   = J(`I', 1, .)
	matrix tau_pv  = J(`I', 1, .)
	matrix tau_N   = J(`I', 2, `N')
	
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

			// 7. Display them cleanly
			*di in yellow "All interaction labels (`n`-way):"
			*foreach x of local interaction_labels {
			*    di "`x'"
			*}

			local lab_list = "`interaction_labels'"
		
		}
		
	}	
	

	
************************************************************************************************************
****** Bandwidth Selection  ****	
************************************************************************************************************

	****** Bwselect
	if ("`h'" == "" & "`h_l'" == "" & "`h_r'" == "") {    
		tempvar h
		qui gen `h' = .
			
		if ("`bwselect'"=="") local bwselect = "mserd"

		if ("`bwjoint'" ~= "" | "`is_factor'" ~= "true") {
			qui rdbwselect `y' `x', c(`c') p(`p') q(`q') vce(`vce_rdbw') kernel(`kernel') bwselect(`bwselect') covs(`covs_eff')
						
			qui replace `h' = e(mat_h)[1,1] if `T'==0			
			qui replace `h' = e(mat_h)[1,2] if `T'==1
	

		if ( "`is_factor'" == "true") {
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
		else {
				matrix tau_h[1, 1] = e(mat_h)[1,1]
				matrix tau_h[1, 2] = e(mat_h)[1,2]
				
				qui count if (abs(`Xc') <= e(mat_h)[1,1] & `T'==0)
				matrix tau_N [1, 1] = r(N)
				qui count if (abs(`Xc') <= e(mat_h)[1,2] & `T'==1)
				matrix tau_N [1, 2] = r(N)			
			}
		}
				
					
        if ("`bwjoint'" == "" & "`is_factor'" == "true") {
			

			local i = 1
			foreach fvar in `f_list' {
				qui rdbwselect `y' `x' if `fvar'==1, c(`c') p(`p') q(`q') vce(`vce_rdbw') kernel(`kernel') bwselect(`bwselect') covs(`covs_eff')
				qui replace `h' = e(mat_h)[1,1] if `fvar'==1 & `T'==0
				qui replace `h' = e(mat_h)[1,2] if `fvar'==1 & `T'==1
								
				matrix tau_h[`i', 1] = e(mat_h)[1,1]
				matrix tau_h[`i', 2] = e(mat_h)[1,2]
					
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,1]  & `T'==0
				matrix tau_N [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= e(mat_h)[1,2] & `T'==1
				matrix tau_N [`i', 2] = r(N)
				
				local ++i
            }
        }			
	}
	else {
		local bwselect = "Manual"		/* Manual bwselect*/
		
		if ("`h_l'"=="" & "`h_r'"=="") { /* same h on each side */
				
		local count: word count `h'
		
		* one common h
		if (`count'==1) {
			matrix tau_h = J(`I', 2, `h')
						
						
			if ( "`is_factor'" == "true") {
				local i = 1
				foreach fvar in `f_list' {
					qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 1]  & `T'==0
					matrix tau_N [`i', 1] = r(N)
					qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 1]  & `T'==1
					matrix tau_N [`i', 2] = r(N)
					local ++i
				}
			}
			else {
					qui count if abs(`Xc') <= tau_h[1, 1]  & `T'==0
					matrix tau_N [1, 1] = r(N)
					qui count if abs(`Xc') <= tau_h[1, 2]  & `T'==1
					matrix tau_N [1, 2] = r(N)				
			}
				
		}
		else {   		/* different h by groups*/
  			
			if (`count'!=`I') {
				di as error "{err}{cmd:h()} incorrectly specified"  
				exit 125
			}
			
			
			local i = 1
			foreach v of local h {
				matrix tau_h[`i', 1] =  `v'
				matrix tau_h[`i', 2] =  `v'
				local ++i
			}
		
			tempvar h
			qui gen `h' = .
		
		
		
		if ( "`is_factor'" == "true") {
			local i = 1
			foreach fvar in `f_list' {
				qui replace `h' = tau_h[`i', 1] if `fvar'==1
				
				qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 1] & `T'==0
				matrix tau_N [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 1] & `T'==1
				matrix tau_N [`i', 2] = r(N)
				
				local ++i
            }
		}
		else {
				qui replace `h' = tau_h[1, 1] 
				
				qui count if abs(`Xc') <= tau_h[1, 1] & `T'==0
				matrix tau_N [1, 1] = r(N)
				qui count if abs(`Xc') <= tau_h[1, 1] & `T'==1
				matrix tau_N [1, 2] = r(N)
		}
			
			
			
			
			
		}	
	}
	else {  /* different h on each side */
		
		local count_l: word count `h_l'
		local count_r: word count `h_r'
				
		if (`count_l'+`count_r' ==2) {
			matrix tau_h_l = J(`I', 1, `h_l')
			matrix tau_h_r = J(`I', 1, `h_r')
			matrix tau_h = tau_h_l, tau_h_r
								
			tempvar h
			qui gen `h' = .
			qui replace `h' = `h_l' if `T'==0
			qui replace `h' = `h_r' if `T'==1
					
			
			
			if ( "`is_factor'" == "true") {
			local i = 1
			foreach fvar in `f_list' {
				qui count if (`fvar'==1 & abs(`Xc') <= tau_h[`i', 1] &  `T'==0)
				matrix tau_N [`i', 1] = r(N)
				qui count if (`fvar'==1 & abs(`Xc') <= tau_h[`i', 2]  & `T'==1)
				matrix tau_N [`i', 2] = r(N)
				local ++i
            }
			}
			else {				
				qui count if (abs(`Xc') <= tau_h[1, 1] &  `T'==0)
				matrix tau_N [1, 1] = r(N)
				qui count if (abs(`Xc') <= tau_h[1, 2]  & `T'==1)
				matrix tau_N [1, 2] = r(N)
				
				
				
				
			}
			
			
			
		}
		
		else {
			if (`count_l'!=`I' | `count_r'!=`I') {
				di as error "{err}{cmd:h_l()} and/or {cmd:h_r()} incorrectly specified"  
				exit 125
			}
			
			local i = 1
			foreach v of local h_l {
				matrix tau_h[`i', 1] =  `v'
				local ++i
			}
			local i = 1
			foreach v of local h_r {
				matrix tau_h[`i', 2] =  `v'
				local ++i
			}

			tempvar h
			qui gen `h' = .
			local i = 1
			foreach fvar in `f_list' {
				qui replace `h' = tau_h[`i', 1] if `fvar'==1 & `T'==0
				qui replace `h' = tau_h[`i', 2] if `fvar'==1 & `T'==1

				qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 1] & `T'==0
				matrix tau_N [`i', 1] = r(N)
				qui count if `fvar'==1 & abs(`Xc') <= tau_h[`i', 2] & `T'==1
				matrix tau_N [`i', 2] = r(N)
				
				local ++i
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
	qui gen `kw' = .

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
		
	* Create the regression formula depending on the presence of W
    if "`w'" != "" {
		if "`is_factor'"=="true" {
			
			*di "Running regression with heterogeneity by factor `w'..."
			
			qui reg `y' c.(`Xp')##`T'##(`w') c.(`covs_eff')##(`w') [aw = `kw'], vce(`vce') level(`level') 
			local i = 1
			foreach fvar in `f_list' {
				local tmp =  _b[1.`T'] + _b[1.`T'#`fvar']
				matrix tau_hat [`i', 1] = `tmp'
				local ++i
			}
			
			qui reg `y' c.(`Xq')##`T'##(`w') c.(`covs_eff')##(`w') [aw = `kw'], vce(`vce') level(`level')
			mat CV = e(V)					
				
			local i = 1
			foreach fvar in `f_list' {								
				local ind1 = "1.`T'"
				local ind2 = "1.`T'#`fvar'"
				
				local tmp1 =  _b[`ind1'] + _b[`ind2']				
				local tmp2 = CV["1.`T'","1.`T'"] + CV["1.`T'#`fvar'","1.`T'#`fvar'"] + 2*CV["1.`T'", "1.`T'#`fvar'"]
				local tmp3 = CV["1.`T'","1.`T'"] + CV["1.`T'", "1.`T'#`fvar'"]
				
				matrix tau_bc[`i', 1] = `tmp1'
				matrix tau_se[`i', 1] = sqrt(`tmp2')
				matrix tau_V[`i', `i'] = `tmp2'
				matrix tau_t[`i', 1]  = `tmp1'/sqrt(`tmp2')
			
				matrix tau_pv[`i', 1] = 2*normal(-abs(tau_t[`i', 1]))

				if (`i'>1 & `i' <= `I') {
					matrix tau_V [1, `i'] = `tmp3'
					matrix tau_V [`i', 1] = `tmp3'
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
					local tmp3 = CV["1.`T'","1.`T'"] + CV["1.`T'", "1.`T'#`lev_i'"] +  CV["1.`T'", "1.`T'#`lev_j'"] + CV["1.`T'#`lev_i'", "1.`T'#`lev_j'"]
					matrix tau_V [`i', `j'] = `tmp3'
					matrix tau_V [`j', `i'] = `tmp3'							
				}				
			}
			}	
						
		}
		else{
			*di "Running regression with heterogeneity by `w'..."
			
			qui reg `y' c.(`Xp')##`T'##c.(`w') c.(`covs_eff')##c.(`w') [aw = `kw'], vce(`vce') level(`level') 			
			
			local i = 1
			foreach fvar in `f_list' {
				local tmp =  _b[`fvar']
				matrix tau_hat [`i', 1] = `tmp'
				local ++i
			}
						
			qui reg `y' c.(`Xq')##`T'##c.(`w') c.(`covs_eff')##c.(`w') [aw = `kw'], vce(`vce') level(`level') 
			mat CV = e(V)
						
			local i = 1
			foreach fvar in `f_list' {
				local tmp1 =   _b[`fvar']
				local tmp2 =  _se[`fvar']
				matrix tau_bc[`i', 1] = `tmp1'
				matrix tau_se[`i', 1] = `tmp2'
				matrix tau_t[`i', 1] = `tmp1'/`tmp2'
				matrix tau_pv[`i', 1] = 2*normal(-abs(tau_t[`i', 1]))
				*matrix tau_V  [`i', `i'] = CV[`fvar',`fvar']
				*local j = 1
				*foreach fvar_j in `f_list' {
				*	matrix tau_V [`i', `j'] = CV[`fvar',`fvar_j']
				*	local ++j
				*}
				local ++i
			}			
		
		
		*** Covariances
			tokenize "`f_list'"
			forvalues i = 1/`I' {
				local k = `i' 
				forvalues j = `k'/`I' {					
					local lev_i : word `i' of `f_list'
					local lev_j : word `j' of `f_list'					
					*local tmp3 =  CV["1.`T'#`lev_i'", "1.`T'#`lev_j'"]
					local tmp3 =  CV["`lev_i'", "`lev_j'"]
					matrix tau_V [`i', `j'] = `tmp3'
					matrix tau_V [`j', `i'] = `tmp3'							
				}				
			}
		}		
    }
    else {
        *di "Running regression without heterogeneity..."
	    
        qui reg `y' c.(`Xp')##`T'  c.(`covs_eff') [aw = `kw'], vce(`vce') level(`level')
				local tmp =  _b[1.`T']
				matrix tau_hat [1, 1] = `tmp'
        qui reg `y' c.(`Xq')##`T'  c.(`covs_eff') [aw = `kw'], vce(`vce') level(`level')
				local tmp1 =   _b[1.`T']
				local tmp2 =  _se[1.`T']
				matrix tau_bc [1, 1] = `tmp1'
				matrix tau_se [1, 1] = `tmp2'
				matrix tau_t  [1, 1] = `tmp1'/`tmp2'
				matrix tau_pv [1, 1] = 2*normal(-abs(tau_t  [1, 1]))
    }
		
	************************************************
	********* OUTPUT TABLE *************************
	************************************************
	mat tau_ci_l = tau_bc - 1.96*tau_se
	mat tau_ci_r = tau_bc + 1.96*tau_se

	
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
	disp in smcl in gr "{ralign 18: Cutoff c = `c'}"        _col(19) " {c |} " _col(21) in gr "Left of " in yellow "c"  _col(33) in gr "Right of " in yellow "c"         _col(55) in gr "Number of obs = "  in yellow %10.0f `N'
	disp in smcl in gr "{hline 19}{c +}{hline 22}"                                                                                                                       _col(55) in gr "BW type       = "  in yellow "{ralign 10:`bwselect'}" 
	disp in smcl in gr "{ralign 18:Number of obs}"          _col(19) " {c |} " _col(21) as result %9.0f `N_l'           _col(34) %9.0f  `N_r'                            _col(55) in gr "Kernel        = "  in yellow "{ralign 10:`kernel_type'}" 
		
	if "`is_factor'"!="true" {                            
		disp in smcl in gr "{ralign 18:Eff. Number of obs}"     _col(19) " {c |} " _col(21) as result %9.0f tau_N[1,1]     _col(34) %9.0f  tau_N[1,2]                    _col(55) in gr "VCE method    = "  in yellow "{ralign 10:`vce_type'}" 
		disp in smcl in gr "{ralign 18:Order est. (p)}"         _col(19) " {c |} " _col(21) as result %9.0f `p'            _col(34) %9.0f  `p'         
		disp in smcl in gr "{ralign 18:Order bias (q)}"         _col(19) " {c |} " _col(21) as result %9.0f `q'            _col(34) %9.0f  `q'         
		disp in smcl in gr "{ralign 18:BW est. (h)}"            _col(19) " {c |} " _col(21) as result %9.3f tau_h[1,1]     _col(34) %9.3f  tau_h[1,2]                 
	}
	else {
		disp in smcl in gr "{ralign 18:Order est. (p)}"         _col(19) " {c |} " _col(21) as result %9.0f `p'             _col(34) %9.0f  `p'          _col(55) in gr "VCE method    = "  in yellow "{ralign 10:`vce_type'}" 
		disp in smcl in gr "{ralign 18:Order bias (q)}"         _col(19) " {c |} " _col(21) as result %9.0f `q'             _col(34) %9.0f  `q'          
	}
	disp ""
		
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
	local ml3b = `ml2b' + 10	/* Z   */
	local ml4b = `ml3b' + 14	/* pv  */
	local ml5b = `ml4b' + 12	/* CIl */
	local ml6b = `ml5b' + 10	/* CIr */
	local ml7b = `ml6b' + 8 	/* Nh- */
	local ml8b = `ml7b' + 8 	/* Nh+ */
	local ml9b = `ml8b' + 14 	/* h-  */
	local ml10b = `ml9b' + 9  	/* h+  */

	disp ""		
	disp           "Outcome: `y'. Running variable: `x'."
	
	if "`is_factor'"=="true" {
		disp in smcl in gr "{hline `ml1'}{c TT}{hline 100}"

		disp in smcl in gr "{ralign `ml':}"      _col(`ml1') " {c |} " _col(`ml2') "Point"    _col(`ml3') " {c |} " "Robust Inference"  
		disp in smcl in gr "{ralign `ml': `w'}"  _col(`ml1') " {c |} " _col(`ml2') "Estimate" _col(`ml3') " {c |} " "z-stat"  _col(`ml4') "P>|z|"  _col(`ml5')  `"[`level'% Conf. Interval]"' _col(`ml6')  "Nh-"   _col(`ml7')  "Nh+" _col(`ml8') "h-" _col(`ml9') "h+" 
						
		disp in smcl in gr "{hline `ml1'}{c +}{hline 100}"
	
		forval i = 1/`I'  {
			local lab : word `i' of `lab_list'
			
			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') in ye %5.3f scalar(tau_hat[`i',1]) _col(`ml3b')  " {c |} " %5.3f scalar(tau_t[`i',1]) _col(`ml4b')  %5.3f scalar(tau_pv[`i',1]) _col(`ml5b') %5.3f scalar(tau_ci_l[`i',1]) _col(`ml6b') %5.3f scalar(tau_ci_r[`i',1]) _col(`ml7b') %8.0f tau_N[`i',1] _col(`ml8b') %8.0f tau_N[`i',2] _col(`ml9b') %5.3f tau_h[`i',1] _col(`ml10b') %5.3f tau_h[`i',2]
		
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
			disp in smcl in gr "{ralign `ml':`lab'}" _col(`ml1b') " {c |} " _col(`ml2b') in ye %5.3f scalar(tau_hat[`i',1]) _col(`ml3b')  " {c |} " %5.3f scalar(tau_t[`i',1]) _col(`ml4b')  %5.3f scalar(tau_pv[`i',1]) _col(`ml5b') %5.3f scalar(tau_ci_l[`i',1]) _col(`ml6b') %5.3f scalar(tau_ci_r[`i',1]) 
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
		
	restore
		
    * Return results
    ereturn clear
	
	matrix b = tau_bc'
	matrix V = tau_V
		
	matrix tau_hat = tau_hat'
	matrix tau_bc = tau_bc'
		
	if "`is_factor'"=="true" {
		matrix colnames b = `f_list'
		matrix rownames V = `f_list'
		matrix colnames V = `f_list'
		
		
		*matrix rownames tau_hat = `f_list'
		*matrix rownames tau_bc  = `f_list'
		matrix rownames tau_se  = `f_list'
		matrix rownames tau_V   = `f_list'
		matrix rownames tau_t   = `f_list'
		matrix rownames tau_pv  = `f_list'
		matrix rownames tau_N   = `f_list'
				
	}
	else {
		matrix colnames b = `f_list2'
		matrix rownames V = `f_list2'
		matrix colnames V = `f_list2'
	}

	cap ereturn post b V, esample(`touse')
	
	ereturn local outcomevar "`y'"
	ereturn local runningvar "`x'"
	ereturn local depvar "`y'"
	ereturn local cmd "rdhte"
	
	*ereturn scalar N = `N'
	*ereturn scalar c = `c'
	*ereturn scalar p = `p'
	*ereturn scalar level   = `level'
	*ereturn local kernel     = "`Kernel'"
	ereturn local bwselect   = "`bwselect'"
	*ereturn local vce_select = "`vce'"
	*if ("`covs_eff'"!="")    ereturn local covs "`covs_list'"
	*if ("`cluster'"!="")     ereturn local clustvar "`clustvar'"	
	
	ereturn matrix h = tau_h
	ereturn matrix tau_hat = tau_hat
	ereturn matrix tau_bc  = tau_bc
	ereturn matrix tau_se  = tau_se
	ereturn matrix tau_V   = tau_V
	ereturn matrix tau_t   = tau_t
	ereturn matrix tau_pv  = tau_pv
	ereturn matrix tau_N   = tau_N			
	ereturn matrix tau_ci_lb = tau_ci_l
	ereturn matrix tau_ci_ub = tau_ci_r			
		
end




