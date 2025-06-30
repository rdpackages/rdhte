************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_lincom
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.1.0  2025-06-30

capture program drop rdhte_lincom
program define rdhte_lincom, eclass
    version 15
    syntax anything(name=expr)

	* Inference
	matrix b = e(tau_bc)	
	ereturn repost b = b
    qui lincom `expr'
	local ci_l = r(lb)  
	local ci_r = r(ub)  
	local pv = r(p)
	local level = r(level) 
	
	* Estimation
	matrix b = e(tau_hat)
	ereturn repost b = b
	qui lincom `expr'
	local est = r(estimate) 
	matrix b = e(tau_bc)	
	ereturn repost b = b
	
	************************************************
	********* OUTPUT TABLE *************************
	************************************************
		
	disp ""
	disp "RD Heterogeneous Treatment Effects: Linear combinations of parameters" 
	disp ""
    
	local maxlen = 18
	local maxlen1 = `maxlen'  + 1
	local maxlen2 = `maxlen1' + 6
	local maxlen3 = `maxlen2' + 14
	local maxlen4 = `maxlen3' + 12
	local maxlen5 = `maxlen4' + 12

	disp in smcl in gr "{ralign `maxlen':( 1)  `expr' = 0}"  
	disp ""

	disp in smcl in gr "{hline `maxlen1'}{c TT}{hline 60}"
	disp in smcl in gr "{ralign `maxlen':Y1}"  _col(`maxlen1') " {c |} " _col(`maxlen2') "Coefficient"  _col(`maxlen3') `"[`level'% Conf. Interval]"'   _col(`maxlen5') "P>|z|"    
	disp in smcl in gr "{hline `maxlen1'}{c +}{hline 60}"
	disp in smcl in gr "{ralign `maxlen': (1)}" _col(`maxlen1') " {c |} " _col(`maxlen2') in ye %7.0g `est' _col(`maxlen3') %5.3f `ci_l' _col(`maxlen4') %5.3f `ci_r'  _col(`maxlen5') %5.3f `pv'
	disp in smcl in gr "{hline `maxlen1'}{c BT}{hline 60}"
	disp ""
	
	* Return results
    *return clear
	ereturn local rdhte_lincom_est `est'
	ereturn local rdhte_lincom_lb `ci_l'
	ereturn local rdhte_lincom_ub `ci_r'
	ereturn local rdhte_lincom_pv `pv' 
		
end
