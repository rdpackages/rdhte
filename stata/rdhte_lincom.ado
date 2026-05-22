************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_lincom
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
* 2026-05-15 (later): switched 5 numeric returns from `ereturn local` to
*                     `ereturn scalar` so callers can do arithmetic on them.
*!version 0.2.0  2026-05-15

capture program drop rdhte_lincom
program define rdhte_lincom, eclass
    version 16.0
    syntax anything(name=expr) [, level(real 95)]

	* Inference (BC: e(tau_bc) reposted as e(b), e(V) unchanged from rdhte)
	matrix b = e(tau_bc)
	ereturn repost b = b
	capture qui lincom `expr', level(`level')
	if _rc != 0 {
		di as error "rdhte_lincom: lincom failed for hypothesis: `expr'"
		exit `_rc'
	}
	local ci_l = r(lb)
	local ci_r = r(ub)
	local pv   = r(p)

	* Estimation (conventional: repost tau_hat as e(b); SE/CI come from BC pass)
	matrix b = e(tau_hat)
	ereturn repost b = b
	capture qui lincom `expr', level(`level')
	if _rc != 0 {
		di as error "rdhte_lincom: lincom failed for hypothesis: `expr'"
		exit `_rc'
	}
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
	
	* Return results. These are NUMERIC values, so stored as scalars
	* (so callers can do `e(rdhte_lincom_est) + 1`); the previous
	* `ereturn local` stored them as strings, breaking arithmetic use.
	ereturn scalar rdhte_lincom_est   = `est'
	ereturn scalar rdhte_lincom_lb    = `ci_l'
	ereturn scalar rdhte_lincom_ub    = `ci_r'
	ereturn scalar rdhte_lincom_pv    = `pv'
	ereturn scalar rdhte_lincom_level = `level'

end
