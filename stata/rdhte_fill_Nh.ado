************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_fill_Nh
* Internal helper. Fills the I x 2 effective-sample-size matrix tau_N
* with counts of observations on each side of the cutoff that fall
* within the row's bandwidths. Used by rdhte's bandwidth-selection
* block to consolidate six near-duplicate count loops into one call.
*
* Args:
*   tau_h(name)     name of an I x 2 matrix of (h_left, h_right) per row
*   tau_N(name)     name of an I x 2 matrix to be filled (left N, right N)
*   is_factor(str)  "true" if the rdhte fit uses categorical covs_hte;
*                   anything else falls back to the single-row continuous case
*   f_list(str)     Stata word-list of binary group indicators, one per row
*                   of tau_h. Only used when is_factor == "true".
*   xc(varname)     centered running variable
*   t(varname)      0/1 treatment indicator (T = X >= 0 internally)
*
* For each row i:
*   tau_N[i, 1] = #{ obs : group_i & |xc| <= tau_h[i, 1] & T == 0 }
*   tau_N[i, 2] = #{ obs : group_i & |xc| <= tau_h[i, 2] & T == 1 }
* Continuous case: all parameters share one bandwidth, so count once
* using row 1's h and replicate across every row of tau_N.
************************************************************************************************
*!version 0.2.0 2026-05-15

capture program drop rdhte_fill_Nh
program define rdhte_fill_Nh
    version 16.0
    syntax , TAU_h(name) TAU_N(name) IS_factor(string) ///
             XC(varname) T(varname) [F_list(string)]

    if ("`is_factor'" == "true") {
        local i = 1
        foreach fvar in `f_list' {
            local _h_l = el(`tau_h', `i', 1)
            local _h_r = el(`tau_h', `i', 2)
            qui count if `fvar'==1 & abs(`xc') <= `_h_l' & `t'==0
            matrix `tau_N'[`i', 1] = r(N)
            qui count if `fvar'==1 & abs(`xc') <= `_h_r' & `t'==1
            matrix `tau_N'[`i', 2] = r(N)
            local ++i
        }
    }
    else {
        local _h_l = el(`tau_h', 1, 1)
        local _h_r = el(`tau_h', 1, 2)
        qui count if abs(`xc') <= `_h_l' & `t'==0
        local _nl = r(N)
        qui count if abs(`xc') <= `_h_r' & `t'==1
        local _nr = r(N)
        local _I = rowsof(`tau_N')
        forvalues i = 1/`_I' {
            matrix `tau_N'[`i', 1] = `_nl'
            matrix `tau_N'[`i', 2] = `_nr'
        }
    }
end
