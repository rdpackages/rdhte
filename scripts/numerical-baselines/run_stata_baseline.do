version 15
clear all
set more off

args repo_root output
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
if `"`output'"' == "" {
    local output "`repo_root'/docs/audit/baselines/stata-current.json"
}
local output : subinstr local output "\" "/", all

capture mkdir "`repo_root'/docs"
capture mkdir "`repo_root'/docs/audit"
capture mkdir "`repo_root'/docs/audit/baselines"

adopath ++ "`repo_root'/stata"

program define write_matrix_array
    args handle field matname comma
    matrix M = `matname'
    local rows = rowsof(M)
    local cols = colsof(M)
    file write `handle' `"        "`field'": ["'
    forvalues i = 1/`rows' {
        if `i' > 1 file write `handle' `","'
        file write `handle' `"["'
        forvalues j = 1/`cols' {
            if `j' > 1 file write `handle' `","'
            if missing(M[`i', `j']) {
                file write `handle' `"null"'
            }
            else {
                file write `handle' %24.16e (M[`i', `j'])
            }
        }
        file write `handle' `"]"'
    }
    file write `handle' `"]`comma'"' _n
end

program define write_rdhte_case
    args handle casename comma
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "rdhte": {"' _n
    write_matrix_array `handle' tau_hat e(tau_hat) ","
    write_matrix_array `handle' tau_bc e(tau_bc) ","
    write_matrix_array `handle' tau_se e(tau_se) ","
    write_matrix_array `handle' tau_t e(tau_t) ","
    write_matrix_array `handle' tau_pv e(tau_pv) ","
    write_matrix_array `handle' tau_ci_lb e(tau_ci_lb) ","
    write_matrix_array `handle' tau_ci_ub e(tau_ci_ub) ","
    write_matrix_array `handle' tau_n e(tau_N) ","
    write_matrix_array `handle' h e(h) ","
    write_matrix_array `handle' tau_v e(tau_V) ","
    file write `handle' `"        "bwselect": "`e(bwselect)'""' _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

program define write_rdbwhte_case
    args handle casename comma
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "rdbwhte": {"' _n
    write_matrix_array `handle' tau_n e(tau_N) ","
    write_matrix_array `handle' h e(h) ","
    file write `handle' `"        "bwselect": "`e(bwselect)'""' _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

program define write_lincom_case
    args handle casename comma
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "lincom": {"' _n
    file write `handle' `"        "estimate": "' %24.16e (`e(rdhte_lincom_est)') `","' _n
    file write `handle' `"        "conf_low": "' %24.16e (`e(rdhte_lincom_lb)') `","' _n
    file write `handle' `"        "conf_high": "' %24.16e (`e(rdhte_lincom_ub)') `","' _n
    file write `handle' `"        "p_value": "' %24.16e (`e(rdhte_lincom_pv)') _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

file open jout using "`output'", write replace
file write jout `"{"' _n
file write jout `"  "schema_version": 1,"' _n
file write jout `"  "package": "rdhte","' _n
file write jout `"  "language": "stata","' _n
file write jout `"  "source": "working-tree","' _n
file write jout `"  "timestamp_utc": null,"' _n
file write jout `"  "environment": {"stata_version": "`c(stata_version)'", "platform": "`c(os)'"}, "' _n
file write jout `"  "cases": {"' _n

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, covs_hte(w_left) vce(hc2 cluster_var)
write_rdhte_case jout binary_left ""

rdhte_lincom 1.w_left - 0.w_left
write_lincom_case jout binary_left_lincom ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, covs_hte(w_left) vce(hc2 cluster_var) bwjoint
write_rdhte_case jout binary_left_joint ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdbwhte y x, covs_hte(w_left) vce(cluster cluster_var)
write_rdbwhte_case jout binary_left_bw ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, covs_hte(i.w_ideology) vce(hc2 cluster_var)
write_rdhte_case jout categorical_ideology ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, covs_hte(w_strength) kernel(uni) vce(hc2 cluster_var)
write_rdhte_case jout continuous_strength ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, covs_hte(i.w_left##c.w_strength) h(0.1) vce(hc2 cluster_var)
write_rdhte_case jout interaction_strength ","

use "`repo_root'/stata/rdhte_dataset.dta", clear
rdhte y x, h(0.1) vce(hc3)
write_rdhte_case jout average_manual ","

file write jout `"  }"' _n
file write jout `"}"' _n
file close jout

display as text "Wrote `output'"
