************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_countvars
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.2.0  2026-05-15


capture program drop rdhte_countvars
program define rdhte_countvars, eclass
    version 16.0
    syntax varlist(fv)

    // Initialize counters
    local n_continuous = 0
    local n_factor = 0
    local n_interaction_factor = 0
    local n_interaction_continuous = 0
    
    // Loop over each term in varlist. Each term is classified into
    // exactly ONE bucket (the original chain double-counted i.W#i.Z and
    // higher-order factor-factor interactions because the catch-all
    // `(^|#)[a-zA-Z0-9_]+(#|$)` regex matched even after the
    // factor-factor branch already incremented).
    foreach v of local varlist {

        if regexm("`v'", "#") {
            // Interaction term: route to ONE of the two interaction buckets.
            if regexm("`v'", "(c\..*#|#c\..*|c\..*#c\..*)") {
                // Any c.var participating in the interaction => continuous
                local ++n_interaction_continuous
            }
            else if regexm("`v'", "i\..*#i\..*") {
                // Pure factor#factor (possibly higher-order)
                local ++n_interaction_factor
            }
            else {
                // Unprefixed terms inside an interaction (e.g. T#1.W) -- treat as continuous
                local ++n_interaction_continuous
            }
        }
        else {
            // Single (non-interaction) term.
            if regexm("`v'", "^(i|ib|ibn)[0-9]*\..*") {
                local ++n_factor
            }
            else if regexm("`v'", "^c\..*") {
                local ++n_continuous
            }
            else {
                // Unprefixed bare varname => continuous by default
                local ++n_continuous
            }
        }
    }

    // Return results to eclass
    ereturn clear
    ereturn local n_c  `n_continuous'
    ereturn local n_f  `n_factor'
    ereturn local n_if `n_interaction_factor' 
    ereturn local n_ic `n_interaction_continuous'
end
