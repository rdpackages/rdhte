************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_countvars
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.1.0  2025-06-27 


capture program drop rdhte_countvars
program define rdhte_countvars, eclass
    syntax varlist(fv)

    // Initialize counters
    local n_continuous = 0
    local n_factor = 0
    local n_interaction_factor = 0
    local n_interaction_continuous = 0
    
    // Loop over each term in varlist
    foreach v of local varlist {

        // Check if it's an interaction
        if regexm("`v'", "#") {
            // Count factor-factor interaction
            if regexm("`v'", "i\..*#i\..*") {
                local ++n_interaction_factor
            }
            // Count interactions involving continuous variables
            if regexm("`v'", "(c\..*#|#c\..*)") | regexm("`v'", "c\..*#c\..*") {
                local ++n_interaction_continuous
            }
            // Catch unprefixed interactions like W3#i.W2 or i.W2#W3
            else if regexm("`v'", "(^|#)[a-zA-Z0-9_]+(#|$)") {
                local ++n_interaction_continuous
            }
        }
        // Not an interaction term
        else {
            if regexm("`v'", "^(i|ib|ibn)[0-9]*\..*") {
                local ++n_factor
            }
            else if regexm("`v'", "^c\..*") {
                local ++n_continuous
            }
            else {
                // unprefixed outside interaction → continuous
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
