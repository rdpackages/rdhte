************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_fvexpand 
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.2.0  2026-05-15

capture program drop rdhte_fvexpand
program define rdhte_fvexpand, rclass
    version 16.0
    local fvbase = c(fvbase)
    local rc 0
    nobreak {
        set fvbase off
        capture noisily fvexpand `0'
        local rc = _rc
        if `rc' == 0 {
            * Capture r() values inside nobreak before any reset clears them.
            local _varlist `"`r(varlist)'"'
            local _fvops   `"`r(fvops)'"'
        }
        set fvbase `fvbase'
    }
    if `rc' exit `rc'
    return local varlist `"`_varlist'"'
    return local fvops   `"`_fvops'"'
end
