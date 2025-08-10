************************************************************************************************
* RDHTE STATA PACKAGE -- rdhte_fvexpand 
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max Farrell, Filippo Palomba, Rocio Tititunik
************************************************************************************************
*!version 0.1.1  2025-08-08

capture program drop rdhte_fvexpand 
program rdhte_fvexpand // , rclass
    local fvbase = c(fvbase)    
    nobreak {        
        set fvbase off        
        capture noisily fvexpand `0'        
        set fvbase `fvbase'        
    }    
end
