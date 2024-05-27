# lp_var_inference

MATLAB code for inference on impulse responses using Local Projections or VARs in local-to-VAR models.

**Reference:**
[Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2024), "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

Tested in: MATLAB R2023b on Macbook Pro 2023 (M3 Pro)

## Contents

**[documents](documents):** Paper and supplement
- [lp_varithmetic.pdf](documents/lp_varithmetic.pdf): Main paper
- [lp_varithmetic_supplement.pdf](documents/lp_varithmetic_supplement.pdf): Online supplement

**[functions](functions):** Auxiliary routines
- [estimation](functions/estimation): Functions for LP and VAR estimation
- [utilities](functions/utilities): Utility functions

**[simulations](simulations):** Simulation studies
- [arma](simulations/arma): Map ARMA into local-to-AR representation ([get_arma.m](simulations/arma/inputs/get_arma.m), [report_arma.m](simulations/arma/inputs/report_arma.m)) + run simulations ([run_arma.m](simulations/arma/run_arma.m)) + generate figures ([figures_arma.m](simulations/arma/figures_arma.m))
- [sw](simulations/sw): Map Smets-Wouters model into local-to-VAR representation ([get_varma_sw.m](simulations/sw/inputs/get_varma_sw.m), [report_varma_sw.m](simulations/sw/inputs/report_varma_sw.m)) + run simulations ([run_varma_sw.m](simulations/sw/run_varma_sw.m)) + generate figures ([figures_varma_sw.m](simulations/sw/figures_varma_sw.m)), for three different variable selections and identification schemes
- [auxiliary_functions](simulations/auxiliary_functions): Auxiliary simulation-specific routines

## Acknowledgements
Our estimation routines build on the replication materials for [Montiel-Olea and Plagborg-Møller (2021)](https://github.com/jm4474/Lag-augmented_LocalProjections).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).
