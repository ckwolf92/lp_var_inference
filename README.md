# lp_var_inference

Matlab code for inference on impulse responses using Local Projections or VARs in local-to-VAR models.

**Reference:**
[Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2024), "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

Tested in: Matlab R2023b on Macbook Pro 2023

## Contents

**[documents](documents):** paper and supplement
- [lp_varithmetic.pdf](documents/lp_varithmetic.pdf): main paper
- [lp_varithmetic_supplement.pdf](documents/lp_varithmetic_supplement.pdf): online supplement

**[functions](functions):** auxiliary routines
- [estimation](functions/estimation): functions for LP and VAR estimation
- [utilities](functions/utilities): utility functions

**[simulations](simulations):** simulation studies
- [arma](simulations/arma): map arma into local-to-AR representation ([get_arma.m](simulations/arma/inputs/get_arma.m), [report_arma.m](simulations/arma/inputs/report_arma.m)) + run simulations ([run_arma.m](simulations/arma/run_arma.m)) + generate figures ([figures_arma.m](simulations/arma/figures_arma.m))
- [sw](simulations/sw_lshock): map Smets-Wouters model into local-to-VAR representation ([get_varma_sw.m](simulations/sw_lshock/inputs/get_varma_sw.m), [report_varma_sw.m](simulations/sw_lshock/inputs/report_varma_sw.m)) + run simulations ([run_varma_sw.m](simulations/sw_lshock/run_varma.m)) + generate figures ([figures_varma_sw.m](simulations/sw_lshock/figures_varma.m)), for three different variable selections and identification schemes
- [aux](simulations/aux): auxiliary simulation-specific routines

## Acknowledgements
Our estimation routines build on the replication materials for [Montiel-Olea and Plagborg-Møller (2021)](https://github.com/jm4474/Lag-augmented_LocalProjections).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).
