# lp_var_inference

MATLAB code for inference on impulse responses using Local Projections or VARs.

**Reference:**
[Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2024), "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

Tested in: MATLAB R2025b on Macbook Pro 2023 (M3 Pro)

## Contents

**[documents](documents):** Paper and supplement
- [lp_varithmetic.pdf](documents/lp_varithmetic.pdf): Main paper
- [lp_varithmetic_supplement.pdf](documents/lp_varithmetic_supplement.pdf): Online supplement

**[analytics](analytics):** Figures for analytical results
- [plot_arbias.m](analytics/plot_arbias.m): Plot the least-favorable MA polynomial for the local-to-AR(1) model
- [plot_worstcase.m](analytics/plot_worstcase.m): Generate the various analytical worst-case figures

**[emp_ses](emp_ses):** Compute empirical standard error ratios based on Ramey (2016)
- [run_ramey_ses.m](emp_ses/run_ramey_ses.m): Compute LP and VAR SEs
- [plot_ramey_ses.m](emp_ses/plot_ramey_ses.m): Generate supplementary figure that shows standard error ratios by horizon

**[estimation](estimation):** Functions for LP and VAR estimation and inference
- [ir_estim.m](estimation/ir_estim.m): Main function for estimates and confidence intervals

**[simulations](simulations):** Simulation studies
- [sw](simulations/sw): Map Smets-Wouters model into local-to-VAR representation ([get_varma_sw.m](simulations/sw/inputs/get_varma_sw.m), [report_varma_sw.m](simulations/sw/inputs/report_varma_sw.m)) + run simulations ([run_varma_sw.m](simulations/sw/run_varma_sw.m)) + generate figures ([figures_varma_sw.m](simulations/sw/figures_varma_sw.m)), for three different variable selections and identification schemes
- [auxiliary_functions](simulations/auxiliary_functions): Auxiliary simulation-specific routines

## Detailed replication instructions

- _Figures 4.1-4.3 and A.2-A.3_. First run the file [run_ramey_ses.m](emp_ses/run_ramey_ses.m). Then execute [plot_worstcase.m](analytics/plot_worstcase.m).
- _Figures 5.1 and D.1_. First run [gen_varma_oil.m](simulations/oil/gen_varma_oil.m). Then execute [simul_oil.m](simulations/oil/simul_oil.m) after setting the option "boot" to "true" for experiments "exp_id=1" and "exp_id=2." The relevant experiment type needs to be selected in the [plot_oil.m](simulations/oil/plot_oil.m) and [plot_oil_rmse.m](simulations/oil/plot_oil_rmse.m).
- _Figure A.1_. Run the file [plot_arbias.m](analytics/plot_arbias.m).



## Acknowledgements
Our estimation routines build on the replication materials for [Montiel-Olea and Plagborg-Møller (2021)](https://github.com/jm4474/Lag-augmented_LocalProjections).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).
