# lp_var_inference

MATLAB code for inference on impulse responses using Local Projections or VARs.

**Reference:**
[Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2024), "Double Robustness of Local Projections and Some Unpleasant VARithmetic"

Tested in: MATLAB R2023b on Macbook Pro 2023 (M3 Pro)

## Contents

**[documents](documents):** Paper and supplement
- [lp_varithmetic.pdf](documents/lp_varithmetic.pdf): Main paper
- [lp_varithmetic_supplement.pdf](documents/lp_varithmetic_supplement.pdf): Online supplement

**[analytics](analytics):** Figures for analytical results
- [plot_worstcase.m](analytics/plot_worstcase.m): Generate figures

**[emp_ses](emp_ses):** Compute empirical standard error ratios based on Ramey (2016)
- [run_ramey_ses.m](emp_ses/run_ramey_ses.m): Compute LP and VAR SEs

**[functions](functions):** Auxiliary routines
- [estimation](functions/estimation): Functions for LP and VAR estimation and inference
- [utilities](functions/utilities): Utility functions

**[simulations](simulations):** Simulation studies
- [arma](simulations/arma): Map ARMA into local-to-AR representation ([get_arma.m](simulations/arma/inputs/get_arma.m), [report_arma.m](simulations/arma/inputs/report_arma.m)) + run simulations ([run_arma.m](simulations/arma/run_arma.m)) + generate figures ([figures_arma.m](simulations/arma/figures_arma.m))
- [sw](simulations/sw): Map Smets-Wouters model into local-to-VAR representation ([get_varma_sw.m](simulations/sw/inputs/get_varma_sw.m), [report_varma_sw.m](simulations/sw/inputs/report_varma_sw.m)) + run simulations ([run_varma_sw.m](simulations/sw/run_varma_sw.m)) + generate figures ([figures_varma_sw.m](simulations/sw/figures_varma_sw.m)), for three different variable selections and identification schemes
- [auxiliary_functions](simulations/auxiliary_functions): Auxiliary simulation-specific routines

## Instructions

To replicate the figures for the analytical results, first run [run_ramey_ses.m](emp_ses/run_ramey_ses.m) to generate the empirical LP and VAR standard errors, and then run [plot_worstcase.m](analytics/plot_worstcase.m) to get the figures.

To replicate the main simulation results reported in the paper it is always necessary to proceed in three steps: first, generate the DGP, in the inputs subfolder (e.g., [get_arma.m](simulations/arma/inputs/get_arma.m)); second, run the simulations (e.g., [run_arma.m](simulations/arma/run_arma.m)); and third, generate the figures (e.g., [figures_arma.m](simulations/arma/figures_arma.m)). Supplementary results for the total amount of misspecification are also generated in the inputs subfolder (e.g., [report_arma.m](simulations/arma/inputs/report_arma.m)).

Detailed instructions for the main text figures and tables follow.

- _Figures 4.1-4.3 and A.2-A.3_. Run the file [plot_worstcase.m](analytics/plot_worstcase.m).
- _Table 5.1_, _Figure 5.1_. Run the file [report_arma.m](simulations/arma/inputs/report_arma.m).
- _Figure 5.2_. First run the file [get_arma.m](simulations/arma/inputs/get_arma.m). Then execute [run_arma.m](simulations/arma/run_arma.m) and finally [figures_arma.m](simulations/arma/figures_arma.m). Note that the option "boot" needs to be set to "true" to replicate the full set of results reported in the paper.
- _Table 5.2_. Run the file [report_varma_sw.m](simulations/sw/inputs/report_varma_sw.m). The shock identification scheme needs to be set to "lshock".
- _Figure 5.3_. First run the file [get_varma_sw.m](simulations/sw/inputs/get_varma_sw.m). Then execute [run_varma_sw.m](simulations/sw/run_varma_sw.m) and finally [figures_varma_sw.m](simulations/sw/figures_varma_sw.m). The shock identification scheme needs to be set to "lshock", the option "boot" in the various specifications needs to be set to "true", and the relevant DGP type needs to be selected in the figures file.

## Acknowledgements
Our estimation routines build on the replication materials for [Montiel-Olea and Plagborg-Møller (2021)](https://github.com/jm4474/Lag-augmented_LocalProjections).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).
