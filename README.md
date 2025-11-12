# lp_var_inference

MATLAB code for inference on impulse responses using Local Projections or VARs.

**Reference:**
[Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2025), "Double Robustness of Local Projections and Some Unpleasant VARithmetic" [(arXiv)](https://arxiv.org/abs/2405.09509)

Tested in: MATLAB R2025b on Macbook Pro 2023 (M3 Pro)

## Contents

**[documents](documents):** Replication materials
- [lp_varithmetic_replication.pdf](documents/lp_varithmetic_replication.pdf): Supplementary information for replication of our empirical literature surveys

**[analytics](analytics):** Figures for analytical results
- [plot_arbias.m](analytics/plot_arbias.m): Plot the least-favorable MA polynomial for the local-to-AR(1) model
- [plot_worstcase.m](analytics/plot_worstcase.m): Generate the various analytical worst-case figures

**[emp_ses](emp_ses):** Compute empirical standard error ratios based on Ramey (2016)
- [run_ramey_ses.m](emp_ses/run_ramey_ses.m): Compute LP and VAR SEs
- [plot_ramey_ses.m](emp_ses/plot_ramey_ses.m): Generate supplementary figure that shows standard error ratios by horizon

**[estimation](estimation):** Functions for LP and VAR estimation and inference
- [ir_estim.m](estimation/ir_estim.m): Main function for estimates and confidence intervals

**[simulations](simulations):** Simulation study
- [auxiliary_functions](simulations/auxiliary_functions): Auxiliary simulation-specific routines
- [data](simulations/data): Data from the oil shock application of [Känzig (2020)](https://www.aeaweb.org/articles?id=10.1257/aer.20190964)
- [oil](simulations/oil): Maps oil VAR with 18 lags into a local-to-VAR representation ([gen_varma_oil.m](simulations/oil/gen_varma_oil.m)) + run simulations ([simul_oil.m](simulations/oil/simul_oil.m)) + generate figures ([plot_oil.m](simulations/oil/plot_oil.m))


## Detailed replication instructions

- _Figures 4.1-4.3 and A.2-A.3_. First run the file [run_ramey_ses.m](emp_ses/run_ramey_ses.m). Then execute [plot_worstcase.m](analytics/plot_worstcase.m).
- _Figures 5.1 and D.1_. First run [gen_varma_oil.m](simulations/oil/gen_varma_oil.m). Then execute [simul_oil.m](simulations/oil/simul_oil.m) after setting the option "boot" to "true" for experiments "exp_id=1" and "exp_id=2." The relevant experiment type needs to be selected in [plot_oil.m](simulations/oil/plot_oil.m).
- _Figure A.1_. Run the file [plot_arbias.m](analytics/plot_arbias.m).



## Acknowledgements
Our estimation routines build on the replication materials for [Montiel-Olea and Plagborg-Møller (2021)](https://github.com/jm4474/Lag-augmented_LocalProjections).

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049), and Wolf does the same for [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).
