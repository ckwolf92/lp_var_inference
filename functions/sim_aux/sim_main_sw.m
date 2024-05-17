% Run specifications for the Smets-Wouters exercises

for jSpec = 1:length(Specifications)

    % Setup
    clearvars -except jSpec Specifications resp_ind innov_ind horzs exercise
    load('inputs/varma_sw_dgps')  % Load dynare results

    spec  = Specifications(jSpec).spec;
    boot  = Specifications(jSpec).boot;
    longT = Specifications(jSpec).longT;

    sim_setup_general  % General setup
    sim_setup_sw       % sw setup
    dgps = dgp.ps;     % dgps considered by sw
    sim_run            % Run simulations
end