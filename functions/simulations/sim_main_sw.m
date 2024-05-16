% sim_main_sw.m   Setup and run specifications for the Smets-Wouters
%                 exercises
% Note: "Specifications" contains the following fields:
%  - spec: worst (worst-case), estp (estimated p), fixp (fixed p)
%  - boot: Run bootstrap. 
%  - longT: true: Run T=2000. false: Run T=240.

%% SETUP

exercise          = 'varma';
fields            = {'spec', 'boot', 'longT'};
Specifications    = cell2struct(cell(length(fields), 1), fields);
Specifications(1) = struct('spec', 'worst', 'boot', false, 'longT', false);
Specifications(2) = struct('spec', 'estp',  'boot', false, 'longT', false);
Specifications(3) = struct('spec', 'fixp',  'boot', false, 'longT', false);
Specifications(4) = struct('spec', 'fixp',  'boot', false, 'longT', true);


%% RUN

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