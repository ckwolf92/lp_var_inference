%% LP vs VAR INFERENCE: LIMIT SIMULATIONS
% this version: 01/16/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath('../auxiliary_functions')

%% DATA-GENERATING PROCESS (DGP)

dgp = struct;

% parameters

dgp.rhos      = [0.7]; % AR(1) parameters to consider
dgp.T         = 1e5;  % sample size
dgp.T_scale   = 240; % artificial sample size to scale mis-specification
dgp.sigma     = 1; % shock volatility
dgp.zeta      = [1/3]; % mis-specification dampening
dgp.alpha_raw = [0 0.6 0.4 0.25 0.1 -0.1 0.05]; % mis-specification MA, to be scaled later

dgp.shares    = [0,0.1,0.25]; % targeted volatility share of mis-specification term

% system size

dgp.n_s   = length(dgp.alpha_raw);
dgp.n_y   = 1;
dgp.n_eps = 1;

%% SETTINGS

settings = struct;

%----------------------------------------------------------------
% Monte Carlo Simulation
%----------------------------------------------------------------

settings.sim.numrep      = 5e3; % no. of repetitions
settings.sim.ng_seed     = 202007251; % random number seed
settings.sim.num_workers = 4; % no. of parallel workers (=0: run serial)

%----------------------------------------------------------------
% Estimation
%----------------------------------------------------------------

settings.est.type = 'limit'; % estimation methods type

settings.est.p         = 1; % lag length used for estimation
settings.est.horzs     = [1 3 6 9 12]; % horizons of interest
settings.est.no_const  = true; % true: omit intercept
settings.est.se_homosk = true; % true: homoskedastic ses
settings.est.alpha     = 0.1; % significance level

%% SPECIFICATIONS

specs = cell(2,1);           % specifications for the simulations

specs{1} = {'estimator','var',...
            'lag_aug',false,...
            'bias_corr',false};
specs{2} = {'estimator','lp',...
            'lag_aug',true,...
            'bias_corr',false};

%% PRELIMINARIES

rng(settings.sim.ng_seed, 'twister'); % set RNG seed

aux1  = repmat(dgp.rhos,size(dgp.shares));
aux2  = repmat(dgp.shares',size(dgp.rhos))';
dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
clear aux1 aux2

numdgp  = size(dgps,2);                 % no. of DGPs
numhorz = length(settings.est.horzs);   % no. of horizons
numspec = length(specs);                % no. of regression specifications
numrep  = settings.sim.numrep;          % no. of repetitions

spec_shared = {'alpha', settings.est.alpha,...
               'no_const', settings.est.no_const,...
               'se_homosk', settings.est.se_homosk};

dgp.irs_true = nan(numdgp,numhorz);

%% RUN SIMULATIONS

%----------------------------------------------------------------
% Placeholder Matrices
%----------------------------------------------------------------

estims = zeros(numdgp, numspec, numhorz, numrep);
ses    = estims;

cis_lower = zeros(numdgp, numspec, numhorz, numrep);
cis_upper = cis_lower;

%----------------------------------------------------------------
% Start Parallel Workers
%----------------------------------------------------------------

if settings.sim.num_workers > 0
    poolobj = parpool(settings.sim.num_workers);
end

%----------------------------------------------------------------
% Simulations
%----------------------------------------------------------------

for i_dgp = 1:numdgp

disp(['I am currently at DGP ' num2str(i_dgp)])

% set up DGP

dgp.rho   = dgps(1,i_dgp);
dgp.share = dgps(2,i_dgp);

dgp.scale = sqrt(dgp.share/(1-dgp.share) * 1./(dgp.T_scale.^(-2*dgp.zeta) .* sum(dgp.alpha_raw.^2)));
dgp.alpha = dgp.scale * dgp.alpha_raw;

set_up_ar1

% seed

i_rand_seeds = randi(2^32-1,1,numrep);

% Monte Carlo iterations

parfor i_rep = 1:numrep

% seed

rng(i_rand_seeds(i_rep), 'twister');

% simulate data

data_y = generate_data(dgp);

% placeholders

i_rep_estims    = zeros(numspec, numhorz);
i_rep_ses       = i_rep_estims;
i_rep_cis_lower = nan(numspec, numhorz);
i_rep_cis_upper = i_rep_cis_lower;

% run specifications

for i_spec = 1:numspec

[i_rep_estims(i_spec,:),i_rep_ses(i_spec,:),i_cis_dm] ...
    = ir_estim(data_y, settings.est.p, settings.est.horzs, spec_shared{:}, specs{i_spec}{:});

i_rep_cis_lower(i_spec,:,1) = i_cis_dm(1,:);
i_rep_cis_upper(i_spec,:,1) = i_cis_dm(2,:);

end

% store results

estims(i_dgp,:,:,i_rep) = i_rep_estims;
ses(i_dgp,:,:,i_rep) = i_rep_ses;
cis_lower(i_dgp,:,:,i_rep) = i_rep_cis_lower;
cis_upper(i_dgp,:,:,i_rep) = i_rep_cis_upper;

end

end

clear data_y i_cis_dm i_dgp i_rand_seeds i_rep i_rep_cis_lower i_rep_cis_upper i_rep_estims i_rep_ses ...
    i_spec num2dgp numhorz numrep numspec numdgp spec_shared

%----------------------------------------------------------------
% Store Results
%----------------------------------------------------------------

dgp.dgps = dgps;
clear dgps

results = struct;

results.estims    = estims;
results.ses       = ses;
results.cis_lower = cis_lower;
results.cis_upper = cis_upper;

clear estims ses cis_lower cis_upper

%----------------------------------------------------------------
% Post-Computations
%----------------------------------------------------------------

% coverage

irs_true_reshape = permute(dgp.irs_true, [1 3 2]);

results.cover_inds = (irs_true_reshape >= results.cis_lower ...
                        & irs_true_reshape <= results.cis_upper) + 0; % coverage indicator

results.coverage_prob = mean(results.cover_inds,4); % coverage probability

clear irs_true_reshape

% median length

results.lengths       = results.cis_upper-results.cis_lower;
results.median_length = median(results.lengths,4);

%% SAVE RESULTS

status = mkdir('results');
results_filename = sprintf('%s%s', 'sim_ar1_', settings.est.type);
save(strcat('results/',results_filename, '.mat'), 'dgp', 'specs', 'settings', 'results');