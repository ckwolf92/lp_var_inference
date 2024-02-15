%% LP vs VAR INFERENCE: SIMULATIONS
% this version: 01/24/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath(genpath('../../functions'))

%% DATA-GENERATING PROCESS (DGP)

load('inputs/varma_sw_dgps')

dgp = struct;

% parameters

dgp.ps         = dgp_settings.ps(1); % lags
dgp.T          = 240;  % sample size
dgp.T_scale    = dgp_settings.T; % artificial sample size to scale mis-specification
dgp.zeta       = dgp_settings.zeta; % mis-specification dampening
dgp.worst_hor  = 4; % horizon for maximal distortion \alpha(L)

% system size

dgp.n_y   = 3;
dgp.n_eps = dgp.n_y;

%% SETTINGS

settings = struct;

%----------------------------------------------------------------
% Monte Carlo Simulation
%----------------------------------------------------------------

settings.sim.numrep      = 2e3; % no. of repetitions
settings.sim.rng_seed    = 202007252; % random number seed
settings.sim.num_workers = 8; % no. of parallel workers (=0: run serial)

%----------------------------------------------------------------
% Estimation
%----------------------------------------------------------------

settings.est.type = 'worst'; % estimation methods type

settings.est.ps        = dgp.ps; % lag length used for estimation
settings.est.horzs     = [1:1:16]; % horizons of interest
settings.est.no_const  = true; % true: omit intercept
settings.est.se_homosk = true; % true: homoskedastic ses
settings.est.alpha     = 0.1; % significance level

settings.est.p_select  = 1; % use information criterion to select lag length? (1 for AIC, 2 for BIC)
settings.p_max         = 10; % maximal lag length to consider

settings.est.resp_ind  = 1; % index of response variable of interest
settings.est.innov_ind = 3; % index of innovation of interest

settings.est.boot      = false; % true: bootstrap
settings.est.boot_num  = 1e2;  % number of bootstrap samples

%% SPECIFICATIONS

specs = cell(2,1);           % specifications for the simulations

if settings.est.boot

    specs{1} = {'estimator', 'var',...
            'bias_corr', false,...
            'bootstrap', 'var'};
    specs{2} = {'estimator', 'lp',...
            'bias_corr', false,...
            'bootstrap', 'var'};

else

    specs{1} = {'estimator', 'var',...
            'bias_corr', false};
    specs{2} = {'estimator', 'lp',...
            'bias_corr', false};

end

%% PRELIMINARIES

rng(settings.sim.rng_seed, 'twister'); % set RNG seed

dgps = dgp.ps;

numdgp  = size(dgps,2);                 % no. of DGPs
numhorz = length(settings.est.horzs);   % no. of horizons
numspec = length(specs);                % no. of regression specifications
numrep  = settings.sim.numrep;          % no. of repetitions

spec_shared = {'resp_ind', settings.est.resp_ind, ...
               'innov_ind', settings.est.innov_ind, ...
               'alpha', settings.est.alpha,...
               'no_const', settings.est.no_const,...
               'se_homosk', settings.est.se_homosk,...
               'boot_num', settings.est.boot_num};

dgp.irs_true = nan(numdgp, numhorz);

%% RUN SIMULATIONS

%----------------------------------------------------------------
% Placeholder Matrices
%----------------------------------------------------------------

estims = zeros(numdgp, numspec, numhorz, numrep);
ses    = estims;

cis_lower = zeros(numdgp, numspec, numhorz, 4, numrep);
cis_upper = cis_lower;

ps = zeros(numdgp, numrep);

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

    indx_p   = find(dgp_settings.ps == dgps(1,i_dgp));

    dgp.p           = dgp.ps(i_dgp);
    dgp.n_yp        = dgp.n_y * dgp.p;

    settings.est.p  = dgp.p;

    dgp.A           = dgp_inputs{indx_p}.A;
    dgp.A_c         = dgp_inputs{indx_p}.A_c;
    dgp.H           = dgp_inputs{indx_p}.H;
    dgp.H_c         = dgp_inputs{indx_p}.H_c;
    dgp.D           = dgp_inputs{indx_p}.D;
    dgp.Sigma       = dgp_inputs{indx_p}.Sigma;
    dgp.alpha_tilde = dgp_inputs{indx_p}.alpha_worst(:,:,:,dgp.worst_hor);

    dgp.alpha_tilde(2:end,:,:) = dgp.T_scale^(dgp.zeta) * dgp.alpha_tilde(2:end,:,:);

    dgp.n_s         = dgp.n_y * dgp.p + (dgp.n_y * size(dgp.alpha_tilde,1)-dgp.n_eps);
    
    set_up_varma
    
    % seed
    
    i_rand_seeds = randi(2^32-1,1,numrep);
    
    % Monte Carlo iterations
    
    parfor (i_rep = 1:numrep, settings.sim.num_workers)
    
        % seed
        
        rng(i_rand_seeds(i_rep), 'twister');
        
        % simulate data
        
        data_y = generate_data(dgp);
        if settings.est.p_select == 0
            p = dgp.p;
        else
            p = ic_var(data_y,settings.p_max,settings.est.p_select);
        end
        
        % placeholders
        
        i_rep_estims    = zeros(numspec, numhorz);
        i_rep_ses       = i_rep_estims;
        i_rep_cis_lower = nan(numspec, numhorz, 4);
        i_rep_cis_upper = i_rep_cis_lower;
        
        % run specifications
        
        for i_spec = 1:numspec
        
            [i_rep_estims(i_spec,:),i_rep_ses(i_spec,:),i_cis_dm,i_cis_boot] ...
                = ir_estim(data_y, p, settings.est.horzs, spec_shared{:}, specs{i_spec}{:});
            
            i_rep_cis_lower(i_spec,:,1) = i_cis_dm(1,:);
            i_rep_cis_upper(i_spec,:,1) = i_cis_dm(2,:);
            
            i_rep_cis_lower(i_spec,:,2:end) = i_cis_boot(1,:,:);
            i_rep_cis_upper(i_spec,:,2:end) = i_cis_boot(2,:,:);
        
        end
        
        % store results
        
        estims(i_dgp,:,:,i_rep) = i_rep_estims;
        ses(i_dgp,:,:,i_rep) = i_rep_ses;
        cis_lower(i_dgp,:,:,:,i_rep) = i_rep_cis_lower;
        cis_upper(i_dgp,:,:,:,i_rep) = i_rep_cis_upper;
        ps(i_dgp,i_rep) = p;
        
        if mod(i_rep, ceil(numrep/10)) == 0
            fprintf('%6d%s\n', round(i_rep/numrep*100), '%');
        end
    
    end

end

clear data_y i_cis_dm i_dgp i_rand_seeds i_rep i_rep_cis_lower i_rep_cis_upper i_rep_estims i_rep_ses ...
    i_spec num2dgp numhorz numrep numspec numdgp spec_shared

if settings.sim.num_workers > 0
    delete(poolobj); % Stop parallel workers
end

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
results.ps        = ps;

clear estims ses cis_lower cis_upper

%----------------------------------------------------------------
% Post-Computations
%----------------------------------------------------------------

% coverage

irs_true_reshape = permute(dgp.irs_true, [1 3 2]);

results.cover_inds = (irs_true_reshape >= results.cis_lower ...
                        & irs_true_reshape <= results.cis_upper) + 0; % coverage indicator

results.coverage_prob = mean(results.cover_inds,5); % coverage probability

clear irs_true_reshape

% median length

results.lengths       = results.cis_upper-results.cis_lower;
results.median_length = median(results.lengths,5);

%% SAVE RESULTS

results_filename = sprintf('%s%s', 'sim_varma_', settings.est.type);
save(strcat('results/',results_filename, '.mat'), 'dgp', 'specs', 'settings', 'results');