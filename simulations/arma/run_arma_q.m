%% LP vs VAR INFERENCE: SIMULATIONS
% this version: 02/14/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath(genpath('../../functions'))

%% DATA-GENERATING PROCESSES (DGPs)

load('inputs/arma_dgps')

dgp = struct;

% parameters

dgp.rhos       = dgp_settings.rhos; % AR parameters
dgp.p          = dgp_settings.VAR_estimlaglength; % lags
dgp.T          = 240;  % sample size
dgp.T_scale    = dgp_settings.T; % artificial sample size to scale mis-specification
dgp.zeta       = dgp_settings.zeta; % mis-specification dampening
dgp.thetas     = dgp_settings.thetas(1:2); % MA term

% system size

dgp.n_y   = 1;
dgp.n_yp  = dgp.n_y * dgp.p;
dgp.n_eps = dgp.n_y;
dgp.n_s   = dgp.n_yp + dgp.n_y * dgp_settings.alpha_lags;

%% SETTINGS

settings = struct;

%----------------------------------------------------------------
% Monte Carlo Simulation
%----------------------------------------------------------------

settings.sim.numrep      = 2e4; % no. of repetitions
settings.sim.rng_seed    = 202007252; % random number seed
settings.sim.num_workers = 8; % no. of parallel workers (=0: run serial)

%----------------------------------------------------------------
% Estimation
%----------------------------------------------------------------

settings.est.type = 'quarterly'; % estimation methods type

settings.est.p         = 1; % lag length used for estimation
settings.est.horzs     = [1:1:4]; % horizons of interest
settings.est.no_const  = true; % true: omit intercept
settings.est.se_homosk = true; % true: homoskedastic ses
settings.est.alpha     = 0.1; % significance level

settings.est.resp_ind  = 1; % index of response variable of interest
settings.est.innov_ind = 1; % index of innovation of interest

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

aux1  = repmat(dgp.rhos,size(dgp.thetas));
aux2  = repmat(dgp.thetas',size(dgp.rhos))';
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

cis_lower = zeros(numdgp, numspec, numhorz, 4, numrep);
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

    indx_rho   = find(dgp_settings.rhos == dgps(1,i_dgp));
    indx_theta = find(dgp_settings.thetas == dgps(2,i_dgp));

    dgp.A           = dgp_inputs{indx_rho,indx_theta}.A;
    dgp.A_c         = dgp_inputs{indx_rho,indx_theta}.A_c;
    dgp.H           = dgp_inputs{indx_rho,indx_theta}.H;
    dgp.H_c         = dgp_inputs{indx_rho,indx_theta}.H_c;
    dgp.D           = dgp_inputs{indx_rho,indx_theta}.D;
    dgp.Sigma       = dgp_inputs{indx_rho,indx_theta}.Sigma;
    dgp.alpha_tilde = dgp_inputs{indx_rho,indx_theta}.alpha_tilde;

    dgp.alpha_tilde(2:end) = dgp.T_scale^(dgp.zeta) * dgp.alpha_tilde(2:end);
    
    set_up_varma
    
    % seed
    
    i_rand_seeds = randi(2^32-1,1,numrep);
    
    % Monte Carlo iterations
    
    parfor (i_rep = 1:numrep, settings.sim.num_workers)
    
        % seed
        
        rng(i_rand_seeds(i_rep), 'twister');
        
        % simulate data
        
        data_y = generate_data(dgp);
        
        % placeholders
        
        i_rep_estims    = zeros(numspec, numhorz);
        i_rep_ses       = i_rep_estims;
        i_rep_cis_lower = nan(numspec, numhorz, 4);
        i_rep_cis_upper = i_rep_cis_lower;
        
        % run specifications
        
        for i_spec = 1:numspec
        
            [i_rep_estims(i_spec,:),i_rep_ses(i_spec,:),i_cis_dm,i_cis_boot] ...
                = ir_estim(data_y, settings.est.p, settings.est.horzs, spec_shared{:}, specs{i_spec}{:});
            
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

results_filename = sprintf('%s%s', 'sim_arma_', settings.est.type);
save(strcat('results/',results_filename, '.mat'), 'dgp', 'specs', 'settings', 'results');