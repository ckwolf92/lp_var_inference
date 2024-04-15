%% GET POPULATION VARMA(p,infty) AS DGP
% this version: 02/07/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath(genpath('../../../functions'))
addpath('_subroutines')

%% SETTINGS

%----------------------------------------------------------------
% Lag Lengths
%----------------------------------------------------------------

settings.ps  = [1 2 4 8];
settings.n_p = length(settings.ps);

dgps = cell(settings.n_p,1);

%----------------------------------------------------------------
% Computational Settings
%----------------------------------------------------------------

settings.VMA_hor            = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength   = 250; % lag length in population VAR
settings.alpha_lags         = 100;

%----------------------------------------------------------------
% DGP Settings
%----------------------------------------------------------------

settings.max_hor_h          = 10; % maximal IRF horizon of interest
settings.max_hor_alpha_l    = settings.VMA_hor; % maximal lag length for worst-case \alpha(L)
settings.resp_ind           = 2; % response variable of interest
settings.innov_ind          = 1; % innovation variable of interest
settings.T                  = 240; % sample size for DGP
settings.zeta               = 1/2; % mis-specification scaling

%% GET DGPs

for i_p = 1:settings.n_p

    settings.VAR_estimlaglength = settings.ps(i_p); % number of lags in mis-specified VAR

    %----------------------------------------------------------------
    % SW Model
    %----------------------------------------------------------------
    
    % run model
    
    dynare SW_Model noclearall
    clean_folder_SW
    
    SW_model.decision = decision(2:end,:);
    
    % ABCD representation
    
    SW_model.obs = [1 19 20 21]; % (shock,pi,w,l)
    
    SW_model.n_y   = size(SW_model.obs,2);
    SW_model.n_eps = M_.exo_nbr;
    SW_model.n_s   = M_.nspred;
    
    SW_model.ABCD = ABCD_fun_SW(SW_model);
    
    % clean-up
    
    clean_workspace_SW
    
    model = SW_model;
    clear SW_model
    
    %----------------------------------------------------------------
    % VARMA(p,\infty)
    %----------------------------------------------------------------
    
    % VAR(infty)
    
    VAR_infty = popVAR(model,settings);
    y_aux     = get2ndmoments_VAR(VAR_infty,model,settings);
    
    % VAR(p)
    
    VAR_p  = popVARp(model,settings,y_aux);
    
    % Residual VMA(infty)
    
    [VMA,VARMA] = getresidVMA(VAR_infty,VAR_p,model,settings);
    
    %----------------------------------------------------------------
    % Simulation DGP
    %----------------------------------------------------------------
    
    dgps{i_p} = dgp_fn(VAR_p,VMA,model,settings);

    clear model VAR_infty VAR_p VARMA VMA y_aux

end

clear i_p

dgp_inputs   = dgps;
dgp_settings = settings;
clear settings dgps

results_filename = 'varma_sw_dgps';
save(strcat(results_filename, '.mat'), 'dgp_inputs', 'dgp_settings');