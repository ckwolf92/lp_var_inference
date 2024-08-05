%% GET POPULATION VARMA(p,infty) AS DGP
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/21/2024

%% HOUSEKEEPING

clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath(fullfile('..','..','auxiliary_functions')))

%% SETTINGS

%----------------------------------------------------------------
% ARMA(1,1) Models
%----------------------------------------------------------------

settings.rhos   = 0.9;
settings.thetas = [0 1/4 1/2];

settings.n_rho   = length(settings.rhos);
settings.n_theta = length(settings.thetas);

dgps = cell(settings.n_rho,settings.n_theta);

%----------------------------------------------------------------
% Computational Settings
%----------------------------------------------------------------

settings.VMA_hor            = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength   = 250; % lag length in population VAR
settings.alpha_lags         = 40;

%----------------------------------------------------------------
% DGP Settings
%----------------------------------------------------------------

settings.VAR_estimlaglength = 1; % number of lags in mis-specified VAR
settings.max_hor_h          = 10; % maximal IRF horizon of interest
settings.max_hor_alpha_l    = settings.VMA_hor; % maximal lag length for worst-case \alpha(L)
settings.resp_ind           = 1; % response variable of interest
settings.innov_ind          = 1; % innovation variable of interest
settings.T                  = 240; % sample size for DGP
settings.zeta               = 1/2; % mis-specification scaling

%% GET DGPs

for i_rho = 1:settings.n_rho

    for i_theta = 1:settings.n_theta

        %----------------------------------------------------------------
        % ABCD Representation
        %----------------------------------------------------------------
        
        % parameters
        
        model.param.rho   = settings.rhos(i_rho);
        model.param.theta = settings.thetas(i_theta);
        model.param.sigma = 1;
        
        % ABCD
        
        model.ABCD.A = [model.param.rho, model.param.sigma * model.param.theta; 0, 0];
        model.ABCD.B = [model.param.sigma; 1];
        model.ABCD.C = [model.param.rho, model.param.sigma * model.param.theta];
        model.ABCD.D = model.param.sigma;
        
        % system size
        
        model.n_s   = size(model.ABCD.A,1);
        model.n_eps = size(model.ABCD.B,2);
        model.n_y   = size(model.ABCD.C,1);
        
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
        
        dgps{i_rho,i_theta} = dgp_fn(VAR_p,VMA,model,settings);

        clear model VAR_infty VAR_p VARMA VMA y_aux

    end

end

clear i_rho i_theta

dgp_inputs   = dgps;
dgp_settings = settings;
clear settings dgps

results_filename = 'arma_dgps';
save(strcat(results_filename, '.mat'), 'dgp_inputs', 'dgp_settings');