%% GET OIL-SHOCK POPULATION VARMA(p,infty) AS DGP
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 09/08/2025

%% HOUSEKEEPING

clc
clear
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = cd;

addpath(genpath('../auxiliary_functions'))
addpath(genpath('../data'))
addpath('../../estimation')

%% SETTINGS

%----------------------------------------------------------------
% Import Data
%----------------------------------------------------------------

% estimation sample

smplStart = '1974M01'; 
smplEnd   = '2017M12'; 

% identification sample

smplStartProxy = '1974M01';  
smplEndProxy   = '2017M12';

% data construction

loadOilData;

data_level = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100];
data_level = data_level(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)),:);

Y = [proxy, data_level];

data_oil = Y - mean(Y);

clearvars -except data_oil path

%----------------------------------------------------------------
% Lag Lengths
%----------------------------------------------------------------

settings.ps  = [12 15 18];
settings.n_p = length(settings.ps);

dgps     = cell(settings.n_p,1);
VARs_pop = cell(settings.n_p,1);
VARs_p   = cell(settings.n_p,1);

%----------------------------------------------------------------
% Computational Settings
%----------------------------------------------------------------

settings.VMA_hor            = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength   = 18; % lag length in population VAR
settings.alpha_lags         = 100;

%----------------------------------------------------------------
% DGP Settings
%----------------------------------------------------------------

settings.max_hor_h          = 50; % maximal IRF horizon of interest
settings.max_hor_alpha_l    = settings.VMA_hor; % maximal lag length for worst-case \alpha(L)
settings.T                  = 720; % sample size for DGP
settings.zeta               = 1/2; % mis-specification scaling

settings.resp_ind  = 7; % response variable
settings.innov_ind = 1; % impulse variable
settings.horzs     = 0:settings.max_hor_h;

%% GET DGPs

for i_p = 1:settings.n_p

    settings.VAR_estimlaglength = settings.ps(i_p); % number of lags in mis-specified VAR
    
    %----------------------------------------------------------------
    % VARMA(p,\infty)
    %----------------------------------------------------------------
    
    % true VAR(p_0)
    
    VAR_pop = var_fn(data_oil,settings.VAR_poplaglength,settings);
     
    % mis-specified VAR(p)
    
    VAR_p  = var_fn(data_oil,settings.VAR_estimlaglength,settings);
    
    % residual VMA(infty)
    
    [VMA,VARMA] = getresidVMA(VAR_pop,VAR_p,settings);
    
    %----------------------------------------------------------------
    % Simulation DGP
    %----------------------------------------------------------------
    
    dgps{i_p}     = dgp_fn(VAR_p,VMA,settings);
    VARs_pop{i_p} = VAR_pop;
    VARs_p{i_p}   = VAR_p;
    
    clear model VARMA VMA VAR_pop VAR_p

end

clear i_p

dgp_inputs   = dgps;
dgp_settings = settings;
clear settings dgps

save oil_dgps