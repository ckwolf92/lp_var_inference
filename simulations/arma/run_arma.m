%% LP vs AR INFERENCE: SIMULATIONS
% this version: 05/16/2024

% HOUSEKEEPING
clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))


% arma-specific settings
resp_ind  = 1;      % index of response variable of interest
innov_ind = 1;      % index of innovation variable
horzs     = 1:6;    % horizons of interest
exercise  = 'arma';
boot      = false;  % boot: Run bootstrap. 
longT     = false;  % longT: true: Run T=2000. false: Run T=240.
spec      = '';     % Empty for arma: Argument used for VARMA exercises.

load('inputs/arma_dgps')


% Setup
sim_setup_general
sim_setup_arma

% Setup arma dgps
aux1  = repmat(dgp.rhos,size(dgp.thetas));
aux2  = repmat(dgp.thetas',size(dgp.rhos))';
dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
clear aux1 aux2

% Run simulations
sim_run