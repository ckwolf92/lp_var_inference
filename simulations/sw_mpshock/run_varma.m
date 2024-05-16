%% LP vs VAR INFERENCE: SIMULATIONS
% this version: 05/16/2024

% HOUSEKEEPING
clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))


% sw_mpshock-specific settings
resp_ind  = 2;    % index of response variable of interest
innov_ind = 1;    % index of innovation variable
horzs     = 0:40; % horizons of interest

sim_main_sw