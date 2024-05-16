%% LP vs VAR INFERENCE: SIMULATIONS
% this version: 05/16/2024

% HOUSEKEEPING
clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))


% sw_mprecursive-specific settings
resp_ind  = 1;    % index of response variable of interest
innov_ind = 3;    % index of innovation variable
horzs     = 1:40; % horizons of interest

sim_main_sw