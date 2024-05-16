%% LP vs VAR INFERENCE: GENERATE FIGURES
% this version: 05/16/2024

%% HOUSEKEEPING

clc
clear all
close all
addpath(genpath('../../functions'))
warning('off','MATLAB:dispatcher:nameConflict')

%% SETTINGS

% DGP type
% either arma_ or arma_boot
dgp_type = 'arma_';  
sim_genfigures

