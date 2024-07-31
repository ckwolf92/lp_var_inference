%% LP vs AR INFERENCE: SIMULATIONS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/21/2024

%% HOUSEKEEPING

clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))
addpath(genpath('../auxiliary_functions'))

% Get DGP
cd inputs/; run('get_arma.m'); cd ..; clear;

%% SETUP

%--------------------------------------------------------------------------
% Settings
%--------------------------------------------------------------------------

resp_ind  = 1;      % index of response variable of interest
innov_ind = 1;      % index of innovation variable
horzs     = 1:6;    % horizons of interest
exercise  = 'arma';
boot      = false;  % boot: Run bootstrap. 
longT     = false;  % longT: true: Run T=2000. false: Run T=240.
spec      = '';     % empty for arma: Argument used for VARMA exercises.
scheme    = '';     % empty for arma: Argument used for VARMA exercises. 

%--------------------------------------------------------------------------
% Inputs
%--------------------------------------------------------------------------

load('inputs/arma_dgps')

%--------------------------------------------------------------------------
% Complete Setup
%--------------------------------------------------------------------------

sim_setup_general
sim_setup_arma

%% RUN SIMULATIONS

sim_run