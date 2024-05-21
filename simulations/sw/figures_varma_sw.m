%% LP vs VAR INFERENCE: GENERATE FIGURES
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/21/2024

%% HOUSEKEEPING

clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))
addpath(genpath('../auxiliary_functions'))

%% SETTINGS

% -------------------------------------------------------------------------
% MODIFY TO SET SW SCHEME
% -------------------------------------------------------------------------
scheme   = 'mprecursive';  % Either mpshock, lshock, or mprecursive

% -------------------------------------------------------------------------
% MODIFY TO SET dgp_type
% -------------------------------------------------------------------------
%    Either 'varma_fixp',  'varma_fixplongT',  'varma_estp',     'varma_worst', 
%           'varma_bfixp', 'varma_bfixplongT', 'varma_bestp', or 'varma_bworst'
dgp_type = 'varma_worst';

%% GENERATE FIGURES

sim_genfigures  