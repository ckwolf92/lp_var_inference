%% LP vs VAR INFERENCE: SIMULATIONS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/21/2024

%% HOUSEKEEPING

clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(fullfile('..','..','estimation'))
addpath(genpath(fullfile('..','auxiliary_functions')))

%% SETUP

%----------------------------------------------------------------
% Set Identification Scheme
%----------------------------------------------------------------

scheme = 'lshock';  % mpshock, lshock, or mprecursive
sim_setscheme_sw
cd inputs/; run('get_varma_sw.m'); cd ..;  % Get DGP 
clearvars -except scheme resp_ind innov_ind horzs SW_model_obs;

%--------------------------------------------------------------------------
% Specifications
%--------------------------------------------------------------------------

% Note: "Specifications" contains the following fields:
%   - spec: worst (worst-case), estp (estimated p), fixp (fixed p)
%   - boot: Run bootstrap. 
%   - longT: true: Run T=2000. false: Run T=240.

exercise          = 'varma';
fields            = {'spec', 'boot', 'longT'};
Specifications    = cell2struct(cell(length(fields), 1), fields);
Specifications(1) = struct('spec', 'worst', 'boot', false, 'longT', false);
Specifications(2) = struct('spec', 'estp',  'boot', false, 'longT', false);
Specifications(3) = struct('spec', 'fixp',  'boot', false, 'longT', false);
Specifications(4) = struct('spec', 'fixp',  'boot', false, 'longT', true);

%% MAIN

sim_main_sw       % Run simulations