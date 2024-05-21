%% LP vs VAR INFERENCE: SIMULATIONS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/21/2024

%% HOUSEKEEPING

clc
clear
close all
warning('off','MATLAB:dispatcher:nameConflict')
addpath(genpath('../../functions'))
addpath(genpath('../auxiliary_functions'))

% -------------------------------------------------------------------------
% MODIFY TO SET SW SCHEME
% -------------------------------------------------------------------------
scheme = 'mprecursive';  % mpshock, lshock, or mprecursive
sim_setscheme_sw     % Settings for the SW shock ID scheme

%% SETUP

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