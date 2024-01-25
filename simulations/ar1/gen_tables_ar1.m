%% LP vs VAR INFERENCE: GENERATE TABLES
% this version: 01/24/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

%% SETTINGS

% DGP type
dgp_type = 'ar1_boot'; % either 'ar1_delta', 'ar1_limit', or 'ar1_boot'

% file names
load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat'));  % load results from this file
save_filename = dgp_type; % file name for saved table

% select DGPs
shares_sel = [0,0.1,0.25];
rhos_sel   = [0.7];

% select CI procedures
switch dgp_type(5:8)
    case 'delt'
        proc_names = {'$\text{AR}$', '$\text{LP}$'};
        procs = [1 1; % First index: inference procedure; second index: type of confidence interval
                 2 1];
    case 'limi'
        proc_names = {'$\text{AR}$', '$\text{LP}$'};
        procs = [1 1; % First index: inference procedure; second index: type of confidence interval
                 2 1];
    case 'boot'
        proc_names = {'$\text{AR}$', '$\text{AR}_b$', '$\text{LP}$', '$\text{LP}_b$'};
        procs = [1 1;
                 1 2;
                 2 1;
                 2 4];
end

%% LOAD RESULTS

load(load_filename);

% pick out indices of selected DGPs
numdgp_sel = length(shares_sel);
dgp_sel = zeros(1,numdgp_sel);
for j=1:numdgp_sel
    dgp_sel(j) = find(dgp.dgps(1,:)==rhos_sel & dgp.dgps(2,:)==shares_sel(j));
end

numproc = length(proc_names); % number of inference procedures

horzs_sel = [1 4 8 12];
numhorzs_sel = length(horzs_sel);
horzs_sel_indic = zeros(1,numhorzs_sel);
for j = 1:numhorzs_sel
    horzs_sel_indic(j) = find(settings.est.horzs==horzs_sel(j));
end


%% WRITE TABLE

status = mkdir('tables');
f = fopen(fullfile('tables', strcat(save_filename, '.tex')), 'w'); % open file for writing

fprintf(f, '%s%s%s%s%s\n', '\begin{tabular}{r|', repmat('c', 1, numproc), '|', repmat('c', 1, numproc), '}');
fprintf(f, '%s%d%s%d%s\n', '& \multicolumn{', numproc, '}{c|}{Coverage} & \multicolumn{', numproc, '}{c}{Median length} \\');
fprintf(f, '%s', '$h$');
for i=1:2
    for j=1:numproc
        fprintf(f, '%s%s', ' & ', proc_names{j});
    end
end
fprintf(f, '%s\n%s\n', ' \\', '\hline');

for d=dgp_sel
    share = dgp.dgps(2,d);
    fprintf(f, '%s%d%s%4.2f%s\n', '\multicolumn{', 1+2*numproc, '}{c}{\emph{var. share $= ', share, '$}} \\');
    for ih=1:numhorzs_sel
        h = horzs_sel(ih);
        fprintf(f, '%3d', h);
        for j=1:numproc
            cp = results.coverage_prob(d,procs(j,1),horzs_sel_indic(ih),procs(j,2));
            fprintf(f, '%s%5.3f', ' & ', cp);
        end
        for j=1:numproc
            ml = results.median_length(d,procs(j,1),horzs_sel_indic(ih),procs(j,2));
            fprintf(f, '%s%5.3f', ' & ', ml);
        end
        fprintf(f, '%s\n', ' \\');
    end
end

fprintf(f, '%s', '\end{tabular}');

fclose(f);