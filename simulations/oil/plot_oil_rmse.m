%% FIGURES FOR OIL SHOCK
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
addpath('results/')


%% LOAD RESULTS

exp_id = 2; % 1 is baseline p, 2 is AIC p

load(['sim_', num2str(exp_id)])

%% CREATE BIAS, VARIANCE, AND MSE

% bias squared

irs_true_reshape = permute(dgp.irs_true, [1 3 2]);

results.bias2 = (mean(results.estims,4) - irs_true_reshape).^2;

clear irs_true_reshape

% variance

results.vce = var(results.estims,1,4);

% MSE

results.mse   = 0.5 * results.bias2 + 0.5 * results.vce;

%% GENERATE FIGURES

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

% plot objects

horzs = settings.est.horzs;

% lines

colors.blue    = [102/255 178/255 255/255];
colors.dblue   = 0.7 * colors.blue;
colors.lblue   = 0.4 * colors.blue + 0.6 * [1 1 1];
colors.red     = [204/255 0/255 0/255];
colors.dred    = 0.7 * colors.red;
colors.lred    = 0.4 * colors.red + 0.6 * [1 1 1];
colors.brown   = [102/255 0/255 0/255];
colors.lbrown  = 0.5 * colors.brown + 0.5 * [1 1 1];
colors.purple  = [153/255 0/255 153/255];
colors.lpurple = 0.5 * colors.purple + 0.5 * [1 1 1];

line_colors = [colors.red; ...
                colors.blue];

line_specs = {'-','-'}; 

line_width = 5 * ones(length(line_specs),1);

%----------------------------------------------------------------
% Estimation: Bias-Variance-MSE Plots
%----------------------------------------------------------------

% preparations

results.target_irf = dgp.irs_true';

dgp_sel = 1:1:size(results.coverage_prob,1);
% dgp_sel = 1;

the_rms_irf = sqrt(mean(results.target_irf.^2)); % root average squared true IRF across horizons
the_objects = {'bias2','vce','mse'}; % objects to plot
the_titles = {'Bias','Standard Deviation','RMSE'};

proc_names = {'VAR', 'LP'};

% figures
cd([path, '/figures'])
fig_count = 1;

for d = dgp_sel

for j = 1:length(the_objects)

    the_result = squeeze(sqrt(results.(the_objects{j})(d,:,:)))';
    the_result = the_result./the_rms_irf(1,d);

    figure(fig_count)

    pos = get(gca, 'Position');
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Position', pos)
    hold on
    for i_proc = 1:2
        plot(horzs, the_result(:,i_proc), ...
            line_specs{i_proc}, ...
            'Color', line_colors(i_proc,:), ...
            'LineWidth', line_width(i_proc));
        hold on
    end

    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    if j ~= 1
        legend(proc_names, 'Location', 'southeast', ...
            'NumColumns', 1, 'interpreter', 'latex', 'FontSize', 18);
    end
    grid on
    title(the_titles(j), 'interpreter', 'latex', 'FontSize', 22);

    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1.2*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    exportgraphics(gcf, fullfile([the_objects{j} '.eps']))

    fig_count = fig_count + 1;

end

end
cd(path)