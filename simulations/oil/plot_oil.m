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

%% GET RESULTS

%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

exp_id = 1; % 1 is fixed p, 2 is AIC p

load(['sim_', num2str(exp_id)])

%----------------------------------------------------------------
% Bias, Variance, MSE
%----------------------------------------------------------------

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
% Inference
%----------------------------------------------------------------

% some auxiliary variables

horzs   = settings.est.horzs;
numhorz = length(horzs); % no. of estimated impulse response horizons

proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
procs = [1 1; % first index: inference procedure; second index: type of confidence interval
    1 2;
    2 1;
    2 4];
line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
    102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
line_specs = {'-', '-.', ':', '--'};
numproc = length(proc_names); % number of inference procedures

% axis limits

ylim_cover         = [0 1]; % y-limits for coverage prob plot
yticks_length      = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

% size

plotwidth     = 0.39;
gapsize       = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos      = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% figure name

save_suffix   = '.eps'; % suffix for saved figures
save_filename = ['oil_ci_', num2str(exp_id)];

% plot figures

cd([path, '/figures'])

dgp_sel = 1:1:size(results.coverage_prob,1);

for d = dgp_sel

    var_asymp_covg = dgp.var_asymp_covg(d, :);

    the_f = figure;

    % coverage probability
    subplot(1,2,1)
    pos = get(gca, 'Position');
    set(gca,'FontSize',18)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(1);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on;

    for j=1:numproc
        plot(horzs, squeeze(results.coverage_prob(d,procs(j,1),:,procs(j,2))), ...
            line_specs{j}, 'Color', line_colors(j,:),'LineWidth',5);
    end
    plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), ...
        'Color', 'k', 'LineStyle', ':', 'LineWidth', 3.5); % Nominal confidence level


    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim(ylim_cover);
    title('coverage probability','interpreter','latex');

    legend(proc_names, 'Location', 'SouthEast','interpreter','latex');
    grid on

    % median length
    subplot(1,2,2)
    pos = get(gca, 'Position');
    set(gca,'FontSize',18)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(2);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on;
    for j=1:numproc
        plot(horzs, squeeze(log10(results.median_length(d,procs(j,1),:,procs(j,2)))), ...
            line_specs{j}, 'Color', line_colors(j,:),'LineWidth',5);
    end
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim([min(yticks_length) max(yticks_length)]);
    yticks(yticks_length);
    yticklabels(yticklabels_length);
    title('median length, log scale','interpreter','latex');
    grid on;

    % size
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');

    if length(dgp_sel) == 1
        graphic_path = [save_filename, save_suffix];
    else
        graphic_path = [[save_filename '_' num2str(d)], save_suffix];
    end

    % save
    exportgraphics(the_f,  graphic_path)

end

cd([path])

%----------------------------------------------------------------
% Point Estimation
%----------------------------------------------------------------

% preparations

results.target_irf = dgp.irs_true';

dgp_sel = 1:1:size(results.coverage_prob,1);

the_rms_irf = sqrt(mean(results.target_irf.^2)); % root average squared true IRF across horizons
the_objects = {'bias2','vce','mse'}; % objects to plot
the_titles = {'Bias','Standard Deviation','RMSE'};

proc_names = {'VAR', 'LP'};

line_colors = line_colors([1 3],:);
line_specs  = line_specs([1,3]);

% plot figures

cd([path, '/figures'])

for d = dgp_sel

for j = 1:length(the_objects)

    the_result = squeeze(sqrt(results.(the_objects{j})(d,:,:)))';
    the_result = the_result./the_rms_irf(1,d);

    the_f = figure;

    pos = get(gca, 'Position');
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Position', pos)
    hold on
    for i_proc = 1:2
        plot(horzs, the_result(:,i_proc), ...
            line_specs{i_proc}, ...
            'Color', line_colors(i_proc,:), ...
            'LineWidth', 5);
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

    % figure name

    if length(dgp_sel) == 1
        graphic_path = [['oil_', num2str(the_objects{j}) '_' num2str(exp_id)], save_suffix];
    else
        graphic_path = [['oil_', num2str(the_objects{j}) '_' num2str(exp_id) '_' num2str(d)], save_suffix];
    end

    % save
    exportgraphics(the_f,  graphic_path)

end

end
cd(path)