%% LP vs VAR INFERENCE: GENERATE FIGURES
% this version: 01/24/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

%% SETTINGS

% DGP type
dgp_type = 'varma_bfixp_longT'; % either 'varma_fixp', 'varma_estp', 'varma_worst', 'varma_bfixp', 'varma_bestp', or 'varma_bworst'

% file names
load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat'));  % load results from this file
save_suffix   = '.eps'; % suffix for saved figures

% select DGPs, CI procedures, and line specs
switch dgp_type(7:8)
    case 'fi'
        p_sel = [1,2,4,8];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'es'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'wo'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'bf'
        p_sel = [1,2,4,8];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 1 4;
                 2 1;
                 2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
    case 'be'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 1 4;
                 2 1;
                 2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
    case 'bw'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 1 4;
                 2 1;
                 2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
end


if strcmp(dgp_type, 'varma_bfixp_longT')
    p_sel = [2,4];
end


% axis limits
ylim_cover = [0 1]; % y-limits for coverage prob plot
yticks_length = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

%% LOAD RESULTS

load(load_filename);

% pick out indices of selected DGPs
numdgp_sel = length(p_sel);
dgp_sel = zeros(1,numdgp_sel);
for i=1:length(p_sel)
    dgp_sel(i) = find(dgp.dgps(1,:)==p_sel(i));
end

numproc = length(proc_names); % number of inference procedures

%% GENERATE FIGURES

% set-up

status = mkdir('figures');
save_filename = fullfile('figures', dgp_type); % first part of file name for saved figures

horzs   = settings.est.horzs;
numhorz = length(horzs); % no. of estimated impulse response horizons

% size

plotwidth = 0.39;
gapsize = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% plot figure
for d=dgp_sel

    compute_var_acovg

    the_f = figure;
    
    % coverage probability
    subplot(1,2,1)
    pos = get(gca, 'Position');
    set(gca,'FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(1);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on;
    
    for j=1:numproc
        plot(horzs, squeeze(results.coverage_prob(d,procs(j,1),:,procs(j,2))), ...
            line_specs{j}, 'Color', line_colors(j,:),'LineWidth',3);
    end
    plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), ...
        'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Nominal confidence level
    plot(horzs, var_asymp_covg, 'Color', 'k', ...
        'LineStyle', ':', 'LineWidth', 2)  % Asymptotic VAR coverage

    
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim(ylim_cover);
    title('coverage probability','interpreter','latex');
    grid on

    % median length
    subplot(1,2,2)
    pos = get(gca, 'Position');
    set(gca,'FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(2);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on;
    for j=1:numproc
        plot(horzs, squeeze(log10(results.median_length(d,procs(j,1),:,procs(j,2)))), ...
            line_specs{j}, 'Color', line_colors(j,:),'LineWidth',3);
    end
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim([min(yticks_length) max(yticks_length)]);
    yticks(yticks_length);
    yticklabels(yticklabels_length);
    title('median length, log scale','interpreter','latex');
    legend(proc_names, 'Location', 'SouthEast','interpreter','latex');
    grid on;

    % size
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.3*1*pos(3) 1.3*0.56*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    % save
%    saveas(the_f,sprintf('%s%s%d%s%d%s', save_filename, '_dgp', d, save_suffix));
    exportgraphics(the_f, sprintf('%s%s%d%s%d%s', save_filename, '_dgp', d, save_suffix))
    close(the_f);

end