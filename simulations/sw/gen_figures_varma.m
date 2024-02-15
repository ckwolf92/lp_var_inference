%% LP vs VAR INFERENCE: GENERATE FIGURES
% this version: 01/24/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

%% SETTINGS

% DGP type
dgp_type = 'varma_estp'; % either 'varma_fixp', 'varma_estp', or 'varma_worst'

% file names
load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat'));  % load results from this file
save_suffix   = '.png'; % suffix for saved figures

% select DGPs, CI procedures, and line specs
switch dgp_type(7)
    case 'f'
        p_sel = [1,2,4];
        rhos_sel   = [0.9];
        proc_names = {'AR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [0 0 0; lines(7); 0.5 0.5 0.5];
        line_specs = {'-o', '--o', ':o', '-.o', '-x', '--x', ':x', '-.x', '-s'};
    case 'e'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'AR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [0 0 0; lines(7); 0.5 0.5 0.5];
        line_specs = {'-o', '--o', ':o', '-.o', '-x', '--x', ':x', '-.x', '-s'};
    case 'w'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'AR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 2 1];
        line_colors = [0 0 0; lines(7); 0.5 0.5 0.5];
        line_specs = {'-o', '--o', ':o', '-.o', '-x', '--x', ':x', '-.x', '-s'};
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

status = mkdir('figures');
save_filename = fullfile('figures', dgp_type); % first part of file name for saved figures

horzs   = settings.est.horzs;
numhorz = length(horzs); % no. of estimated impulse response horizons

for d=dgp_sel

    the_f = figure;

    % coverage probability
    subplot(1,2,1);
    set(gca,'TickLabelInterpreter','latex')
    hold on;
    for j=1:numproc
        plot(horzs, squeeze(results.coverage_prob(d,procs(j,1),:,procs(j,2))), line_specs{j}, 'Color', line_colors(j,:));
    end
    plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), 'Color', 'k', 'LineStyle', ':'); % Nominal confidence level
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim(ylim_cover);
    title('coverage probability','interpreter','latex');

    % median length
    subplot(1,2,2);
    set(gca,'TickLabelInterpreter','latex')
    hold on;
    for j=1:numproc
        plot(horzs, squeeze(log10(results.median_length(d,procs(j,1),:,procs(j,2)))), line_specs{j}, 'Color', line_colors(j,:));
    end
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim([min(yticks_length) max(yticks_length)]);
    yticks(yticks_length);
    yticklabels(yticklabels_length);
    title('median length, log scale','interpreter','latex');
    legend(proc_names, 'Location', 'SouthEast','interpreter','latex');

    % title
    if dgp_type(7) == 'f'
%         sgtitle(sprintf('%s %4.2f%s %4.2f%s', '$p$ = ', dgp.dgps(1,d)),'interpreter','latex');
        sgtitle(sprintf('%s %4.0f%s %4.0f%s', '$p$ = ', dgp.dgps(1,d)),'interpreter','latex');
    else
        sgtitle(sprintf('%s %4.2f%s %4.2f%s', '$p$ = AIC'),'interpreter','latex');
    end

    % size
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1*pos(3) 0.56*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    % save
    saveas(the_f,sprintf('%s%s%d%s%d%s', save_filename, '_dgp', d, save_suffix));
    close(the_f);

end