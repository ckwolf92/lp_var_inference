% Generate figures from simulation results

% file names
save_suffix   = '.eps'; % suffix for saved figures

if contains(dgp_type, 'varma')
    load_filename = fullfile('results', strcat('sim_', dgp_type, '_', scheme, '.mat'));  % load results from this file
    sim_genfigures_varmaplot;

elseif contains(dgp_type, 'arma') 
    sim_genfigures_armaplot    
    load_filename = fullfile('results', strcat('sim_', dgp_type, '.mat'));  % load results from this file
end

if contains(dgp_type, 'longT')
    p_sel = [2,4];
end


% axis limits
ylim_cover         = [0 1]; % y-limits for coverage prob plot
yticks_length      = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

%% LOAD RESULTS

sim_genfigures_loadres
numproc = length(proc_names); % number of inference procedures

%% GENERATE FIGURES

% set-up
if contains(dgp_type, 'varma')
    status        = mkdir(['figures/', scheme]);
    save_filename = fullfile('figures/',scheme, '/', dgp_type); % first part of file name for saved figures
elseif contains(dgp_type, 'arma') 
    status        = mkdir(['figures/']);
    save_filename = fullfile('figures/', dgp_type); % first part of file name for saved figures
end


horzs   = settings.est.horzs;
numhorz = length(horzs); % no. of estimated impulse response horizons

% size

plotwidth     = 0.39;
gapsize       = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos      = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% plot figure
for d=dgp_sel

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
    if dgp.T == 2000
        plot(horzs, var_asymp_covg, 'Color', [255/255 165/255 0/255], ...
            'LineStyle', ':', 'LineWidth', 3.5)  % Asymptotic VAR coverage
    end

    
    hold off;
    xlim([min(horzs) max(horzs)])
    xlabel('horizon','interpreter','latex');
    ylim(ylim_cover);
    title('coverage probability','interpreter','latex');
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
    legend(proc_names, 'Location', 'SouthEast','interpreter','latex');
    grid on;

    % size
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
 
    graphic_path = [save_filename, '_dgp', num2str(d), save_suffix];
    
    % save
    exportgraphics(the_f,  graphic_path)
    close(the_f);

end




