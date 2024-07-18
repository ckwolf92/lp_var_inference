function plot_seratio(appl, ses_type)

arguments
    appl
    ses_type = 'delta'
end
    numappl = length(appl);
    f = figure('Visible', 'off', 'Units','inches');
    f.Position(3:4) = [6,4];
    ax              = gca;
    ax.FontSize     = 12;

    Style = {'-', '--', ':', '-.'};

    for i_appl = 1:numappl
        switch ses_type
            case 'delta'
                ses = appl(i_appl).results.ses;
            case 'boot'
                ses = appl(i_appl).results.ses_boot;
        end

        se_ratio = squeeze(ses(1,:,:) ./ ses(2,:,:))';
        horzs = appl(i_appl).est.horzs;
        subplot(2,2,i_appl)
        plot(horzs([1,end]), [1,1], 'Color', 'k', 'LineWidth', 1, ...
            'LineStyle', '-')
        hold on
        
        p = [];

        for i_plot = 1:size(se_ratio, 2)
            p_i = plot(horzs, se_ratio(:, i_plot), 'LineWidth', 1, ...
                'Color', 'b','LineStyle', Style{i_plot});
            p = [p; p_i];
            ylim([0, 1.1*max(se_ratio(:))])
        end
        

        % Visuals
        legend(p, appl(i_appl).data.yname, 'Location', 'southoutside', ...
            'Interpreter','latex', 'Orientation','horizontal', ...
            'FontSize', 9, 'NumColumns',2)
        title(appl(i_appl).name, 'Interpreter','latex')
        box on; grid on;
        xlim([horzs(1) horzs(end)])

        % Set ticks
        ticksize = 12/(month(appl(i_appl).data.date(2)) - ...
                       month(appl(i_appl).data.date(1)));
        xticks(horzs(1):ticksize:horzs(end))
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        sgtitle('SE ratio: VAR to LP', 'Interpreter', 'latex')
    end

end