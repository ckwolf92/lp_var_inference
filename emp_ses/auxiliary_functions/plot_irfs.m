function plot_irfs(Appl_i, scale)

% Unpack
horzs     = Appl_i.est.horzs;
estims    = Appl_i.results.estims;
yname     = Appl_i.data.yname;
name      = Appl_i.name;
cis_lower = Appl_i.results.cis_lower;
cis_upper = Appl_i.results.cis_upper;
numresp   = size(Appl_i.results.estims, 2);


% Get tick size (frequency-dependent)
ticksize = 12/(month(Appl_i.data.date(2)) - month(Appl_i.data.date(1)));


f                       = figure('Visible', 'off', 'Units', 'inches');
f.Position(3:4)         = [6, 3];
ax                      = gca;
ax.FontSize             = 12;
ax.TickLabelInterpreter = 'latex';

for i_resp = 1:numresp  % By response variable
    subplot(2,2, i_resp)

    % Plot estimates
    plot(horzs, squeeze(estims(1, i_resp, :))*scale(i_resp), 'LineWidth', 1.5)
    hold on
    plot(horzs, squeeze(estims(2, i_resp, :))*scale(i_resp), 'LineWidth', 1.5, 'Color', 'r')

    % Plot LP bands
    plot(horzs, squeeze(cis_lower(2, i_resp, :, 1))*scale(i_resp), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', ':')
    plot(horzs, squeeze(cis_upper(2, i_resp, :, 1))*scale(i_resp), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', ':')

    % Visuals
    box on; grid on;
    title(yname{i_resp}, 'Interpreter', 'latex')
    xlim([horzs(1), horzs(end)])
    xticks(horzs(1):ticksize:horzs(end))

end

% Subplot title
sgtitle([name, ' (blue=VAR and red=LP)'], 'Interpreter', 'latex')

end
