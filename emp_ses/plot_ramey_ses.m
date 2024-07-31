%% PLOT EMPIRICAL STANDARD ERRORS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 07/18/2024

%% HOUSEKEEPING

clear
clc
close all

addpath('auxiliary_functions/')
addpath(genpath('../functions'))
addpath('results')

%% PREPARATIONS

%----------------------------------------------------------------
% Import Results
%----------------------------------------------------------------

load res_application

%----------------------------------------------------------------
% Organize
%----------------------------------------------------------------

QQ      = [.1, .5, .9];  % Quantiles to be plotted
numappl = length(appl);

% organize quarterly applications

se_ratio_q = [];
for i_appl = 2:numappl
    se_ratio_q = [se_ratio_q;
        squeeze(appl(i_appl).results.ses_boot(1,:,:) ./ appl(i_appl).results.ses_boot(2,:,:))];
end
se_ratio_q = se_ratio_q';

% organize monetary application (monthly)

se_ratio_m = squeeze(appl(1).results.ses_boot(1,:,:) ./appl(1).results.ses_boot(2,:,:))';

% define groups (by half-year)

groups_q = {1:2,3:4,5:6,7:8,9:10,11:12,13:14,15:21};  
groups_m = {1:6,7:12,13:18,19:24,25:30,31:36,37:42,43:49};
numgroups = length(groups_q);

se_ratio_grouped = nan(numgroups, length(QQ));

% group IRFs

for igroup = 1:numgroups
    se_ratio_grouped_i = [reshape(se_ratio_q(groups_q{igroup}, :),[], 1);
                          reshape(se_ratio_m(groups_m{igroup}, :), [], 1)];
    se_ratio_grouped(igroup, :) = quantile(se_ratio_grouped_i, QQ);
end

%% PLOT

figure('Units','inches','Position',[2, 2, 8, 4]);
set(gca,'TickLabelInterpreter','latex', 'FontSize', 12);
grid on
hold on
scatter(4*[0.5:0.5:4], se_ratio_grouped(:,2), 'k', 'filled', 'Marker','square')  % Plot center
hold on
scatter(4*[0.5:0.5:4], se_ratio_grouped(:,[1,3]), 'k', 'Marker', '_')  % Plot upper/lower
for j = 1:numgroups
    plot([j*2,j*2], se_ratio_grouped(j, [1,3]), 'Color', 'k', 'LineWidth', 1)  % Vertical line
end
xlim([1 17])
xticks([2:2:16])
xticklabels({'2', '4', '6', '8', '10', '12', '14', '16+'})
xlabel('Horizon (quarters)', 'Interpreter', 'latex')
ylim([0 1.1])
print('results/ses_emp','-depsc');