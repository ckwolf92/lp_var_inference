%% LAG LENGTHS IN VAR LITERATURE
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 01/06/2025

%% HOUSEKEEPING

clc
clear all
close all

%% DATA

% import

data_table = readtable('lit_varlags_raw.csv',detectImportOptions('lit_varlags_raw.csv'));
data = table2array(data_table);

categories = data_table.Properties.VariableNames;

% save individual series

lags           = data(:,1); % lag length
maxhor         = data(:,2); % maximal horizon
freq           = data(:,3); % data frequency (1,4,12,365)
lag_crit_indic = data(:,4); % use of any lag length criterion?
bayes_indic    = data(:,5); % bayesian?

% compute some statistics

lags_freq_rel = lags ./ freq;
lags_maxhor_rel = lags ./ maxhor;

%% REPORT RESULTS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = 0.6 * [1 1 1];
settings.colors.green  = [37/255 152/255 14/255];

plotwidth = 0.41;
gapsize = 0.05;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

%----------------------------------------------------------------
% Lag Length
%----------------------------------------------------------------

figure(1)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
histogram(lags_freq_rel,10,'normalization','probability',...
    'FaceColor',settings.colors.grey,'EdgeColor',settings.colors.grey)
hold on
plot(mean(lags_freq_rel) * [1 1], [0 0.6],'linewidth',3,'linestyle','-','color',settings.colors.black)
hold on
plot(sum(lags_freq_rel .* lag_crit_indic) / sum(lag_crit_indic) * [1 1], ...
    [0 0.6],'linewidth',3,'linestyle','--','color',settings.colors.green)
hold on
set(gcf,'color','w')
xlabel('Lag Length/Frequency','interpreter','latex','fontsize',22)
xticks([0 1 2 3 4])
yticks([0:0.1:0.6])
ylim([0 0.6])
text(1.05,0.5,'Avg.: $0.96$','interpreter','latex','FontSize',20,'color',settings.colors.black)
text(0.1,0.5,'IC: $0.83$','interpreter','latex','FontSize',20,'color',settings.colors.green)
grid on
hold off

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
histogram(lags_maxhor_rel,10,'normalization','probability',...
    'FaceColor',settings.colors.grey,'EdgeColor',settings.colors.grey)
hold on
plot(mean(lags_maxhor_rel) * [1 1], [0 0.6],'linewidth',3,'linestyle','-','color',settings.colors.black)
hold on
plot(sum(lags_maxhor_rel .* lag_crit_indic) / sum(lag_crit_indic) * [1 1], ...
    [0 0.6],'linewidth',3,'linestyle','--','color',settings.colors.green)
hold on
set(gcf,'color','w')
xlabel('Lag Length/Max.\ Horizon','interpreter','latex','fontsize',22)
xticks([0:0.25:1.5])
yticks([0:0.1:0.6])
ylim([0 0.6])
text(0.3,0.5,'Avg.: $0.28$','interpreter','latex','FontSize',20,'color',settings.colors.black)
text(-0.02,0.5,'IC: $0.21$','interpreter','latex','FontSize',20,'color',settings.colors.green)
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.2*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('_figures/figure_readme_1','-depsc');

%----------------------------------------------------------------
% Shrinkage
%----------------------------------------------------------------

disp(['A fraction ' num2str(mean(lag_crit_indic)) ...
    ' of specifications select the lag length using some kind of information criterion.'])

disp(['A fraction ' num2str(mean(bayes_indic)) ...
    ' emply additional Bayesian shrinkage.'])