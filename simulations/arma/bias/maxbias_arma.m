%% WORST-CASE ALPHA(L), AR
% this version: 02/07/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath(genpath('../../../functions'))

%% COMPUTATIONS

% model parameters

rho_list = [0.3, 0.6, 0.95];
n_rho    = length(rho_list);

% settings

max_hor_bias = 30;
max_hor_irf  = 20;

h_list = [1:1:max_hor_irf];

% bias

bias_all = NaN(max_hor_bias,max_hor_irf,n_rho);
for i_rho = 1:n_rho
    rho = rho_list(i_rho);
    for i = 1:max_hor_irf
        for j = 1:max_hor_bias
            h = h_list(i);
            l = j;
            bias_all(j,i,i_rho) = h * rho^(h-1) * (1-rho^2) * rho^(l-1) - (l <= h) * rho^(h-l);
        end
    end
end

bias_direct_all = NaN(max_hor_bias,max_hor_irf,n_rho);
for i_rho = 1:n_rho
    rho = rho_list(i_rho);
    for i = 1:max_hor_irf
        for j = 1:max_hor_bias
            h = h_list(i);
            l = j;
            bias_direct_all(j,i,i_rho) = - (l <= h) * rho^(h-l);
        end
    end
end

bias_indirect_all = NaN(max_hor_bias,max_hor_irf,n_rho);
for i_rho = 1:n_rho
    rho = rho_list(i_rho);
    for i = 1:max_hor_irf
        for j = 1:max_hor_bias
            h = h_list(i);
            l = j;
            bias_indirect_all(j,i,i_rho) = h * rho^(h-1) * (1-rho^2) * rho^(l-1);
        end
    end
end

%% PLOT RESULTS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [150/255 150/255 150/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];
settings.colors.hor = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255; ... % blue
                            37/255 152/255 14/255]; % green

% figure spacing

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% figure selection

h_plot_indx = [1 5 10];
n_h         = length(h_plot_indx);

%----------------------------------------------------------------
% Figures
%----------------------------------------------------------------

status = mkdir('figures');

figure(1)

for i_rho = 1:n_rho

    subplot(1,3,i_rho)
    pos = get(gca, 'Position');
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(i_rho);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on
    for i_h = 1:n_h
        plot(1:1:max_hor_bias,bias_all(:,h_plot_indx(i_h),i_rho),'-','Color',settings.colors.hor(i_h,:),'LineWidth',4)
    end
    hold on
    xlim([1,20])
    ylim([-1 0.5])
    set(gcf,'color','w')
    title([' $\rho$ = ' num2str(rho_list(i_rho))],'interpreter','latex','fontsize',24)
    xlabel('$\ell$','interpreter','latex','FontSize',20)
    % ylabel('\% Deviation','interpreter','latex','FontSize',20)
    if i_rho == 2
        legend({'$h = 1$','$h = 5$', '$h = 10$'},'Location','Northeast','fontsize',18,'interpreter','latex','Orientation','horizontal')
    end
    grid on
    hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 0.95*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figures/ar_bias','-dpng');

figure(2)

for i_rho = 1:n_rho

    subplot(1,3,i_rho)
    pos = get(gca, 'Position');
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(i_rho);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on
    for i_h = 1:n_h
        plot(1:1:max_hor_bias,bias_direct_all(:,h_plot_indx(i_h),i_rho),'-','Color',settings.colors.hor(i_h,:),'LineWidth',4)
    end
    hold on
    xlim([1,20])
    ylim([-1 0.5])
    set(gcf,'color','w')
    title([' $\rho$ = ' num2str(rho_list(i_rho))],'interpreter','latex','fontsize',24)
    xlabel('$\ell$','interpreter','latex','FontSize',20)
    % ylabel('\% Deviation','interpreter','latex','FontSize',20)
    if i_rho == 2
        legend({'$h = 1$','$h = 5$', '$h = 10$'},'Location','Northeast','fontsize',18,'interpreter','latex','Orientation','horizontal')
    end
    grid on
    hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 0.95*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figures/ar_bias_direct','-dpng');

figure(3)

for i_rho = 1:n_rho

    subplot(1,3,i_rho)
    pos = get(gca, 'Position');
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(i_rho);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on
    for i_h = 1:n_h
        plot(1:1:max_hor_bias,bias_indirect_all(:,h_plot_indx(i_h),i_rho),'-','Color',settings.colors.hor(i_h,:),'LineWidth',4)
    end
    hold on
    xlim([1,20])
    ylim([0 1])
    set(gcf,'color','w')
    title([' $\rho$ = ' num2str(rho_list(i_rho))],'interpreter','latex','fontsize',24)
    xlabel('$\ell$','interpreter','latex','FontSize',20)
    % ylabel('\% Deviation','interpreter','latex','FontSize',20)
    if i_rho == 2
        legend({'$h = 1$','$h = 5$', '$h = 10$'},'Location','Northeast','fontsize',18,'interpreter','latex','Orientation','horizontal')
    end
    grid on
    hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 0.95*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figures/ar_bias_indirect','-dpng');