%% WORST-CASE ALPHA(L), AR(1)
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 07/18/2024

%% HOUSEKEEPING

clear
clc
close all

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

% linestyle

settings.linestyle = {'-','--','-.',':','.'};

% figure spacing

plotwidth = 0.265;
gapsize = 0.055;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% figure selection

h_plot_indx = [1 5 10];
n_h         = length(h_plot_indx);

%----------------------------------------------------------------
% Figures
%----------------------------------------------------------------

figure(1)

for i_rho = 1:n_rho

    subplot(1,3,i_rho)
    pos = get(gca, 'Position');
    set(gca,'FontSize',18)
    set(gca,'TickLabelInterpreter','latex')
    pos(1) = left_pos(i_rho);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    hold on
    for i_h = 1:n_h
        plot(1:1:max_hor_bias,bias_all(:,h_plot_indx(i_h),i_rho),settings.linestyle{i_h},'LineWidth',4)
    end
    hold on
    xlim([1,20])
    ylim([-1 0.5])
    set(gcf,'color','w')
    title([' $\rho$ = ' num2str(rho_list(i_rho))],'interpreter','latex','fontsize',25)
    xlabel('$\ell$','interpreter','latex','FontSize',20)
    if i_rho == 3
        legend({'$h = 1$','$h = 5$', '$h = 10$'},'Location','Southeast','fontsize',18,'interpreter','latex','Orientation','horizontal','NumColumns',2)
    end
    grid on
    hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.4*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
mkdir('results');
print('results/ar_bias','-depsc');