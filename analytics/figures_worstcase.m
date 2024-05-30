clear all;

% Worst-case coverage of conventional VAR CI,
% comparison of bias-aware VAR CI to LP CI,
% and optimal bias-aware CI


%% Settings

signif = 0.1;                           % Significance level
se_ratio = linspace(0.001,0.999,100);   % Values of SE ratio to plot (sqrt{aVar(AR)/aVar(LP)})
M = [0.1 1 1.5 2 3];                    % Values of M to plot
linestyles = {'-','--','-.',':','.'};   % Line style for each value of M
plot_pos = [0 0 8 4];                   % Figure position (in inches)
plot_ylim_length = [0 2];               % Vertical axis limits for relative length plot


%% Compute

n_M = length(M);
n_r = length(se_ratio);

z = norminv(1-signif/2); % Normal critical value
r = @(b,c) normcdf(b-c)+normcdf(-b-c); % Rejection probability
cv = @(b) fzero(@(c) r(b,c)-signif,z+b); % Bias-aware critical value
prob_fct = @(k,b) log(r(b,z))+log(1-r(b/k,z)); % Log joint probability Hausman rejects and VAR CI doesn't cover

% Worst-case coverage of conventional CI as function of relative length
wc_cov = 1-r(M(:).*sqrt(1./(se_ratio(:).^2)'-1),z);

cv_var_aware = nan(n_M,n_r);
weight_opt_aware = nan(n_M,n_r);
length_opt_aware = nan(n_M,n_r);
wc_prob = nan(1,n_r);
wc_bias = nan(1,n_r);
opts = optimoptions('fmincon','Display','notify');

for i_r=1:n_r

    the_aux = sqrt(1./se_ratio(i_r)^2-1);

    for i_M=1:n_M

        % Critical value for bias-aware VAR CI
        cv_var_aware(i_M,i_r) = cv(M(i_M)*the_aux);

        % Optimal bias-aware CI
        the_obj = @(w) cv((1-w)*M(i_M)*the_aux/sqrt(1+w^2*the_aux^2))*sqrt(1+w^2*the_aux^2);
        [the_w,the_val] = fmincon(the_obj, M(i_M)^2/(1+M(i_M)^2), [], [], [], [], 0, 1, [], opts);
        weight_opt_aware(i_M,i_r) = the_w;
        length_opt_aware(i_M,i_r) = the_val*se_ratio(i_r)/z;

    end

    % Worst-case asy. probability of the joint event:
    % (a) Hausman test fails to reject correct VAR specification
    % AND (b) conventional VAR t-test rejects true parameter
    [the_b,the_val] = fmincon(@(b) -prob_fct(the_aux,b), the_aux, [], [], [], [], 1e-9, Inf, [], opts);
    wc_prob(i_r) = exp(-the_val);
    wc_bias(i_r) = the_b;

end

length_var_aware = (cv_var_aware/z).*se_ratio; % Relative length of bias-aware VAR CI


%% Plot

% Worst-case coverage of VAR CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,wc_cov(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
end
hold off;
the_xlim = xlim;
line(the_xlim,[1 1]*(1-signif),'Color','k','LineStyle','-'); % Mark nominal coverage level on vertical axis
xlim(the_xlim);
ylim([0 1]);
label_legend(M);

% Relative length of bias-aware VAR CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,length_var_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
end
hold off;
the_xlim = xlim;
line(the_xlim,[1 1],'Color','k','LineStyle','-'); % Mark 1 on vertical axis
xlim(the_xlim);
ylim(plot_ylim_length);
label_legend(M);

% Weight on LP in optimal bias-aware CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,weight_opt_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
end
hold off;
ylim([0 1]);
label_legend(M);

% Relative length of optimal bias-aware CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,length_opt_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
end
hold off;
ylim([0 1]);
label_legend(M);

% Worst-case joint probability
figure('Units','inches','Position',plot_pos);
plot(se_ratio,wc_prob,'Color','k','LineWidth',2);
the_xlim = xlim;
line(the_xlim,[1 1]*signif,'Color','k','LineStyle',':'); % Mark nominal significance level on vertical axis
xlim(the_xlim);
ylim([0 1]);
label_legend();


%% Plotting function

function label_legend(varargin)
    xlabel('$\sqrt{\mathrm{aVar}(\hat{\delta}_h)/\mathrm{aVar}(\hat{\beta}_h)}$','Interpreter','Latex');
    set(gca,'FontSize',12);
    if ~isempty(varargin)
        legend(strcat('M=', arrayfun(@num2str, varargin{1}, 'UniformOutput', 0)),'Location','SouthEast','FontSize',12);
    end
end