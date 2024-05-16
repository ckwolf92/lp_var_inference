% sim_genfigures_armaplot.m  ARMA: Select DGPs, CI procedures, and line specs for
%                            sim_genfigures.m

switch dgp_type(5:end)
    case '_b'
        thetas_sel = [0,0.25];
        rhos_sel   = [0.9];
        proc_names = {'AR', 'AR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
                 1 4;
                 2 1;
                 2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
    case '_'
        thetas_sel = [0,0.25];
        rhos_sel   = [0.9];
        proc_names = {'AR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    otherwise
        error('Invalid dgp_type...')
end