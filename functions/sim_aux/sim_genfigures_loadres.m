% Load results files for sim_genfigures.m

load(load_filename);

if contains(dgp_type, 'varma')

    % pick out indices of selected DGPs
    p_sel      = intersect(p_sel, dgp.dgps(1,:));
    numdgp_sel = length(p_sel);
    dgp_sel    = zeros(1,numdgp_sel);

    for i=1:length(p_sel)
        dgp_sel(i) = find(dgp.dgps(1,:)==p_sel(i));
    end

elseif contains(dgp_type, 'arma')
    
    % pick out indices of selected DGPs
    numdgp_sel = length(thetas_sel) * length(rhos_sel);
    dgp_sel    = zeros(1,numdgp_sel);
    
    for i=1:length(rhos_sel)
        for j=1:length(thetas_sel)
            dgp_indx = i + (j-1) * length(rhos_sel);
            dgp_sel(dgp_indx) = find(dgp.dgps(1,:)==rhos_sel(i) & dgp.dgps(2,:)==thetas_sel(j));
        end
    end

end