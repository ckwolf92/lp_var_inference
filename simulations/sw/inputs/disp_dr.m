function disp_dr(M_,options_,dr,order,var_list)
%disp_dr(M_,options_,dr,order,var_list)
% Display the decision rules
%
% INPUTS
%  - M         [structure]   storing the model information
%  - options   [structure]   storing the options
% - dr         [struct]      decision rules.
% - order      [integer]     order of approximation.
% - var_list   [cell]        list of endogenous variables for which the decision rules should be printed.
%
% OUTPUTS
% none

% Copyright © 2001-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

if order~=1 && M_.hessian_eq_zero
    order = 1;
    warning('disp_dr: using order = 1 because Hessian is equal to zero');
end

nx =size(dr.ghx,2);
nu =size(dr.ghu,2);
k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
klag = dr.kstate(k,[1 2]);
k1 = dr.order_var;

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

nvar = length(var_list);

ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list{i}, M_.endo_names(k1), 'exact');
    if isempty(i_tmp)
        disp(var_list{i});
        error ('One of the variable specified does not exist') ;
    else
        ivar(i) = i_tmp;
    end
end

% get length of display strings
header_label_length=16; %default
for ii=1:length(ivar)
    header_label_length = max(header_label_length,length(M_.endo_names{k1(ivar(ii))})+2);
end
if options_.loglinear
    header_label_length=header_label_length+5;
end
header_label_format  = sprintf('%%%ds',header_label_length);
value_format_float  = sprintf('%%%d.6f',header_label_length);
value_format_overflow  = sprintf('%%%d.6e',header_label_length);
value_format_zero  = sprintf('%%%dd',header_label_length);

% account for additional characters introduced by auxiliary variables
if ~isempty(M_.aux_vars)
    aux_vars_type = [M_.aux_vars.type];
    if any(aux_vars_type==4)
        aux_var_additional_characters=17;
    else
        aux_var_additional_characters=3;
    end
else
    aux_var_additional_characters=0;
end

if options_.loglinear
    var_name_width = max([cellofchararraymaxlength(M_.endo_names(k1(ivar)))+5, cellofchararraymaxlength(M_.exo_names)]);
else
    var_name_width = max([cellofchararraymaxlength(M_.endo_names(k1(ivar))), cellofchararraymaxlength(M_.exo_names)]);
end
%deal with covariances
if order > 1
    var_name_width=max(2*(var_name_width+aux_var_additional_characters)+2,20); %account for covariances, separated by comma
else
    var_name_width=max(var_name_width+aux_var_additional_characters,20);
end
label_format  = sprintf('%%-%ds',var_name_width);


%% start displayimg
disp('POLICY AND TRANSITION FUNCTIONS')
% variable names
str = char(32*ones(1,var_name_width));
for i=1:nvar
    if options_.loglinear
        str = [str sprintf(header_label_format, ['log(',M_.endo_names{k1(ivar(i))},')'])];
    else
        str = [str sprintf(header_label_format, M_.endo_names{k1(ivar(i))})];
    end
end
disp(str);
%
% constant
%
str=sprintf(label_format,'Constant');

decision = [];

flag = 0;
for i=1:nvar
    x = dr.ys(k1(ivar(i)));
    if order > 1
        x = x + dr.ghs2(ivar(i))/2;
    end
    [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,value_format_overflow,header_label_length,flag,options_);

    decision = [decision x];
end
if flag
    disp(str)
end
if order > 1
    decisiontemp = [];
    str = sprintf(label_format,'(correction)');
    flag = 0;
    for i=1:nvar
        x = dr.ghs2(ivar(i))/2;
        [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,value_format_overflow,header_label_length,flag,options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];

end
%
% ghx
%
for k=1:nx
    decisiontemp = [];
    flag = 0;
    if isfield(dr,'state_var')
        str1 = subst_auxvar(dr.state_var(k),-1,M_);
    else
        str1 = subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2,M_);
    end
    if options_.loglinear
        str = sprintf(label_format,['log(',str1,')']);
    else
        str = sprintf(label_format,str1);
    end
    for i=1:nvar
        x = dr.ghx(ivar(i),k);
        [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,value_format_overflow,header_label_length,flag,options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];
end
%
% ghu
%
for k=1:nu
    decisiontemp = [];
    flag = 0;
    str = sprintf(label_format, M_.exo_names{k});
    for i=1:nvar
        x = dr.ghu(ivar(i),k);
        [str,flag] = get_print_string(str, x, value_format_zero, value_format_float, value_format_overflow, header_label_length, flag, options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];
end

if order > 1
    % ghxx
    for k = 1:nx
        for j = 1:k
            decisiontemp = [];
            flag = 0;
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2,M_), ...
                           subst_auxvar(k1(klag(j,1)),klag(j,2)-M_.maximum_lag-2,M_));
            str = sprintf(label_format, str1);
            for i=1:nvar
                if k == j
                    x = dr.ghxx(ivar(i),(k-1)*nx+j)/2;
                else
                    x = dr.ghxx(ivar(i),(k-1)*nx+j);
                end
                [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,value_format_overflow,header_label_length,flag,options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
            decision = [decision;decisiontemp];
        end
    end
    %
    % ghuu
    %
    for k = 1:nu
        for j = 1:k
            decisiontemp = [];
            flag = 0;
            str = sprintf(label_format, [M_.exo_names{k} ',' M_.exo_names{j}]);
            for i=1:nvar
                if k == j
                    x = dr.ghuu(ivar(i),(k-1)*nu+j)/2;
                else
                    x = dr.ghuu(ivar(i),(k-1)*nu+j);
                end
                [str,flag]=get_print_string(str, x, value_format_zero, value_format_float, value_format_overflow, header_label_length, flag, options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
            decision = [decision;decisiontemp];
        end
    end
    %
    % ghxu
    %
    for k = 1:nx
        for j = 1:nu
            decisiontemp = [];
            flag = 0;
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2,M_), M_.exo_names{j});
            str = sprintf(label_format,str1);
            for i=1:nvar
                x = dr.ghxu(ivar(i),(k-1)*nu+j);
                [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,value_format_overflow,header_label_length,flag,options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
            decision = [decision; decisiontemp];
        end
    end
end
save polfunction decision
end


function [str,flag]=get_print_string(str, x, value_format_zero, value_format_float, value_format_overflow, max_length, flag, options_)
if abs(x) >= options_.dr_display_tol
    flag = 1;
    temp=sprintf(value_format_float, x);
    if length(temp)>max_length
        temp=sprintf(value_format_overflow, x);
    end
    str = [str temp];
else
    str = [str sprintf(value_format_zero, 0)];
end
end
