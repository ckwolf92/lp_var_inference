% clean_tax()
%   Contains a MATLAB translation of jordatax.do from the replication package 
% % of Ramey (2016)
% Input:
%   - path [string] 
% Output:
%   - 

function [y, shock, date, yname, shockname] = clean_tax(data_raw)
df       = data_raw;


% Read data
df.date = clean_date(df.quarter);
df      = df(year(df.date) >= 1950 & year(df.date) <= 2007, :);  % Set date


df.rfedtaxrev = df.nfedtaxrev./df.pgdp;       % real federal tax revenue 
df.rgov       = df.ngov./df.pgdp;             % real government spending, using GDP deflator
df.rfedgov    = df.nfed./df.pgdp;             % real federal spending, using GDP deflator
df.taxy       = 100*df.nfedtaxrev./df.ngdp;   % average tax rate


% ------------------------------------------------------------------------- 
% CREATE PER CAPITA LOG VARIABLES;
% -------------------------------------------------------------------------
vars_lpc = {'rgdp', 'rcons', 'rinv', 'tothours', ...
            'rfedtaxrev', 'rgov', 'rfedgov'};
for j = 1:length(vars_lpc)
    df.(['l', vars_lpc{j}]) = log_percapita(df.(vars_lpc{j}), df.pop);
end


% RENAME TO SHORTER NAMES
df = renamevars(df, strcat('l', {'rgdp', 'rcons', 'rinv',...
                             'rfedtaxrev', 'rgov', 'rfedgov'}), ...
                strcat('l', {'y', 'c', 'i', 'h', 'tax', 'g'}));


% DEFINE QUADRATIC TREND;
df.t = (1:height(df))';
df.t2 = df.t.^2;

% BLANCHARD-PEROTTI DUMMY VARIABLE
df.dum75q2 = zeros(height(df), 1);
df.dum75q2(year(df.date) == 1975 & quarter(df.date) == 2) = 1;


% DEMEAN THE TAX SHOCK AS MERTENS-RAVEN DO;
df.rrtaxu_dm = zeros(height(df), 1);
df.rrtaxu_dm = df.rrtaxu - mean(df.rrtaxu(df.rrtaxu ~= 0), 'omitmissing');

df.dtaxy = [NaN; df.taxy(1:end-1)];
df.dltax = [NaN; df.ltax(1:end-1)];

% ------------------------------------------------------------------------- 
% Residualize and prepare output
% ------------------------------------------------------------------------- 

% Residualize a constant, time trend, and dummies
df = df(2:end, :);

control = [ones(size(df,1),1), df.t, df.t2, df.dum75q2];  % Control variables
M       = eye(size(control,1)) - control*inv(control'*control)*control';  % annihilator
y       = M*[df.lg df.ly, df.ltax];
shock   = M*df.rrtaxu; 

% Save
date      = df.date;
yname      = {'lg', 'ly', 'ltax'};
shockname = 'rrtaxu';

end
