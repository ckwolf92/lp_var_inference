function [y, shock, date, yname, shockname] = clean_monetary(data_raw)

shockname = 'FF4_TC';                       % Shock variable name
yname      = {'LIP', 'LCPI', 'GS1', 'EBP'};  % Outcome variable names

df      = data_raw;
df.date = clean_date(data_raw.DATES);

% Estimation sample
keep = find(year(df.date) == 1990 & month(df.date)==1):...
       find(year(df.date) == 2012 & month(df.date)==6);
df      = df(keep, :);  % Set date
y       = df{:, yname};
shock   = df{:, shockname}; 
date    = df.date;


end