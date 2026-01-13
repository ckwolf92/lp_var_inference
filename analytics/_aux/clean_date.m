function date = clean_date(date_raw)

% -------------------------------------------------------------------------
% Get increment 
% -------------------------------------------------------------------------

diff_dates = date_raw(2) - date_raw(1);

if abs(diff_dates - 1/12) < 1e-10
    month_shift = 1;
elseif abs(diff_dates - 1/4) < 1e-10
    month_shift = 3;
else
    error('Incorrect date input...')
end

% -------------------------------------------------------------------------
% Decimal to datetime: Robust to truncation error
% -------------------------------------------------------------------------

date = datetime(floor(date_raw), round((date_raw-floor(date_raw))*12)+month_shift, 1);

end