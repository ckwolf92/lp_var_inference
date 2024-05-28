% delete SW_Model.log
% delete SW_Model.m
% delete SW_Model_dynamic.m
% delete SW_Model_results.mat
% delete SW_Model_set_auxiliary_variables.m
% delete SW_Model_static.m
% rmdir SW_Model/Output
% rmdir('SW_Model','s')
% load polfunction
% clc

try
    delete SW_Model.log
    delete SW_Model/Output/SW_Model_results.mat
    rmdir SW_Model/Output
    rmdir('SW_Model','s')
    rmdir('+SW_Model','s')
    delete +SW_Model
catch
end
load polfunction
clc