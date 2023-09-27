%% Fitting results supposed to be saved in the folder 'model/out/' (set by models_exe_script.m); the results have been moved to '/results_fits' 

clear all;

path_tmp = split(pwd, '/model'); path_home = path_tmp{1}; 
addpath([path_home, '/model/func']); addpath([path_home, '/model/lib_ver101']); 

mode_lapse  = 1; 

% imod = 1; 'BMBU'
% imod = 2; 'RLVU'
% imod = 3; 'HYBR'
% imod = 4; 'CFIX',
% imod = 5; 'ZFIX'};
for imod = 5
    models_exe_script(imod, mode_lapse, path_home); 
end

