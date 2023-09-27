%% RC results supposed to be saved in the folder 'model/out_rc/' (set by models_RCexe_script.m); the results have then been  moved to '/results_model_rc' 
clear all;

path_tmp = split(pwd, '/model'); path_home = path_tmp{1}; 
addpath([path_home, '/model/func']); addpath([path_home, '/model/lib_ver101']); 

mode_lapse  = 1; 
% imod = 1; 'BMBU'
% imod = 2; 'RLVU'
% imod = 3; 'HYBR'
% imod = 5; 'ZFIX'};

for imod = 5
    models_RCexe_script(imod, mode_lapse, path_home); 

end

