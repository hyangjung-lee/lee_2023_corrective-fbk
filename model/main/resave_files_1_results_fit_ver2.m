
clear all; 
save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/res_fits/BMBU/';
load_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_fits/BMBU/';


%% model BMBU
mkdir([save_path])
for iSub = 1 : 30
load([load_path, 'fit_BMBU_subj', num2str(iSub), '.mat'])
fitResults = rmfield(fitResults, {'exitflag', 'output'})
% 
% fitResults.fitted_params = vbmcResults(iSub).fitted_params;
% fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,funInputFurther.designM,options_ibs,funInputFurther)';

save([save_path, 'fit_BMBU_subj', num2str(iSub), '.mat'], 'fitResults');

end



%% model HYBR

save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/res_fits/HYBR/';
load_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_fits/HYBR/';
mkdir([save_path])

for iSub = 1 : 30
load([load_path, 'fit_HYBR_subj', num2str(iSub), '.mat'])
  fitResults = rmfield(fitResults, {'exitflag', 'output'})
% 
% fitResults.fitted_params = vbmcResults(iSub).fitted_params;
% fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_HYBR,theta,R,funInputFurther.designM,options_ibs,funInputFurther)';

save([save_path, 'fit_HYBR_subj', num2str(iSub), '.mat'], 'fitResults');

  
  
end

%% model RLVU

save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/res_fits/RLVU/';
load_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_fits/RLVU/';
mkdir([save_path])


for iSub = 1 : 30
load([load_path, 'fit_RLVU_subj', num2str(iSub), '.mat'])
    fitResults = rmfield(fitResults, {'exitflag', 'output'})
% 
% fitResults.fitted_params = vbmcResults(iSub).fitted_params;
% fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_RLVU,theta,R,funInputFurther.designM,options_ibs,funInputFurther)';

save([save_path, 'fit_RLVU_subj', num2str(iSub), '.mat'], 'fitResults');

  
  
end

%% model CFIX

save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/res_fits/CFIX/';
load_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_fits/CFIX/';
mkdir([save_path])



for iSub = 1 : 30
load([load_path, 'fit_CFIX_subj', num2str(iSub), '.mat']); 
  fitResults = rmfield(fitResults, {'exitflag', 'output'})

% fitResults.fitted_params = vbmcResults(iSub).fitted_params;
% fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_CFIX,theta,R,funInputFurther.designM,options_ibs,funInputFurther)';

save([save_path, 'fit_CFIX_subj', num2str(iSub), '.mat'], 'fitResults');

  
  
end

%% model ZFIX
save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/res_fits/ZFIX/';
load_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_fits/ZFIX/';
mkdir([save_path])


for iSub = 1 : 30
load([load_path, 'fit_ZFIX_subj', num2str(iSub), '.mat']); 
fitResults = rmfield(fitResults, {'exitflag', 'output'})

% fitResults.fitted_params = vbmcResults(iSub).fitted_params;
% fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_ZFIX,theta,R,funInputFurther.designM,options_ibs,funInputFurther)';

save([save_path, 'fit_ZFIX_subj', num2str(iSub), '.mat'], 'fitResults');

  
  
end
