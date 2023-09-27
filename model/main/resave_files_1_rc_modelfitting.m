
clear all

%% model BMBU

save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_rc/bads-fitResults_randstrt30/mainData_numBMBU';
mkdir(save_path)

% mainData_numBMBU -> (1) FitBMBU 
info_fitWwhich = 'FitBMBU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numBMBU'; 
mkdir(final_path)
for iSub = 1 : 30


load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_FitBMBU_toSimBMBU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'], 'fitResults');

end



% mainData_numBMBU -> (2) FitRLVU 
info_fitWwhich = 'FitRLVU'; 
final_path = [save_path, '/', info_fitWwhich , '/'] 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numBMBU'; 
mkdir(final_path)
for iSub = 1 : 30


load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_RLVU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'], 'fitResults');

end



% mainData_numBMBU -> (3) FitHYBR 
info_fitWwhich = 'FitHYBR'; 
final_path = [save_path, '/', info_fitWwhich , '/'] 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numBMBU'; 
mkdir(final_path)
for iSub = 1 : 30


load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_HYBR,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'], 'fitResults');

end


% mainData_numBMBU -> (3) FitZFIX
info_fitWwhich = 'FitZFIX'; 
final_path = [save_path, '/', info_fitWwhich , '/'] 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numBMBU'; 
mkdir(final_path)
for iSub = 1 : 30


load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_ZFIX,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimBMBU', num2str(iSub), '.mat'], 'fitResults');

end



%% model HYBR
save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_rc/bads-fitResults_randstrt30/mainData_numHYBR';
mkdir(save_path)

% mainData_numHYBR -> (1) FitBMBU 
info_fitWwhich = 'FitBMBU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numHYBR'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_FitBMBU_toSimHYBR', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'], 'fitResults');

  
  
end




% mainData_numHYBR -> (2) FitRLVU 
info_fitWwhich = 'FitRLVU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numHYBR'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'], 'fitResults');

  
  
end





% mainData_numHYBR -> (3) FitHYBR
info_fitWwhich = 'FitHYBR'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numHYBR'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'], 'fitResults');

  
  
end






% mainData_numHYBR -> (4) FitZFIX
info_fitWwhich = 'FitZFIX'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numHYBR'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimHYBR', num2str(iSub), '.mat'], 'fitResults');

  
  
end


%% model RLVU
save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_rc/bads-fitResults_randstrt30/mainData_numRLVU';
mkdir(save_path)

% mainData_numRLVU -> (1) FitBMBU 
info_fitWwhich = 'FitBMBU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numRLVU'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_FitBMBU_toSimRLVU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'], 'fitResults');

  
  
end



% mainData_numRLVU -> (2) FitRLVU
info_fitWwhich = 'FitRLVU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numRLVU'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'], 'fitResults');

  
  
end


% mainData_numRLVU -> (3) FitHYBR
info_fitWwhich = 'FitHYBR'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numRLVU'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'], 'fitResults');

  
  
end



% mainData_numRLVU -> (4) FitZFIX
info_fitWwhich = 'FitZFIX'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numRLVU'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimRLVU', num2str(iSub), '.mat'], 'fitResults');

  
  
end
%% model ZFIX
save_path = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_rc/bads-fitResults_randstrt30/mainData_numZFIX';
mkdir(save_path)

% mainData_numZFIX -> (1) FitBMBU 
info_fitWwhich = 'FitBMBU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numZFIX'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_FitBMBU_toSimZFIX', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'], 'fitResults');

  
  
end




% mainData_numZFIX -> (2) FitRLVU 
info_fitWwhich = 'FitRLVU'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numZFIX'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'], 'fitResults');

  
  
end




% mainData_numZFIX -> (3) FitHYBR 
info_fitWwhich = 'FitHYBR'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numZFIX'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'], 'fitResults');

  
  
end




% mainData_numZFIX -> (4) FitZFIX 
info_fitWwhich = 'FitZFIX'; 
final_path = [save_path, '/', info_fitWwhich , '/']; 
orig_frompath = '/Users/eva/Dropbox/Lee_2023_plosb_corrective-fbk/results_model_recovery/RCresults_final/FF/bads-fitResults_randstrt30/mainData_numZFIX'; 
mkdir(final_path)




for iSub = 1 : 30
    
    
load([orig_frompath, '/', info_fitWwhich, '/numerical_interim_', info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'])


fitResults.fitted_params = vbmcResults(iSub).fitted_params;
fitResults.fval = vbmcResults(iSub).fval;
% fitResults.exitflag = vbmcResults(iSub).exitflag;
% fitResults.output = vbmcResults(iSub).output;
% fitResults.output.function = '@(theta)ibslike(@fitting_BMBU,theta,R,fInput.designM,options_ibs,fInput)';

save([final_path, info_fitWwhich, '_toSimZFIX', num2str(iSub), '.mat'], 'fitResults');

  
  
end

