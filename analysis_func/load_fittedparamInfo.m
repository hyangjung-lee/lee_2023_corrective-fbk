

function [fParam] =  load_fittedparamInfo(imod, numSub, path_loc )
%% FF best-fitting parameters

for iSub = 1 : numSub
    load([path_loc, '/results_fits/BMBU/fit_BMBU_subj', num2str(iSub), '.mat'])
fit_Param.BMBU(iSub,:) = fitResults.fitted_params{1};
end

%% RL best-fitting parameters  

for iSub = 1 : numSub
    load([path_loc, '/results_fits/RLVU/fit_RLVU_subj', num2str(iSub), '.mat'])
fit_Param.RLVU(iSub,:) = fitResults.fitted_params{1};
end
%% HYBRID best-fitting parameters 


for iSub = 1 : numSub
    load([path_loc, '/results_fits/HYBR/fit_HYBR_subj', num2str(iSub), '.mat'])
fit_Param.HYBR(iSub,:) = fitResults.fitted_params{1};
end



%% Fixed best-fitting parameters 


for iSub = 1 : numSub
    load([path_loc, '/results_fits/CFIX/fit_CFIX_subj', num2str(iSub), '.mat'])
fit_Param.Fixed(iSub,:) = fitResults.fitted_params{1};
end


%% Base best-fitting parameters 


for iSub = 1 : numSub
    load([path_loc, '/results_fits/ZFIX/fit_ZFIX_subj', num2str(iSub), '.mat'])
fit_Param.Base(iSub,:) = fitResults.fitted_params{1};
end




if (imod == 1)
    fParam = fit_Param.BMBU;
elseif (imod == 2)
     fParam = fit_Param.RLVU;
    
elseif (imod == 3)
     fParam = fit_Param.HYBR;
     
elseif (imod == 4)
     fParam = fit_Param.Fixed;
     
     
elseif (imod == 5)
    fParam = fit_Param.Base;
    
end

end