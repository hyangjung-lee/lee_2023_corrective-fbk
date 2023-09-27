
clear all; close all; 

cd ../
path_home = pwd; 
load([path_home, '/data/raw/exp_data.mat'],'subjCmat'); 
[nSub, ndim.cond] = size(subjCmat);
[ndim.run, ndim.trial] = size(subjCmat{1,1});
all_dataP   = ndim.trial*ndim.cond*ndim.run; 
addpath([path_home, '/analysis_func']); 

clear raw_pop pop bpop numpK fit_p
nmodels = 4; 
% recvry_perf = raw_pop = recvry_perf; pop = recvry_perf; 
raw_pop =  zeros(nSub, nmodels); pop = raw_pop;
name_fold = [path_home, '/results_model_rc'];

imod = 1; 
for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([name_fold, '/mainData_numBMBU/FitBMBU/FitBMBU_toSimBMBU', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;  
end

imod = 2; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numBMBU/FitRLVU/FitRLVU_toSimBMBU', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
       
end


imod = 3; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numBMBU/FitHYBR/FitHYBR_toSimBMBU', num2str(iSub),'.mat']) 
   
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
               
end



imod = 4; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numBMBU/FitZFIX/FitZFIX_toSimBMBU', num2str(iSub),'.mat']) 

    % eliminate stored extra param.
    fitResults.fitted_params{1} =  fitResults.fitted_params{1}(1); 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
        
end



addpath(genpath([path_home, '/lib/VBA-toolbox-master']))
 
model2compare = [1,2,3, 4]; 
lml_AIC = pop*(-0.5);  % Log model evidence, Cao et al., 

[poster, out] = VBA_groupBMC(lml_AIC(:, model2compare)');


row_1 = out.ep;

for iSub = 1 : nSub
   [~, pA(iSub,:)] = sort(pop(iSub,:) );
end

frow_1  = [length(find(pA(:,1) == 1))/size(pop,1), length(find(pA(:,1) == 2))/size(pop,1), length(find(pA(:,1) == 3))/size(pop,1)]; 

%% fitting models to RLVu data



imod = 1; 
for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numRLVU/FitBMBU/FitBMBU_toSimRLVU', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
    
end
imod = 2; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numRLVU/FitRLVU/FitRLVU_toSimRLVU', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
    
end


imod = 3; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numRLVU/FitHYBR/FitHYBR_toSimRLVU', num2str(iSub),'.mat']) 
   
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;   
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
    
end





imod = 4; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numRLVU/FitZFIX/FitZFIX_toSimRLVU', num2str(iSub),'.mat']) 

    % eliminate stored extra param.
    fitResults.fitted_params{1} =  fitResults.fitted_params{1}(1); 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
        
        
end




 
model2compare = [1,2,3,4]; 
lml_AIC = pop*(-0.5); % Log model evidence, Cao et al., 

[poster, out] = VBA_groupBMC(lml_AIC(:, model2compare)');


row_2 = out.ep;



for iSub = 1 : nSub
   [~, pA(iSub,:)] = sort(pop(iSub,:) );
end
frow_2  = [length(find(pA(:,1) == 1))/size(pop,1), length(find(pA(:,1) == 2))/size(pop,1), length(find(pA(:,1) == 3))/size(pop,1)]; 


%% fitting models to HYBR data



imod = 1; 
for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numHYBR/FitBMBU/FitBMBU_toSimHYBR', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;

end
imod = 2; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numHYBR/FitRLVU/FitRLVU_toSimHYBR', num2str(iSub),'.mat']) 
  
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;

end


imod = 3; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;
    load([ name_fold, '/mainData_numHYBR/FitHYBR/FitHYBR_toSimHYBR', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;

end




imod = 4; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numHYBR/FitZFIX/FitZFIX_toSimHYBR', num2str(iSub),'.mat']) 

     % eliminate stored extra param.
    fitResults.fitted_params{1} =  fitResults.fitted_params{1}(1); 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;   
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;


end



model2compare = [1,2,3, 4]; 
lml_AIC = pop*(-0.5);  % Log model evidence, Cao et al., 

[poster, out] = VBA_groupBMC(lml_AIC(:, model2compare)');

row_3 = out.ep;


for iSub = 1 : nSub
   [~, pA(iSub,:)] = sort(pop(iSub,:) );
end

frow_3  = [length(find(pA(:,1) == 1))/size(pop,1), length(find(pA(:,1) == 2))/size(pop,1), length(find(pA(:,1) == 3))/size(pop,1)]; 


%% fitting ZFIT Data

 clear raw_pop pop bpop numpK fit_p


imod = 1; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numZFIX/FitBMBU/FitBMBU_toSimZFIX', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
        
end



imod = 2; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numZFIX/FitRLVU/FitRLVU_toSimZFIX', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;

        
end

imod = 3; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([ name_fold, '/mainData_numZFIX/FitHYBR/FitHYBR_toSimZFIX', num2str(iSub),'.mat']) 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
     
end










imod = 4; 

for iSub = 1 : nSub
    nan_dataP = 0; num_dataP = all_dataP - nan_dataP;

    load([name_fold, '/mainData_numZFIX/FitZFIX/FitZFIX_toSimZFIX', num2str(iSub),'.mat']) 

    % eliminate stored extra param.
    fitResults.fitted_params{1} =  fitResults.fitted_params{1}(1); 

    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
    raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;     
    numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
        
end


model2compare = [1,2,3, 4]; 
lml_AIC = pop*(-0.5); % Log model evidence, Cao et al., 

[poster, out] = VBA_groupBMC(lml_AIC(:, model2compare)');




row_4 = out.ep;
clear pA pB

for iSub = 1 : nSub
   [~, pA(iSub,:)] = sort(pop(iSub,:) );
end

frow_4  = [length(find(pA(:,1) == 1))/size(pop,1), length(find(pA(:,1) == 2))/size(pop,1), length(find(pA(:,1) == 3))/size(pop,1)]; 



%%% Confusion matrix
ConfusM = [row_1;row_2;row_3; row_4]; 

%%% Inversion matrix
% InvM = ConfusM./sum(ConfusM, 1); 


%% S3 Fig
figure(131); clf
subplot(1,2,1);
a = imagesc(ConfusM); 
xlabel('fit model'); ylabel('simulated model'); 
colormap(pink); colorbar;caxis([0 1])
xlabel('fit model'); ylabel('simulated model');
set(gca,'fontsize', 20, 'ytick',1:4,'yticklabel',{},'ticklength', [ 0.02 0.02],'tickdir','out'); axis square

round(ConfusM, 3)
