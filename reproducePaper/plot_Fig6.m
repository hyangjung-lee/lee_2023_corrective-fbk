
clear all; close all; 
cd ../
path_str.currpath   = pwd;
path_str.functions  = [path_str.currpath, '/analysis_func'];
path_str.results    = [path_str.currpath, '/results_fits'];
path_str.data       = [path_str.currpath, '/data/raw']; 
addpath(path_str.functions);

str_model   = {'World', 'Hybrid', 'Value', 'Fixed','Base'}; nmod = length(str_model);
name_model  = {'BMBU', 'HYBR', 'RLVU', 'CFIX','ZFIX'}

% raw choice data
load([path_str.data , '/exp_data.mat'], 'subjCmat')
[nSub, ndim.cond] = size(subjCmat);
[ndim.run, ndim.trial] = size(subjCmat{1,1}); % iSub, iCond

raw_pop     = NaN(nSub, nmod); pop = NaN(nSub, nmod); bpop = NaN(nSub, nmod); numpK = NaN(1, nmod);
all_dataP   = ndim.trial*ndim.cond*ndim.run; 

imod = 1; str_model{imod}

fit_p = cell(size(str_model));
for iSub = 1:nSub
    nan_dataP = length(find(subjCmat{iSub, 1}(:) == 0))+ length(find(subjCmat{iSub, 2}(:) == 0)) + length(find(subjCmat{iSub, 3}(:) == 0)); num_dataP = all_dataP - nan_dataP;
    load([path_str.results, '/', name_model{imod}, '/fit_', name_model{imod} ,'_subj', num2str(iSub), '.mat'], 'fitResults');
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
        raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
        numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
end


imod = 2; str_model{imod}

fit_p = cell(size(str_model));
for iSub = 1:nSub
    nan_dataP = length(find(subjCmat{iSub, 1}(:) == 0))+ length(find(subjCmat{iSub, 2}(:) == 0)) + length(find(subjCmat{iSub, 3}(:) == 0)); num_dataP = all_dataP - nan_dataP;
    load([path_str.results, '/', name_model{imod}, '/fit_', name_model{imod} ,'_subj', num2str(iSub), '.mat'], 'fitResults');
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
        raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
        numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
end


imod = 3; str_model{imod}

fit_p = cell(size(str_model));
for iSub = 1:nSub
    nan_dataP = length(find(subjCmat{iSub, 1}(:) == 0))+ length(find(subjCmat{iSub, 2}(:) == 0)) + length(find(subjCmat{iSub, 3}(:) == 0)); num_dataP = all_dataP - nan_dataP;
    load([path_str.results, '/', name_model{imod}, '/fit_', name_model{imod} ,'_subj', num2str(iSub), '.mat'], 'fitResults');
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
        raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
        numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
end


imod = 4; str_model{imod}

fit_p = cell(size(str_model));
for iSub = 1:nSub
    nan_dataP = length(find(subjCmat{iSub, 1}(:) == 0))+ length(find(subjCmat{iSub, 2}(:) == 0)) + length(find(subjCmat{iSub, 3}(:) == 0)); num_dataP = all_dataP - nan_dataP;
    load([path_str.results, '/', name_model{imod}, '/fit_', name_model{imod} ,'_subj', num2str(iSub), '.mat'], 'fitResults');
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
        raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
        numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
end

imod = 5; str_model{imod}

fit_p = cell(size(str_model));
for iSub = 1:nSub
    nan_dataP = length(find(subjCmat{iSub, 1}(:) == 0))+ length(find(subjCmat{iSub, 2}(:) == 0)) + length(find(subjCmat{iSub, 3}(:) == 0)); num_dataP = all_dataP - nan_dataP;
    load([path_str.results, '/', name_model{imod}, '/fit_', name_model{imod} ,'_subj', num2str(iSub), '.mat'], 'fitResults');
    [gof.llkd, gof.aicc, gof.bic, gof.prmbest] = load_gof_models(fitResults, num_dataP); 
        raw_pop(iSub, imod) = gof.llkd;                 pop(iSub, imod)     = gof.aicc;    
        numpK(:, imod)      = length( gof.prmbest );    fit_p{imod}(iSub,:) = gof.prmbest;
end




model2compare = [1,2,3,4,5];
subjlist = 1:30;



raw_pop = raw_pop(subjlist,model2compare);
pop     = pop(subjlist,model2compare);




%% plotting Fig 6 B, C
%%%% compared to ZFIX 
wdlength = 0.5; 
rel_aic =  pop(:,length(model2compare)) - pop;
rel_bic =  bpop(:,length(model2compare)) - bpop;
rel_ll =   ( - raw_pop(:,length(model2compare)) ) - (- raw_pop);


clear rev_cmap_model rev_model_axis rev_str_model; 
rev_model_axis = [ 4, 3, 1, 2]; 
rev_str_model = {'Fixed', 'Value', 'World', 'Hybrid'} ;

rel_aic = rel_aic( :, rev_model_axis); 
rel_bic = rel_bic( :, rev_model_axis); 
rel_ll = rel_ll( :, rev_model_axis); 

rev_cmap_model(1,:) = [205 205 205]/255;
rev_cmap_model(2,:) = [198 178 191]/255;
rev_cmap_model(3,:) = [191 200 179]/255;
rev_cmap_model(4,:) = [197 179 234]/255;

figure(61);clf;  % set(gcf,'position', [ 2004         524         382         210]); 
subplot(1,2,1); 
for imodel = 1 : length(rev_model_axis)
    h=bar(imodel, mean(rel_ll(:, imodel )  ), 'EdgeColor','None', 'FaceColor',rev_cmap_model(imodel,:) ); hold on; h.BarWidth= wdlength;
    errorbar( imodel, mean(rel_ll(:, imodel )  ), std(rel_ll(:, imodel )  )/sqrt(nSub), 'k-', 'linew', 0.75); 
end
% Lines for guide
imodel = 1; 
plot( [0 size(rel_ll, 2) + 1], [ 1 1] * mean(rel_ll(:, imodel )  ), '--','color', rev_cmap_model(imodel,:), 'linew',2); 
imodel = 4; 
plot( [0 size(rel_ll, 2) + 1], [ 1 1] * mean(rel_ll(:, imodel )  ), '--','color', rev_cmap_model(imodel,:), 'linew',2); 

set(gca, 'fontsize',20, 'xtick', 1: length(rev_model_axis) , 'xticklabel', rev_str_model); xlim([0 length(rev_model_axis) + 1])
box off; set(gca,'XColor',[ 1 1 1], 'xtick',[],'linew', 2,'ticklength',[0.02 0.02])
ylabel({'\delta (LL)'}); 

subplot(1,2,2); 

for imodel = 1 :  length(rev_model_axis)
    h=bar(imodel, mean(rel_aic(:, imodel )  ), 'EdgeColor','None', 'FaceColor',rev_cmap_model(imodel,:) ); hold on; h.BarWidth= wdlength;
    errorbar( imodel, mean(rel_aic(:, imodel )  ), std(rel_aic(:, imodel )  )/sqrt(nSub), 'k-', 'linew', 0.75); 
end
% Lines for guide
imodel = 1; 
plot( [0 size(rel_aic, 2) + 1], [ 1 1] * mean(rel_aic(:, imodel )  ), '--','color', rev_cmap_model(imodel,:), 'linew',2); 
imodel = 4; 
plot( [0 size(rel_aic, 2) + 1], [ 1 1] * mean(rel_aic(:, imodel )  ), '--','color', rev_cmap_model(imodel,:), 'linew',2); 


set(gca, 'fontsize',20, 'xtick', 1: length(rev_model_axis) , 'xticklabel', rev_str_model); xlim([0 length(rev_model_axis) + 1])
box off; set(gca,'XColor',[ 1 1 1], 'xtick',[],'linew', 2,'ticklength',[0.02 0.02])
ylabel({'\Delta (AICc)'}); 



%% plotting Fig 6 D
addpath(genpath([path_str.currpath, '/lib/VBA-toolbox-master']))
 
lml_AIC = pop*(-0.5);  % Log model evidence, Cao et al., 

[poster, out] = VBA_groupBMC(lml_AIC(:, model2compare)');

wdlength = 0.5;
clear rev_cmap_model rev_model_axis rev_str_model; 
rev_cmap_model(1,:) = [1 1 1];
rev_cmap_model(2,:) = [205 205 205]/255;
rev_cmap_model(3,:) = [198 178 191]/255;
rev_cmap_model(4,:) = [191 200 179]/255;
rev_cmap_model(5,:) = [197 179 234]/255;


rev_model_axis = [5,  4, 3, 1, 2]; 
chance_lev = 1/length(rev_model_axis); 

rev_str_model = {'Base', 'Fixed', 'Value', 'World', 'Hybrid'} ;

out.Ef = out.Ef( rev_model_axis, :); out.Vf = out.Vf( rev_model_axis, :); 
out.pxp = out.pxp(:,  rev_model_axis); 


figure(62);clf; 
subplot(1,2,1); 
    plot([0, length(rev_model_axis)+1], [chance_lev chance_lev], 'k--','linew',1) ; hold on; box off; set(gca,'XColor',[ 1 1 1], 'xtick',[]);
    set(gca,'linew', 2,'ticklength',[0.02 0.02])
    for imodel = 1 : length(rev_model_axis)
        if (imodel == 1)
            h=bar(imodel, out.Ef(imodel )  , 'EdgeColor','k', 'FaceColor',rev_cmap_model(imodel,:) ); hold on; 
            h.BarWidth = wdlength; 
        else
            h=bar(imodel, out.Ef(imodel )  , 'EdgeColor','None', 'FaceColor',rev_cmap_model(imodel,:) ); 
            h.BarWidth = wdlength;
        end
        errorbar( imodel, out.Ef(imodel ) , sqrt(out.Vf(imodel )   ), 'k-', 'linew', 0.75); 
    end
    if (out.bor < 0.05)
        for imodel = 1 : length(rev_model_axis)
            plot( imodel, out.pxp(imodel), 'o','markersize', 9, 'markerfacecolor',rev_cmap_model(imodel,:),'markeredgecolor','k'); 
            plot( imodel+[-0.2 0.2], [ 1 1 ]*out.pxp(imodel), 'k-','linew', 1.5); 
        end
    end
ylabel('Expected posterior (p)'); ylim([0 1]); set(gca, 'fontsize', 20); set(gcf,'position', [ 2004         468         382         266]); 
out.pxp

%% t-tests

% paired tests with only the two models (2)

% Value (vs. Fixed)
[h, p,ci, stats] = ttest(pop(:,3), pop(:,4),'tail','left')

% Value (vs. Hybrid)
[h, p,ci, stats] = ttest(pop(:,3), pop(:,2),'tail','right')


% Value (vs. World)
[h, p, ci, stats] = ttest(pop(:,3), pop(:,1),'tail','right')


% World (vs. Fixed)
[h, p, ci, stats] = ttest(pop(:,1), pop(:,4),'tail','left')


% World (vs. Hybrid)
[h, p, ci, stats] = ttest(pop(:,1), pop(:,2),'tail','left')

