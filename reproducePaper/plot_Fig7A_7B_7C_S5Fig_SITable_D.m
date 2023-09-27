clear all; close all;
cd ../

path_str.currpath   = pwd;

color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];
xxd = [-2:2]; Xxd = [ones(size(xxd')), xxd']; linealph = 1; 


%% obs PSE: plotting Fig 7A  & models' PSE: Fig 7B, C and S5 Fig  
%% %% plotting observed data (Fig 7A)

path_str.data       = [path_str.currpath, '/data']; 
load([path_str.data, '/behav_retro'],'obs_psepre'); load([path_str.data, '/behav_prosp'],'obs_psepost');
    r_sSam = obs_psepost.small_corr - obs_psepre.small_corr; 
    r_lSam = obs_psepost.large_corr - obs_psepre.large_corr;   
    r_sDiff = obs_psepost.small_incor - obs_psepre.small_incor; 
    r_lDiff = obs_psepost.large_incor - obs_psepre.large_incor;


mk_type = 'o';
stp = 0; 



for behav_mode = 1:3

if (behav_mode == 1)
str_att = 'prosp'; figure(71);  clf; 
    sim_sSam = obs_psepost.small_corr; 
    sim_lSam = obs_psepost.large_corr; 
    sim_sDiff = obs_psepost.small_incor; 
    sim_lDiff = obs_psepost.large_incor;

elseif (behav_mode == 2)
str_att = 'retro'; figure(72);  clf; 
    sim_sSam = obs_psepre.small_corr; 
    sim_lSam = obs_psepre.large_corr; 
    sim_sDiff = obs_psepre.small_incor; 
    sim_lDiff = obs_psepre.large_incor;

elseif (behav_mode == 3)
str_att = 'subtract'; figure(73);  clf; 
    sim_sSam = r_sSam; 
    sim_lSam = r_lSam; 
    sim_sDiff = r_sDiff; 
    sim_lDiff = r_lDiff;
    
end

    len_sim = size(sim_sSam, 1); 
    b_simAggsSam = nanmean(sim_sSam,1);
    b_simAgglSam = nanmean(sim_lSam,1); 
    b_simAggsDiff = nanmean(sim_sDiff,1); 
    b_simAgglDiff = nanmean(sim_lDiff,1); 

    e_sSam = nanstd(sim_sSam,1)/sqrt(len_sim);
    e_lSam = nanstd(sim_lSam,1)/sqrt(len_sim); 
    e_sDiff = nanstd(sim_sDiff,1)/sqrt(len_sim); 
    e_lDiff = nanstd(sim_lDiff,1)/sqrt(len_sim); 


axis_small = 1:3; axis_large = 3:5;
    
subplot(1,2,1);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]);xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),'-', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAggsSam(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAggsSam(:, axis_small),  e_sSam(:, axis_small),'color', choiceColor, 'linewidth', 2 ,'linestyle','none'); 


choiceColor = color_yellow;
h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),'-', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAgglSam(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)' + stp, b_simAgglSam(:, axis_large),  e_lSam(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 
set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'},'linewidth', 2,'ticklength', [0.03 003]);

subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),'-', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAggsDiff(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)'+ stp, b_simAggsDiff(:, axis_large),  e_sDiff(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

choiceColor = color_yellow; 
h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), '-', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAgglDiff(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAgglDiff(:, axis_small),  e_lDiff(:, axis_small),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 
set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'},'linewidth', 2,'ticklength', [0.03 003]);


% Non-veridical trials (dotted lines)
axis_small = 4:5; axis_large = 1:2;
 
subplot(1,2,1);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),':', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAggsSam(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAggsSam(:, axis_small),  e_sSam(:, axis_small),'color', choiceColor, 'linewidth', 2 ,'linestyle','none'); 

choiceColor = color_yellow;
h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),':', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAgglSam(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)' + stp, b_simAgglSam(:, axis_large),  e_lSam(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 


subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),':', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAggsDiff(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)'+ stp, b_simAggsDiff(:, axis_large),  e_sDiff(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

choiceColor = color_yellow; 
h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), ':', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAgglDiff(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAgglDiff(:, axis_small),  e_lDiff(:, axis_small),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 



set(gcf,'position', [ 2887         334         761         339]); 


yl = 1.2; 


if strcmp(str_att, 'prosp')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
elseif strcmp(str_att, 'retro')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
elseif strcmp(str_att, 'subtract')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}');  ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
end


end

obs_r_sSam = r_sSam;
obs_r_lSam = r_lSam;
obs_r_sDiff = r_sDiff;
obs_r_lDiff = r_lDiff;

%% %% plotting ex-post simulation (Fig 7B, C and S5 Fig)

path_str.data       = [path_str.currpath, '/results_expost']; 
str_model   = {'BMBU', 'HYBR', 'RLVU'}

stp = 0; 

for imodel = 1:3
    if (imodel == 1)        % World
        mk_type = 's';
    elseif (imodel == 2)    % Hybrid
        mk_type = 'd';
    elseif (imodel == 3)    % Value
             mk_type = 'v';   
    end
    
    load([path_str.data, '/', str_model{imodel}, '/expost_', str_model{imodel}, '_retro'],'psepre'); 
    load([path_str.data, '/', str_model{imodel}, '/expost_', str_model{imodel}, '_prosp'],'psepost'); 
    
    r_sSam = psepost.small_corr - psepre.small_corr; 
    r_lSam = psepost.large_corr - psepre.large_corr;   
    r_sDiff = psepost.small_incor - psepre.small_incor; 
    r_lDiff = psepost.large_incor - psepre.large_incor;


    
    

for behav_mode = 1 : 3

if (behav_mode == 1)
str_att = 'prosp'; figure(10*imodel + 71);  clf; 
    sim_sSam = psepost.small_corr; 
    sim_lSam = psepost.large_corr; 
    sim_sDiff = psepost.small_incor; 
    sim_lDiff = psepost.large_incor;

elseif (behav_mode == 2)
str_att = 'retro'; figure(10*imodel + 72);  clf; 
    sim_sSam = psepre.small_corr; 
    sim_lSam = psepre.large_corr; 
    sim_sDiff = psepre.small_incor; 
    sim_lDiff = psepre.large_incor;

elseif (behav_mode == 3)
str_att = 'subtract'; figure(10*imodel + 73);  clf; 
    sim_sSam = r_sSam; 
    sim_lSam = r_lSam; 
    sim_sDiff = r_sDiff; 
    sim_lDiff = r_lDiff;
    
end

    len_sim = size(sim_sSam, 1); 
    b_simAggsSam = nanmean(sim_sSam,1);
    b_simAgglSam = nanmean(sim_lSam,1); 
    b_simAggsDiff = nanmean(sim_sDiff,1); 
    b_simAgglDiff = nanmean(sim_lDiff,1); 

    e_sSam = nanstd(sim_sSam,1)/sqrt(len_sim);
    e_lSam = nanstd(sim_lSam,1)/sqrt(len_sim); 
    e_sDiff = nanstd(sim_sDiff,1)/sqrt(len_sim); 
    e_lDiff = nanstd(sim_lDiff,1)/sqrt(len_sim); 



axis_small = 1:3; axis_large = 3:5;
    
subplot(1,2,1);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]);xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),'-', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAggsSam(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAggsSam(:, axis_small),  e_sSam(:, axis_small),'color', choiceColor, 'linewidth', 2 ,'linestyle','none'); 


choiceColor = color_yellow;
h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),'-', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAgglSam(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)' + stp, b_simAgglSam(:, axis_large),  e_lSam(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 
set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'},'linewidth', 2,'ticklength', [0.03 003]);


subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),'-', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAggsDiff(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)'+ stp, b_simAggsDiff(:, axis_large),  e_sDiff(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

choiceColor = color_yellow; 
h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), '-', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAgglDiff(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAgglDiff(:, axis_small),  e_lDiff(:, axis_small),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 
set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'},'linewidth', 2,'ticklength', [0.03 003]);



% Non-veridical trials (dotted lines)
axis_small = 4:5; axis_large = 1:2;
 
subplot(1,2,1);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),':', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAggsSam(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAggsSam(:, axis_small),  e_sSam(:, axis_small),'color', choiceColor, 'linewidth', 2 ,'linestyle','none'); 

choiceColor = color_yellow;
h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),':', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAgglSam(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)' + stp, b_simAgglSam(:, axis_large),  e_lSam(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 


subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); xlim([-3 3]);

choiceColor = color_blue;
h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),':', 'color', choiceColor , 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_large,2)' , b_simAggsDiff(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_large,2)'+ stp, b_simAggsDiff(:, axis_large),  e_sDiff(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

choiceColor = color_yellow; 
h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), ':', 'color', choiceColor, 'linewidth', 4); h.Color(4) = linealph; 
plot(Xxd(axis_small,2)' , b_simAgglDiff(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
errorbar(Xxd(axis_small,2)'+ stp, b_simAgglDiff(:, axis_small),  e_lDiff(:, axis_small),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 




set(gcf,'position', [ 2887         334         761         339]); 



yl = 1.2; 


if strcmp(str_att, 'prosp')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
elseif strcmp(str_att, 'retro')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
elseif strcmp(str_att, 'subtract')
    subplot(1,2,1); xlabel('Stimulus_{\ittoi}'); ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
    subplot(1,2,2); xlabel('Stimulus_{\ittoi}');  ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}'); box off; set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});ylim([-yl yl]); 
end


end


end

%% %% statistical tests: model ex-post behavior vs. observed behavior (Table D, PSE comparisons)


bonferroni_corrected = 1; 
totR = NaN(1,3);  
for behav_mode = 1:3 

clear tmpResults_stats

test_type = 2 % 1: wilcoxon sign-rank test / 2: t-test

crit_stat(1) = 0.05;
crit_stat(2) = 0.01;
crit_stat(3) = 0.001;
if (bonferroni_corrected == 1)
    crit_stat = crit_stat./(3*( size(r_sSam, 2) + size(r_lSam, 2) +  size(r_sDiff, 2) + size(r_lDiff, 2))); % 60 tests =  20 datapoints x 3 cases (prospective + retrospective + subtracted)
end
if (behav_mode == 1)
    mod_sSam = psepost.small_corr;
    mod_lSam = psepost.large_corr;
    mod_sDiff = psepost.small_incor;
    mod_lDiff = psepost.large_incor;
     
    obs_sSam = obs_psepost.small_corr;
    obs_lSam = obs_psepost.large_corr;
    obs_sDiff = obs_psepost.small_incor;
    obs_lDiff = obs_psepost.large_incor;

elseif (behav_mode == 2)
    
    mod_sSam = psepre.small_corr;
    mod_lSam = psepre.large_corr;
    mod_sDiff = psepre.small_incor;
    mod_lDiff = psepre.large_incor;
    
    
    obs_sSam = obs_psepre.small_corr;
    obs_lSam = obs_psepre.large_corr;
    obs_sDiff = obs_psepre.small_incor;
    obs_lDiff = obs_psepre.large_incor;
        
elseif (behav_mode == 3)
    mod_sSam = r_sSam;
    mod_lSam = r_lSam;
    mod_sDiff = r_sDiff;
    mod_lDiff = r_lDiff;
     
    obs_sSam = obs_r_sSam;
    obs_lSam = obs_r_lSam;
    obs_sDiff = obs_r_sDiff;
    obs_lDiff = obs_r_lDiff;

end
clear wil wh
for j = 1 : 5
    if (test_type == 1)
    [wil(j), wh(j)] = signrank( mod_sSam(:,j) , obs_sSam(:,j));
    
    elseif (test_type == 2)
   [wh(j), wil(j)] = ttest( mod_sSam(:,j) , obs_sSam(:,j));
    end

end

pwil = wil;
pwil(pwil< crit_stat(3)) = 3;
pwil(pwil< crit_stat(2)) = 2;
pwil(pwil<crit_stat(1)) = 1; 
pwil(pwil>crit_stat(1) & pwil<1) = 0;

tmpResults_stats(1,:) = round( pwil );

%
for j = 1 : 5
      if (test_type == 1)
            [wil(j), wh(j)] = signrank( mod_lSam(:,j) , obs_lSam(:,j));
        elseif (test_type == 2)
            [wh(j), wil(j)] = ttest( mod_lSam(:,j) , obs_lSam(:,j));
      end
end

pwil = wil;
pwil(pwil< crit_stat(3)) = 3;
pwil(pwil< crit_stat(2)) = 2;
pwil(pwil<crit_stat(1)) = 1; 
pwil(pwil>crit_stat(1) & pwil<1) = 0;

tmpResults_stats(2,:) = round( pwil );

for j = 1 : 5
    if (test_type == 1)
[wil(j), wh(j)] = signrank( mod_sDiff(:,j) , obs_sDiff(:,j));
    elseif (test_type == 2) 
   [wh(j), wil(j)] = ttest( mod_sDiff(:,j) , obs_sDiff(:,j));
    end

end



pwil = wil;
pwil(pwil< crit_stat(3)) = 3;
pwil(pwil< crit_stat(2)) = 2;
pwil(pwil<crit_stat(1)) = 1; 
pwil(pwil>crit_stat(1) & pwil<1) = 0;


tmpResults_stats(3,:) = round( pwil );

for j = 1 : 5
     if (test_type == 1)
[wil(j), wh(j)] = signrank( mod_lDiff(:,j) , obs_lDiff(:,j));
     elseif (test_type == 2)
   [wh(j), wil(j)] = ttest( mod_lDiff(:,j) , obs_lDiff(:,j));
     end

end




pwil = wil;
pwil(pwil< crit_stat(3)) = 3;
pwil(pwil< crit_stat(2)) = 2;
pwil(pwil<crit_stat(1)) = 1; 
pwil(pwil>crit_stat(1) & pwil<1) = 0;


tmpResults_stats(4,:) = round( pwil );


Results_stats = [ length(find( tmpResults_stats == 0)), length(find( tmpResults_stats == 1)), length(find( tmpResults_stats == 2)), length(find( tmpResults_stats == 3))]
totR(behav_mode) =  Results_stats(1);


end
sum(totR)
