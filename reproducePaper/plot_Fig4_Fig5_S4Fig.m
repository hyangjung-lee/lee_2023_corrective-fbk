clear all ; close all; 
cd ../


color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];

color_bp1 = [154, 208, 229]/255;
color_bp2 = [201, 221, 230]/255;

color_yp2 =  [237, 232, 196]/255;
color_yp1 =  [234, 223, 131]/255;

cmap_curves = [color_blue; color_bp1; color_bp2; color_yp2; color_yp1; color_yellow];


%% Psychometric curves (conceptual)

figure(41); clf; set(gcf,'position', [ 224        1019         159         139]); 
plot([-5:0.05:5], normcdf([-5:0.05:5], 0, 1),'-','color', 'k','linew', 2); hold on; 
plot([-5:0.05:5], normcdf([-5:0.05:5], 1, 1),'-','color', color_blue,'linew', 2);
plot([-5:0.05:5], normcdf([-5:0.05:5], -1, 1),'-','color', color_yellow,'linew', 2);

set_str{1} = 'retro'; 
set_str{2} = []; 
set_str{3} = 'subtract'; 

%% BMBU (Fig 4 or 5 E, G , H)
xxd = [-2:2]; Xxd = [ones(size(xxd')), xxd'];
linealph = 1; 
mk_type = 's';
for iplot  = 1 : length(set_str)
    str_att = set_str{iplot}; 
    if (iplot == 1)
        load('./results_exante/BMBU/exante_BMBU_retro.mat')
    elseif (iplot  == 2)
        load('./results_exante/BMBU/exante_BMBU_prosp.mat')
    elseif (iplot == 3)
        load('./results_exante/BMBU/exante_BMBU_subtract.mat')
    end
  
axis_small = 1:3; axis_large = 3:5;
len_sim = size(sim_sSam, 1); 

b_simAggsSam = nanmean(sim_sSam,1);
b_simAgglSam = nanmean(sim_lSam,1); 
b_simAggsDiff = nanmean(sim_sDiff,1); 
b_simAgglDiff = nanmean(sim_lDiff,1); 

e_sSam = nanstd(sim_sSam,1)/sqrt(len_sim);
e_lSam = nanstd(sim_lSam,1)/sqrt(len_sim); 
e_sDiff = nanstd(sim_sDiff,1)/sqrt(len_sim); 
e_lDiff = nanstd(sim_lDiff,1)/sqrt(len_sim); 


figure(41+iplot);

stp = 0; 
subplot(1,2,1);
    plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
    plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); hold on; 

    choiceColor = cmap_curves(1:3,:); col_regime = color_blue; 
    h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),'-', 'color', col_regime , 'linewidth', 4); h.Color(4) = linealph; 

    for ic = 1 : 3
        h=plot(Xxd(axis_small(ic),2)' , b_simAggsSam(:, axis_small(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) ,'markeredgecolor', choiceColor(ic,:), 'linewidth', 1, 'markersize', 12); 
        errorbar(Xxd(axis_small(ic),2)'+ stp, b_simAggsSam(:, axis_small(ic)),  e_sSam(:, axis_small(ic)),'color', 'k', 'linewidth', 1 ,'linestyle','none'); 
    end
xlim([-3 3]); set(gca,'fontsize', 20,'xtick', -2:2, 'ticklength', [0.03 0.03]); xlabel('Stimulus_{\itt}_{\it-}_{\it1}'); ylabel('PSE_{\itt}')


choiceColor = cmap_curves(4:6,:); col_regime = color_yellow; 
    h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),'-', 'color', col_regime, 'linewidth', 4); h.Color(4) = linealph; 
    for ic = 1 : 3
        h=plot(Xxd(axis_large(ic),2)' , b_simAgglSam(:, axis_large(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) , 'markeredgecolor', choiceColor(ic,:) , 'linewidth', 1, 'markersize', 12); 
        errorbar(Xxd(axis_large(ic),2)' + stp, b_simAgglSam(:, axis_large(ic)),  e_lSam(:, axis_large(ic)),'color', 'k', 'linewidth', 1 ,'linestyle','none'); 
    end


set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'}); xlabel('Stimulus_{\itt}_{\it-}_{\it1}') 
xlim([-3 3]);  set(gca,'linewidth', 2,'ticklength', [0.03 003])
set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
xlabel('Stimulus_{\ittoi}') 

subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); hold on; 



choiceColor = flipud(cmap_curves(1:3,:)); col_regime = color_blue; 


h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),'-', 'color', col_regime , 'linewidth', 4); 
h.Color(4) = linealph; 


for ic = 1 : 3
h=plot(Xxd(axis_large(ic),2)' , b_simAggsDiff(:, axis_large(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) , 'markeredgecolor', choiceColor(ic,:) ,'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_large(ic),2)'+ stp, b_simAggsDiff(:, axis_large(ic)),  e_sDiff(:, axis_large(ic)),'color', 'k', 'linewidth', 1,'linestyle','none'); 
end


choiceColor = flipud( cmap_curves(4:6,:) ); col_regime = color_yellow; 


h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), '-', 'color', col_regime, 'linewidth', 4); 
h.Color(4) = linealph; 


for ic = 1 : 3
h=plot(Xxd(axis_small(ic),2)' , b_simAgglDiff(:, axis_small(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:)  , 'markeredgecolor', choiceColor(ic,:)  ,'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_small(ic),2)'+ stp, b_simAgglDiff(:, axis_small(ic)),  e_lDiff(:, axis_small(ic)),'color', 'k', 'linewidth', 1,'linestyle','none'); 
end

xlim([-3 3]); set(gca,'fontsize', 25,'xtick', -2:2, 'ticklength', [0.03 0.03]); 
set(gca,'linewidth', 2,'ticklength', [0.03 003]); set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
xlabel('Stimulus_{\ittoi}') 



set(gcf,'position', [ 2887         334         761         339]); 

yl = 1;
subplot(1,2,1); 
set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});
ylim([-yl yl]); 


if (isempty(str_att))    
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off;
    subplot(1,2,2); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; 
elseif strcmp(str_att, 'retro')
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
    subplot(1,2,2);ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
elseif strcmp(str_att, 'subtract')
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');box off;
    subplot(1,2,2);  ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');box off;
end
subplot(1,2,2); 
set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});
ylim([-yl yl]); 


end


%% RLVU (Fig 4 or 5 B, D, H)

mk_type = 'v';
for iplot  = 1 : 3
    figure(45+iplot);clf; 

    str_att = set_str{iplot}; 
    if (iplot == 1)
load('./results_exante/RLVU/exante_RLVU_retro.mat')

    elseif (iplot  == 2)
load('./results_exante/RLVU/exante_RLVU_prosp.mat')

    elseif (iplot == 3)
load('./results_exante/RLVU/exante_RLVU_subtract.mat')


    end





  
axis_small = 1:3; axis_large = 3:5;
len_sim = size(sim_sSam, 1); 

b_simAggsSam = nanmean(sim_sSam,1);
b_simAgglSam = nanmean(sim_lSam,1); 
b_simAggsDiff = nanmean(sim_sDiff,1); 
b_simAgglDiff = nanmean(sim_lDiff,1); 

e_sSam = nanstd(sim_sSam,1)/sqrt(len_sim);
e_lSam = nanstd(sim_lSam,1)/sqrt(len_sim); 
e_sDiff = nanstd(sim_sDiff,1)/sqrt(len_sim); 
e_lDiff = nanstd(sim_lDiff,1)/sqrt(len_sim); 


figure(45+iplot);


 stp = 0; 
subplot(1,2,1);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); hold on; 

choiceColor = cmap_curves(1:3,:); col_regime = color_blue; 
h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),'-', 'color', col_regime , 'linewidth', 4); 
h.Color(4) = linealph; 

for ic = 1 : 3
h=plot(Xxd(axis_small(ic),2)' , b_simAggsSam(:, axis_small(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) ,'markeredgecolor', choiceColor(ic,:), 'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_small(ic),2)'+ stp, b_simAggsSam(:, axis_small(ic)),  e_sSam(:, axis_small(ic)),'color', 'k', 'linewidth', 1 ,'linestyle','none'); 
end

xlim([-3 3]); set(gca,'fontsize', 20,'xtick', -2:2, 'ticklength', [0.03 0.03]); xlabel('Stimulus_{\itt}_{\it-}_{\it1}'); ylabel('PSE_{\itt}')


choiceColor = cmap_curves(4:6,:); col_regime = color_yellow; 
h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),'-', 'color', col_regime, 'linewidth', 4); 
h.Color(4) = linealph; 
for ic = 1 : 3
h=plot(Xxd(axis_large(ic),2)' , b_simAgglSam(:, axis_large(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) , 'markeredgecolor', choiceColor(ic,:) , 'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_large(ic),2)' + stp, b_simAgglSam(:, axis_large(ic)),  e_lSam(:, axis_large(ic)),'color', 'k', 'linewidth', 1 ,'linestyle','none'); 


end


set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'}); xlabel('Stimulus_{\itt}_{\it-}_{\it1}') 
xlim([-3 3]);  set(gca,'linewidth', 2,'ticklength', [0.03 003])
set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
xlabel('Stimulus_{\ittoi}') 

subplot(1,2,2);
plot( [-3 3], [0 0 ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
plot([0 0], [-4 4], '-','color',[0.5 0.5 0.5]); hold on; 



choiceColor = flipud(cmap_curves(1:3,:)); col_regime = color_blue; 


h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),'-', 'color', col_regime , 'linewidth', 4); 
h.Color(4) = linealph; 


for ic = 1 : 3
h=plot(Xxd(axis_large(ic),2)' , b_simAggsDiff(:, axis_large(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:) , 'markeredgecolor', choiceColor(ic,:) ,'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_large(ic),2)'+ stp, b_simAggsDiff(:, axis_large(ic)),  e_sDiff(:, axis_large(ic)),'color', 'k', 'linewidth', 1,'linestyle','none'); 
end


choiceColor = flipud( cmap_curves(4:6,:) ); col_regime = color_yellow; 


h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), '-', 'color', col_regime, 'linewidth', 4); 
h.Color(4) = linealph; 


for ic = 1 : 3
h=plot(Xxd(axis_small(ic),2)' , b_simAgglDiff(:, axis_small(ic)),mk_type, 'markerfacecolor', choiceColor(ic,:)  , 'markeredgecolor', choiceColor(ic,:)  ,'linewidth', 1, 'markersize', 12); 
errorbar(Xxd(axis_small(ic),2)'+ stp, b_simAgglDiff(:, axis_small(ic)),  e_lDiff(:, axis_small(ic)),'color', 'k', 'linewidth', 1,'linestyle','none'); 
end

xlim([-3 3]); set(gca,'fontsize', 25,'xtick', -2:2, 'ticklength', [0.03 0.03]); 
set(gca,'linewidth', 2,'ticklength', [0.03 003]); set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
xlabel('Stimulus_{\ittoi}') 



set(gcf,'position', [ 2887         334         761         339]); 

yl = 1;
subplot(1,2,1); 
set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});
ylim([-yl yl]); 


if (isempty(str_att))    
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off;
    subplot(1,2,2); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; 
elseif strcmp(str_att, 'retro')
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
    subplot(1,2,2);ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
elseif strcmp(str_att, 'subtract')
    subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');box off;
    subplot(1,2,2);  ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');box off;
end
subplot(1,2,2); 
set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)});
ylim([-yl yl]); 


end


%% S4 Fig A, B


figure(141);clf; 


load('./data/raw/exp_data.mat'); 
nSub = size(subjCmat,1); perf_sub = zeros(nSub, 1); 
nTrial = size(subjCmat{1,1},2); 
for iSub = 1: nSub
    
    cd1 = mean(sum(subjCmat{iSub,1} == subjCLmat{iSub,1},2)/nTrial) ;
    cd2 = mean(sum(subjCmat{iSub,2} == subjCLmat{iSub,2},2)/nTrial) ;
    cd3 = mean(sum(subjCmat{iSub,3} == subjCLmat{iSub,3},2)/nTrial) ;
    perf_sub(iSub,:) =   mean([cd1, cd2, cd3] );
end

    
subplot(1,2,1); 

which_model = 'RLVU'; 

if (strcmp(which_model,'RLVU'))
    load('./results_exante/RLVU/exante_RLVU_batchInfo.mat', 'simperf')
    h=histogram(simperf*100,[ 0.5:0.04:0.8] *100); hold on;
    h.FaceColor = [160 63 116]/255;
end
if (strcmp(which_model,'BMBU'))
    load('./results_exante/BMBU/exante_BMBU_batchInfo.mat', 'simperf')
    h=histogram(simperf*100,[ 0.5:0.04:0.8] *100); hold on;
    h.FaceColor = [86 137 56]/255;
end

h.EdgeColor = 'none'; h.FaceAlpha = 0.3; set(gca,'fontsize', 20); xlabel('Task performance (%)'); YA = ylim; 
human_min = min(perf_sub)*100; human_max = max(perf_sub)*100; 
plot([human_min human_min], [0 15],'k--','linew',2 )
plot([human_max human_max], [0 15],'k--','linew',2 )


h=histogram(perf_sub*100,[ 0.5:0.04:0.8] *100, 'linew',3);hold on; 
h.FaceColor = 'none';h.EdgeColor = 'k';



subplot(1,2,2); 

which_model = 'BMBU'; 

if (strcmp(which_model,'RLVU'))
    load('./results_exante/RLVU/exante_RLVU_batchInfo.mat', 'simperf')
    h=histogram(simperf*100,[ 0.5:0.04:0.8] *100); hold on;
    h.FaceColor = [160 63 116]/255;
end
if (strcmp(which_model,'BMBU'))
    load('./results_exante/BMBU/exante_BMBU_batchInfo.mat', 'simperf')
    h=histogram(simperf*100,[ 0.5:0.04:0.8] *100); hold on;
    h.FaceColor = [86 137 56]/255;
end

h.EdgeColor = 'none'; h.FaceAlpha = 0.3; set(gca,'fontsize', 20); xlabel('Task performance (%)'); YA = ylim; 
human_min = min(perf_sub)*100; human_max = max(perf_sub)*100; 
plot([human_min human_min], [0 15],'k--','linew',2 )
plot([human_max human_max], [0 15],'k--','linew',2 )


h=histogram(perf_sub*100,[ 0.5:0.04:0.8] *100, 'linew',3);hold on; 
h.FaceColor = 'none';h.EdgeColor = 'k';

% Check out these values
% mean(perf_sub) = 0.6785
 % min[0.6065] max[0.7394] 
 
 % World -updating [max - min] = [0.7456 - 0.5948]= 0.1508
 % Value -updating  [ 0.65 learning rate]            = [0.7801 - 0.5636] = 0.2165
 %                  [ 0.4 learning rate]                [0.7798 - 0.5811 ] = 0.1987
 %                  [ 0.35 learning rate]              [0.7783 - 0.5863 ] = 0.1920
 %                 [ 0.3 learning rate]                [0.7795 - 0.5895] = 0.1900
  %                 
