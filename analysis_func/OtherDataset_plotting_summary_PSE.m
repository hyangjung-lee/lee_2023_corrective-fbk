%% Plotting PSE figures of other datsets 
function OtherDataset_plotting_summary_PSE(acSsimIndvsSam, acSsimIndvlSam, acSsimIndvsDiff, acSsimIndvlDiff, input_info)
color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];
save_fold= input_info.save_fold ;

mod_str = input_info.mod_str;

paramsetID = input_info.paramsetID ;
numSim = input_info.numSim ;

str_att = input_info.prefix;

if (input_info.dataind == 11)
  xcent = 50;  ycent = 50;  
   
   if (strcmp(str_att, 'subtract'))
       ycent = 0;  
   end
else
    xcent = 0; 
    ycent = 0 ;    
end

linealph = 1; 

xxd = input_info.S; % Urai et al., 2017 S =[    20    35    45    48    52    55    65    80];
Xxd = [ones(size(xxd')), xxd'];


sim_sSam = acSsimIndvsSam;
sim_lSam = acSsimIndvlSam; 
sim_sDiff = acSsimIndvsDiff; 
sim_lDiff = acSsimIndvlDiff; 


b_simAggsSam = nanmean(sim_sSam,1);
b_simAgglSam = nanmean(sim_lSam,1);  
  
axis_small = find(~isnan(b_simAggsSam)); axis_large = find(~isnan(b_simAgglSam)); len_sim = size(sim_sSam, 1); 

% b_simAggsSam = nanmean(sim_sSam,1);
% b_simAgglSam = nanmean(sim_lSam,1); 
b_simAggsDiff = nanmean(sim_sDiff,1); 
b_simAgglDiff = nanmean(sim_lDiff,1); 

e_sSam = nanstd(sim_sSam,1)/sqrt(len_sim);
e_lSam = nanstd(sim_lSam,1)/sqrt(len_sim); 
e_sDiff = nanstd(sim_sDiff,1)/sqrt(len_sim); 
e_lDiff = nanstd(sim_lDiff,1)/sqrt(len_sim); 


mk_type = 'o'; stp = 0;

if (mk_type == 'o')
    figure(33); clf; 
else
    figure(33);
end
 
    subplot(1,2,1);
            plot( [-116 116], [ycent ycent ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
            plot([xcent xcent], [-116 116], '-','color',[0.5 0.5 0.5]); hold on; 

            choiceColor = color_blue;
            h=plot(Xxd(axis_small,2)', b_simAggsSam(:, axis_small),'-', 'color', choiceColor , 'linewidth', 4); 
            h.Color(4) = linealph; 


            h=plot(Xxd(axis_small,2)' , b_simAggsSam(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
            errorbar(Xxd(axis_small,2)'+ stp, b_simAggsSam(:, axis_small),  e_sSam(:, axis_small),'color', choiceColor, 'linewidth', 2 ,'linestyle','none'); 

            choiceColor = color_yellow;
            h=plot(Xxd(axis_large,2)', b_simAgglSam(:, axis_large),'-', 'color', choiceColor, 'linewidth', 4); 
            h.Color(4) = linealph; 

            h=plot(Xxd(axis_large,2)' , b_simAgglSam(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
            errorbar(Xxd(axis_large,2)' + stp, b_simAgglSam(:, axis_large),  e_lSam(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

            set(gca,'fontsize', 25, 'xtick', -2:2, 'xticklabel',{'-2','-1','0','1','2'}, 'ticklength', [0.03 0.03]); xlabel('Stimulus_{\itt}_{\it-}_{\it1}') 


            xlim([-3 3]); ylim([-0.5 0.5]); set(gca,'linewidth', 2,'ticklength', [0.03 003]); set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
            xlabel('Stimulus_{\ittoi}') ; 

    subplot(1,2,2);
            plot( [-116 116], [ycent ycent ],'-','color', [0.7 0.7 0.7], 'linewidth', 0.8);  hold on;
            plot([xcent xcent], [-116 116], '-','color',[0.5 0.5 0.5]); hold on; 

            choiceColor = color_blue;

            h=plot(Xxd(axis_large,2)', b_simAggsDiff(:, axis_large),'-', 'color', choiceColor , 'linewidth', 4); 
            h.Color(4) = linealph; 
            h=plot(Xxd(axis_large,2)' , b_simAggsDiff(:, axis_large),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
            errorbar(Xxd(axis_large,2)'+ stp, b_simAggsDiff(:, axis_large),  e_sDiff(:, axis_large),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 

            choiceColor = color_yellow; 

            h=plot(Xxd(axis_small,2)',b_simAgglDiff(:, axis_small), '-', 'color', choiceColor, 'linewidth', 4); 
            h.Color(4) = linealph; 

            h=plot(Xxd(axis_small,2)' , b_simAgglDiff(:, axis_small),mk_type, 'color', choiceColor , 'linewidth', 2, 'markersize', 12); 
            errorbar(Xxd(axis_small,2)'+ stp, b_simAgglDiff(:, axis_small),  e_lDiff(:, axis_small),'color', choiceColor, 'linewidth', 2,'linestyle','none'); 


            xlim([-3 3]); set(gca,'fontsize', 25,'xtick', -2:2, 'ticklength', [0.03 0.03]); 
            ylim([-0.5 0.5]); set(gca,'linewidth', 2,'ticklength', [0.03 003]); set(gca, 'ytick', [-0.5 0 0.5], 'yticklabel', {'-','0','+'}) ;
            xlabel('Stimulus_{\ittoi}') 



    ylim([-0.7 0.7]);set(gcf,'position', [ 2887         334         761         339]); 

    %%%% polish for visual
    subplot(1,2,1); 
    yl = 1.2;
    set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)}); ylim([-yl yl]); 


    if (isempty(str_att))    
        subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off;
        subplot(1,2,2); ylabel('PSE_{\ittoi}_{\it +}_{\it1}'); box off; 
    elseif strcmp(str_att, 'SlowDrift_')
        subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
        subplot(1,2,2);ylabel('PSE_{\ittoi}_{\it -}_{\it1}'); box off;
    elseif strcmp(str_att, 'subtract')
        subplot(1,2,1); ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');
        subplot(1,2,2);  ylabel('PSE_{\ittoi}_{\it +}_{\it1} - PSE_{\ittoi}_{\it -}_{\it1}');
    end
    subplot(1,2,2); 
    set(gca,'ytick', [-yl,0, yl], 'yticklabel', {['-',num2str(yl)], '0', num2str(yl)}); ylim([-yl yl]); 


    if (strcmp( input_info.mod_str , 'urai') )
        subplot(1,2,1); xlim([10 90]); ylim([40 60]); set(gca,'xtick', [20 35 45 48 50 52 55 65 80], 'ytick', [40,50, 60], 'xticklabel', {'20','35','', '', '50', '', '','65', '80'}, 'yticklabel', {'40','50','60'}); 
        subplot(1,2,2); xlim([10 90]); ylim([40 60]);  set(gca,'xtick',[20 35 45 48 50 52 55 65 80], 'ytick', [40,50, 60], 'xticklabel', {'20','35','', '', '50', '', '','65', '80'}, 'yticklabel', {'40','50','60'});
       if (strcmp(str_att, 'subtract'))
            subplot(1,2,1); xlim([10 90]); ylim([-6 6]); set(gca,'xtick', [20 35 45 48 50 52 55 65 80], 'ytick', [-6,0, 6], 'xticklabel', {'20','35','', '', '50', '', '','65', '80'}, 'yticklabel', {'-6','0','6'}); 
            subplot(1,2,2); xlim([10 90]); ylim([-6 6]);  set(gca,'xtick',[20 35 45 48 50 52 55 65 80], 'ytick', [-6,0, 6], 'xticklabel', {'20','35','', '', '50', '', '','65', '80'}, 'yticklabel', {'-6','0','6'});

       end

    elseif (strcmp(input_info.mod_str , 'hachen')  || strcmp(input_info.mod_str , 'hachen_rat'))
        subplot(1,2,1); xlim([-5 5]); set(gca,'fontsize', 25, 'xtick', -4:4, 'xticklabel', {'-4','','-2','','0','','2','','4'});  xlim([-5 5])
        subplot(1,2,2); xlim([-5 5]);set(gca,'fontsize', 25, 'xtick', -4:4, 'xticklabel', {'-4','','-2','','0','','2','','4'});   xlim([-5 5])


    end
    fig_fold = [save_fold,'/figures']; mkdir(fig_fold)

    if exist( [fig_fold, '/', str_att, 'PSE_', mod_str,'_setID', num2str(input_info.iPpair),   '_', num2str(numSim),'.pdf'], 'file') ~= 2 
        save_figureToPDF(gcf,[fig_fold, '/', str_att, 'PSE_', mod_str,'_setID', num2str(input_info.iPpair),   '_', num2str(numSim)]);
    end

end