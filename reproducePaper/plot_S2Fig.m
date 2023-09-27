%% S2A Fig
%% the mock-BMBU-simulation ;  figure for 6-trial examples
clear all; close all; 
cd ../
load('./results_demo/v2_allF1_conceptfig_sim_expostBMBU1.mat'); % load mock simulation examples
color_yellow = [0.9290 0.6940, 0.1250];

cmap_code{1} = color_yellow; 
cmap_code{2} = [255, 147, 0]/255; 
cmap_code{3} = [148, 33, 146]/255; 



iRun = 1;iSim = 1;
trialseq = 1: 5;



figure(121); clf; 
xLi = [0 trialseq(end)+1]; 
yLi = [-3 3]; 

plot([xLi], [0 0 ], 'k--') ; hold on; 
plot(xLi, 2.6*[1 1], 'k:', 'linew',2);
plot(xLi, -2.6*[1 1], 'k:', 'linew',2);

for iTrial = 1: length(trialseq)
    plot(iTrial*[1 1 ], [-2.5 2.5], 'k--')

    step_x = 0.1;
    md_x =  raw_Simdata{2}.S(iRun,iTrial,1); 
    plot( iTrial, md_x , 'x',  'markerfacecolor', 'k', 'markeredgecolor',0*[1 1 1 ], 'markersize', 12, 'linew',2 ) 

    step_m = 0.25;
    md_m = raw_Simdata{2}.m_rnd(iRun,iTrial,iSim); 
    plot( iTrial+ step_m, md_m , 'o',  'markerfacecolor', 'none', 'markeredgecolor',0*[1 1 1 ], 'markersize', 10) 

    step_ch = 0.25; ext_heigh = 2.6; 
    md_ch = raw_Simdata{2}.C(iRun,iTrial,iSim);

    plot( iTrial+ step_ch, ext_heigh*md_ch , 's',   'color','k', 'markersize', 13) 


    ch  = md_ch; 
    fb = raw_Simdata{2}.F(iRun,iTrial, 1);

    corincor = (ch == fb); 
    c_idx = find(corincor==1 & fb ~= 0); 
    in_idx = find(corincor==0 & fb ~= 0);
    none_idx = find(fb == 0); 


    if (~isempty(c_idx))

        plot( iTrial+ step_ch, ext_heigh*md_ch , 's',   'markerfacecolor','g', 'markeredgecolor','none', 'markersize', 13) 

    elseif (~isempty(in_idx))
          plot( iTrial+ step_ch, ext_heigh*md_ch , 's',   'markerfacecolor','r', 'markeredgecolor','none', 'markersize', 13) 

    elseif (~isempty(none_idx))


    end
    step_mm = 0.6;
    md_mm = raw_Simdata{2}.km_rnd(iRun,iTrial,iSim);


    color_km = [ cmap_code{1} ; cmap_code{1}; cmap_code{3}; cmap_code{1}; cmap_code{3}];
    plot(iTrial+ step_mm, md_mm , 'o', 'markerfacecolor', color_km(iTrial,:), 'markeredgecolor','none','markersize', 13)



    step_bd = 0.3;

    md_bd =  raw_Simdata{2}.pri(iRun, iTrial,1); 

    if (iTrial == 1)

        plot([iTrial-step_bd, iTrial+step_bd],  [md_bd, md_bd] , 'k-', 'linew',5)



    else
          plot([iTrial-step_bd, iTrial+step_bd],  [md_bd, md_bd] , 'k-', 'linew',5)

    end

end

md_bd =  raw_Simdata{2}.pri(iRun, iTrial+1,iSim) ; 
plot([iTrial+1-step_bd, iTrial+1+step_bd],  [md_bd, md_bd] , 'k-', 'linew',5)

    
whole_tc = raw_Simdata{2}.pri(iRun, 1: length(trialseq)+1,iSim);
    
    
    for iTrial = 1 : length(trialseq)
        its = iTrial +1 ;   
    end


box off; ylim(yLi); set(gca,'fontsize', 20, 'xtick', 1: 6);


%% Fig S2B
%% repeating the mock-BMBU-simulation using the same sequence 

figure(122); clf; 

color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];


iRun = 1;
trialseq = 1: 5;
xLi = [0 trialseq(end)+1]; 
yLi = [-3 3]; 

plot([xLi], [0 0 ], 'k--') ; hold on; 
plot(xLi, 2.6*[1 1], 'k:', 'linew',2);
plot(xLi, -2.6*[1 1], 'k:', 'linew',2);

for iTrial = 1: length(trialseq)
    plot(iTrial*[1 1 ], [-2.5 2.5], 'k--')

end
for iSim = 1:100
%     whole_tc = nanmean( raw_Simdata{2}.pri(iRun, 1: length(trialseq)+1,iSim), 3);
    whole_tc =  raw_Simdata{2}.pri(iRun, 1: length(trialseq)+1,iSim);
    h= plot( [1:length(trialseq)+1]-step_bd, whole_tc , '-k');
    h.Color(4) = 0.3;
    
    Z(iSim,:) = whole_tc; 
end


for iTrial = 1: length(trialseq)

    
    step_bd = 0.25;

    md_bd =  nanmean( raw_Simdata{2}.pri(iRun, iTrial, :), 3); 

    if (iTrial == 1)
        plot(iTrial- step_bd, md_bd , '.', 'color','k' , 'linew',2, 'markersize', 12)

    else

        plot(iTrial- step_bd, md_bd , '.', 'color','k' , 'linew',2, 'markersize', 12)

    end

    Q(:, iTrial) = md_bd; 
end
meansimm = nanmean( raw_Simdata{2}.pri(iRun, 1:length(trialseq)+1, :), 3);
plot([1:length(trialseq)+1]-step_bd, meansimm , '-', 'color','k' , 'linew',2);box off; 

ylim(yLi); xlim([0 length(trialseq)+1+0.5]);set(gca,'fontsize', 20);

