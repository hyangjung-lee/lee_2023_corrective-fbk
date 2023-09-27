
%% Covering Figure 2,4,5
% Figure 2D
% Figure 2C
% Figure 4A
% Figure 5A

%% Figure 2D; example retrospective and prospective psychometric curves conditioned to toi=[M, small, correct] and [M->LARGE -> correct]
clear all; close all; 
cd ../
S=[-2:2]; 
numBoot = 1;
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment

options.fixedPars = NaN(1,5);
options.fixedPars(3) = 0; options.fixedPars(4) = 0;
options.fixedPars(5) = 0; % Eta;

path_str.currpath = pwd; 
path_str.data       = [path_str.currpath, '/data']; 
addpath(genpath([path_str.currpath, '/analysis_func']))

%%% load data sorted by episode types
load([path_str.data, '/Episode_freq_OBS.mat'])
        

%%% How data is stored (20 toi episodes defined for conditioning)
%%% toiplusone_Pcub.EpisodeCode.SMALL(LARGE) made at toi
%%% : psychometric data points (toi+1) when conditioned to a trial of interest (toi) 
%%%%%% for example, toiplusone_Pcub.extPure{2}.SMALL(:,:,:,iSub): 5-by-2-by-2  (S_{toi+1}) -by- (C_{toi+1}=SM,C_{toi+1}=LG) -by- (F_{toi+1}=Incor, F_{toi+1}=Corr)

%%%%%%%%%%% extPure{2}.SMALL: (toi+1) behavior when conditioned to toi=[XS -> SMALL -> correct]
%%%%%%%%%%% norPure{2}.SMALL: (toi+1) behavior when conditioned to toi=[S -> SMALL -> correct]
%%%%%%%%%%% extPure{6}.SMALL: (toi+1) behavior when conditioned to toi=[M -> SMALL -> correct]
%%%%%%%%%%% norPure{4}.SMALL: (toi+1) behavior when conditioned to toi=[L -> SMALL -> correct]
%%%%%%%%%%% extPure{4}.SMALL: (toi+1) behavior when conditioned to toi=[XL -> SMALL -> correct]

%%%%%%%%%%% extPure{3}.SMALL: (toi+1) behavior when conditioned to toi=[XS -> SMALL -> incorrect]
%%%%%%%%%%% norPure{3}.SMALL: (toi+1) behavior when conditioned to toi=[S -> SMALL -> incorrect]
%%%%%%%%%%% extPure{7}.SMALL: (toi+1) behavior when conditioned to toi=[M -> SMALL -> incorrect]
%%%%%%%%%%% norPure{5}.SMALL: (toi+1) behavior when conditioned to toi=[L -> SMALL -> incorrect]
%%%%%%%%%%% extPure{5}.SMALL: (toi+1) behavior when conditioned to toi=[XL -> SMALL -> incorrect]


%%%%%%%%%%% extPure{2}.LARGE: (toi+1) behavior when conditioned to toi=[XL -> LARGE -> correct]
%%%%%%%%%%% norPure{2}.LARGE: (toi+1) behavior when conditioned to toi=[L -> LARGE -> correct]
%%%%%%%%%%% extPure{6}.LARGE: (toi+1) behavior when conditioned to toi=[M -> LARGE -> correct]
%%%%%%%%%%% norPure{4}.LARGE: (toi+1) behavior when conditioned to toi=[S -> LARGE -> correct]
%%%%%%%%%%% extPure{4}.LARGE: (toi+1) behavior when conditioned to toi=[XS -> LARGE -> correct]

%%%%%%%%%%% extPure{3}.LARGE: (toi+1) behavior when conditioned to toi=[XL -> LARGE -> incorrect]
%%%%%%%%%%% norPure{3}.LARGE: (toi+1) behavior when conditioned to toi=[L -> LARGE -> incorrect]
%%%%%%%%%%% extPure{7}.LARGE: (toi+1) behavior when conditioned to toi=[M -> LARGE -> incorrect]
%%%%%%%%%%% norPure{5}.LARGE: (toi+1) behavior when conditioned to toi=[S -> LARGE -> incorrect]
%%%%%%%%%%% extPure{5}.LARGE: (toi+1) behavior when conditioned to toi=[XS -> LARGE -> incorrect]


indxChoiceL = 2; 
numSub = size(toiplusone_Pcub.extPure{2}.SMALL, 4); 
for which_toi = [-1, 1]        
        
    if (which_toi == -1)
        CUB = toiminusone_Pcub;
    elseif (which_toi == +1)
        CUB = toiplusone_Pcub;
    end       
    nT = NaN(numSub, length(S)); 
    nC = nT; psychmet = nT; LG_nT = nT; LGpsychmet = nT; 
    for iSub = 1 : numSub
        numTrials = sum(sum( CUB.extPure{6}.SMALL(:,:,:,iSub)      , 3), 2);
        numChoiceL = sum(CUB.extPure{6}.SMALL(:,indxChoiceL,:,iSub), 3);


        nT(iSub,:) = numTrials; 
        nC(iSub,:) = numChoiceL; 
        psychmet(iSub,:) = numChoiceL./numTrials; 



        numTrials = sum(sum( CUB.extPure{6}.LARGE(:,:,:,iSub)      , 3), 2);
        numChoiceL = sum(CUB.extPure{6}.LARGE(:,indxChoiceL,:,iSub), 3);

        LG_nT(iSub,:) = numTrials; 

        LGpsychmet(iSub,:) = numChoiceL./numTrials; 

    end
    
    
figure(22 + which_toi); clf; 
    plot(-2:2,   mean(psychmet),'o','markersize', 13, 'color', [0.7 0.7 0.7], 'linew',2); hold on; 
    errorbar(-2:2,  mean(psychmet), std(psychmet)/sqrt(30),'o','markersize', 13, 'color', [0.7 0.7 0.7], 'linew',2);

    plot(-2:2,   mean(LGpsychmet),'o','markersize', 13, 'color', [0 0 0], 'linew',2);
    errorbar(-2:2,   mean(LGpsychmet), std(LGpsychmet)/sqrt(30), 'o','markersize', 13, 'color', [0 0 0], 'linew',2);


    xlim([-3 3]); ylim([0 1]);
    set(gca,'fontsize', 25 , 'linew', 2, 'xtick', -2:2,'ytick', [0 0.5 1], 'tickLength', [0.02 0.02])
    plot([-3 3], [0.5 0.5], 'k--'); hold on;



    data(:,1) = S;
    data(:,2) = sum(nT,1).*mean(psychmet);
    data(:,3) = sum(nT,1);

    [~,stand_aftFit] = fit_psych_PSE(data,  options , numBoot ) ; 
    plot(-2.1:0.01:2.1, normcdf(-2.1:0.01:2.1, stand_aftFit(1), stand_aftFit(2)),'-','linew',2,'color', [0.7 0.7 0.7])


     data(:,1) = S;
     data(:,2) = sum(LG_nT,1).*mean(LGpsychmet);
     data(:,3) = sum(LG_nT,1);

    [~,stand_aftFit] = fit_psych_PSE(data,  options , numBoot ) ; 
    plot(-2.1:0.01:2.1, normcdf(-2.1:0.01:2.1, stand_aftFit(1), stand_aftFit(2)),'-','linew',2,'color', [0 0 0])



    if (which_toi == -1)
       xlabel('Stimulus_{\itt-1}');  
    elseif (which_toi == +1)
         xlabel('Stimulus_{\itt+1}');  
    end
end




%% Figure 2C cube
w=0.85;
u=[0 0;0 0;w w;w w];
v=[0 w;0 w;0 w;0 w];
z=[w w;0 0;0 0;w w];
s=[nan nan];
x=[u;s;v];
y=[v;s;u];
z=[z;s;w-z];

figure(25); clf; 
    axis([0 6 0 6 0 6]);axis off; axis square; view(142,30)

for i=0:4,for j=0:1,for k=0:1
surface(i+x,j+y,k+z,'facec','w', 'facealpha', 0.1)


if (i == 2 & j == 1  & k == 1 ) % (M) - Small choice - correct ; i = 0 (XL) - 2(M) - 4 (XS), j = 0 (d=small) , k = 1 (F=correct)
   surface(i+x,j+y,k+z,'facec','k','edgec','k', 'linew', 3) 
   
   
end

end,end,end

%% Figure 4A cube 
color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];
figure(41);clf; 
    axis([0 6 0 6 0 6]);axis off; axis square; view(142,30)

for i=0:4,for j=0:1,for k=0:1
surface(i+x,j+y,k+z,'facec','w', 'facealpha', 0.1)


if (i <= 2 & j == 1 & k == 1 )
   surface(i+x,j+y,k+z,'facec',color_yellow,'edgec',color_yellow, 'linew', 3) 
   
   
end

if (i >= 2 & j == 0 & k == 1 )
   surface(i+x,j+y,k+z,'facec',color_blue,'edgec',color_blue, 'linew', 3) 
   
   
end

end,end,end

%% Figure 5A cube

figure(51);clf;  
    axis([0 6 0 6 0 6]);axis off; axis square; view(142,30)

for i=0:4,for j=0:1,for k=0:1
surface(i+x,j+y,k+z,'facec','w', 'facealpha', 0.1)


if (i <= 2 & j == 0 & k == 0 )
   surface(i+x,j+y,k+z,'facec',color_blue,'edgec',color_blue, 'linew', 3) 
   
   
end

if (i >= 2 & j == 1 & k == 0 )
   surface(i+x,j+y,k+z,'facec',color_yellow,'edgec',color_yellow, 'linew', 3) 
   
   
end

end,end,end