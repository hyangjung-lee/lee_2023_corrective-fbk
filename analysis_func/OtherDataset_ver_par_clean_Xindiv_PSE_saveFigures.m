%% Fit PSE for other datasets
function OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures( psychD, LG_psychD, condEpiRes, input_info)

WhichQuad = input_info.WhichQuad; 
numSim = input_info.numSim; 
str_att = input_info.str_att; 

SUB = WhichQuad;
save_fold = input_info.save_fold; 
iPpair = input_info.iPpair; 

mod_str = input_info.mod_str;
dataind = input_info.dataind;


paramsetID =   input_info.paramsetID  ;
     
     
     
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment

options.fixedPars = NaN(1,5);
options.fixedPars(3) = 0; options.fixedPars(4) = 0;
options.fixedPars(5) = 0; % Eta;

fit_mode = 0; 
if (fit_mode == 0)
   iBoot = [];  
end



S = input_info.S; 


for iSub = 1 : length(WhichQuad)
   subIDtoFit = WhichQuad(iSub);
   
   if (ismember(subIDtoFit, SUB )) 
   
   
for cond_order = 2 : 7
  
    insertAvg = sum(condEpiRes.extPure{cond_order}.SMALL.numTrials(:,:, 1:numSim), 3); nTri = insertAvg(iSub,:) ; 
    insertAvg = sum(condEpiRes.extPure{cond_order}.SMALL.numChoices(:,:, 1:numSim), 3); nCho = insertAvg(iSub,:) ;

    yyA = nCho./nTri;

    temp_idx = find(isnan(yyA) == 0);

    data = []; 

    data(:,1) = S(temp_idx);
    data(:,2) = nCho(temp_idx);
    data(:,3) = nTri(temp_idx);


    [aftFit, stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
    muEst =  aftFit(1);
    sigEst =  stand_aftFit(2) ;

    sigEstFitted.extPure(cond_order) = sigEst;
    pseFitted.extPure(cond_order) = muEst;

end




for cond_order = 2 : 7

    insertAvg = sum(condEpiRes.extPure{cond_order}.LARGE.numTrials(:,:, 1:numSim), 3); nTri = insertAvg(iSub,:)  ; 
    insertAvg = sum(condEpiRes.extPure{cond_order}.LARGE.numChoices(:,:, 1:numSim), 3); nCho = insertAvg(iSub,:) ;
    yyA = nCho./nTri;
    temp_idx = find(isnan(yyA) == 0);

    data = []; 

    data(:,1) = S(temp_idx);
    data(:,2) = nCho(temp_idx);
    data(:,3) = nTri(temp_idx);


    [aftFit, stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
    muEst =  aftFit(1);
    sigEst =  stand_aftFit(2) ;


    LG_sigEstFitted.extPure(cond_order) = sigEst;
    LG_pseFitted.extPure(cond_order) = muEst;
end



for cond_order = 2 : 7
   
    insertAvg = sum(condEpiRes.norPure{cond_order}.SMALL.numTrials(:,:, 1:numSim), 3); nTri = insertAvg(iSub,:) ;
    insertAvg = sum(condEpiRes.norPure{cond_order}.SMALL.numChoices(:,:, 1:numSim), 3); nCho = insertAvg(iSub,:) ;
   
    yyA = nCho./nTri;
    temp_idx = find(isnan(yyA) == 0);

    data = []; 

    data(:,1) = S(temp_idx);
    data(:,2) = nCho(temp_idx);
    data(:,3) = nTri(temp_idx);

    [aftFit, stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
    muEst =  aftFit(1);
    sigEst =  stand_aftFit(2) ;
       
    sigEstFitted.norPure(cond_order) = sigEst;
    pseFitted.norPure(cond_order) = muEst;

end




for cond_order = 2 : 7

    insertAvg = sum(condEpiRes.norPure{cond_order}.LARGE.numTrials(:,:, 1:numSim), 3); nTri = insertAvg(iSub,:)  ;  
    insertAvg = sum(condEpiRes.norPure{cond_order}.LARGE.numChoices(:,:, 1:numSim), 3); nCho = insertAvg(iSub,:) ;
    yyA = nCho./nTri;
    temp_idx = find(isnan(yyA) == 0);

    data = []; 
    data(:,1) = S(temp_idx);
    data(:,2) = nCho(temp_idx);
    data(:,3) = nTri(temp_idx);

    [aftFit, stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
    muEst =  aftFit(1);
    sigEst =  stand_aftFit(2) ;
   
    LG_sigEstFitted.norPure(cond_order) = sigEst;
    LG_pseFitted.norPure(cond_order) = muEst;
end

%%%% To acquire real PSE 

realPSE.extPure = pseFitted.extPure;
realPSE.norPure = pseFitted.norPure;

ydatOBJ1 = realPSE.extPure;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = realPSE.norPure ;  

if (dataind == 3 || dataind == 4)
    pseAbsouSame = [ydatOBJ1(:,4), ydatOBJ2(:,4),  ydatOBJ1(:,2),ydatOBJ2(:,2), ydatOBJ1(:,6) ,    NaN,             NaN,         NaN,              NaN];
    pseAbsouDiff = [NaN,                NaN,            NaN,            NaN,    ydatOBJ1(:,7) , ydatOBJ2(:,5), ydatOBJ1(:,5) , ydatOBJ2(:,3), ydatOBJ1(:,3) ];
elseif (dataind == 11)  
    pseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ2(:,6), ydatOBJ1(:,6) ,ydatOBJ1(:,4),ydatOBJ2(:,4), ydatOBJ2(:,4),ydatOBJ1(:,4) ];
    pseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ2(:,3),ydatOBJ1(:,3), ydatOBJ1(:,7),ydatOBJ2(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ];
else
    pseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ1(:,6) ,ydatOBJ2(:,4),ydatOBJ1(:,4) ];
    pseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ];  
end



ydatOBJ1 = sigEstFitted.extPure;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = sigEstFitted.norPure ;  

if (dataind == 3 || dataind == 4)
    SDpseAbsouSame = [ydatOBJ1(:,4), ydatOBJ2(:,4),  ydatOBJ1(:,2),ydatOBJ2(:,2), ydatOBJ1(:,6) ,    NaN,             NaN,         NaN,              NaN];
    SDpseAbsouDiff = [NaN,                NaN,            NaN,            NaN,    ydatOBJ1(:,7) , ydatOBJ2(:,5), ydatOBJ1(:,5) , ydatOBJ2(:,3), ydatOBJ1(:,3) ];
elseif (dataind == 11)
    SDpseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ2(:,6), ydatOBJ1(:,6) ,ydatOBJ1(:,4),ydatOBJ2(:,4), ydatOBJ2(:,4),ydatOBJ1(:,4) ];
    SDpseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ2(:,3),ydatOBJ1(:,3), ydatOBJ1(:,7),ydatOBJ2(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ];
else
    SDpseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ1(:,6) ,ydatOBJ2(:,4),ydatOBJ1(:,4) ];
    SDpseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ]; 
end



LG_realPSE.extPure = LG_pseFitted.extPure;
LG_realPSE.norPure = LG_pseFitted.norPure;

ydatOBJ1 = LG_realPSE.extPure ;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = LG_realPSE.norPure ;  


if (dataind == 3 || dataind == 4)
    LG_pseAbsouSame = [NaN,                NaN,            NaN,            NaN, ydatOBJ1(:,6), ydatOBJ2(:,2), ydatOBJ1(:,2) ,ydatOBJ2(:,4), ydatOBJ1(:,4) ];
    LG_pseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,5) ,ydatOBJ2(:,5),ydatOBJ1(:,7), NaN,                NaN,            NaN,            NaN  ] ;
elseif (dataind == 11)
    LG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ2(:,4),ydatOBJ1(:,4), ydatOBJ1(:,6) ,ydatOBJ2(:,6), ydatOBJ2(:,2),ydatOBJ1(:,2) ];
    LG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ2(:,7),ydatOBJ1(:,7) ,ydatOBJ1(:,3) , ydatOBJ2(:,3), ydatOBJ2(:,3),ydatOBJ1(:,3) ];
else
    LG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ1(:,6) ,ydatOBJ2(:,2),ydatOBJ1(:,2) ];
    LG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ1(:,7) ,ydatOBJ2(:,3),ydatOBJ1(:,3) ] ;  
end



ydatOBJ1 = LG_sigEstFitted.extPure ;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = LG_sigEstFitted.norPure ;  


if (dataind == 3 || dataind == 4)
    SDLG_pseAbsouSame = [NaN,                NaN,            NaN,            NaN, ydatOBJ1(:,6), ydatOBJ2(:,2), ydatOBJ1(:,2) ,ydatOBJ2(:,4), ydatOBJ1(:,4) ];
    SDLG_pseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,5) ,ydatOBJ2(:,5),ydatOBJ1(:,7), NaN,                NaN,            NaN,            NaN  ] ;
elseif (dataind == 11)
    SDLG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ2(:,4),ydatOBJ1(:,4), ydatOBJ1(:,6) ,ydatOBJ2(:,6), ydatOBJ2(:,2),ydatOBJ1(:,2) ];
    SDLG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ2(:,7),ydatOBJ1(:,7) ,ydatOBJ1(:,3) , ydatOBJ2(:,3), ydatOBJ2(:,3),ydatOBJ1(:,3) ];
else
    SDLG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ1(:,6) ,ydatOBJ2(:,2),ydatOBJ1(:,2) ];
    SDLG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ1(:,7) ,ydatOBJ2(:,3),ydatOBJ1(:,3) ] ;
end
%%%

acSsimIndvsSam(iSub,:) = pseAbsouSame;
acSsimIndvsDiff(iSub,:)= pseAbsouDiff;
acSsimIndvlSam(iSub,:) = LG_pseAbsouSame;
acSsimIndvlDiff(iSub,:) = LG_pseAbsouDiff;



SDacSsimIndvsSam(iSub,:) = SDpseAbsouSame;
SDacSsimIndvsDiff(iSub,:)= SDpseAbsouDiff;
SDacSsimIndvlSam(iSub,:) = SDLG_pseAbsouSame;
SDacSsimIndvlDiff(iSub,:) = SDLG_pseAbsouDiff;


end
end



%%
if (~isempty(str_att))
    file_to_save = strcat( save_fold,'/backward_pse', '/',str_att, 'AcSsim_stats_', mod_str,'_setID', num2str(iPpair),   '_', num2str(numSim), '.mat'); 
else
    file_to_save = strcat( save_fold,'/forward_pse', '/',str_att, 'AcSsim_stats_', mod_str,'_setID', num2str(iPpair),   '_', num2str(numSim), '.mat'); 
end
save(file_to_save, 'acSsimIndvsSam','acSsimIndvsDiff','acSsimIndvlSam','acSsimIndvlDiff', 'SDacSsimIndvsSam', 'SDacSsimIndvsDiff','SDacSsimIndvlSam', 'SDacSsimIndvlDiff')

input_info.save_fold = save_fold;
input_info.mod_str = mod_str;
input_info.paramsetID = paramsetID;
input_info.numSim = numSim;
input_info.prefix = str_att;
input_info.iPpair = iPpair; 

OtherDataset_plotting_summary_PSE(acSsimIndvsSam, acSsimIndvlSam, acSsimIndvsDiff, acSsimIndvlDiff, input_info)

if (~isempty(str_att))
    file_to_savePre = strcat( save_fold, '/forward_pse' , '/', 'AcSsim_stats_', mod_str,'_setID', num2str(iPpair),   '_', num2str(numSim), '.mat'); 
load(file_to_savePre)
    raw_acSsimIndvsSam = acSsimIndvsSam;
    raw_acSsimIndvlSam = acSsimIndvlSam; 
    raw_acSsimIndvsDiff = acSsimIndvsDiff; 
    raw_acSsimIndvlDiff = acSsimIndvlDiff; 



    SDraw_acSsimIndvsSam = SDacSsimIndvsSam;
    SDraw_acSsimIndvlSam = SDacSsimIndvlSam; 
    SDraw_acSsimIndvsDiff = SDacSsimIndvsDiff; 
    SDraw_acSsimIndvlDiff = SDacSsimIndvlDiff; 



    file_to_subtract = strcat( save_fold,'/backward_pse', '/', str_att, 'AcSsim_stats_', mod_str,'_setID', num2str(iPpair),   '_', num2str(numSim), '.mat'); 
load(file_to_subtract)
    sSam = acSsimIndvsSam;
    lSam = acSsimIndvlSam; 
    sDiff = acSsimIndvsDiff; 
    lDiff = acSsimIndvlDiff; 

    SDsSam = SDacSsimIndvsSam;
    SDlSam = SDacSsimIndvlSam; 
    SDsDiff = SDacSsimIndvsDiff; 
    SDlDiff = SDacSsimIndvlDiff; 

    acSsimIndvsSam = raw_acSsimIndvsSam -  sSam   ; 
    acSsimIndvlSam = raw_acSsimIndvlSam -   lSam ; 
    acSsimIndvsDiff = raw_acSsimIndvsDiff -  sDiff ; 
    acSsimIndvlDiff = raw_acSsimIndvlDiff -  lDiff ; 


    SDacSsimIndvsSam = (SDraw_acSsimIndvsSam + SDsSam ) /2;
    SDacSsimIndvlSam = (SDraw_acSsimIndvlSam + SDlSam ) /2;
    SDacSsimIndvsDiff = (SDraw_acSsimIndvsDiff + SDsDiff ) /2;
    SDacSsimIndvlDiff = (SDraw_acSsimIndvlDiff + SDlDiff ) /2;


    input_info.prefix = 'subtract';
    OtherDataset_plotting_summary_PSE(acSsimIndvsSam, acSsimIndvlSam, acSsimIndvsDiff, acSsimIndvlDiff, input_info)


%%% the folder 'subtracted' should be created 

    save( strcat(save_fold, '/', 'subtracted/', 'paramsetID',  num2str(iPpair),   '_', num2str(numSim), '.mat'), 'acSsimIndvsSam','acSsimIndvsDiff','acSsimIndvlSam','acSsimIndvlDiff', 'SDacSsimIndvsSam', 'SDacSsimIndvsDiff','SDacSsimIndvlSam', 'SDacSsimIndvlDiff')

 
end
end
