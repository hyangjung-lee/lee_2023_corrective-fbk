%% Outputs: bootstrapped PSE
%
clear all; 
cd ../
path_home = pwd; 
addpath([path_home,'/lib']); addpath(genpath([path_home,'/analysis_func'])); 

input_info.histdir = 'prosp';
input_info.histdir = 'retro';


if strcmp(input_info.histdir , 'prosp')
    load([path_home,'/analysis_func/pock/','humanEpisodeInfo_prosp.mat'])
    str_savefile =      [path_home,'/analysis_func/pock/tmp/forward'];
elseif strcmp(input_info.histdir , 'retro')
    
    load([path_home,'/analysis_func/pock/','humanEpisodeInfo_retro.mat'])
    str_savefile =      [path_home,'/analysis_func/pock/tmp/backward'];
end
lapgue_on = 0; %  mu, sigma only as free

 
numBoot = 5000;
numBoot = 4; % for test


BootData = [];
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment

options.fixedPars = NaN(1,5); options.fixedPars(3) = 0; options.fixedPars(4) = 0; options.fixedPars(5) = 0; % lapse, guess, Eta 0

if (lapgue_on == 1)
    options.expType = 'YesNo';
    options.fixedPars(5) = 0; % Eta;
end


fit_mode = 1; 
if (fit_mode == 0)
   iBoot = [];  
end



S=[-2:2]; 

for iSub = 1 : 30
    iSub
    bootdata = []; 
      for iBoot = 1 : numBoot 
          tic; 
            curveYfit  =[]; pseFitted = []; sigEstFitted = []; estPSE_d_small =[]; 
            realPSE = []; LG_realPSE =[]; 
            LG_sigEstFitted = []; LG_pseFitted =[]; LG_curveYfit = []; 
            for cond_order = 2 : 7
                data = []; 
                yyA = psychD.extPure{cond_order}(iSub,:);

                temp_idx = find(isnan(yyA) == 0);
                       nTri = condEpiRes.extPure{cond_order}.SMALL.numTrials(iSub,:);
                       nCho = condEpiRes.extPure{cond_order}.SMALL.numChoices(iSub,:);

                data(:,1) = S(temp_idx);
                data(:,2) = nCho(temp_idx);
                data(:,3) = nTri(temp_idx);

                [aftFit,stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
                muEst =  aftFit(1);
                sigEst =  stand_aftFit(2) ;

                sigEstFitted.extPure(cond_order) = sigEst;
                pseFitted.extPure(cond_order) = muEst;
            end

            for cond_order = 2 : 5
                data = []; 
                yyA = psychD.norPure{cond_order}(iSub,:);

                temp_idx = find(isnan(yyA) == 0);

                nTri = condEpiRes.norPure{cond_order}.SMALL.numTrials(iSub,:);
                nCho = condEpiRes.norPure{cond_order}.SMALL.numChoices(iSub,:);
                data(:,1) = S(temp_idx);
                data(:,2) = nCho(temp_idx);
                data(:,3) = nTri(temp_idx);

                [aftFit,stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
                muEst =  aftFit(1);
                sigEst =  stand_aftFit(2) ;

                sigEstFitted.norPure(cond_order) = sigEst;

                pseFitted.norPure(cond_order) = muEst;
            end

            realPSE.extPure = pseFitted.extPure;
            realPSE.norPure = pseFitted.norPure;

            ydatOBJ1 = realPSE.extPure;   % this would generate corrected episode effects for the effect of M-episode. 
            ydatOBJ2 = realPSE.norPure ;  

            pseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ1(:,6) ,ydatOBJ2(:,4),ydatOBJ1(:,4) ];

            pseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ];


            bootdata.sSam(1,:,iBoot) = pseAbsouSame;
            bootdata.sDiff(1,:,iBoot) = pseAbsouDiff;


            ydatOBJ1 = sigEstFitted.extPure;   % this would generate corrected episode effects for the effect of M-episode. 
            ydatOBJ2 = sigEstFitted.norPure ;  

            SDpseAbsouSame = [ydatOBJ1(:,2),ydatOBJ2(:,2),ydatOBJ1(:,6) ,ydatOBJ2(:,4),ydatOBJ1(:,4) ];

            SDpseAbsouDiff = [ydatOBJ1(:,3),ydatOBJ2(:,3),ydatOBJ1(:,7) ,ydatOBJ2(:,5),ydatOBJ1(:,5) ];

            bootdata.SD_sSam(1,:,iBoot) = SDpseAbsouSame;
            bootdata.SD_sDiff(1,:,iBoot) = SDpseAbsouDiff;





            for cond_order = 2 : 7
                data = []; 

                yyA = LG_psychD.extPure{cond_order}(iSub,:);

                temp_idx = find(isnan(yyA) == 0);


                nTri = condEpiRes.extPure{cond_order}.LARGE.numTrials(iSub,:);
                nCho = condEpiRes.extPure{cond_order}.LARGE.numChoices(iSub,:);
                data(:,1) = S(temp_idx);
                data(:,2) = nCho(temp_idx);
                data(:,3) = nTri(temp_idx);


                [aftFit,stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
                muEst =  aftFit(1);
                sigEst =  stand_aftFit(2) ;


                LG_sigEstFitted.extPure(cond_order) = sigEst;
                LG_pseFitted.extPure(cond_order) = muEst;
            end


            for cond_order = 2 : 5
                data = []; 
                yyA = LG_psychD.norPure{cond_order}(iSub,:);

                temp_idx = find(isnan(yyA) == 0);

                nTri = condEpiRes.norPure{cond_order}.LARGE.numTrials(iSub,:);
                nCho = condEpiRes.norPure{cond_order}.LARGE.numChoices(iSub,:);

                data(:,1) = S(temp_idx);
                data(:,2) = nCho(temp_idx);
                data(:,3) = nTri(temp_idx);

                [aftFit,stand_aftFit] = fit_psych_PSE(data,  options , iBoot ) ; 
                muEst =  aftFit(1);
                sigEst =  stand_aftFit(2) ;


                LG_sigEstFitted.norPure(cond_order) = sigEst;
                LG_pseFitted.norPure(cond_order) = muEst;
            end

LG_realPSE.extPure = LG_pseFitted.extPure;
LG_realPSE.norPure = LG_pseFitted.norPure;


ydatOBJ1 = LG_realPSE.extPure;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = LG_realPSE.norPure ;  

LG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ1(:,6) ,ydatOBJ2(:,2),ydatOBJ1(:,2) ];

LG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ1(:,7) ,ydatOBJ2(:,3),ydatOBJ1(:,3) ];

bootdata.lSam(1,:, iBoot) = LG_pseAbsouSame;
bootdata.lDiff(1,:, iBoot) = LG_pseAbsouDiff;


ydatOBJ1 = LG_sigEstFitted.extPure ;   % this would generate corrected episode effects for the effect of M-episode. 
ydatOBJ2 = LG_sigEstFitted.norPure ;  


SDLG_pseAbsouSame = [ydatOBJ1(:,4),ydatOBJ2(:,4),ydatOBJ1(:,6) ,ydatOBJ2(:,2),ydatOBJ1(:,2) ];

SDLG_pseAbsouDiff = [ydatOBJ1(:,5),ydatOBJ2(:,5),ydatOBJ1(:,7) ,ydatOBJ2(:,3),ydatOBJ1(:,3) ] ;


bootdata.SD_lSam(1,:, iBoot) = SDLG_pseAbsouSame;
bootdata.SD_lDiff(1,:, iBoot) = SDLG_pseAbsouDiff;



elapsed = toc; 
fprintf( '<< Sub# = %d, iBoot= %d', ...
    iSub, iBoot);
fprintf( ' : %.3f sec executed\n', elapsed );
      end

     parsave_boot(str_savefile, bootdata, iSub)
end

store_allsubj(str_savefile)
