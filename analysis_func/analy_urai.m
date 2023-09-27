%% Prepare Urai et al's data 
%{
clear all; 
cd ../

filelist = dir(['./data/other_datasets/Urai_2017_Ncomm/raw_downloaded/CSV/' , '*_sj*.csv']); 
s1_ref = 70; 

for iSub = 1 : length(filelist)
    filename =  [filelist(iSub).folder,'/', filelist(iSub).name]; % '/Users/eva/Downloads/CSV/2ifc_data_sj01.csv'

    opts = detectImportOptions(filename); 
    opts.SelectedVariableNames  = {'stim', 'coherence', 'resp', 'correct', 'rt', 'trialnr', 'blocknr','sessionnr'};
    T = readtable(filename, opts); T = table2array(T);
    if (iSub == 15)
       A = T;  
    end

    stimIdent = T(:,1); % column 'stim'
    unsignedCoh = T(:,2)*100 ;  % column'coherence'
    seqS = NaN(size(stimIdent)); 
    seqS(stimIdent ==  -1) =  s1_ref - unsignedCoh(stimIdent ==  -1); % when s2 < s1 (i.e. s2 is smaller than the ref, which is s1). 
    seqS(stimIdent ==  1) = s1_ref + unsignedCoh(stimIdent ==  1); % when s2 > s1
    seqC = T(:,3); % column 'resp'
    seqF = stimIdent;  % column 'stimulus identity' - class variable
    seqRT = T(:,5); % column rt

    seq_trialnr = T(:, 6);
    seq_blocknr = T(:, 7);
    seq_sessionnr = T(:, 8);

    mainsess_type = unique(seq_sessionnr); 

    nCond = length(unique(seq_sessionnr)); % which is block in Urai et al (2017): this can be  different across subjects

    for iFb = 1 : nCond
    block_type =  unique( seq_blocknr(  (  seq_sessionnr == mainsess_type(iFb)  )  ) ) ; 
    numRun = length( block_type ) ; 
        for iRun = 1 : numRun 
            idx =  find ( seq_sessionnr  ==  mainsess_type(iFb)  & seq_blocknr == block_type( iRun) ); % split trials for each block , and for each session 
            subjSmat{iSub, iFb}{iRun} = seqS ( idx, : );
            subjCLmat{iSub, iFb}{iRun} = seqF  ( idx, : ); % equiv. class variable
            subjCmat{iSub, iFb}{iRun} = seqC  (idx, : );
            subjRTmat{iSub, iFb}{iRun} = seqRT  ( idx, : );

            subjNtrmat{iSub, iFb}(:, iRun ) =  length( idx); 
            subjNsesmat(iSub) = nCond; 
        end
    end
end

% save('ExpRDK_Data_Urai.mat', 'subjSmat', 'subjCLmat','subjCmat','subjRTmat','subjNtrmat','subjNsesmat')

%}
%% Episode Analyses on Urai et al's data
clear all; 
cd ../

path_home = pwd; 
addpath([path_home,'/lib']); addpath(genpath([path_home,'/analysis_func'])); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_info.histdir = 'prosp'; % run first
% input_info.histdir = 'retro'; % run subsequently 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_fold = [path_home, '/data/temp']; % Results to be saved in this folder 


if strcmp(input_info.histdir , 'prosp')
    str_att = []; 
elseif strcmp(input_info.histdir , 'retro')
    str_att = 'SlowDrift_'; 
end

load([path_home, '/data/other_datasets/Urai_2017_Ncomm/Urai_Data_processed/ExpRDK_Data_Urai.mat']) ;

filelist = dir([path_home, '/data/other_datasets/Urai_2017_Ncomm/raw_downloaded/CSV/' , '*_sj*.csv']);

iSub = 1;
filename =  [filelist(iSub).folder,'/', filelist(iSub).name]; 
s1_ref = 70;
opts = detectImportOptions(filename); 
opts.SelectedVariableNames  = {'stim', 'coherence', 'resp', 'correct', 'rt', 'trialnr', 'blocknr','sessionnr'};
T = readtable(filename, opts); T = table2array(T);

stimIdent = T(:,1); % column 'stim'
unsignedCoh = T(:,2)*100 ;  % column'coherence'
seqS = NaN(size(stimIdent)); 
seqS(stimIdent ==  -1) =  s1_ref - unsignedCoh(stimIdent ==  -1); % when s2 < s1 (i.e. s2 is smaller than the ref, which is s1. 
seqS(stimIdent ==  1) = s1_ref + unsignedCoh(stimIdent ==  1); % when s2 > s1
S = [20, 35, 45, 48, 52, 55, 65, 80];
mapS = [ 20, 35,35, 45, 45, 48, 48, 52,52 ,55,55, 65,65, 80];



tauBackConsider = 5; 

imode = 0;
if imode == 0
    numSim = 1;
end

iSim = 1; modeGroup = 1; 

flagRT =  3;  nback = 1; 

numSub = size(subjSmat, 1); 
maxTrial= 50; 
WhichQuad = 1: numSub; 
for iG = 1 : modeGroup
    for subIDtoFit = 1:numSub

        nCond = subjNsesmat(subIDtoFit); 
        iSub = subIDtoFit; 
        for iFb = 1 :nCond
            numRun = length(subjSmat{iSub ,iFb}); 
            for iRun = 1 : numRun 
                MssMap_pre{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 
                MssMap_post{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 

                MccMap_pre{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 
                MccMap_post{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 

                MffMap_pre{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 

                McfMap_pre{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 
                McfMap_post{iSub}{iRun, iFb} = NaN(1, maxTrial - 1, iSim) ; 
               
                if imode == 0 
                    ss = subjSmat{iSub, iFb}{iRun}';  
                    ff = subjCLmat{iSub, iFb}{iRun}'; 
                    cc = subjCmat{iSub, iFb}{iRun}'; 
                end

                cf_tmp = (cc == ff); cf = double(cf_tmp);  
                idx = find(cf == 0);
                cf(idx) = -1; cf(isnan(ff)) = NaN ; 

if strcmp(input_info.histdir , 'prosp')
%%%% sorting (Prospective)
        pre = 1: length(subjCmat{iSub, iFb}{iRun}) - 1;
        post = 2:length(subjCmat{iSub, iFb}{iRun});
elseif strcmp(input_info.histdir , 'retro')
%%%%  sorting (Retrospective) - slow drift       
        pre = 2:length(subjCmat{iSub, iFb}{iRun});   
        post = 1: length(subjCmat{iSub, iFb}{iRun}) - 1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%
                 iback = 1 ;

                    MssMap_pre{iSub}{iRun, iFb}(:, pre, iSim) = ss(:, pre(:, 1:end-iback+1));
                    MssMap_post{iSub}{iRun, iFb}(:, pre, iSim) = ss(:, post(:, iback:end));      

                    MccMap_pre{iSub}{iRun,iFb}(:,pre,iSim)=cc(:,pre(:, 1:end-iback+1));
                    MccMap_post{iSub}{iRun,iFb}(:,pre,iSim)=cc(:,post(:, iback:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added
                    MffMap_pre{iSub}{iRun, iFb}(:,pre,iSim)=ff(:, pre(:, 1:end-iback+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    McfMap_pre{iSub}{iRun, iFb}(:,pre,iSim)=cf(:, pre(:, 1:end-iback+1));
                    McfMap_post{iSub}{iRun,iFb}(:,pre,iSim)=cf(:, post(:, iback:end));
                            
                  
                            
                
            end
        end
    end

    %%
    
   nStim = length(S); 
   for cond_order = 2 : 7
        psychD.extPure{cond_order} = NaN( numSub, nStim, numSim) ;
        psychD.norPure{cond_order} = NaN( numSub, nStim, numSim) ;
        LG_psychD.extPure{cond_order} = NaN( numSub, nStim, numSim) ;
        LG_psychD.norPure{cond_order} = NaN( numSub, nStim, numSim) ;
  
        condEpiRes.extPure{cond_order}.SMALL.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.SMALL.numTrials = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numTrials = NaN( numSub, nStim, numSim) ;
        
        condEpiRes.norPure{cond_order}.SMALL.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.SMALL.numTrials = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numTrials = NaN( numSub, nStim, numSim) ;

        
        
   end
%%

    for subIDtoFit = 1:numSub

        iSub = subIDtoFit
        for iback = 1:1

            ssMap_pre = []; ssMap_post = []; ccMap_pre=[]; ccMap_post=[]; cfMap_pre=[]; cfMap_post = []; crit_pre =[]; crit_post=[];

                     for iSim = 1 : numSim
                         
                         DD.ext  = cell(7,  1); LG_DD.ext = cell(7,1);  
                         DD.nor  = cell(7,  1); LG_DD.nor = cell(7,1);  
                            
                            for ii = 2 : 7                         
                                DD.ext{ii} = cell(1,nStim);  LG_DD.ext{ii} = cell(1,nStim);
                                DD.nor{ii} = cell(1,nStim);  LG_DD.nor{ii} = cell(1,nStim);
                            end
                            nCond = subjNsesmat(subIDtoFit); 
                            for iFb = 1 : nCond
                                numRun = length(subjSmat{iSub ,iFb});
                                     for iRun = 1 : numRun 
                                        ssMap_pre= MssMap_pre{iSub}{iRun,iFb}(:,:,iSim )';
                                        ssMap_post = MssMap_post{iSub}{iRun,iFb}(:,:,iSim )';

                                        ccMap_pre  = MccMap_pre{iSub}{iRun,iFb}(:,:,iSim)';
                                        ccMap_post = MccMap_post{iSub}{iRun,iFb}(:,:,iSim)'; 
                                        cfMap_pre  = McfMap_pre{iSub}{iRun,iFb}(:,:,iSim)';
                                        cfMap_post = McfMap_post{iSub}{iRun,iFb}(:,:,iSim)';
                                        
   
                                   if (flagRT == 3) 
                                        Sprev = ssMap_pre ;
                                        Chprev = ccMap_pre ; 
                                        Fbprev = cfMap_pre ; 


                                        range = fliplr(unique(seqS)'); 
                                        Chpost = ccMap_post ;  Chpost = (-1)*Chpost; %%%% correction ! to bring into the Lak et al (2020) %s1
                                        Spost = ssMap_post ; %%%%% correction ! to bring into the Lak et al (2020) %s
                                        Spost = Spost + 1000;
                                        for irg = 1 : length(range)
                                            Spost(Spost == (1000+ range(irg))) = mapS(irg);
                                        end
                                        %%%%%  
                                   end
                                                      
                                  

                                            %%% [1] b2 * strong-Switch-Condition 

                                            SMALLstrSwM.ext =[];
                                            LARGEstrSwM.ext = []; 
                                            SMALLstrSwM.nor =[];
                                            LARGEstrSwM.nor= []; 
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEstrSwM.extL = [];     LARGEstrSwM.extS = []; 
                                            LARGEstrSwM.norL = [];     LARGEstrSwM.norS = []; 
                                            SMALLstrSwM.extL = [];      SMALLstrSwM.extS = []; 
                                            SMALLstrSwM.norL = [];      SMALLstrSwM.norS = []; 

                                            % (1) XL(larege s1), Lch(%s1: s1 is larger than s2), incorrect -> S more 
                                            idx = find( Sprev == range(14) & Chprev == -1 &  Fbprev == -1); 
 
                                            LARGEstrSwM.ext = [LARGEstrSwM.ext; idx];   
                                            LARGEstrSwM.extL = idx(Chpost(LARGEstrSwM.ext) == 1);
                                            LARGEstrSwM.extS = idx(Chpost(LARGEstrSwM.ext) == -1);
                                            % (2) L, Lch, incorrect -> S more 
                                            idx = find( (Sprev == range(12) | Sprev == range(13) ) & Chprev == -1 &  Fbprev == -1); 
   
                                            LARGEstrSwM.nor = [LARGEstrSwM.nor; idx];         
                                            LARGEstrSwM.norL = idx(Chpost(LARGEstrSwM.nor) == 1);
                                            LARGEstrSwM.norS = idx(Chpost(LARGEstrSwM.nor) == -1);


                                            % (3) S, Sch, incorrect -> L more 
                                            idx = find( ( Sprev == range(2) | Sprev == range(3) ) & Chprev == +1 &  Fbprev == -1);
  
                                            SMALLstrSwM.nor = [SMALLstrSwM.nor; idx];      
                                            SMALLstrSwM.norL = idx(Chpost(SMALLstrSwM.nor) == 1);
                                            SMALLstrSwM.norS = idx(Chpost(SMALLstrSwM.nor) == -1);

                                            % (4) XS, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(1) & Chprev == +1 &  Fbprev == -1); 
      
                                            SMALLstrSwM.ext = [SMALLstrSwM.ext; idx];           
                                            SMALLstrSwM.extL = idx(Chpost(SMALLstrSwM.ext) == 1);
                                            SMALLstrSwM.extS = idx(Chpost(SMALLstrSwM.ext) == -1);
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 3; 

                                            for iS = 1 : nStim
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.ext); tmp_sensIDX = find(Spost(SMALLstrSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                 
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.nor); tmp_sensIDX = find(Spost(SMALLstrSwM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.ext); tmp_sensIDX = find(Spost(LARGEstrSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                 tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.nor); tmp_sensIDX = find(Spost(LARGEstrSwM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                            end

                                            %%% [2] b3 * strong-Stay-Condition 
                                            SMALLstrStyM.ext =[];
                                            LARGEstrStyM.ext =[];

                                            SMALLstrStyM.nor =[];
                                            LARGEstrStyM.nor =[];

                                            LARGEstrStyM.extL = [];     LARGEstrStyM.extS = []; 
                                            LARGEstrStyM.norL = [];     LARGEstrStyM.norS = []; 
                                            SMALLstrStyM.extL = [];      SMALLstrStyM.extS = []; 
                                            SMALLstrStyM.norL = [];      SMALLstrStyM.norS = []; 

                                            % (1) XL, Sch, correct -> S more 
                                            idx = find( Sprev == range(14) & Chprev == +1 &  Fbprev == 1);
     
                                            SMALLstrStyM.ext = [SMALLstrStyM.ext; idx];     
                                            SMALLstrStyM.extL = idx(Chpost(SMALLstrStyM.ext) == 1);
                                            SMALLstrStyM.extS = idx(Chpost(SMALLstrStyM.ext) == -1);

                                            % (2) L, Sch, correct -> S more 
                                            idx = find(  (Sprev == range(12) | Sprev == range(13) ) & Chprev == +1 &  Fbprev == 1); 
     
                                            SMALLstrStyM.nor = [SMALLstrStyM.nor; idx];    
                                            SMALLstrStyM.norL = idx(Chpost(SMALLstrStyM.nor) == 1);
                                            SMALLstrStyM.norS = idx(Chpost(SMALLstrStyM.nor) == -1);

                                            % (3) S, Lch, correct -> L more 
                                            idx = find( ( Sprev == range(2) | Sprev == range(3) ) & Chprev == -1 &  Fbprev == 1); 
   
                                            LARGEstrStyM.nor = [LARGEstrStyM.nor; idx];   
                                            LARGEstrStyM.norL = idx(Chpost(LARGEstrStyM.nor) == 1);
                                            LARGEstrStyM.norS = idx(Chpost(LARGEstrStyM.nor) == -1);

                                            % (4) XS, Lch, correct -> L more 
                                            idx = find( Sprev == range(1) & Chprev == -1 &  Fbprev == 1); 
    
                                            LARGEstrStyM.ext = [LARGEstrStyM.ext; idx];   
                                            LARGEstrStyM.extL = idx(Chpost(LARGEstrStyM.ext) == 1);
                                            LARGEstrStyM.extS = idx(Chpost(LARGEstrStyM.ext) == -1);
                                                                
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 4; 


                                            for iS = 1 : nStim                    
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.ext); tmp_sensIDX = find(Spost(SMALLstrStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.nor); tmp_sensIDX = find(Spost(SMALLstrStyM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];


                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.ext); tmp_sensIDX = find(Spost(LARGEstrStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.nor); tmp_sensIDX = find(Spost(LARGEstrStyM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                            end
                                            

                                            %%% [3] b4 * weak-Switch-Condition 
                                       
                                            SMALLweakSwM.ext=[];
                                            LARGEweakSwM.ext=[];

                                            SMALLweakSwM.nor=[];
                                            LARGEweakSwM.nor=[];

                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEweakSwM.extL = [];     LARGEweakSwM.extS = []; 
                                            LARGEweakSwM.norL = [];     LARGEweakSwM.norS = []; 
                                            SMALLweakSwM.extL = [];      SMALLweakSwM.extS = []; 
                                            SMALLweakSwM.norL = [];      SMALLweakSwM.norS = []; 

                                            % (1) XL, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(14) & Chprev == +1 &  Fbprev == -1); 
                                           
                                            SMALLweakSwM.ext = [SMALLweakSwM.ext; idx];      
                                            SMALLweakSwM.extL = idx(Chpost(SMALLweakSwM.ext) == 1);
                                            SMALLweakSwM.extS = idx(Chpost(SMALLweakSwM.ext) == -1);

                                            % (2) L, Sch, incorrect -> L more 
                                            idx = find(  (Sprev == range(12) | Sprev == range(13) ) & Chprev == +1 &  Fbprev == -1); 
                                    
                                            SMALLweakSwM.nor = [SMALLweakSwM.nor;idx];      
                                            SMALLweakSwM.norL = idx(Chpost(SMALLweakSwM.nor) == 1);
                                            SMALLweakSwM.norS = idx(Chpost(SMALLweakSwM.nor) == -1);



                                            % (3) S, Lch, incorrect -> S more 
                                            idx = find( ( Sprev == range(2) | Sprev == range(3) ) & Chprev == -1 &  Fbprev == -1);

                                            LARGEweakSwM.nor = [LARGEweakSwM.nor; idx]; 
                                            LARGEweakSwM.norL = idx(Chpost(LARGEweakSwM.nor) == 1);
                                            LARGEweakSwM.norS = idx(Chpost(LARGEweakSwM.nor) == -1);

                                            % (4) XS, Lch, incorrect -> S more 
                                            idx = find( Sprev == range(1) & Chprev == -1 &  Fbprev == -1);

                                            LARGEweakSwM.ext = [LARGEweakSwM.ext; idx];
                                            LARGEweakSwM.extL = idx(Chpost(LARGEweakSwM.ext) == 1);
                                            LARGEweakSwM.extS = idx(Chpost(LARGEweakSwM.ext) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 5; 


                                            for iS = 1 : nStim                                                                               
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.ext); tmp_sensIDX = find(Spost(SMALLweakSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.nor); tmp_sensIDX = find(Spost(SMALLweakSwM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.ext); tmp_sensIDX = find(Spost(LARGEweakSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];


                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.nor); tmp_sensIDX = find(Spost(LARGEweakSwM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];                     
                                            end


                                            %%% [4] b5 * weak-Stay-Condition 
                                            SMALLweakStyM.ext=[];
                                            LARGEweakStyM.ext=[];
                                            SMALLweakStyM.nor=[];
                                            LARGEweakStyM.nor=[];
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEweakStyM.extL = [];     LARGEweakStyM.extS = []; 
                                            LARGEweakStyM.norL = [];     LARGEweakStyM.norS = []; 
                                            SMALLweakStyM.extL = [];      SMALLweakStyM.extS = []; 
                                            SMALLweakStyM.norL = [];      SMALLweakStyM.norS = []; 

                                            % (1) XL, Lch, correct -> L more 
                                            idx = find( Sprev == range(14) & Chprev == -1 &  Fbprev == 1); 
        
                                            LARGEweakStyM.ext = [LARGEweakStyM.ext; idx]; 
                                            LARGEweakStyM.extL = idx(Chpost(LARGEweakStyM.ext) == 1);
                                            LARGEweakStyM.extS = idx(Chpost(LARGEweakStyM.ext) == -1);

                                            % (2) L, Lch, correct -> L more 
                                            idx = find( (Sprev == range(12) | Sprev == range(13) ) & Chprev == -1 &  Fbprev == 1); 

                                            LARGEweakStyM.nor = [LARGEweakStyM.nor; idx];  
                                            LARGEweakStyM.norL = idx(Chpost(LARGEweakStyM.nor) == 1);
                                            LARGEweakStyM.norS = idx(Chpost(LARGEweakStyM.nor) == -1);

                                            % (3) S, Sch, correct -> S more 
                                            idx = find( ( Sprev == range(2) | Sprev == range(3) ) & Chprev == +1 &  Fbprev == 1);
 
                                            SMALLweakStyM.nor = [SMALLweakStyM.nor; idx];      
                                            SMALLweakStyM.norL = idx(Chpost(SMALLweakStyM.nor) == 1);
                                            SMALLweakStyM.norS = idx(Chpost(SMALLweakStyM.nor) == -1);

                                            % (4) XS, Sch, correct -> S more 
                                            idx = find( Sprev == range(1) & Chprev == +1 &  Fbprev == 1);
   
                                            SMALLweakStyM.ext = [SMALLweakStyM.ext; idx];
                                            SMALLweakStyM.extL = idx(Chpost(SMALLweakStyM.ext) == 1);
                                            SMALLweakStyM.extS = idx(Chpost(SMALLweakStyM.ext) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 2; 


                                            for iS = 1 : nStim
                                                

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.ext); tmp_sensIDX = find(Spost(SMALLweakStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.nor); tmp_sensIDX = find(Spost(SMALLweakStyM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                   
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.ext); tmp_sensIDX = find(Spost(LARGEweakStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                                                             
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.nor); tmp_sensIDX = find(Spost(LARGEweakStyM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
      
                                            end
                                         
                                            %%% [5] b6 * medium-Switch-Condition 
  
                                            SMALLmediSwM.ext=[];     SMALLmediSwM.nor=[];
                                            LARGEmediSwM.ext=[];   LARGEmediSwM.nor=[];
                                            LARGEmediSwM.extL = [];     LARGEmediSwM.extS = []; 
                                            SMALLmediSwM.extL = [];      SMALLmediSwM.extS = []; 

                                            LARGEmediSwM.norL = [];     LARGEmediSwM.norS = []; 
                                            SMALLmediSwM.norL = [];      SMALLmediSwM.norS = []; 

                                            % (ext) LM, Sch, incorrect -> L more 
                                            idx = find(   ( Sprev == range(8) | Sprev == range(9) )  & Chprev == +1 &  Fbprev == -1);
      
                                            SMALLmediSwM.ext = [SMALLmediSwM.ext; idx];    
                                            SMALLmediSwM.extL = idx(Chpost(SMALLmediSwM.ext) == 1);
                                            SMALLmediSwM.extS = idx(Chpost(SMALLmediSwM.ext) == -1);

                                            % (nor) LLM, Sch, incorrect -> L more 
                                            idx = find(   ( Sprev == range(10) | Sprev == range(11) )  & Chprev == +1 &  Fbprev == -1);

                                            SMALLmediSwM.nor = [SMALLmediSwM.nor; idx];    
                                            SMALLmediSwM.norL = idx(Chpost(SMALLmediSwM.nor) == 1);
                                            SMALLmediSwM.norS = idx(Chpost(SMALLmediSwM.nor) == -1);   


                                            % (ext) SM, Lch, incorrect -> Smore 
                                            idx = find(    ( Sprev == range(6) | Sprev == range(7) )  & Chprev == -1 &  Fbprev == -1);
    
                                            LARGEmediSwM.ext = [LARGEmediSwM.ext; idx];    
                                            LARGEmediSwM.extL = idx(Chpost(LARGEmediSwM.ext) == 1);
                                            LARGEmediSwM.extS = idx(Chpost(LARGEmediSwM.ext) == -1);

                                            % (nor) SSM, Lch, incorrect -> Smore 
                                            idx = find(    ( Sprev == range(4) | Sprev == range(5) )  & Chprev == -1 &  Fbprev == -1);
  
                                            LARGEmediSwM.nor = [LARGEmediSwM.nor; idx];    
                                            LARGEmediSwM.norL = idx(Chpost(LARGEmediSwM.nor) == 1);
                                            LARGEmediSwM.norS = idx(Chpost(LARGEmediSwM.nor) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 7; 


                                            for iS = 1 : nStim
                                                
                                              
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediSwM.ext); tmp_sensIDX = find(Spost(SMALLmediSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediSwM.nor); tmp_sensIDX = find(Spost(SMALLmediSwM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediSwM.ext); tmp_sensIDX = find(Spost(LARGEmediSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
           
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediSwM.nor); tmp_sensIDX = find(Spost(LARGEmediSwM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                            end

                                            %%% [6] b7 * medium-Stay-Condition 
                                  
                                            SMALLmediStyM.ext=[];  LARGEmediStyM.ext=[];    SMALLmediStyM.nor=[];  LARGEmediStyM.nor=[];
                                            LARGEmediStyM.extL = [];     LARGEmediStyM.extS = []; 
                                            SMALLmediStyM.extL = [];      SMALLmediStyM.extS = []; 

                                            LARGEmediStyM.norL = [];     LARGEmediStyM.norS = []; 
                                            SMALLmediStyM.norL = [];      SMALLmediStyM.norS = []; 

                                            % (ext) LM, Lch, correct -> L more 
                                            idx = find(  ( Sprev == range(8) | Sprev == range(9) ) & Chprev == -1 &  Fbprev == 1);
                                            LARGEmediStyM.ext = [LARGEmediStyM.ext; idx];    
                                            LARGEmediStyM.extL = idx(Chpost(LARGEmediStyM.ext) == 1);
                                            LARGEmediStyM.extS = idx(Chpost(LARGEmediStyM.ext) == -1);

                                            % (ext) LLM, Lch, correct -> L more 
                                            idx = find(  ( Sprev == range(10) | Sprev == range(11) ) & Chprev == -1 &  Fbprev == 1);
                                            LARGEmediStyM.nor = [LARGEmediStyM.nor; idx];    
                                            LARGEmediStyM.norL = idx(Chpost(LARGEmediStyM.nor) == 1);
                                            LARGEmediStyM.norS = idx(Chpost(LARGEmediStyM.nor) == -1);

                                            % (ext) SM, Sch, correct -> Smore 
                                            idx = find(   ( Sprev == range(6) | Sprev == range(7) ) & Chprev == +1 &  Fbprev == 1);
                                            SMALLmediStyM.ext = [SMALLmediStyM.ext; idx];    
                                            SMALLmediStyM.extL = idx(Chpost(SMALLmediStyM.ext) == 1);
                                            SMALLmediStyM.extS = idx(Chpost(SMALLmediStyM.ext) == -1);


                                            % (nor) SSM, Sch, correct -> Smore 
                                            idx = find(   ( Sprev == range(4) | Sprev == range(5) ) & Chprev == +1 &  Fbprev == 1);
                                            SMALLmediStyM.nor = [SMALLmediStyM.nor; idx];    
                                            SMALLmediStyM.norL = idx(Chpost(SMALLmediStyM.nor) == 1);
                                            SMALLmediStyM.norS = idx(Chpost(SMALLmediStyM.nor) == -1);
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 6; 

                                            for iS = 1 : nStim 
                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediStyM.ext); tmp_sensIDX = find(Spost(SMALLmediStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediStyM.nor); tmp_sensIDX = find(Spost(SMALLmediStyM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediStyM.ext); tmp_sensIDX = find(Spost(LARGEmediStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediStyM.nor); tmp_sensIDX = find(Spost(LARGEmediStyM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                            end

                                     end
                            end
                            
                            
                                      
                                     for cond_order = 2 : 7
                                         for iS = 1 : nStim
                                            psychD.extPure{iG, cond_order}(iSub, iS,  iSim) = sum(DD.ext{cond_order}{iS} == 1)/length(DD.ext{cond_order}{iS}); % DD: episode variable
                                            condEpiRes.extPure{iG, cond_order}.SMALL.numTrials(iSub, iS,  iSim) = length(DD.ext{cond_order}{iS} );
                                            condEpiRes.extPure{iG, cond_order}.SMALL.numChoices(iSub, iS,  iSim) = sum(DD.ext{cond_order}{iS}  == 1);

                                            LG_psychD.extPure{iG, cond_order}(iSub, iS,  iSim) = sum(LG_DD.ext{cond_order}{iS} == 1)/length(LG_DD.ext{cond_order}{iS}); % DD: episode variable                                          
                                            condEpiRes.extPure{iG, cond_order}.LARGE.numTrials(iSub, iS, iSim) = length(LG_DD.ext{cond_order}{iS});
                                            condEpiRes.extPure{iG, cond_order}.LARGE.numChoices(iSub, iS,  iSim) = sum(LG_DD.ext{cond_order}{iS} == 1);


                                            psychD.norPure{iG, cond_order}(iSub, iS,  iSim) = sum(DD.nor{cond_order}{iS} == 1)/length(DD.nor{cond_order}{iS}); % DD: episode variable
                                            condEpiRes.norPure{iG, cond_order}.SMALL.numTrials(iSub, iS,  iSim) = length(DD.nor{cond_order}{iS});
                                            condEpiRes.norPure{iG, cond_order}.SMALL.numChoices(iSub, iS, iSim) = sum(DD.nor{cond_order}{iS} == 1);


                                            LG_psychD.norPure{iG, cond_order}(iSub, iS,  iSim) = sum(LG_DD.nor{cond_order}{iS} == 1)/length(LG_DD.nor{cond_order}{iS}); % DD: episode variable
                                            condEpiRes.norPure{iG, cond_order}.LARGE.numTrials(iSub, iS, iSim) = length(LG_DD.nor{cond_order}{iS});
                                            condEpiRes.norPure{iG, cond_order}.LARGE.numChoices(iSub, iS, iSim) = sum(LG_DD.nor{cond_order}{iS} == 1);

                                         end
                                     end

                     end



        end


                   

                
    end

    
end

% % save('human_ExpRDK_Urai_shouldtest', 'psychD','LG_psychD','condEpiRes','datRTExt', 'datRTNor')



dataind = 11; 
SUB = WhichQuad; 



% addpath( '/Users/eva/Desktop/FF_simul_psinifit/lib');
% addpath( '/Users/eva/Desktop/FF_simul_psinifit/psignifit-master');
mod_str = 'urai'
iPpair = 0; paramsetID = 0; 

mkdir(save_fold)
mkdir( [save_fold, '/', 'obsPsych/' ] ) 

mkdir( [save_fold, '/', 'working/' ] ) 
mkdir( [save_fold, '/', 'subtracted/' ] ) 
mkdir( [save_fold, '/', 'figures/' ] ) 
mkdir( [save_fold, '/', 'backward_pse/' ] ) 
mkdir( [save_fold, '/', 'forward_pse/' ] ) 

                  
    
input_info.WhichQuad = SUB; 
input_info.dataind = dataind; 

input_info.save_fold = save_fold;
input_info.numSim = numSim;   

input_info.mod_str = mod_str; 
input_info.iPpair = iPpair; 
input_info.paramsetID = paramsetID;

input_info.S = S;

if ( isempty( str_att ) )   
       epi_file_to_save = strcat( save_fold,'/', 'obsPsych/', 'sim_stats_Ingred', mod_str, '_setID', num2str(iPpair), '_', num2str(numSim),'.mat');
variable_type = 2; 
        parsave_files(epi_file_to_save, variable_type, psychD, LG_psychD, condEpiRes)
    
        input_info.str_att = str_att; 
        input_info.epi_file_to_save = epi_file_to_save; 

        OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures( psychD, LG_psychD, condEpiRes, input_info)
end

if (strcmp( str_att, 'SlowDrift_'))
        epi_file_to_save2 = strcat( save_fold,'/',  'obsPsych/', 'SlowDrift_sim_stats_Ingred', mod_str,  '_setID', num2str(iPpair), '_', num2str(numSim),'.mat');
variable_type = 2; 
        parsave_files(epi_file_to_save2,  variable_type, psychD, LG_psychD, condEpiRes)

        epi_file_to_save = epi_file_to_save2;
        
        input_info.str_att = str_att; 
        input_info.epi_file_to_save = epi_file_to_save;
          
        OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures( psychD, LG_psychD, condEpiRes, input_info)

end



